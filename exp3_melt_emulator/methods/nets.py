"""Neural per-cell emulators: deep ensemble M6 and physics + ML hybrid M7 (PyTorch).

Samples are shelf cells: the inputs are ``cell_features`` of every training run on its
own ``shelf_mask``, and the target is that run's melt there (m a-1, negative values kept).
Both methods predict cell by cell, so they work on any run's mask.

Training (shared by M6 and M7)
- Inner validation BY RUN: the fit's runs are permuted with ``split_seed``; the last run
  of the permutation is the inner-validation run. It is used only for early stopping and
  never for gradients, the input standardiser or the M7 initialisation.
- Input standardisation (mean, std; std 0 -> 1) is fitted on the remaining training cells.
  Targets are divided by ``y_scale`` (std of the training melt); M6 also subtracts the mean,
  M7 does not (that would break m = K TF+^n g). This only rescales the loss.
- Loss: Gaussian negative log-likelihood, Lakshminarayanan et al. 2017 Eq. 1,
  0.5 log sigma^2 + (y - mu)^2 / (2 sigma^2), with sigma^2 = softplus(.) + 1e-6
  (their footnote 2: softplus and a minimum variance of 1e-6), in scaled units.
- Adam (lr 1e-3), minibatches of ``batch_size`` (512) cells reshuffled every epoch by a torch
  Generator seeded with the member's seed. Each epoch, the area-weighted RMSE of the mean
  on the inner-validation run (``metrics.rmse_map``) is computed; training stops after
  ``patience`` epochs without improvement (or ``max_epochs``) and the best weights are
  restored. Model selection thus uses the area-weighted error, never the NLL.
- Ensemble: M members, random initialisation (``torch.manual_seed(seed)``) per member,
  all trained on the same cells, no bagging (Lakshminarayanan 2017 Sec. 2.4: "We observed
  that bagging deteriorated performance"). The prediction is the uniformly weighted
  Gaussian mixture, summarised by its mean and variance (Lakshminarayanan 2017 Sec. 2.4;
  these two formulas are unnumbered in the paper, whose only numbered equation is Eq. 1):
  mu* = M^-1 sum_m mu_m,  sigma*^2 = M^-1 sum_m (sigma_m^2 + mu_m^2) - mu*^2.
- Defaults (to be frozen in PD0): hidden (32, 64, 32) ReLU, Adam lr 1e-3, batch 512,
  max_epochs 1000, patience 250. The patience is long because on the real data the
  inner-validation RMSE has an early shallow minimum (epoch ~5-80, n ~ 1.8) and, for some
  members, a later deeper one (epoch ~300-650, n ~ 2.6); patience 30 stopped at the first.
  This was seen on the training runs of the final fold only (no test or reference run).
- Determinism: CPU, float64, ``torch.use_deterministic_algorithms(True)`` and a fixed
  thread count inside ``fit``/``predict`` (the global torch settings are restored after).

M7 (hybrid): m = K * b^n * g_theta(x_g), g_theta = softplus(MLP(x_g) + log(e - 1)) (so g = 1
for a zero MLP output), x_g = HYBRID_G_FEATURES (geometry and position relative to the
thermocline; no temperature variable). b = TF+ = max(TF, 0) ("local", the TF of M2,
Burgard 2022 Eq. 14) or b = (TF+ |<TF>|)^(1/2) ("nonlocal", ruling R3: the TF |<TF>| of M3,
Jourdain 2020 Eq. 1, raised to n/2 so both variants scale as alpha^n). Cells with TF <= 0
(refreezing) therefore have a mean of exactly 0. log K and log n are trainable scalars per
member, initialised with n0 = 2 and K0 from the closed-form M2 fit (local) or M3 fit
(nonlocal) on the training runs: at n = 2 and g = 1 the hybrid is K0 TF+^2 (resp. K0 TF+ |<TF>|),
the positive part of M2 (M3), so K0 = coef * c^2. Note K is degenerate with the overall scale of g (only K * g is
identified); n is identified by the cross-run variation of TF at fixed geometry.
The M7 variance head is the second output of the same g-network: sigma^2 = softplus(.) +
1e-6 of geometry only, so the variance never sees temperature either. (The alternative
sigma ~ mu would grow with TF in extrapolation, but with the NLL it biases the mean low,
by ~8 % for a 30 % relative spread, and leaves sigma ~ 0 in the TF+ = 0 cells.)
The choice local / nonlocal (``tf_kind="auto"``) is made by inner cross-validation by run
on the fit's runs: ``tf_kind_folds`` groups (default 5; None = one per run, i.e. full
leave-one-run-out, ~5x the cost) from a seeded permutation (seed split_seed + 1); criterion = mean ``rmse_map`` over the held-out runs; ties ->
local. Each inner fit is a complete Hybrid fit with the same hyperparameters.
"""
from __future__ import annotations

import copy
import math
from contextlib import contextmanager

import numpy as np
import torch
from torch import nn

from exp3_melt_emulator.data import FEATURES, Run, cell_features
from exp3_melt_emulator.methods.physics import C, QuadLocal, QuadNonLocal
from exp3_melt_emulator.metrics import rmse_map

NET_FEATURES = FEATURES
HYBRID_G_FEATURES = ("draft", "draft_minus_zc", "dist_gl", "dist_front", "slope", "s")
_TEMPERATURE = {"T", "S", "TF", "TF_mean", "Tbot"}
assert not _TEMPERATURE & set(HYBRID_G_FEATURES)

VAR_FLOOR = 1e-6           # Lakshminarayanan 2017 footnote 2 (minimum variance)
G_OFFSET = math.log(math.e - 1.0)   # softplus(G_OFFSET) = 1
N0 = 2.0                   # initial exponent (M2 / M3 form)
DTYPE = torch.float64


@contextmanager
def _deterministic(threads):
    prev_threads = torch.get_num_threads()
    prev_det = torch.are_deterministic_algorithms_enabled()
    torch.set_num_threads(threads)
    torch.use_deterministic_algorithms(True)
    try:
        yield
    finally:
        torch.use_deterministic_algorithms(prev_det)
        torch.set_num_threads(prev_threads)


class _Scaler:
    def __init__(self, X):
        self.mean = X.mean(0)
        sd = X.std(0)
        self.std = np.where(sd > 0, sd, 1.0)

    def __call__(self, X):
        return (X - self.mean) / self.std


def _mlp(n_in, hidden, n_out):
    layers, d = [], n_in
    for h in hidden:
        layers += [nn.Linear(d, h), nn.ReLU()]
        d = h
    layers.append(nn.Linear(d, n_out))
    return nn.Sequential(*layers)


class _MeanVarNet(nn.Module):
    """M6 member: MLP -> (mu, sigma^2 = softplus + floor), in scaled target units."""

    def __init__(self, n_in, hidden):
        super().__init__()
        self.mlp = _mlp(n_in, hidden, 2)

    def forward(self, x):
        out = self.mlp(x)
        return out[:, 0], nn.functional.softplus(out[:, 1]) + VAR_FLOOR


class _HybridNet(nn.Module):
    """M7 member: mu = K b^n g(x_g) / y_scale, sigma^2 = softplus(h(x_g)) + floor."""

    def __init__(self, n_in, hidden, K0, y_scale):
        super().__init__()
        self.mlp = _mlp(n_in, hidden, 2)
        self.log_K = nn.Parameter(torch.tensor(math.log(K0), dtype=DTYPE))
        self.log_n = nn.Parameter(torch.tensor(math.log(N0), dtype=DTYPE))
        self.y_scale = y_scale

    def forward(self, xg, b):
        out = self.mlp(xg)
        g = nn.functional.softplus(out[:, 0] + G_OFFSET)
        pos = b > 0
        logb = torch.log(torch.where(pos, b, torch.ones_like(b)))
        power = torch.where(pos, torch.exp(self.log_K + torch.exp(self.log_n) * logb),
                            torch.zeros_like(b))
        return power * g / self.y_scale, nn.functional.softplus(out[:, 1]) + VAR_FLOOR


def _nll(mu, var, y):
    # Lakshminarayanan 2017 Eq. 1 (constant dropped), mean over the batch
    return torch.mean(0.5 * torch.log(var) + 0.5 * (y - mu) ** 2 / var)


class _Ensemble:
    """Shared fit/predict machinery; subclasses define inputs and the member network."""
    name = ""
    input_features: tuple = ()

    def __init__(self, n_members=5, seeds=(0, 1, 2, 3, 4), hidden=(32, 64, 32),
                 max_epochs=1000, patience=250, lr=1e-3, batch_size=512, split_seed=0,
                 threads=1):
        seeds = tuple(int(s) for s in seeds)
        if len(seeds) != n_members:
            raise ValueError(f"need one seed per member: {n_members} members, {len(seeds)} seeds")
        self.n_members, self.seeds, self.hidden = n_members, seeds, tuple(hidden)
        self.max_epochs, self.patience, self.lr = int(max_epochs), int(patience), float(lr)
        self.batch_size, self.split_seed, self.threads = int(batch_size), int(split_seed), int(threads)
        self.nets: list[nn.Module] = []
        self.info: dict = {}
        self._fitted = False

    def hparams(self):
        return dict(n_members=self.n_members, seeds=list(self.seeds), hidden=list(self.hidden),
                    max_epochs=self.max_epochs, patience=self.patience, lr=self.lr,
                    batch_size=self.batch_size, split_seed=self.split_seed, threads=self.threads,
                    optimiser="Adam", var_floor=VAR_FLOOR, dtype="float64")

    # --- hooks -------------------------------------------------------------------------
    def _prepare(self, train, pooled, y):     # sets y_shift, y_scale and any physics init
        raise NotImplementedError

    def _tensors(self, cells):                 # dict of 1-D cell arrays -> tuple of tensors
        raise NotImplementedError

    def _new_net(self):
        raise NotImplementedError

    def _member_info(self, net):
        return {}

    # --- data --------------------------------------------------------------------------
    @staticmethod
    def _cells(run):
        f = cell_features(run)
        return {k: v[run.shelf_mask] for k, v in f.items()}

    def _X(self, cells):
        return np.stack([np.asarray(cells[k], float) for k in self.input_features], 1)

    # --- fit ---------------------------------------------------------------------------
    def fit(self, runs: list[Run]) -> None:
        runs = list(runs)
        if len(runs) < 2:
            raise ValueError(f"{self.name}: fit needs at least 2 runs (one is the inner-validation run)")
        self._fitted = False
        perm = np.random.default_rng(self.split_seed).permutation(len(runs))
        val = runs[perm[-1]]
        train = [r for i, r in enumerate(runs) if i != perm[-1]]
        with _deterministic(self.threads):
            self._fit_split(train, val)
        self._fitted = True

    def _fit_split(self, train, val):
        cells = [self._cells(r) for r in train]
        pooled = {k: np.concatenate([c[k] for c in cells]) for k in cells[0]}
        y = np.concatenate([r.melt[r.shelf_mask] for r in train])
        if not np.isfinite(y).all():
            raise ValueError(f"{self.name}: non-finite melt on a training shelf cell")
        self.scaler = _Scaler(self._X(pooled))
        self._prepare(train, pooled, y)
        inputs = self._tensors(pooled)
        yt = torch.as_tensor((y - self.y_shift) / self.y_scale, dtype=DTYPE)
        val_inputs = self._tensors(self._cells(val))
        members = []
        self.nets = []
        for seed in self.seeds:
            torch.manual_seed(seed)
            net = self._new_net().to(DTYPE)
            rec = self._train(net, inputs, yt, val_inputs, val, seed)
            rec.update(self._member_info(net))
            self.nets.append(net)
            members.append(rec)
        self.info = {
            "device": "cpu", "torch_version": torch.__version__, "hparams": self.hparams(),
            "input_features": list(self.input_features), "val_run": val.case,
            "train_runs": [r.case for r in train], "n_train_cells": int(len(y)),
            "y_shift": float(self.y_shift), "y_scale": float(self.y_scale), "members": members,
        }

    def _train(self, net, inputs, yt, val_inputs, val, seed):
        opt = torch.optim.Adam(net.parameters(), lr=self.lr)
        gen = torch.Generator().manual_seed(seed)
        n = len(yt)
        best, best_state, best_epoch, epoch, stale = np.inf, None, 0, 0, 0
        for epoch in range(1, self.max_epochs + 1):
            net.train()
            order = torch.randperm(n, generator=gen)
            for i in range(0, n, self.batch_size):
                idx = order[i:i + self.batch_size]
                mu, var = net(*(t[idx] for t in inputs))
                loss = _nll(mu, var, yt[idx])
                opt.zero_grad()
                loss.backward()
                opt.step()
            err = rmse_map(self._map(self._mean_phys(net, val_inputs), val), val)
            if err < best:
                best, best_epoch, stale = err, epoch, 0
                best_state = copy.deepcopy(net.state_dict())
            else:
                stale += 1
                if stale >= self.patience:
                    break
        net.load_state_dict(best_state)
        net.eval()
        return {"seed": seed, "epochs_run": epoch, "best_epoch": best_epoch,
                "best_val_rmse": float(best)}

    # --- predict -----------------------------------------------------------------------
    def _phys(self, net, inputs):
        with torch.no_grad():
            mu, var = net(*inputs)
        return (mu.numpy() * self.y_scale + self.y_shift, var.numpy() * self.y_scale**2)

    def _mean_phys(self, net, inputs):
        net.eval()
        return self._phys(net, inputs)[0]

    @staticmethod
    def _map(v, run):
        out = np.full(run.shelf_mask.shape, np.nan)
        out[run.shelf_mask] = v
        return out

    def member_predictions(self, cells: dict) -> tuple[np.ndarray, np.ndarray]:
        """Per-member (mu, sigma^2) in m a-1 for cells given as a dict of 1-D feature
        arrays (as ``cell_features`` restricted to a mask); each output is (M, n_cells)."""
        if not self._fitted:
            raise RuntimeError(f"{self.name}: predict called before fit")
        with _deterministic(self.threads):
            inputs = self._tensors(cells)
            out = [self._phys(net, inputs) for net in self.nets]
        return np.stack([o[0] for o in out]), np.stack([o[1] for o in out])

    def predict_dist(self, run: Run) -> tuple[np.ndarray, np.ndarray]:
        if not self._fitted:
            raise RuntimeError(f"{self.name}: predict called before fit")
        mus, vars_ = self.member_predictions(self._cells(run))
        # Lakshminarayanan 2017 Sec. 2.4 (unnumbered): mixture mean and variance.
        # sigma*^2 = mean(sigma_m^2 + mu_m^2) - mu*^2 = mean(sigma_m^2) + var_pop(mu_m)
        mu = mus.mean(0)
        var = vars_.mean(0) + mus.var(0)
        return self._map(mu, run), self._map(np.sqrt(var), run)

    def predict(self, run: Run) -> np.ndarray:
        return self.predict_dist(run)[0]


class DeepEnsemble(_Ensemble):
    """M6: deep ensemble of per-cell MLPs on all NET_FEATURES (Lakshminarayanan 2017)."""
    name = "M6"
    input_features = NET_FEATURES

    def _prepare(self, train, pooled, y):
        self.y_shift = float(y.mean())
        sd = float(y.std())
        self.y_scale = sd if sd > 0 else 1.0

    def _tensors(self, cells):
        return (torch.as_tensor(self.scaler(self._X(cells)), dtype=DTYPE),)

    def _new_net(self):
        return _MeanVarNet(len(self.input_features), self.hidden)

    def _member_info(self, net):
        return {"K": None, "n": None}


TF_KINDS = ("local", "nonlocal")


class Hybrid(_Ensemble):
    """M7: m = K b^n g_theta(x_g), b = TF+ or (TF+ |<TF>|)^(1/2) (see module docstring)."""
    name = "M7"
    input_features = HYBRID_G_FEATURES

    def __init__(self, n_members=5, seeds=(0, 1, 2, 3, 4), tf_kind="auto", tf_kind_folds=5,
                 **kw):
        super().__init__(n_members=n_members, seeds=seeds, **kw)
        if tf_kind not in TF_KINDS + ("auto",):
            raise ValueError(f"tf_kind must be 'local', 'nonlocal' or 'auto', got {tf_kind!r}")
        if tf_kind_folds is not None and int(tf_kind_folds) < 2:
            raise ValueError(f"tf_kind_folds must be >= 2 or None, got {tf_kind_folds!r}")
        self.tf_kind = tf_kind
        self.tf_kind_folds = tf_kind_folds
        self.tf_kind_: str | None = None
        self._kw = kw

    def hparams(self):
        return dict(super().hparams(), tf_kind=self.tf_kind, tf_kind_folds=self.tf_kind_folds)

    @staticmethod
    def basis(cells, kind):
        """b with m = K b^n g: TF+ (local) or (TF+ |<TF>|)^(1/2) (nonlocal, ruling R3)."""
        tfp = np.maximum(np.asarray(cells["TF"], float), 0.0)
        if kind == "local":
            return tfp
        return np.sqrt(tfp * np.abs(np.asarray(cells["TF_mean"], float)))

    def fit(self, runs):
        runs = list(runs)
        if len(runs) < 2:
            raise ValueError(f"{self.name}: fit needs at least 2 runs (one is the inner-validation run)")
        self._fitted = False
        selection = None
        if self.tf_kind == "auto":
            kind, selection = self._select_kind(runs)
        else:
            kind = self.tf_kind
        self.tf_kind_ = kind
        super().fit(runs)
        self.info["tf_kind"] = kind
        self.info["tf_kind_selection"] = selection

    def _select_kind(self, runs):
        """Inner CV by run over the fit's runs; returns (kind, record)."""
        m = len(runs)
        if m < 3:      # each inner fit needs >= 2 runs
            return "local", {"folds": 0, "rmse": None, "note": "fewer than 3 runs: local by default"}
        k = m if self.tf_kind_folds is None else min(int(self.tf_kind_folds), m)
        perm = np.random.default_rng(self.split_seed + 1).permutation(m)
        groups = np.array_split(perm, k)
        rmse = {}
        per_run = {}
        for kind in TF_KINDS:
            errs = {}
            for grp in groups:
                held = set(int(i) for i in grp)
                inner = Hybrid(n_members=self.n_members, seeds=self.seeds, tf_kind=kind, **self._kw)
                inner.fit([r for i, r in enumerate(runs) if i not in held])
                for i in sorted(held):
                    errs[runs[i].case] = rmse_map(inner.predict(runs[i]), runs[i])
            per_run[kind] = errs
            rmse[kind] = float(np.mean(list(errs.values())))
        kind = "nonlocal" if rmse["nonlocal"] < rmse["local"] else "local"
        return kind, {"folds": k, "rmse": rmse, "per_run_rmse": per_run}

    def _prepare(self, train, pooled, y):
        self.y_shift = 0.0
        sd = float(y.std())
        self.y_scale = sd if sd > 0 else 1.0
        m2 = QuadLocal() if self.tf_kind_ == "local" else QuadNonLocal()
        m2.fit(train)
        # M2: m = coef c^2 |TF| TF ; M3: m = coef c^2 TF |<TF>|; hybrid at n = 2, g = 1:
        # m = K b^2 with b^2 = TF+^2 resp. TF+ |<TF>|, so K0 = coef * c^2 (C from physics.py).
        self.K0 = float(m2.coef * C**2)
        if not self.K0 > 0:
            raise ValueError(f"{self.name}: M2/M3 initial K0 = {self.K0} is not positive")

    def _tensors(self, cells):
        xg = torch.as_tensor(self.scaler(self._X(cells)), dtype=DTYPE)
        b = torch.as_tensor(self.basis(cells, self.tf_kind_), dtype=DTYPE)
        return xg, b

    def _new_net(self):
        return _HybridNet(len(self.input_features), self.hidden, self.K0, self.y_scale)

    def _member_info(self, net):
        return {"K": float(torch.exp(net.log_K.detach())), "n": float(torch.exp(net.log_n.detach())),
                "K0": self.K0, "n0": N0}
