"""Statistical emulators M0 (mean map), M4 (per-cell regression), M5 (GP on EOFs).

All three work on the union grid of the runs passed to ``fit`` (cells that belong to
the shelf of at least one training run). Everything data-dependent (per-cell means and
standard deviations, input scaling, EOFs, kernels, the linear/quadratic choice) is
computed inside ``fit`` from that fit's runs only (no leakage). ``predict`` never reads
``run.melt``. A cell of the evaluated run that no training run ever had on its shelf
is predicted with the pooled training mean (M5: sigma = pooled training std).

Inputs of M4/M5 are the run-level parameters (Tbot, z_c, s), standardised with the
fit's runs; a zero standard deviation is replaced by 1 (no division by zero).
"""
from __future__ import annotations

import warnings
from itertools import combinations_with_replacement

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

from exp3_melt_emulator.data import Run, union_mask
from exp3_melt_emulator.metrics import rmse_map

SEED = 0
VAR_FRACTION = 0.99   # Myren & Lawrence 2021 Sec. 3.1: keep EOFs explaining 99 % of variance
# M5 GP kernel per EOF coefficient: ConstantKernel * RBF(ARD) + WhiteKernel (init, bounds)
CONST_INIT, CONST_BOUNDS = 1.0, (1e-3, 1e3)
RBF_INIT, RBF_BOUNDS = 1.0, (1e-2, 1e3)
# White-noise nugget (normalised units, normalize_y=True), PD0 decision 4 (Esteban,
# 2026-11-10): Myren & Lawrence 2021 Eq. 6-7 has no nugget; Berdahl et al. 2021 Sec. 2.3 use
# "a small nugget for numerical conditioning" (lower bound following Ranjan et al. 2011).
# Our runs carry a small window-to-window drift (noise floor 0.157 m/a, about 1e-3 in
# normalised coefficient units), so the nugget is small and bounded, unable to absorb signal.
WHITE_INIT, WHITE_BOUNDS = 1e-3, (1e-8, 1e-2)
N_RESTARTS = 5
M4_DEGREES = (1, 2)


def _check(runs):
    if not len(runs):
        raise ValueError("fit requires at least one run")


def _stack(runs, union):
    """(m, n_union) matrix of melt on union cells, NaN where the run's mask is False."""
    return np.stack([np.where(r.shelf_mask, r.melt, np.nan)[union] for r in runs])


def _safe_std(a, axis=0):
    s = np.nanstd(a, axis=axis, ddof=1) if a.shape[axis] > 1 else np.zeros(a.shape[1 - axis])
    return np.where(np.isfinite(s) & (s > 0), s, 1.0)


def _params(runs):
    return np.array([[r.Tbot, r.z_c, r.s] for r in runs], float)


class _InputScaler:
    def __init__(self, runs):
        P = _params(runs)
        self.mean = P.mean(0)
        sd = P.std(0)
        self.std = np.where(sd > 0, sd, 1.0)

    def __call__(self, runs):
        return (_params(runs) - self.mean) / self.std


class _UnionMethod:
    """Shared state: union grid + scatter of union-cell vectors back to maps."""
    name = ""

    def __init__(self):
        self._fitted = False

    def _need_fit(self):
        if not self._fitted:
            raise RuntimeError(f"{self.name}: predict called before fit")

    def _begin(self, runs):
        _check(runs)
        self.union = union_mask(runs)
        self.Y = _stack(runs, self.union)                  # (m, n)
        self.cell_mean = np.nanmean(self.Y, axis=0)
        self.pooled_mean = float(np.nanmean(self.Y))
        self.pooled_std = float(np.nanstd(self.Y))

    def _to_map(self, v, run, fill):
        out = np.full(run.shelf_mask.shape, np.nan)
        full = np.full(run.shelf_mask.shape, fill, float)
        full[self.union] = v
        out[run.shelf_mask] = full[run.shelf_mask]
        return out


class MeanMap(_UnionMethod):
    """M0: per-cell mean of the training maps; absent cells contribute nothing."""
    name = "M0"

    def fit(self, runs):
        self._begin(runs)
        self._fitted = True

    def hparams(self):
        return dict(free_parameters=0, input_features=[], unseen_cells="pooled training mean")

    def predict(self, run):
        self._need_fit()
        return self._to_map(self.cell_mean, run, self.pooled_mean)


def _design(Xs, degree):
    cols = [np.ones(len(Xs))] + [Xs[:, i] for i in range(Xs.shape[1])]
    if degree == 2:
        cols += [Xs[:, i] * Xs[:, j] for i, j in combinations_with_replacement(range(Xs.shape[1]), 2)]
    return np.stack(cols, axis=1)


def _fit_coefs(Xs, Y, degree):
    """Least-squares coefficients per union cell, using only runs present in the cell.

    Underdetermined cells get the minimum-norm solution (np.linalg.lstsq)."""
    A = _design(Xs, degree)
    coef = np.zeros((Y.shape[1], A.shape[1]))
    for c in range(Y.shape[1]):
        ok = np.isfinite(Y[:, c])
        if ok.any():
            coef[c] = np.linalg.lstsq(A[ok], Y[ok, c], rcond=None)[0]
    return coef


class CellRegression(_UnionMethod):
    """M4: per-cell least squares on standardised (Tbot, z_c, s).

    Degree (1 = linear, 2 = full quadratic with 10 terms) is chosen by inner
    leave-one-run-out over the fit's runs, criterion = mean ``metrics.rmse_map`` of the
    held-out run (ties go to linear). With fewer than 3 runs the linear model is used.
    """
    name = "M4"

    def hparams(self):
        return dict(input_features=["Tbot", "z_c", "s"], standardised=True, degrees=list(M4_DEGREES),
                    selection="inner leave-one-run-out over the fit's runs, mean rmse_map, "
                              "ties -> linear", min_runs_for_selection=3,
                    unseen_cells="pooled training mean")

    def _fit_degree(self, runs, degree):
        self._begin(runs)
        self.scaler = _InputScaler(runs)
        self.degree = degree
        self.coef = _fit_coefs(self.scaler(runs), self.Y, degree)
        self._fitted = True

    def fit(self, runs):
        _check(runs)
        runs = list(runs)
        self._fitted = False
        degree = 1
        self.loro_rmse = {}
        if len(runs) >= 3:
            for deg in M4_DEGREES:
                errs = []
                for j in range(len(runs)):
                    inner = CellRegression()
                    inner._fit_degree(runs[:j] + runs[j + 1:], deg)
                    errs.append(rmse_map(inner.predict(runs[j]), runs[j]))
                self.loro_rmse[deg] = float(np.mean(errs))
            degree = 2 if self.loro_rmse[2] < self.loro_rmse[1] else 1
        self._fit_degree(runs, degree)

    def predict(self, run):
        self._need_fit()
        a = _design(self.scaler([run]), self.degree)[0]
        return self._to_map(self.coef @ a, run, self.pooled_mean)


class GPonEOF(_UnionMethod):
    """M5: Myren & Lawrence 2021 Sec. 3.1 (Eqs. 2-6), as a ProbMethod.

    Union-grid maps are standardised per cell (mean/std over the runs having the cell;
    Eq. 2-3), entries outside a run's mask are set to 0 after standardising, an SVD
    Z = U S V' (Eq. 4) gives the EOFs V; the first k explaining >= 99 % of the variance,
    k <= m - 1, are kept. Coefficients w = Z V_k are emulated by independent GPs
    (Eq. 6; ARD squared-exponential + white noise, normalize_y=True, 5 restarts,
    random_state=0) on standardised (Tbot, z_c, s). Prediction (Eq. 5):
    y = mu + delta * (V_k w).

    Sigma: per-coefficient predictive variances propagated through the EOFs assuming
    independent coefficients, var(cell) = sum_h v_h V[cell,h]^2, times delta^2. The GP
    predictive variance includes the fitted white-noise level. EOF truncation error is
    NOT added to sigma (it is a deterministic in-sample residual, not a predictive
    variance, and is <= 1 % of the standardised variance by construction); it is
    reported in ``info`` instead: ``k``, ``explained``, ``recon_rmse`` (in-sample RMSE
    of the truncated reconstruction vs the training maps, m/a, POOLED over all runs'
    mask cells), ``recon_rmse_per_run`` and ``recon_rmse_mean_run`` (mean of the
    per-run RMSEs, comparable with per-run metrics). ``info["gp"]`` has, per coefficient
    h, the fitted kernel string, log-theta, bounds, ``at_bound`` flags (theta within
    1e-6 of a bound) and whether the optimiser raised a ConvergenceWarning.
    """
    name = "M5"

    def __init__(self):
        super().__init__()
        self.info: dict = {}

    def hparams(self):
        return dict(input_features=["Tbot", "z_c", "s"], standardised=True, var_fraction=VAR_FRACTION,
                    k_cap="m - 1", kernel="ConstantKernel * RBF(ARD) + WhiteKernel",
                    const_init=CONST_INIT, const_bounds=list(CONST_BOUNDS),
                    rbf_init=RBF_INIT, rbf_bounds=list(RBF_BOUNDS),
                    white_init=WHITE_INIT, white_bounds=list(WHITE_BOUNDS),
                    normalize_y=True, n_restarts=N_RESTARTS, random_state=SEED,
                    sigma="EOF propagation of independent coefficient variances; "
                          "truncation error not added")

    def fit(self, runs):
        self._fitted = False
        self._begin(runs)
        runs = list(runs)
        m = len(runs)
        self.delta = _safe_std(self.Y)
        present = np.isfinite(self.Y)
        Z = np.where(present, (self.Y - self.cell_mean) / self.delta, 0.0)
        _, S, Vt = np.linalg.svd(Z, full_matrices=False)
        var = S**2
        total = var.sum()
        cum = np.cumsum(var) / total if total > 0 else np.ones_like(var)
        k = int(np.searchsorted(cum, VAR_FRACTION) + 1)
        k = max(1, min(k, m - 1, len(S)))
        self.V = Vt[:k].T                                   # (n, k)
        W = Z @ self.V                                      # (m, k)
        self.scaler = _InputScaler(runs)
        X = self.scaler(runs)
        self.gps = []
        gp_diag = []
        for h in range(k):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                kern = (ConstantKernel(CONST_INIT, CONST_BOUNDS)
                        * RBF(np.full(X.shape[1], RBF_INIT), RBF_BOUNDS)
                        + WhiteKernel(WHITE_INIT, WHITE_BOUNDS))
                gp = GaussianProcessRegressor(kern, normalize_y=True,
                                              n_restarts_optimizer=N_RESTARTS, random_state=SEED)
                self.gps.append(gp.fit(X, W[:, h]))
            th, bd = gp.kernel_.theta, gp.kernel_.bounds
            gp_diag.append({
                "kernel": str(gp.kernel_),
                "theta": [float(t) for t in th],
                "bounds": [[float(a), float(b)] for a, b in bd],
                "at_bound": [bool(t - a < 1e-6 or b - t < 1e-6) for t, (a, b) in zip(th, bd)],
                "convergence_warning": any(issubclass(w.category, ConvergenceWarning)
                                           for w in caught),
            })
        recon = (Z @ self.V @ self.V.T) * self.delta + self.cell_mean
        sq = np.where(present, (recon - self.Y) ** 2, np.nan)
        self.info = {
            "k": k, "n_runs": m, "explained": float(cum[k - 1]),
            "recon_rmse": float(np.sqrt(np.nanmean(sq))),
            "recon_rmse_mean_run": float(np.mean([np.sqrt(np.nanmean(sq[i])) for i in range(m)])),
            "gp": gp_diag,
            "recon_rmse_per_run": [float(np.sqrt(np.nanmean(sq[i]))) for i in range(m)],
        }
        self._fitted = True

    def predict_dist(self, run):
        self._need_fit()
        x = self.scaler([run])
        w = np.empty(len(self.gps))
        v = np.empty(len(self.gps))
        for h, gp in enumerate(self.gps):
            mu, sd = gp.predict(x, return_std=True)
            w[h], v[h] = mu[0], sd[0] ** 2
        mean = self.cell_mean + self.delta * (self.V @ w)
        sigma = self.delta * np.sqrt((self.V**2) @ v)
        return (self._to_map(mean, run, self.pooled_mean),
                self._to_map(sigma, run, self.pooled_std))

    def predict(self, run):
        return self.predict_dist(run)[0]
