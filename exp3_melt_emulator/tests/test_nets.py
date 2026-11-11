import dataclasses
from types import SimpleNamespace

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from exp3_melt_emulator.data import FEATURES, cell_features  # noqa: E402
from exp3_melt_emulator.metrics import rmse_map  # noqa: E402
from exp3_melt_emulator.methods import REGISTRY, Method, ProbMethod  # noqa: E402
from exp3_melt_emulator.methods.physics import C, QuadLocal, QuadNonLocal  # noqa: E402
from exp3_melt_emulator.methods.nets import (  # noqa: E402
    HYBRID_G_FEATURES, NET_FEATURES, N_VAL, DeepEnsemble, Hybrid, select_validation_runs,
)
from exp3_melt_emulator.tests.synth import synthetic_runs  # noqa: E402

TEMPERATURE = {"T", "S", "TF", "TF_mean", "Tbot"}
FAST = dict(n_members=2, seeds=(0, 1), max_epochs=40, patience=10)


def _factorial_runs(fn):
    """Synthetic runs on a common mask with Tbot and z_c varied INDEPENDENTLY (3 x 3),
    so that geometry/thermocline features cannot encode Tbot; melt = fn(features)."""
    base = synthetic_runs(n=6)[0]
    out = []
    for i, (Tbot, z_c) in enumerate([(t, z) for t in (-0.5, 0.5, 1.5) for z in (250.0, 300.0, 350.0)]):
        d = base.d
        r = np.clip((-base.z - (z_c - d / 2)) / d, 0, 1)
        run = dataclasses.replace(base, case=f"fact_{i}", split="train", Tbot=Tbot, z_c=z_c,
                                  T_profile=-1.9 + (Tbot + 1.9) * r, S_profile=34.0 + 0.7 * r)
        melt = np.where(run.shelf_mask, fn(cell_features(run), run), np.nan)
        out.append(dataclasses.replace(run, melt=melt))
    return out


def _power_law(f, run):
    dg = f["dist_gl"]
    return 2e-5 * np.maximum(f["TF"], 0.0) ** 2.6 * (1 + 0.3 * dg / np.nanmax(dg[run.shelf_mask]))


def _cells(run):
    f = cell_features(run)
    return {k: v[run.shelf_mask] for k, v in f.items()}


def test_interface():
    runs = synthetic_runs()
    assert isinstance(REGISTRY["M6"](), DeepEnsemble) and isinstance(REGISTRY["M7"](), Hybrid)
    assert NET_FEATURES == FEATURES
    for cls in (DeepEnsemble, Hybrid):
        m = cls(**FAST) if cls is DeepEnsemble else cls(tf_kind="local", **FAST)
        assert isinstance(m, Method) and isinstance(m, ProbMethod)
        assert m.name == ("M6" if cls is DeepEnsemble else "M7")
        with pytest.raises(RuntimeError):
            m.predict(runs[0])
        with pytest.raises(ValueError):
            m.fit([])
        with pytest.raises(ValueError):
            m.fit(runs[:1])
    with pytest.raises(ValueError):
        Hybrid(tf_kind="bogus")
    with pytest.raises(ValueError):
        Hybrid(tf_kind_folds=1)
    with pytest.raises(ValueError):
        DeepEnsemble(n_members=3, seeds=(0, 1))


def test_hybrid_never_sees_temperature():
    # Review Focus 5, part 1
    assert not TEMPERATURE & set(HYBRID_G_FEATURES)
    assert HYBRID_G_FEATURES == ("draft", "draft_minus_zc", "dist_gl", "dist_front", "slope", "s")
    runs = synthetic_runs()
    m = Hybrid(tf_kind="nonlocal", **FAST); m.fit(runs)
    for net in m.nets:
        assert net.mlp[0].in_features == len(HYBRID_G_FEATURES)
    # T and S (and Tbot) do not enter at all; TF and TF_mean enter only through the power law.
    # The variance heads s0, s1 (R11) are unchanged when every temperature variable changes;
    # the variance itself follows sigma^2 = s0^2 + (s1 mu)^2 + floor.
    c = _cells(runs[0])
    mu0, var0 = m.member_predictions(c)
    hot = dict(c, T=c["T"] + 5.0, S=c["S"] - 3.0, TF=3.0 * c["TF"] + 1.0, TF_mean=2.0 * c["TF_mean"])
    mu1, var1 = m.member_predictions(hot)
    a0, b0, floor = m.variance_heads(c)
    a1, b1, _ = m.variance_heads(hot)
    np.testing.assert_array_equal(a0, a1)
    np.testing.assert_array_equal(b0, b1)
    for mu, var in ((mu0, var0), (mu1, var1)):
        np.testing.assert_allclose(var, a0**2 + (b0 * mu) ** 2 + floor, rtol=1e-10)
    assert not np.allclose(mu0, mu1)                   # the mean does respond (test is sensitive)
    only_ts = dict(c, T=c["T"] + 5.0, S=c["S"] - 3.0)
    np.testing.assert_array_equal(m.member_predictions(only_ts)[0], mu0)


@pytest.mark.parametrize("kind", ["local", "nonlocal"])
def test_hybrid_power_law_scaling(kind):
    # Review Focus 5, part 2 (R3: nonlocal scales TF and <TF> together)
    runs = synthetic_runs()
    m = Hybrid(tf_kind=kind, **FAST); m.fit(runs)
    for alpha in (0.6, 1.7):
        for r in runs:
            c = _cells(r)
            mu0, var0 = m.member_predictions(c)
            mu1, var1 = m.member_predictions(dict(c, TF=alpha * c["TF"], TF_mean=alpha * c["TF_mean"]))
            pos = c["TF"] > 0
            assert pos.any() and mu0[:, pos].min() > 0
            for k, mem in enumerate(m.info["members"]):
                np.testing.assert_allclose(mu1[k, pos], alpha ** mem["n"] * mu0[k, pos], rtol=1e-6)
            assert (mu0[:, ~pos] == 0).all() and (mu1[:, ~pos] == 0).all()
            a, b, floor = m.variance_heads(c)
            np.testing.assert_allclose(var1, a**2 + (b * mu1) ** 2 + floor, rtol=1e-10)


def test_hybrid_refreezing_cells_predict_zero():
    runs = synthetic_runs()
    r = runs[0]
    m = Hybrid(tf_kind="local", **FAST); m.fit(runs)
    c = _cells(r)
    c = dict(c, TF=np.where(np.arange(c["TF"].size) % 2 == 0, -0.3, c["TF"]))
    mu, _ = m.member_predictions(c)
    assert (mu[:, c["TF"] <= 0] == 0).all() and (mu[:, c["TF"] > 0] > 0).all()


def test_ensemble_deterministic_cpu():
    runs = synthetic_runs()
    for make in (lambda s: DeepEnsemble(n_members=2, seeds=s, max_epochs=30, patience=10),
                 lambda s: Hybrid(tf_kind="local", n_members=2, seeds=s, max_epochs=30, patience=10)):
        a, b, c = make((0, 1)), make((0, 1)), make((5, 6))
        a.fit(runs); b.fit(runs); c.fit(runs)
        for r in runs:
            ma, sa = a.predict_dist(r)
            mb, sb = b.predict_dist(r)
            assert np.array_equal(ma, mb, equal_nan=True) and np.array_equal(sa, sb, equal_nan=True)
        assert not np.array_equal(a.predict(runs[0]), c.predict(runs[0]), equal_nan=True)
        assert a.info["device"] == "cpu" and a.info["torch_version"] == torch.__version__


def test_fit_restores_global_torch_state():
    threads = torch.get_num_threads()
    det = torch.are_deterministic_algorithms_enabled()
    DeepEnsemble(**FAST).fit(synthetic_runs())
    assert torch.get_num_threads() == threads
    assert torch.are_deterministic_algorithms_enabled() == det


@pytest.mark.parametrize("split_seed", [0, 1])
def test_hybrid_recovers_n(split_seed):
    # Default stopping rule (3 stratified validation runs, max_epochs and patience at their
    # defaults). Only batch_size is reduced: 6 training runs x 20 cells = 120 cells, so the
    # default 512 would be a single full-batch optimiser step per epoch on this tiny set.
    runs = _factorial_runs(_power_law)
    m = Hybrid(tf_kind="local", batch_size=16, split_seed=split_seed)
    m.fit(runs)
    assert all(len(mem["val_runs"]) == N_VAL for mem in m.info["members"])
    ns = [mem["n"] for mem in m.info["members"]]
    assert len(ns) == 5
    assert 2.4 <= float(np.median(ns)) <= 2.8, ns


@pytest.mark.parametrize("make", [lambda: DeepEnsemble(**FAST), lambda: Hybrid(tf_kind="local", **FAST)])
def test_sigma_positive(make):
    runs = synthetic_runs()           # run 1 has a smaller mask, run 3 has negative melt
    m = make(); m.fit(runs)
    for r in runs:
        mu, sigma = m.predict_dist(r)
        assert mu.shape == sigma.shape == r.shelf_mask.shape
        assert np.isnan(mu[~r.shelf_mask]).all() and np.isnan(sigma[~r.shelf_mask]).all()
        assert np.isfinite(mu[r.shelf_mask]).all()
        s = sigma[r.shelf_mask]
        assert np.isfinite(s).all() and (s > 0).all()


@pytest.mark.parametrize("make", [lambda: DeepEnsemble(**FAST), lambda: Hybrid(tf_kind="nonlocal", **FAST)])
def test_mixture_moments(make):
    # Lakshminarayanan 2017 Sec. 2.4: mu* = mean mu_m, sigma*^2 = mean(sigma_m^2 + mu_m^2) - mu*^2
    runs = synthetic_runs()
    m = make(); m.fit(runs)
    r = runs[2]
    mus, vars_ = m.member_predictions(_cells(r))
    mu, sigma = m.predict_dist(r)
    mstar = mus.mean(0)
    np.testing.assert_allclose(mu[r.shelf_mask], mstar, rtol=1e-12)
    np.testing.assert_allclose(sigma[r.shelf_mask] ** 2, (vars_ + mus**2).mean(0) - mstar**2,
                               rtol=1e-8, atol=1e-12)
    np.testing.assert_array_equal(m.predict(r), mu)


@pytest.mark.parametrize("make", [lambda: DeepEnsemble(**FAST), lambda: Hybrid(tf_kind="local", **FAST)])
def test_inner_validation_by_run_and_standardiser(make):
    # per member (R12): validation runs, standardiser, y_scale and restored weights
    runs = synthetic_runs()
    m = make(); m.fit(runs)
    names = m.input_features
    for k, (st, mem) in enumerate(zip(m.members_, m.info["members"])):
        vi = select_validation_runs(runs, st.seed, m.split_seed)
        val = [runs[i] for i in vi]
        train = [r for i, r in enumerate(runs) if i not in vi]
        assert len(val) == N_VAL
        assert mem["val_runs"] == [r.case for r in val]
        assert sorted(mem["train_runs"]) == sorted(r.case for r in train)
        X = np.concatenate([np.stack([_cells(r)[f] for f in names], 1) for r in train])
        sd = X.std(0)
        np.testing.assert_allclose(st.scaler.mean, X.mean(0), rtol=1e-12)
        np.testing.assert_allclose(st.scaler.std, np.where(sd > 0, sd, 1.0), rtol=1e-12)
        y = np.concatenate([r.melt[r.shelf_mask] for r in train])
        assert mem["y_scale"] == pytest.approx(y.std(), rel=1e-12)
        errs = []
        for r in val:
            pred = np.full(r.shelf_mask.shape, np.nan)
            pred[r.shelf_mask] = m.member_predictions(_cells(r))[0][k]
            errs.append(rmse_map(pred, r))
        assert np.mean(errs) == pytest.approx(mem["best_val_rmse"], rel=1e-10)
        assert 1 <= mem["best_epoch"] <= mem["epochs_run"] <= FAST["max_epochs"]
        assert mem["hit_max_epochs"] == (mem["epochs_run"] - mem["best_epoch"] < FAST["patience"])
    assert m.info["n_hit_max_epochs"] == sum(mem["hit_max_epochs"] for mem in m.info["members"])


def test_epoch_budget_bookkeeping():
    runs = synthetic_runs()
    m = DeepEnsemble(n_members=2, seeds=(0, 1), max_epochs=5, patience=1000); m.fit(runs)
    assert all(mem["hit_max_epochs"] and mem["epochs_run"] == 5 for mem in m.info["members"])
    assert m.info["n_hit_max_epochs"] == 2
    m = DeepEnsemble(n_members=2, seeds=(0, 1), max_epochs=10_000, patience=3); m.fit(runs)
    for mem in m.info["members"]:
        assert not mem["hit_max_epochs"] and mem["epochs_run"] == mem["best_epoch"] + 3
        assert mem["val_gain_last_window"] == 0.0 and mem["val_gain_window_before_best"] >= 0.0
    assert DeepEnsemble().max_epochs == Hybrid().max_epochs >= 3000     # shared cap (R13)


def test_hybrid_K0_from_member_training_runs():
    runs = synthetic_runs()
    for kind, cls in (("local", QuadLocal), ("nonlocal", QuadNonLocal)):
        m = Hybrid(tf_kind=kind, **FAST); m.fit(runs)
        for st, mem in zip(m.members_, m.info["members"]):
            ref = cls(); ref.fit(st.train)
            assert mem["K0"] == pytest.approx(ref.coef * C**2, rel=1e-12)
            allr = cls(); allr.fit(runs)
            assert mem["K0"] != pytest.approx(allr.coef * C**2, rel=1e-6)   # validation runs excluded


def _tbot_runs(tbots):
    base = synthetic_runs(n=6)
    return [dataclasses.replace(base[i % 6], case=f"r{i}", Tbot=float(t)) for i, t in enumerate(tbots)]


# (case, Tbot) of the 25 training runs of the 'final' fold (exp3 dataset, sha 3d2eceb1...)
REAL_TRAIN = [("ens_T-0.5_zc250_s1.0", -0.5), ("ens_T-0.5_zc400_s1.0", -0.5), ("ens_T-0.5_zc550_s1.0", -0.5),
              ("ens_T0.0_zc250_s1.0", 0.0), ("ens_T0.0_zc400_s1.0", 0.0), ("ens_T0.0_zc550_s1.0", 0.0),
              ("ens_T0.5_zc250_s1.0", 0.5), ("ens_T0.5_zc400_s1.0", 0.5), ("ens_T0.5_zc550_s1.0", 0.5),
              ("ens_T1.0_zc250_s1.0", 1.0), ("ens_T1.0_zc400_s1.0", 1.0), ("ens_T1.0_zc550_s1.0", 1.0),
              ("ens_T1.5_zc250_s1.0", 1.5), ("ens_T1.5_zc400_s1.0", 1.5), ("ens_T1.5_zc550_s1.0", 1.5),
              ("ens_T0.75_zc527_s0.8", 0.75), ("ens_T1.43_zc295_s0.8", 1.43), ("ens_T-0.05_zc328_s0.8", -0.05),
              ("ens_T-0.42_zc402_s0.8", -0.42), ("ens_T0.34_zc444_s0.8", 0.34), ("ens_T1.42_zc481_s1.2", 1.42),
              ("ens_T0.02_zc291_s1.2", 0.02), ("ens_T-0.39_zc503_s1.2", -0.39), ("ens_T0.97_zc386_s1.2", 0.97),
              ("ens_T0.66_zc311_s1.2", 0.66)]


def test_select_validation_runs_stratified():
    runs = synthetic_runs()                       # Tbot = -0.5, -0.1, 0.3, 0.7, 1.1, 1.5 (one run each)
    for seed in range(5):
        assert select_validation_runs(runs, seed) == [5, 2, 0]   # highest, lower-median, lowest
    runs = _tbot_runs([0, 0, 1, 1, 2, 2, 2])
    picks = []
    for seed in range(20):
        vi = select_validation_runs(runs, seed)
        assert vi == select_validation_runs(runs, seed)              # deterministic
        assert len(set(vi)) == 3
        assert runs[vi[0]].Tbot == 2 and runs[vi[1]].Tbot == 1 and runs[vi[2]].Tbot == 0
        picks.append(vi[0])
    assert set(picks) == {4, 5, 6}               # the top-level pick is random over that level
    assert select_validation_runs(runs, 0, split_seed=1) != select_validation_runs(runs, 0) or \
        select_validation_runs(runs, 1, split_seed=1) != select_validation_runs(runs, 1)
    # a single Tbot level: three distinct runs; < N_VAL + 1 runs: one run at the highest level
    assert len(set(select_validation_runs(_tbot_runs([1.0] * 5), 0))) == 3
    assert select_validation_runs(_tbot_runs([0, 2, 1]), 0) == [1]
    assert select_validation_runs(_tbot_runs([0, 2]), 3) == [1]


def test_every_run_trains_in_some_member_real_fold():
    # R12 on the real final-fold training list: with the default member seeds no run is held
    # out by every member, and the top-Tbot validation run differs between members
    runs = [SimpleNamespace(case=c, Tbot=t) for c, t in REAL_TRAIN]
    sets = [set(select_validation_runs(runs, seed)) for seed in (0, 1, 2, 3, 4)]
    assert not set.intersection(*sets)
    assert len({select_validation_runs(runs, seed)[0] for seed in (0, 1, 2, 3, 4)}) > 1
    for seed in (0, 1, 2, 3, 4):
        vi = select_validation_runs(runs, seed)
        assert [runs[i].Tbot for i in vi] == [1.5, 0.5, -0.5]


def test_every_run_trains_in_some_member_synthetic():
    runs = _tbot_runs([0, 0, 1, 1, 2, 2, 2])
    m = DeepEnsemble(**dict(FAST, n_members=5, seeds=(0, 1, 2, 3, 4))); m.fit(runs)
    assert m.info["never_trained"] == []
    trained = set().union(*[set(mem["train_runs"]) for mem in m.info["members"]])
    assert trained == {r.case for r in runs}


@pytest.mark.parametrize("n", [2, 3])
def test_fit_with_few_runs_uses_one_val_run(n):
    runs = synthetic_runs()[:n]
    for m in (DeepEnsemble(**FAST), Hybrid(tf_kind="local", **FAST)):
        m.fit(runs)
        for mem in m.info["members"]:
            assert mem["val_runs"] == [runs[-1].case]          # the only run at the highest Tbot
            assert len(mem["train_runs"]) == n - 1


def test_hybrid_variance_mean_detached():
    # R11: the variance term cannot move the physics mean: d(var)/d(log n) = d(var)/d(log K) = 0
    runs = synthetic_runs()
    m = Hybrid(tf_kind="local", **FAST); m.fit(runs)
    net = m.nets[0]
    xg, b = m._tensors(_cells(runs[0]), m.members_[0])
    mu, var = net(xg, b)
    net.zero_grad()
    var.sum().backward(retain_graph=True)
    assert net.log_n.grad is None or float(net.log_n.grad) == 0.0
    assert net.log_K.grad is None or float(net.log_K.grad) == 0.0
    assert any(p.grad is not None and p.grad.abs().sum() > 0 for p in net.mlp.parameters())
    net.zero_grad()
    mu.sum().backward()
    assert float(net.log_n.grad) != 0.0                      # the mean does depend on n
    assert net.mlp[0].in_features == len(HYBRID_G_FEATURES)


@pytest.mark.parametrize("make", [lambda: DeepEnsemble(**FAST), lambda: Hybrid(tf_kind="nonlocal", **FAST)])
def test_predict_members_decomposition(make):
    runs = synthetic_runs()
    m = make(); m.fit(runs)
    r = runs[1]                                   # smaller mask
    d = m.predict_members(r)
    mu, sigma = m.predict_dist(r)
    M = m.n_members
    assert d["mu"].shape == d["var"].shape == (M,) + r.shelf_mask.shape
    assert np.isnan(d["mu"][:, ~r.shelf_mask]).all()
    np.testing.assert_allclose(d["aleatoric"] + d["epistemic"], sigma**2, rtol=1e-10)
    sh = r.shelf_mask
    np.testing.assert_allclose(d["mu"][:, sh].mean(0), mu[sh], rtol=1e-12)
    assert d["mean_aleatoric"] == pytest.approx(np.nanmean(d["aleatoric"]))
    assert d["mean_epistemic"] == pytest.approx(np.nanmean(d["epistemic"]))
    assert d["mean_epistemic"] > 0
    np.testing.assert_allclose(d["integrated"], [np.nansum(v) for v in d["mu"]], rtol=1e-12)
    assert d["integrated_mixture"] == pytest.approx(np.nansum(mu))


def test_fit_restores_rng_state():
    torch.manual_seed(123)
    state = torch.get_rng_state()
    Hybrid(tf_kind="local", **FAST).fit(synthetic_runs())
    assert torch.equal(torch.get_rng_state(), state)


def test_non_finite_validation_raises():
    class Broken(DeepEnsemble):
        def _val_error(self, st, val_inputs):
            return float("nan")
    with pytest.raises(RuntimeError, match="not finite"):
        Broken(**FAST).fit(synthetic_runs())


def test_no_leakage_from_unused_run():
    runs = synthetic_runs()
    train, held = runs[:-1], runs[-1]
    a = DeepEnsemble(**FAST); a.fit(train)
    altered = dataclasses.replace(held, melt=np.where(held.shelf_mask, 99.0, np.nan))
    np.testing.assert_array_equal(a.predict(held), a.predict(altered))


def test_tf_kind_selection_by_inner_cv():
    runs = _factorial_runs(_power_law)
    m = Hybrid(tf_kind="auto", tf_kind_folds=3, **FAST); m.fit(runs)
    sel = m.info["tf_kind_selection"]
    assert sel["folds"] == 3 and set(sel["rmse"]) == {"local", "nonlocal"}
    best = min(sel["rmse"], key=lambda k: (sel["rmse"][k], k != "local"))
    assert m.tf_kind_ == best == m.info["tf_kind"]
    for mem in m.info["members"]:
        assert mem["K"] > 0 and mem["n"] > 0
