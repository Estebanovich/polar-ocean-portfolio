import dataclasses

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from exp3_melt_emulator.data import FEATURES, cell_features  # noqa: E402
from exp3_melt_emulator.metrics import rmse_map  # noqa: E402
from exp3_melt_emulator.methods import REGISTRY, Method, ProbMethod  # noqa: E402
from exp3_melt_emulator.methods.nets import (  # noqa: E402
    HYBRID_G_FEATURES, NET_FEATURES, DeepEnsemble, Hybrid,
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
    # T and S (and Tbot) do not enter at all; TF and TF_mean enter only through the power law,
    # so the variance (geometry-only head) is unchanged when every temperature variable changes.
    c = _cells(runs[0])
    mu0, var0 = m.member_predictions(c)
    hot = dict(c, T=c["T"] + 5.0, S=c["S"] - 3.0, TF=3.0 * c["TF"] + 1.0, TF_mean=2.0 * c["TF_mean"])
    mu1, var1 = m.member_predictions(hot)
    np.testing.assert_array_equal(var0, var1)
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
            np.testing.assert_array_equal(var0, var1)


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
    # 9 runs x 20 cells: batch_size=16 gives several optimiser steps per epoch (the default
    # 512 would be one full-batch step per epoch on this tiny set)
    runs = _factorial_runs(_power_law)
    m = Hybrid(tf_kind="local", batch_size=16, max_epochs=300, patience=60, split_seed=split_seed)
    m.fit(runs)
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
    runs = synthetic_runs()
    m = make(); m.fit(runs)
    perm = np.random.default_rng(m.split_seed).permutation(len(runs))
    val = runs[perm[-1]]
    assert m.info["val_run"] == val.case
    assert sorted(m.info["train_runs"]) == sorted(r.case for r in runs if r is not val)
    # standardiser fitted on the training cells only (the inner-validation run excluded)
    names = m.input_features
    train = [r for r in runs if r is not val]
    X = np.concatenate([np.stack([_cells(r)[k] for k in names], 1) for r in train])
    sd = X.std(0)
    np.testing.assert_allclose(m.scaler.mean, X.mean(0), rtol=1e-12)
    np.testing.assert_allclose(m.scaler.std, np.where(sd > 0, sd, 1.0), rtol=1e-12)
    # early stopping restored the best weights: each member's val RMSE equals the recorded best
    mus, _ = m.member_predictions(_cells(val))
    for k, mem in enumerate(m.info["members"]):
        pred = np.full(val.shelf_mask.shape, np.nan)
        pred[val.shelf_mask] = mus[k]
        assert rmse_map(pred, val) == pytest.approx(mem["best_val_rmse"], rel=1e-10)
        assert 1 <= mem["best_epoch"] <= mem["epochs_run"] <= FAST["max_epochs"]


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
