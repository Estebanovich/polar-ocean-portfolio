import dataclasses
import warnings

import numpy as np
import pytest

pytest.importorskip("sklearn")

from exp3_melt_emulator.metrics import rmse_map  # noqa: E402
from exp3_melt_emulator.methods import REGISTRY, Method, ProbMethod  # noqa: E402
from exp3_melt_emulator.methods.stats import CellRegression, GPonEOF, MeanMap  # noqa: E402
from exp3_melt_emulator.tests.synth import synthetic_runs  # noqa: E402


def _smooth_runs(n=8, fn=None):
    """Synthetic runs (common mask) with noise-free melt, smooth in (Tbot, z_c) and column."""
    runs = synthetic_runs(n=n)
    rng = np.random.default_rng(3)
    runs = [dataclasses.replace(r, shelf_mask=runs[0].shelf_mask, draft=runs[0].draft,
                                Tbot=float(rng.uniform(-1, 1)), z_c=float(rng.uniform(300, 500)))
            for r in runs]
    fn = fn or (lambda r, col: 1.0 + 1.5 * r.Tbot + 0.004 * (r.z_c - 300) + 0.3 * col)
    out = []
    for r in runs:
        col = np.broadcast_to(np.arange(r.melt.shape[1])[None, :] - 3.0, r.melt.shape)
        melt = np.where(r.shelf_mask, fn(r, col), np.nan)
        out.append(dataclasses.replace(r, melt=melt, split="train"))
    return out


@pytest.mark.parametrize("name,cls", [("M0", MeanMap), ("M4", CellRegression), ("M5", GPonEOF)])
def test_interface(name, cls):
    runs = synthetic_runs()
    assert isinstance(REGISTRY[name](), cls)
    m = cls()
    assert isinstance(m, Method) and m.name == name
    with pytest.raises(RuntimeError):
        m.predict(runs[0])
    with pytest.raises(ValueError):
        m.fit([])
    m.fit(runs)
    for r in runs:
        p = m.predict(r)
        assert p.shape == r.shelf_mask.shape
        assert np.isnan(p[~r.shelf_mask]).all() and np.isfinite(p[r.shelf_mask]).all()
    assert isinstance(GPonEOF(), ProbMethod)


def test_M0_ignores_absent_cells():
    runs = synthetic_runs()           # run 1 has a smaller mask
    m = MeanMap(); m.fit(runs)
    cell = tuple(np.argwhere(runs[0].shelf_mask & ~runs[1].shelf_mask)[0])
    expected = np.mean([r.melt[cell] for r in runs if r.shelf_mask[cell]])
    assert m.predict(runs[0])[cell] == pytest.approx(expected)
    zero = np.mean([r.melt[cell] if r.shelf_mask[cell] else 0.0 for r in runs])
    assert abs(expected - zero) > 1e-3     # the test would detect zero-filling


@pytest.mark.parametrize("cls", [CellRegression, GPonEOF, MeanMap])
def test_no_leakage(cls):
    runs = synthetic_runs()
    train, held = runs[:-1], runs[-1]
    a = cls(); a.fit(train)
    altered = dataclasses.replace(held, melt=np.where(held.shelf_mask, 99.0, np.nan))
    b = cls(); b.fit(train)                       # fit again, held-out run's melt changed
    for r in train + [held, altered]:
        np.testing.assert_array_equal(a.predict(r), b.predict(r))
    # predict never reads the evaluated run's melt
    np.testing.assert_array_equal(a.predict(held), a.predict(altered))
    # and adding the held-out run to training does change the fit (test is sensitive)
    c = cls(); c.fit(train + [altered])
    assert not np.allclose(a.predict(train[0]), c.predict(train[0]))


def test_deterministic():
    runs = _smooth_runs()
    a, b = GPonEOF(), GPonEOF()
    a.fit(runs); b.fit(runs)
    for r in runs:
        np.testing.assert_array_equal(a.predict_dist(r)[0], b.predict_dist(r)[0])
        np.testing.assert_array_equal(a.predict_dist(r)[1], b.predict_dist(r)[1])


def test_M5_exact_on_training_low_noise():
    runs = _smooth_runs()
    m = GPonEOF(); m.fit(runs)
    for r in runs:
        rms = np.sqrt(np.mean(r.melt[r.shelf_mask] ** 2))
        assert rmse_map(m.predict(r), r) < 0.05 * rms


def test_M5_info_records_k():
    runs = _smooth_runs()
    m = GPonEOF(); m.fit(runs)
    k = m.info["k"]
    assert isinstance(k, int) and 1 <= k <= len(runs) - 1
    assert m.info["explained"] >= 0.99 or k == len(runs) - 1
    assert m.info["recon_rmse"] >= 0.0
    mu, sigma = m.predict_dist(runs[0])
    assert np.isfinite(sigma[runs[0].shelf_mask]).all() and (sigma[runs[0].shelf_mask] >= 0).all()
    assert np.isnan(sigma[~runs[0].shelf_mask]).all()


def test_zero_variance_cells_no_nan():
    runs = _smooth_runs(fn=lambda r, col: np.full(col.shape, 2.0))   # constant melt everywhere
    for cls in (CellRegression, GPonEOF, MeanMap):
        m = cls(); m.fit(runs)
        p = m.predict(runs[0])
        assert np.isfinite(p[runs[0].shelf_mask]).all()
        np.testing.assert_allclose(p[runs[0].shelf_mask], 2.0, atol=1e-6)


def test_M4_chooses_quadratic_on_quadratic_data():
    rng = np.random.default_rng(1)
    runs = synthetic_runs(n=14)
    runs = [dataclasses.replace(r, Tbot=float(rng.uniform(-1, 2)), z_c=float(rng.uniform(250, 550)),
                                s=float(rng.uniform(0.8, 1.2))) for r in runs]
    q = lambda r, col: (1.0 + 2.0 * r.Tbot**2 + 0.5 * (r.z_c - 400) / 100 * r.Tbot + 0.2 * col)
    runs = [dataclasses.replace(r, melt=np.where(r.shelf_mask, q(r, np.broadcast_to(
        np.arange(r.melt.shape[1])[None, :] - 3.0, r.melt.shape)), np.nan)) for r in runs]
    m = CellRegression(); m.fit(runs)
    assert m.degree == 2 and m.loro_rmse[2] < m.loro_rmse[1]


def test_M4_chooses_linear_on_linear_data():
    runs = _smooth_runs(n=10)
    m = CellRegression(); m.fit(runs)
    assert m.degree == 1


def test_M5_info_gp_diagnostics_and_recon():
    runs = _smooth_runs()
    m = GPonEOF()
    with warnings.catch_warnings():
        warnings.simplefilter("error")      # no warning may leak out of fit
        m.fit(runs)
    gp = m.info["gp"]
    assert len(gp) == m.info["k"]
    for d in gp:
        n = len(d["theta"])
        assert isinstance(d["kernel"], str) and "RBF" in d["kernel"]
        assert len(d["bounds"]) == n and len(d["at_bound"]) == n
        assert all(isinstance(b, bool) for b in d["at_bound"])
        assert isinstance(d["convergence_warning"], bool)
        for t, (a, b), ab in zip(d["theta"], d["bounds"], d["at_bound"]):
            assert ab == (t - a < 1e-6 or b - t < 1e-6)
    assert m.info["recon_rmse_mean_run"] == pytest.approx(np.mean(m.info["recon_rmse_per_run"]))
