import dataclasses

import numpy as np
import pytest

from exp3_melt_emulator import metrics as M
from exp3_melt_emulator.data import cell_features


def test_rmse_ignores_outside_mask(synth_runs):
    run = synth_runs[1]                      # smaller mask
    pred = np.where(run.shelf_mask, run.melt, 0.0)
    assert M.rmse_map(pred, run) == pytest.approx(0.0, abs=1e-12)
    wild = np.where(run.shelf_mask, run.melt, 1e6)
    assert M.rmse_map(wild, run) == pytest.approx(0.0, abs=1e-12)
    nan_out = np.where(run.shelf_mask, run.melt, np.nan)
    assert M.rmse_map(nan_out, run) == pytest.approx(0.0, abs=1e-12)
    assert M.integrated_rel_error(wild, run) == pytest.approx(0.0, abs=1e-12)


def test_nan_inside_mask_raises(synth_runs):
    run = synth_runs[0]
    pred = run.melt.copy()
    pred[tuple(np.argwhere(run.shelf_mask)[0])] = np.nan
    for f in (M.rmse_map, M.integrated_rel_error, M.rmse_gl):
        with pytest.raises(ValueError):
            f(pred, run)


def test_rmse_zero_pred_negative_cells(synth_runs):
    run = synth_runs[3]                      # has refreezing cells
    m = run.melt[run.shelf_mask]
    assert (m < 0).any()
    assert M.rmse_map(np.zeros_like(run.melt), run) == pytest.approx(np.sqrt(np.mean(m**2)))


def test_integrated_error_sign(synth_runs):
    run = synth_runs[0]
    assert M.integrated_rel_error(1.1 * run.melt, run) == pytest.approx(0.10)
    assert M.integrated_rel_error(0.9 * run.melt, run) == pytest.approx(-0.10)


def test_gl_region_fraction(synth_runs):
    run = synth_runs[0]
    reg = M.gl_region(run)
    assert reg.dtype == bool and reg.shape == run.melt.shape
    assert reg.sum() == round(0.10 * run.shelf_mask.sum())
    assert not (reg & ~run.shelf_mask).any()
    d = cell_features(run)["dist_gl"]
    assert d[reg].max() <= d[run.shelf_mask & ~reg].min()
    assert np.array_equal(reg, M.gl_region(run))          # deterministic
    assert M.gl_region(run, 0.5).sum() == round(0.5 * run.shelf_mask.sum())


def test_rmse_gl_is_gl_mean_error(synth_runs):
    run = synth_runs[0]
    reg = M.gl_region(run)
    pred = run.melt + 1.0
    assert M.rmse_gl(pred, run) == pytest.approx(1.0)
    pred = run.melt.copy()
    pred[reg] += 2.0
    pred[run.shelf_mask & ~reg] -= 5.0
    assert M.rmse_gl(pred, run) == pytest.approx(2.0)


def test_skill():
    assert M.skill(0.5, 2.0) == pytest.approx(0.75)
    assert M.skill(2.0, 2.0) == 0.0


def _big_run(synth_runs, n, rng, sigma):
    mu = rng.uniform(0, 5, n)
    melt = (mu + sigma * rng.standard_normal(n))[None, :]
    run = dataclasses.replace(synth_runs[0], melt=melt, shelf_mask=np.ones((1, n), bool))
    return run, mu[None, :]


def test_coverage_gaussian(synth_runs):
    rng = np.random.default_rng(1)
    sigma = np.full((1, 100_000), 0.7)
    run, mu = _big_run(synth_runs, 100_000, rng, 0.7)
    for level in (0.5, 0.8, 0.95):
        assert abs(M.coverage(mu, sigma, run, level) - level) < 0.01


def test_crps_zero_sigma_is_mae(synth_runs):
    run = synth_runs[3]
    mu = run.melt + 0.3 * np.cos(np.arange(run.melt.size)).reshape(run.melt.shape)
    mae = np.mean(np.abs(mu - run.melt)[run.shelf_mask])
    assert M.crps_gaussian(mu, np.full_like(mu, 1e-9), run) == pytest.approx(mae, abs=1e-6)


def test_crps_known_value(synth_runs):
    run = dataclasses.replace(synth_runs[0], melt=np.zeros((1, 1)), shelf_mask=np.ones((1, 1), bool))
    # y = mu: CRPS = sigma (2 phi(0) - 1/sqrt(pi)) = sigma (sqrt(2/pi) - 1/sqrt(pi))
    got = M.crps_gaussian(np.zeros((1, 1)), np.full((1, 1), 2.0), run)
    assert got == pytest.approx(2.0 * (np.sqrt(2 / np.pi) - 1 / np.sqrt(np.pi)))


def test_lambert_quadratic():
    # melt(T) = 2T^2 + T at T = 0, 1, 2 -> M = 0, 3, 10
    f = lambda t: 2 * t**2 + t
    out = M.lambert_sensitivities({t: f(t) for t in (2.0, 0.0, 1.0)})
    assert out["dMdT_ref_p1"] == pytest.approx(3.0)       # Eq. 1
    assert out["dMdT_p1_p2"] == pytest.approx(7.0)
    assert out["d2MdT2"] == pytest.approx(4.0)            # Eq. 2: (7-3)/(1.5-0.5)
    # offset reference temperature and non-unit spacing
    out = M.lambert_sensitivities({t: f(t) for t in (-1.0, -0.5, 0.0)})
    assert out["d2MdT2"] == pytest.approx(4.0)
    with pytest.raises(ValueError):
        M.lambert_sensitivities({0.0: 1.0, 1.0: 2.0})


def test_deep_amplification(synth_runs):
    a, b = synth_runs[0], synth_runs[0]
    n10 = round(0.10 * a.shelf_mask.sum())
    # uniform warming response: amplification 1
    assert M.deep_amplification(a.melt, a.melt + 1.0, a, b) == pytest.approx(1.0)
    # extra response only in the deepest cells
    order = np.argsort(np.where(a.shelf_mask, a.draft, np.inf), axis=None, kind="stable")[:n10]
    dm = np.zeros(a.melt.shape)
    dm.flat[order] = 10.0
    pb = a.melt + 1.0 + dm
    m_all = (pb - a.melt)[a.shelf_mask].mean()
    m_deep = (pb - a.melt).flat[order].mean()
    assert M.deep_amplification(a.melt, pb, a, b) == pytest.approx(m_deep / m_all)
    assert m_deep / m_all > 1
