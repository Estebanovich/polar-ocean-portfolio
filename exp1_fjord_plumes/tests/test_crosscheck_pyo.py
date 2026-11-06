"""Tests of the pyOCEAN cross-check of the (c) open boundary (crosscheck/pyocean_c): the exterior shelf series
and the evaluator of the four pre-registered criteria on synthetic series.  Tiny and fast: no model is built or
run and LaCie is never touched (only boundaries.obcs and forcing.boundary_time are imported from pyOCEAN)."""
from __future__ import annotations

import json

import numpy as np
import pytest

from exp1_fjord_plumes.crosscheck.pyocean_c import eval_pyo_c as ev
from exp1_fjord_plumes.crosscheck.pyocean_c import run_pyo_c as rp

DAY = 86400.0


@pytest.fixture(scope="module")
def ext():
    inp = rp.c_inputs(0.0)
    return inp, rp.exterior(inp, 2 * DAY)


def test_exterior_pycnocline_follows_shelf_series(ext):
    """Exterior pycnocline (T = 1.25 degC isotherm, linear in z) = 200 - 30 sin(2 pi t / 6 d) at record times,
    between records and at the periodic wrap (t < dt_rec/2)."""
    inp, f = ext
    zc = inp["zc"]
    for t in (0.0, 900.0, 1800.0, 5400.0, 0.25 * 6 * DAY, 1.37 * DAY, 1.5 * DAY):
        ob = f(t)["east"]
        assert ob.T.shape == (inp["nz"], inp["p"]["ny"]) and ob.zw[-1, 0] == -inp["p"]["H"]
        z_pyc = ev.pycnocline_depth(ob.T[:, 5], zc)
        target = 200.0 - 30.0 * np.sin(2 * np.pi * t / (6 * DAY))
        # 10 m levels and hourly records: < 0.3 m (vertical sampling of the tanh, linear time interpolation)
        assert abs(z_pyc - target) < 0.3, (t, z_pyc, target)
        assert np.all(ob.u_normal == 0.0) and ob.eta is None and ob.ubar is None


def test_exterior_uses_make_case_records(ext):
    """At a record centre the exterior equals the OBEt/OBEs record make_case writes for MITgcm."""
    inp, f = ext
    ob = f(10.5 * 3600.0)["east"]
    np.testing.assert_allclose(ob.T, inp["obT"][10], rtol=0, atol=1e-12)
    np.testing.assert_allclose(ob.S, inp["obS"][10], rtol=0, atol=1e-12)


def test_pycnocline_depth_linear():
    z = -np.array([5.0, 15.0, 25.0])
    assert ev.pycnocline_depth(np.array([0.0, 1.0, 2.0]), z, T_mid=1.5) == pytest.approx(20.0)
    assert np.isnan(ev.pycnocline_depth(np.full(3, 3.0), z, T_mid=1.5))          # no crossing
    assert ev.pycnocline_depth(np.array([2.0, 1.0, 0.0]), z[::-1], T_mid=1.5) == pytest.approx(20.0)


def _series(Q, eta_fn, qout_fn, days=30.0):
    t = np.arange(1, int(days * 24) + 1) * 3600.0
    eta = eta_fn(t)
    q = qout_fn(t)
    cum = np.cumsum(q * 3600.0)
    return t, eta, np.abs(eta), cum


def test_c2_q0_pass_and_fail():
    t, e, ea, c = _series(0, lambda t: 1e-4 * np.sin(t / DAY), lambda t: 0 * t)
    assert ev.criterion2(t, e, ea, c, 0.0)["passed"]
    ea2 = ea.copy()
    ea2[100] = 0.011
    r = ev.criterion2(t, e, ea2, c, 0.0)
    assert not r["passed"] and r["max_abs_eta_mean_m"] == pytest.approx(0.011)


def test_c2_q100_transport_and_drift():
    ok = ev.criterion2(*_series(100, lambda t: 1e-3 + 0 * t, lambda t: 101.5 + 0 * t), 100.0)
    assert ok["passed"] and ok["Q_out_mean_d20_30"] == pytest.approx(101.5)
    bad_q = ev.criterion2(*_series(100, lambda t: 0 * t, lambda t: 97.9 + 0 * t), 100.0)
    assert not bad_q["passed"] and not bad_q["q_passed"] and bad_q["drift_passed"]
    # drift of 0.012 m over days 20-30 (linear rise), transport exact
    drift = ev.criterion2(*_series(100, lambda t: 1.2e-3 * np.maximum(t / DAY - 20, 0), lambda t: 100 + 0 * t), 100.0)
    assert not drift["passed"] and drift["q_passed"] and drift["eta_drift_d20_30_m"] == pytest.approx(0.012)
    short = ev.criterion2(*_series(100, lambda t: 0 * t, lambda t: 100 + 0 * t, days=25.0), 100.0)
    assert not short["passed"]


def test_c3_ratio():
    d = np.arange(1, 31)
    ok = ev.criterion3(d, np.full(30, 0.15), np.full(30, 0.1))
    assert ok["passed"] and ok["max_ratio"] == pytest.approx(1.5)
    ob = np.full(30, 0.1)
    ob[12] = 0.4
    bad = ev.criterion3(d, ob, np.full(30, 0.1))
    assert not bad["passed"] and bad["day_of_max"] == 13 and bad["n_days_failed"] == 1


def test_c4_sine_fit_amplitude_and_lag():
    t = np.arange(0, 30 * 24 + 1) * 3600.0
    w = 2 * np.pi / (6 * DAY)
    z = 195.0 - 12.0 * np.sin(w * t - np.radians(40.0))
    r = ev.criterion4(t, z)
    assert r["amplitude"] == pytest.approx(12.0) and r["phase_deg"] == pytest.approx(40.0)
    assert r["amplitude_ratio"] == pytest.approx(0.4) and r["mean"] == pytest.approx(195.0)
    assert r["lag_days"] == pytest.approx(40.0 / 360.0 * 6.0)


def test_c1_and_evaluate_run(tmp_path):
    """End to end on a synthetic run directory (status, series, daily u)."""
    run = tmp_path / "pyo_shelf"
    run.mkdir()
    (run / "status.json").write_text(json.dumps(dict(config=dict(case="pyo_shelf", Q_runoff=0.0), state="OK",
                                                     finite=True, completed=True, t_s=30 * DAY)))
    t = np.arange(1, 721) * 3600.0
    w = 2 * np.pi / (6 * DAY)
    cols = {c: np.zeros_like(t) for c in rp.SERIES_COLS}
    cols["t"] = t
    cols["pyc_mid"] = 200.0 - 10.0 * np.sin(w * t)
    np.savez(run / "series.npz", data=np.column_stack([cols[c] for c in rp.SERIES_COLS]),
             columns=np.array(rp.SERIES_COLS))
    for d in range(0, 31):
        u = np.full((2, 3, 201), 0.1)
        u[:, :, -3:] = 0.12
        np.savez(run / f"day_{d:02d}.npz", u=u)
    r = ev.evaluate_run(run)
    assert r["c1"]["passed"] and r["c2"]["passed"] and r["c3"]["passed"] and r["c3"]["n_days"] == 30
    assert r["c4"]["amplitude"] == pytest.approx(10.0) and abs(r["c4"]["phase_deg"]) < 1e-6
    assert not ev.criterion1(dict(completed=True, finite=False, t_s=30 * DAY))["passed"]
    assert not ev.criterion1(dict(completed=False, finite=True, t_s=12.8 * DAY))["passed"]
