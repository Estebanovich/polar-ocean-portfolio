"""Tests for exp2_isomip_cavity/analysis/calibrate.py: the K1 Gamma_T calibration (Task 4).

pd0 K1 (binding): deep (< -300 m) Ocean0 melt over the last 6 months of a 2-year run,
interpolated piecewise linearly in (ln Gamma_T, ln <m>) to 30 m/yr; no extrapolation; runs
refused by require_ok (N3) stop the calibration. The synthetic runs come from the
test_metrics writer, so every run goes through metrics.load_run.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import calibrate as cal
from test_metrics import default_geometry, make_run

YEAR = 31536000.0


def deep_draft():
    """Shelf columns (x < 640 km) with base at -400 m: the whole shelf is 'deep300'."""
    d = default_geometry()[0].copy()
    d[d < 0] = -400.0
    return d


def melt_run(root, gammaT, melt_m_per_a, *, nrec=24, status="OK", case_extra=None):
    """A run whose deep melt is ``melt_m_per_a`` [m/yr] in every record (SHIfwFlx < 0)."""
    fw = -melt_m_per_a * 1000.0 / YEAR
    extra = {"params": {"gammaT": gammaT, "gammaS": gammaT / 35.0, "dx_km": 4.0, "s": 1.0}}
    if case_extra:
        extra["params"].update(case_extra)
    return make_run(root, nrec=nrec, shi=lambda k: {"SHIfwFlx": fw}, draft=deep_draft(),
                    status=status, case_extra=extra)


def power_law_runs(tmp_path, A=500.0, p=0.8, gammas=(0.015, 0.03, 0.06)):
    return [melt_run(tmp_path / f"r{g}", g, A * g ** p) for g in gammas]


# ---- deep melt (K1 quantity) -------------------------------------------------------------

def test_deep_melt_last6_reads_window(tmp_path):
    """Records 18-23 only: months 1-18 carry another melt and must not enter."""
    fw = lambda k: {"SHIfwFlx": -(10.0 if k < 18 else 40.0) * 1000.0 / YEAR}
    r = make_run(tmp_path / "r", shi=fw, draft=deep_draft(),
                 case_extra={"params": {"gammaT": 0.03, "gammaS": 0.03 / 35}})
    pt = cal.calibration_point(r)
    assert pt["deep_melt_last6_m_per_a"] == pytest.approx(40.0, rel=1e-12)
    assert pt["gammaT"] == 0.03


def test_deep_mask_excludes_shallow_shelf(tmp_path):
    """Only ice base < -300 m (AD16 Eq. 37) enters K1; shallow shelf melt does not."""
    draft = deep_draft()
    draft[:, 0] = -100.0                          # shallow shelf column
    xmelt = np.where(draft < -300.0, 40.0, 5.0)
    xmelt[draft >= 0] = 0.0
    fw = -xmelt * 1000.0 / YEAR
    r = make_run(tmp_path / "r", shi=lambda k: {"SHIfwFlx": fw}, draft=draft,
                 case_extra={"params": {"gammaT": 0.03, "gammaS": 0.03 / 35}})
    pt = cal.calibration_point(r)
    assert pt["deep_melt_last6_m_per_a"] == pytest.approx(40.0, rel=1e-12)
    assert pt["shelf_melt_last6_m_per_a"] < 40.0


# ---- interpolation ------------------------------------------------------------------------

def test_exact_recovery_power_law(tmp_path):
    """m = A Gamma^p is a straight line in log-log: Gamma_cal = (30 / A)^(1/p) exactly."""
    A, p = 500.0, 0.8
    res = cal.calibrate(power_law_runs(tmp_path, A, p), target=30.0)
    assert res["status"] == "OK"
    assert res["gammaT_cal_raw"] == pytest.approx((30.0 / A) ** (1 / p), rel=1e-12)
    assert res["segment"] == [0.015, 0.03]


def test_piecewise_linear_not_power_law():
    """Non power law: the result uses only the bracketing segment, linear in ln-ln."""
    g = [0.015, 0.03, 0.06]
    m = [12.0, 24.0, 60.0]           # slope 1 on the first segment, ln(2.5)/ln 2 on the second
    gc, seg = cal.interpolate_loglog(g, m, 40.0)
    s2 = math.log(60.0 / 24.0) / math.log(2.0)
    want = math.exp(math.log(0.03) + (math.log(40.0) - math.log(24.0)) / s2)
    assert gc == pytest.approx(want, rel=1e-13)
    assert seg == (0.03, 0.06)
    # differs from a single log-log fit through all points (proves piecewise behaviour)
    k, b = np.polyfit(np.log(g), np.log(m), 1)
    assert abs(math.exp((math.log(40.0) - b) / k) - want) / want > 1e-3
    # a node is returned exactly, and the input order does not matter
    assert cal.interpolate_loglog([0.06, 0.015, 0.03], [60.0, 12.0, 24.0], 24.0)[0] == \
        pytest.approx(0.03, rel=1e-14)


def test_four_point_second_iteration():
    """K1 second iteration: 4 points (3 + the verification run), same rule."""
    g = [0.015, 0.03, 0.06, 0.04]
    m = [12.0, 24.0, 60.0, 33.0]
    gc, seg = cal.interpolate_loglog(g, m, 30.0)
    assert seg == (0.03, 0.04)
    s = math.log(33.0 / 24.0) / math.log(0.04 / 0.03)
    assert gc == pytest.approx(math.exp(math.log(0.03) + math.log(30.0 / 24.0) / s), rel=1e-13)


@pytest.mark.parametrize("m", [[31.0, 40.0, 50.0], [10.0, 20.0, 29.0]])
def test_not_bracketed_is_refused(m):
    """No extrapolation (pd0 K1): 30 outside [min m, max m] -> NotBracketed."""
    with pytest.raises(cal.NotBracketed, match="PD2"):
        cal.interpolate_loglog([0.015, 0.03, 0.06], m, 30.0)


def test_not_bracketed_report_written(tmp_path):
    """calibrate() reports NOT_BRACKETED (gammaT_cal None) and the CLI exits non-zero."""
    runs = [melt_run(tmp_path / f"r{g}", g, v) for g, v in [(0.015, 40.0), (0.03, 50.0),
                                                            (0.06, 60.0)]]
    res = cal.calibrate(runs, target=30.0)
    assert res["status"] == "NOT_BRACKETED" and res["gammaT_cal"] is None
    assert "PD2" in res["message"]
    out = tmp_path / "cal.json"
    rc = cal.main(["--res", "4km", "--out", str(out)] + [str(r) for r in runs])
    assert rc == 3
    assert json.loads(out.read_text())["status"] == "NOT_BRACKETED"


def test_non_monotonic_is_refused():
    with pytest.raises(cal.NotBracketed, match="monotonic"):
        cal.interpolate_loglog([0.015, 0.03, 0.06], [20.0, 35.0, 25.0], 30.0)


def test_duplicate_or_bad_inputs_refused():
    with pytest.raises(ValueError):
        cal.interpolate_loglog([0.03, 0.03, 0.06], [20.0, 25.0, 40.0], 30.0)
    with pytest.raises(ValueError):
        cal.interpolate_loglog([0.015, 0.03], [20.0, -1.0], 30.0)   # ln of a non-positive melt


# ---- N3 and consistency guards ------------------------------------------------------------

def test_failed_run_is_refused(tmp_path):
    """A run that fails N3 (require_ok) stops the calibration; no partial result."""
    runs = power_law_runs(tmp_path)
    bad = melt_run(tmp_path / "bad", 0.045, 20.0, status="FAILED")
    with pytest.raises(RuntimeError):
        cal.calibrate(runs + [bad], target=30.0)


def test_nan_run_is_refused(tmp_path):
    runs = power_law_runs(tmp_path)
    (runs[1] / "STDOUT.0000").write_text("(PID.TID 0000.0001) %MON dynstat_theta_max = NaN\n")
    with pytest.raises(RuntimeError):
        cal.calibrate(runs, target=30.0)


def test_short_run_is_refused(tmp_path):
    """Not a 2-year run of 24 monthly records -> the 'last6' window does not exist."""
    runs = power_law_runs(tmp_path)[:2] + [melt_run(tmp_path / "short", 0.06, 60.0, nrec=12)]
    with pytest.raises(ValueError, match="24"):
        cal.calibrate(runs, target=30.0)


def test_mixed_resolution_or_gammaS_refused(tmp_path):
    runs = power_law_runs(tmp_path)[:2]
    other = melt_run(tmp_path / "o", 0.06, 60.0, case_extra={"dx_km": 2.0})
    with pytest.raises(ValueError, match="resolution"):
        cal.calibrate(runs + [other], target=30.0)
    bad_s = melt_run(tmp_path / "s", 0.06, 60.0, case_extra={"gammaS": 0.06 / 30})
    with pytest.raises(ValueError, match="gammaS"):
        cal.calibrate(runs + [bad_s], target=30.0)


def test_non_ocean0_run_refused(tmp_path):
    runs = power_law_runs(tmp_path)[:2]
    other = melt_run(tmp_path / "o", 0.06, 60.0, case_extra={"init_profile": "ramp"})
    with pytest.raises(ValueError, match="Ocean0"):
        cal.calibrate(runs + [other], target=30.0)


# ---- K1 verdict and CLI -------------------------------------------------------------------

@pytest.mark.parametrize("m,ok", [(28.0, True), (32.0, True), (30.5, True), (27.99, False),
                                  (32.01, False)])
def test_k1_verdict_closed_interval(m, ok):
    assert cal.k1_verdict(m)["pass"] is ok


def test_cli_writes_json_and_prints(tmp_path, capsys):
    runs = power_law_runs(tmp_path)
    out = tmp_path / "figures" / "calibration_4km.json"
    rc = cal.main(["--res", "4km", "--out", str(out)] + [str(r) for r in runs])
    assert rc == 0
    d = json.loads(out.read_text())
    assert d["resolution"] == "4km" and d["status"] == "OK"
    assert len(d["points"]) == 3 and d["points"][0]["gammaT"] == 0.015
    assert all("N1_drift" in p for p in d["points"])
    printed = capsys.readouterr().out.strip().splitlines()[-1]
    assert printed == f"{d['gammaT_cal']:.5g}"
    assert float(printed) == d["gammaT_cal"]
    assert d["gammaT_cal"] == pytest.approx(d["gammaT_cal_raw"], rel=1e-4)
