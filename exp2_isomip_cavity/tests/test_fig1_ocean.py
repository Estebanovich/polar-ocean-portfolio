"""Tests for analysis/fig1_ocean0_validation.py: pd0 V3 band, V2 ranges, Ocean1 verdicts and
the Figure 1 data assembly (synthetic inputs, no model output)."""
from __future__ import annotations

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import fig1_ocean0_validation as f1


# ---- V3 band (pd0 table: m = 5 -> bt 0.17-0.30, moc 0.10-0.18; m = 10 -> bt 0.25-0.45,
# moc 0.19-0.28), hand-computed values -------------------------------------------------------

def test_v3_band_at_the_anchors():
    assert f1.v3_band(5.0, "bt") == pytest.approx((0.17, 0.30), abs=1e-12)
    assert f1.v3_band(10.0, "bt") == pytest.approx((0.25, 0.45), abs=1e-12)
    assert f1.v3_band(5.0, "moc") == pytest.approx((0.10, 0.18), abs=1e-12)
    assert f1.v3_band(10.0, "moc") == pytest.approx((0.19, 0.28), abs=1e-12)


def test_v3_band_linear_between_anchors():
    # m = 7.5: bt lo = 0.17 + 0.5 * 0.08 = 0.21, hi = 0.30 + 0.5 * 0.15 = 0.375
    assert f1.v3_band(7.5, "bt") == pytest.approx((0.21, 0.375), abs=1e-12)
    # m = 6: moc lo = 0.10 + 0.2 * 0.09 = 0.118, hi = 0.18 + 0.2 * 0.10 = 0.20
    assert f1.v3_band(6.0, "moc") == pytest.approx((0.118, 0.20), abs=1e-12)


def test_v3_band_extends_the_anchor_lines_to_the_evaluable_range():
    # m = 12: bt lo = 0.17 + 1.4 * 0.08 = 0.282, hi = 0.30 + 1.4 * 0.15 = 0.51
    assert f1.v3_band(12.0, "bt") == pytest.approx((0.282, 0.51), abs=1e-12)
    # m = 4: moc lo = 0.10 - 0.2 * 0.09 = 0.082, hi = 0.18 - 0.2 * 0.10 = 0.16
    assert f1.v3_band(4.0, "moc") == pytest.approx((0.082, 0.16), abs=1e-12)


def test_v3_not_evaluable_outside_4_12_and_boundaries_included():
    assert f1.v3_verdict(3.999, 0.3, "bt")["verdict"] == "NOT_EVALUABLE"
    assert f1.v3_verdict(12.001, 0.3, "bt")["verdict"] == "NOT_EVALUABLE"
    assert f1.v3_verdict(4.0, 0.12, "moc")["verdict"] == "PASS"         # band 0.082-0.16
    assert f1.v3_verdict(12.0, 0.30, "bt")["verdict"] == "PASS"         # band 0.282-0.51
    assert f1.v3_band(3.0, "bt") is None


def test_v3_pass_fail_inside_and_outside_the_band():
    assert f1.v3_verdict(10.0, 0.25, "bt")["verdict"] == "PASS"         # lower edge included
    assert f1.v3_verdict(10.0, 0.2499, "bt")["verdict"] == "FAIL"
    assert f1.v3_verdict(10.0, 0.4501, "bt")["verdict"] == "FAIL"
    v = f1.v3_verdict(7.5, 0.30, "bt")
    assert v["verdict"] == "PASS" and v["band"] == pytest.approx([0.21, 0.375])
    with pytest.raises(ValueError):
        f1.v3_verdict(7.5, 0.3, "psi")


def test_v3_strict_interpolation_variant_clamps_to_the_anchors():
    """Reported alongside: the band held at the anchor values outside [5, 10]."""
    assert f1.v3_band(12.0, "bt", extend=False) == pytest.approx((0.25, 0.45))
    assert f1.v3_band(4.0, "moc", extend=False) == pytest.approx((0.10, 0.18))


# ---- V2 ranges and fits ---------------------------------------------------------------

def test_v2_ranges_closed_intervals():
    assert f1.range_verdict(4.0, (4.0, 12.0)) == "PASS"
    assert f1.range_verdict(12.0, (4.0, 12.0)) == "PASS"
    assert f1.range_verdict(0.1999, f1.V2_RANGES["V2b"]["range"]) == "FAIL"
    assert f1.V2_RANGES["V2c"]["range"] == (0.09, 0.43)
    assert f1.V2_RANGES["V2d"]["range"] == (0.14, 0.68)


def test_linear_fit_recovers_slope_intercept():
    x = np.linspace(1, 12, 50)
    fit = f1.linear_fit(x, 0.024 * x + 0.071)
    assert fit["slope"] == pytest.approx(0.024) and fit["intercept"] == pytest.approx(0.071)
    assert fit["r"] == pytest.approx(1.0)


def test_cfl_window_only_counts_monitor_blocks_after_t_min():
    log = ""
    for t, w in ((0.0, 0.9), (100.0, 0.2), (200.0, 0.3)):
        log += f"(PID.TID 0000.0001) %MON time_secondsf = {t:.13E}\n"
        for k in f1.CFL_KEYS:
            log += f"(PID.TID 0000.0001) %MON {k} = {w if 'wvel' in k else 0.1:.13E}\n"
    out = f1.cfl_maxima_window(log, t_min=100.0)
    assert out["max"] == pytest.approx(0.3) and out["n_monitor_blocks"] == 2
    assert not out["above_0.5"]
    assert f1.cfl_maxima_window(log, t_min=0.0)["above_0.5"]


# ---- verdicts and figure data from a synthetic summary -----------------------------------

def _summary(melt20=10.3, bt=0.26, moc=0.186, v2b=0.185):
    n = 240
    melt = list(np.linspace(0.2, melt20, n))
    melt[228:] = [melt20] * 12
    return {"n_records": n, "window": {"name": "year20", "records_1based": list(range(229, 241))},
            "mean_melt_shelf": {"value": melt20}, "mean_melt_deep300": {"value": 30.0},
            "overturning_max": {"value": v2b}, "psi_bt_max_cavity": {"value": 0.26},
            "psi_bt_max_front": {"value": 0.49},
            "circulation_strength": {"window_mean": {"melt_shelf": melt20, "bt": bt, "moc": moc},
                                     "monthly": {"melt_shelf": melt,
                                                 "bt": [0.024 * m + 0.071 for m in melt],
                                                 "moc": [0.02 * m + 0.019 for m in melt]}},
            "N1": {"value": 0.001, "pass": True}, "N2": {"value": 1e-6, "pass": True}}


def test_ocean1_verdicts_mechanical():
    v = f1.ocean1_verdicts(_summary())
    assert v["V2a"]["verdict"] == "PASS" and v["V2b"]["verdict"] == "FAIL"
    assert v["V2c"]["verdict"] == "PASS" and v["V2d"]["verdict"] == "PASS"
    assert v["V3"]["bt"]["verdict"] == "PASS" and v["V3"]["moc"]["verdict"] == "FAIL"
    assert v["V3"]["bt"]["fit_240_months"]["slope"] == pytest.approx(0.024)
    assert v["V3"]["moc"]["fit_published"] == {"slope": 0.02, "intercept": 0.019, "r": 0.97}
    assert v["N1"]["pass"] and v["N2"]["pass"]


def test_ocean1_verdicts_v3_not_evaluable_when_melt_outside():
    v = f1.ocean1_verdicts(_summary(melt20=13.0))
    assert v["V2a"]["verdict"] == "FAIL"
    assert v["V3"]["bt"]["verdict"] == "NOT_EVALUABLE"
    assert v["V3"]["moc"]["verdict"] == "NOT_EVALUABLE"


def test_figure_data_assembly():
    s = _summary()
    ny, nx = 5, 8
    melt_map = np.arange(ny * nx, dtype=float).reshape(ny, nx)
    shelf = np.zeros((ny, nx), bool)
    shelf[1:4, 2:6] = True
    wet = shelf.copy()
    wet[:, 6:] = True
    cal = {"2km": {"points": [{"gammaT": 0.015, "deep_melt_last6_m_per_a": 13.0},
                              {"gammaT": 0.03, "deep_melt_last6_m_per_a": 29.2}],
                   "gammaT_cal": 0.030959},
           "4km": {"points": [{"gammaT": 0.015, "deep_melt_last6_m_per_a": 14.2}],
                   "gammaT_cal": 0.029356}}
    ver = {"2km": {"point": {"gammaT": 0.030959, "deep_melt_last6_m_per_a": 30.0}},
           "4km": {"point": {"gammaT": 0.029356, "deep_melt_last6_m_per_a": 30.02}}}
    d = f1.figure_data(s, f1.ocean1_verdicts(s), melt_map, shelf, wet,
                       np.arange(nx) * 2e3, np.arange(ny) * 2e3, cal, ver)
    assert np.isnan(d["a"]["melt"][0][0]) and d["a"]["melt"][1][2] == 10.0   # off-shelf masked
    assert d["b"]["years"][0] == pytest.approx(1 / 24) and len(d["b"]["melt"]) == 240
    assert d["b"]["year20"] == [19.0, 20.0]
    assert [r["name"] for r in d["c"]] == ["V2a", "V2b", "V2c", "V2d"]
    assert d["c"][1]["verdict"] == "FAIL"
    assert d["c"][0]["position"] == pytest.approx((10.3 - 4) / 8)
    assert len(d["d"]["bt"]["monthly_melt"]) == 240 and d["d"]["bt"]["verdict"] == "PASS"
    assert d["d"]["band_label"] == "read from Yung et al. 2026 Fig. 13 (approx.)"
    assert d["d"]["bt"]["band_m"][0] == 4.0 and d["d"]["bt"]["band_m"][-1] == 12.0
    assert d["e"]["2km"]["gammaT"] == [0.015, 0.03]
    assert d["e"]["4km"]["verification"] == [0.029356, 30.02]
    assert d["e"]["k1_band"] == [28.0, 32.0]


# ---- PD2 (2026-11-10): V2b post hoc row, V3 margins, extrapolated band segments ----------

def test_v2b_full_domain_posthoc_row():
    xg = np.array([600e3, 700e3, 766e3])
    rf = np.array([0.0, -200.0, -220.0, -240.0])
    psi = np.zeros((4, 3))
    psi[1, 0] = 0.185                                   # cavity peak (x < 640 km)
    psi[2, 2] = -0.2224                                 # northern-boundary peak, |psi| larger
    r = f1.v2b_full_domain_posthoc(psi, xg, rf, jpl_value=0.234)
    assert r["value"] == pytest.approx(0.2224) and r["verdict_against_range"] == "PASS"
    assert r["x_km"] == pytest.approx(766.0) and r["z_m"] == pytest.approx(-220.0)
    assert r["post_hoc"] is True and "Yung" in r["label"] and r["jpl_value"] == 0.234
    assert r["range"] == [0.20, 0.25]


def test_v3_margins_reported_for_both_band_variants():
    # moc at m = 10.34: extended lo = 0.10 + 1.068 * 0.09 = 0.19612; held lo = 0.19
    v = f1.v3_verdict(10.34, 0.1865, "moc")
    assert v["verdict"] == "FAIL"
    assert v["margin_extended"]["abs_Sv"] == pytest.approx(0.19612 - 0.1865, abs=1e-9)
    assert v["margin_extended"]["rel"] == pytest.approx((0.19612 - 0.1865) / 0.19612, rel=1e-9)
    assert v["margin_held_at_anchors"]["abs_Sv"] == pytest.approx(0.0035, abs=1e-9)
    assert v["margin_held_at_anchors"]["rel"] == pytest.approx(0.0035 / 0.19, rel=1e-9)
    inside = f1.v3_verdict(10.0, 0.30, "bt")
    assert inside["margin_extended"]["abs_Sv"] == 0.0                # inside: no margin


def test_band_segments_outside_anchors_flagged_extrapolated():
    s = _summary()
    d = f1.figure_data(s, f1.ocean1_verdicts(s), np.zeros((2, 2)), np.ones((2, 2), bool),
                       np.ones((2, 2), bool), np.arange(2) * 2e3, np.arange(2) * 2e3,
                       {r: {"points": [], "gammaT_cal": 0.03} for r in ("2km", "4km")},
                       {r: {"point": {"gammaT": 0.03, "deep_melt_last6_m_per_a": 30.0}}
                        for r in ("2km", "4km")})
    m = np.array(d["d"]["bt"]["band_m"])
    ex = np.array(d["d"]["bt"]["band_extrapolated"])
    assert 5.0 in m and 10.0 in m
    assert (ex == ((m < 5.0) | (m > 10.0))).all()
    assert d["d"]["bt"]["m_bar_in_extrapolated_band"] is True          # m_bar = 10.3
