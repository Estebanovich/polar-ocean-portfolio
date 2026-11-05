"""Tests of the mechanical PD2 evaluator on synthetic inputs (no real runs)."""
from __future__ import annotations

import math

import pytest

from exp1_fjord_plumes.analysis import pd2_convergence_eval as ev

OK_R = {"conclusive": True, "R_ext": 1.0}
BAD_R = {"conclusive": False, "R_ext": None}
U, D, P = "sube", "baja", "plano"


# ---------------------------------------------------------------- classify
@pytest.mark.parametrize("args,expected", [
    # V first: invalid beats everything, even a width-fixed pattern
    ((False, (U, U), (D, D), True, OK_R), "V"),
    # G before S: width fixed beats a clean supporting pattern
    ((True, (U, U), (D, D), True, OK_R), "G"),
    ((True, (D, D), (P, P), True, BAD_R), "G"),
    # S / P1 / P2 (R_Q up twice, b down twice)
    ((True, (U, U), (D, D), False, {"conclusive": True, "R_ext": 0.75}), "S"),
    ((True, (U, U), (D, D), False, {"conclusive": True, "R_ext": 1.25}), "S"),
    ((True, (U, U), (D, D), False, {"conclusive": True, "R_ext": 0.7499}), "P1"),
    ((True, (U, U), (D, D), False, {"conclusive": True, "R_ext": 1.2501}), "P2"),
    ((True, (U, U), (D, D), False, BAD_R), "P2"),
    # W: R_Q up twice but b does not go down in both
    ((True, (U, U), (D, P), False, OK_R), "W"),
    ((True, (U, U), (U, D), False, OK_R), "W"),
    # D, F, M
    ((True, (D, D), (D, D), False, OK_R), "D"),
    ((True, (P, P), (D, D), False, OK_R), "F"),
    ((True, (U, D), (D, D), False, OK_R), "M"),
    ((True, (U, P), (D, D), False, OK_R), "M"),
    ((True, (P, D), (D, D), False, OK_R), "M"),
    ((True, (D, U), (D, D), False, OK_R), "M"),
])
def test_classify_every_class_and_order(args, expected):
    assert ev.classify(*args) == expected


def test_every_class_of_the_note_is_covered():
    assert set(ev.CLASS_TEXT) == {"V", "G", "S", "P1", "P2", "W", "D", "F", "M"}


# ------------------------------------------------------- thresholds / trend
def test_threshold_floors_and_trend_labels():
    assert ev.threshold(0.0, 0.0, 0.03) == 0.03
    assert ev.threshold(0.03, 0.04, 0.03) == pytest.approx(2 * 0.05)
    assert ev.threshold(0.1, 0.1, 2.5) == 2.5  # floor wins over 2 sigma = 0.283
    assert ev.trend(0.5, 0.6, 0.05) == "sube"
    assert ev.trend(0.5, 0.4, 0.05) == "baja"
    assert ev.trend(0.5, 0.549, 0.05) == "plano"
    assert ev.trend(0.5, 0.5, 0.05) == "plano"
    assert ev.trend(10.0, 6.0, 2.5) == "baja"


def test_stationarity():
    assert ev.is_stationary(0.50, 0.505, 0.5025)
    assert not ev.is_stationary(0.50, 0.52, 0.51)


def test_width_fixed_G_flag():
    b5 = 10.0
    # b proportional to the local grid: both ratios equal r_eff exactly
    assert ev.width_fixed_flag(b5 * ev.R32, b5, b5 / ev.R21)[0]
    # within 20 % in both steps still counts as fixed (ratio/r_eff = 0.85 and 1.15)
    assert ev.width_fixed_flag(b5 * ev.R32 * 0.85, b5, b5 / ev.R21 / 1.15)[0]
    # b independent of the grid (ratios 1): |1/1.9 - 1| > 0.2
    assert not ev.width_fixed_flag(b5, b5, b5)[0]
    # only one of the two steps grid-proportional: not G
    assert not ev.width_fixed_flag(b5 * ev.R32, b5, b5)[0]


# ---------------------------------------------------------------- Richardson
def _phi(p, c=0.1, ext=1.0, r21=ev.R21, r32=ev.R32):
    h1, h2, h3 = 1.0, r21, r21 * r32
    return ext + c * h1**p, ext + c * h2**p, ext + c * h3**p


@pytest.mark.parametrize("p", [0.6, 1.0, 1.5, 2.0, 2.9])
def test_richardson_recovers_known_order(p):
    r = ev.richardson(*_phi(p))
    assert r["conclusive"]
    assert r["p"] == pytest.approx(p, rel=1e-6)
    assert r["R_ext"] == pytest.approx(1.0, abs=1e-6)


def test_richardson_equal_ratios_reduces_to_log_formula():
    phi1, phi2, phi3 = _phi(1.5, r21=2.0, r32=2.0)
    r = ev.richardson(phi1, phi2, phi3, r21=2.0, r32=2.0)
    assert r["p"] == pytest.approx(math.log((phi3 - phi2) / (phi2 - phi1)) / math.log(2.0))


@pytest.mark.parametrize("phis,reason", [
    ((0.5, 0.6, 0.55), "non-monotone"),      # s < 0
    ((0.5, 0.5, 0.6), "epsilon"),            # eps21 = 0
    ((0.5, 0.6, 0.6), "epsilon"),            # eps32 = 0
])
def test_richardson_inconclusive_branches(phis, reason):
    r = ev.richardson(*phis)
    assert not r["conclusive"] and r["R_ext"] is None
    assert reason in r["reason"]


def test_richardson_p_out_of_range_is_inconclusive():
    low = ev.richardson(*_phi(0.3))   # p < 0.5
    high = ev.richardson(*_phi(3.6))  # p > 3
    for r in (low, high):
        assert not r["conclusive"] and r["R_ext"] is None and "outside" in r["reason"]
    assert low["p"] == pytest.approx(0.3, rel=1e-6)
    assert high["p"] == pytest.approx(3.6, rel=1e-6)


def test_richardson_p_r2_is_diagnostic_only():
    r = ev.richardson(*_phi(1.5))
    assert r["p_r2"] == pytest.approx(
        abs(math.log((_phi(1.5)[2] - _phi(1.5)[1]) / (_phi(1.5)[1] - _phi(1.5)[0]))) / math.log(2))


# ------------------------------------------------------------ validity checks
def _log(n, cg=500, skip=None, cfl=0.3, cap_at=None):
    lines = []
    for i in range(1, n + 1):
        if i == skip:
            continue
        it = 2000 if i == cap_at else cg
        lines.append(f" %MON time_tsnumber = {i}\n %MON advcfl_uvel_max = {cfl:.6E}\n"
                     f" %MON advcfl_wvel_max = {cfl:.6E}\n      cg3d_iters (last) =   {it}\n")
    return "# cg3d_iters is printed\n" + "".join(lines)


def test_monitor_valid_log():
    r = ev.check_monitor(ev.parse_monitor(_log(10)), 10)
    assert r["valid"] and r["per_step_verified"] and r["cg3d_max"] == 500


def test_monitor_step_at_cap_is_invalid():
    r = ev.check_monitor(ev.parse_monitor(_log(10, cap_at=7)), 10)
    assert not r["valid"] and any("cap" in x for x in r["reasons"])


def test_monitor_missing_iteration_line_is_invalid():
    r = ev.check_monitor(ev.parse_monitor(_log(10, skip=4)), 10)
    assert not r["valid"] and any("found 9" in x for x in r["reasons"])


def test_monitor_cfl_over_limit_is_invalid():
    r = ev.check_monitor(ev.parse_monitor(_log(10, cfl=0.51)), 10)
    assert not r["valid"] and any("CFL" in x for x in r["reasons"])
    assert ev.check_monitor(ev.parse_monitor(_log(10, cfl=0.5)), 10)["valid"]  # <= 0.5


def test_monitor_not_per_step_reference_is_flagged_not_invalid():
    r = ev.check_monitor(ev.parse_monitor(_log(10, skip=4)), 10, require_per_step=False)
    assert r["valid"] and not r["per_step_verified"] and "not per-step verified" in r["note"]


# ------------------------------------------------------------------ assembly
def _m(rq, b, valid=True, sR=0.005, sb=0.3, reasons=()):
    return {"valid": valid, "reasons": list(reasons), "R_Q": rq, "b": b, "sigma_R_Q": sR,
            "sigma_b": sb, "C2_pass": True, "C2pp_pass": True}


def test_assemble_supported_and_sensitivities():
    metrics = {"conv_dx10": _m(0.40, 38.0), "conv_dx5": _m(0.55, 25.0), "conv_dx2p5": _m(0.68, 14.0),
               ev.REF: _m(0.56, 24.0), "sens_freeslip": _m(0.70, 24.0),
               "sens_xu2012": _m(0.56, 24.0)}
    # note: dx10->dx5 b ratio 1.52/1.912 = 0.79 (|.-1| = 0.21 >= 0.2) so G is not triggered
    res = ev.assemble(metrics)
    c = res["convergence"]
    assert c["class"] in {"S", "P1", "P2"}
    assert c["steps"]["dx10_to_dx5"]["R_Q_trend"] == "sube"
    s = res["sensitivity"]
    assert s["wall_free_slip"]["R_Q_trend"] == "sube" and s["wall_free_slip"]["wall_drag_contributes"]
    assert s["xu2012_closure_joint_with_wall_drag"]["R_Q_trend"] == "plano"
    assert s["dz_sensitivity"]["R_Q_trend"] == "plano"


def test_assemble_invalid_run_gives_V_and_missing_gives_incomplete():
    base = {"conv_dx10": _m(0.4, 38.0), "conv_dx5": _m(0.5, 25.0, valid=False, reasons=["cg3d cap"]),
            "conv_dx2p5": _m(0.6, 14.0)}
    res = ev.assemble(base)
    assert res["convergence"]["class"] == "V" and "conv_dx5" in res["convergence"]["invalid_runs"]
    res2 = ev.assemble({"conv_dx10": _m(0.4, 38.0)})
    assert res2["convergence"]["class"] == "INCOMPLETE"
