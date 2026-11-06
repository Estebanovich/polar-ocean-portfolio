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


# ------------------------------------------------------------ fix round 1 (R27)
def test_richardson_monotone_divergence_is_inconclusive_R27():
    # R_Q 0.40 -> 0.50 -> 0.70 under refinement: increments GROW (|eps32/eps21| = 0.5 < 1)
    r = ev.richardson(0.70, 0.50, 0.40)
    assert not r["conclusive"] and r["R_ext"] is None
    assert "monotone divergence" in r["reason"]
    assert ev.classify(True, (U, U), (D, D), False, r) == "P2"


def test_richardson_ratio_exactly_one_is_divergence_R27():
    r = ev.richardson(0.75, 0.5, 0.25)  # eps21 = eps32 = -0.25 exactly
    assert not r["conclusive"] and "monotone divergence" in r["reason"]


def test_richardson_nonconvergence_branch():
    r = ev.richardson(*_phi(1.5), max_iter=1)
    assert not r["conclusive"] and r["p"] is None and "did not converge" in r["reason"]


def test_stationarity_near_two_percent():
    assert not ev.is_stationary(0.50, 0.5105, 0.5)   # 0.021
    assert ev.is_stationary(0.50, 0.5095, 0.5)       # 0.019


def test_width_fixed_guard_zero_width():
    fixed, v = ev.width_fixed_flag(10.0, 0.0, 5.0)
    assert fixed is False and not all(map(math.isfinite, v))
    fixed, v = ev.width_fixed_flag(float("nan"), 5.0, 2.0)
    assert fixed is False and not all(map(math.isfinite, v))


def test_assemble_zero_width_is_V():
    res = ev.assemble({"conv_dx10": _m(0.4, 38.0), "conv_dx5": _m(0.5, 0.0), "conv_dx2p5": _m(0.6, 14.0)})
    assert res["convergence"]["class"] == "V"
    assert any("width" in x for r in res["convergence"]["invalid_runs"].values() for x in r)


def test_assemble_failed_run_without_metrics_gives_V_and_json(tmp_path):
    bad = {"run": "conv_dx5", "valid": False, "reasons": ["run state is FAILED"]}
    res = ev.assemble({"conv_dx10": _m(0.4, 38.0), "conv_dx5": bad, "conv_dx2p5": _m(0.6, 14.0)})
    assert res["convergence"]["class"] == "V"
    assert res["convergence"]["invalid_runs"]["conv_dx5"] == ["run state is FAILED"]
    # end to end: three FAILED run dirs -> V and a written JSON
    import json

    dirs = {}
    for k in ev.CONV:
        d = tmp_path / k
        d.mkdir()
        (d / "run_status.json").write_text(json.dumps(
            {"state": "FAILED", "last_iter": 1, "n_timesteps": 10, "wall_s": 1}))
        dirs[k] = d
    out = tmp_path / "out"
    r = ev.evaluate_study(dirs, out)
    assert r["convergence"]["class"] == "V"
    assert json.loads((out / "pd2_convergence.json").read_text())["convergence"]["class"] == "V"


def test_assemble_exact_class_and_floor_pin():
    metrics = {"conv_dx10": _m(0.40, 38.0), "conv_dx5": _m(0.55, 25.0), "conv_dx2p5": _m(0.68, 14.0)}
    assert ev.assemble(metrics)["convergence"]["class"] == "P2"  # Richardson inconclusive (p < 0.5)
    # b threshold: 2 sigma_AB = 0.85 m < 2.5 m floor -> the floor decides
    a, b = _m(0.5, 24.0), _m(0.5, 26.4)
    assert ev._diff(a, b)["b_trend"] == "plano"
    assert ev._diff(a, _m(0.5, 26.6))["b_trend"] == "sube"
    assert ev._diff(a, b)["thr_b"] == 2.5


def test_p1_carries_dx_dz_caveat():
    metrics = {"conv_dx10": _m(0.30, 38.0), "conv_dx5": _m(0.40, 25.0), "conv_dx2p5": _m(0.45, 12.0)}
    c = ev.assemble(metrics)["convergence"]
    assert c["class"] == "P1"
    assert c["numerical_part_R_ext_minus_RQ_dx5"] == pytest.approx(c["richardson"]["R_ext"] - 0.40)
    assert "mixes dx and dz refinement" in c["numerical_part_caveat"]


def test_make_figure_smoke(tmp_path):
    metrics = {"conv_dx10": _m(0.30, 38.0), "conv_dx5": _m(0.40, 25.0), "conv_dx2p5": _m(0.45, 12.0),
               ev.REF: _m(0.41, 24.0), "sens_freeslip": _m(0.5, 22.0)}
    for m in metrics.values():
        m["theory_b_lower_half"] = 7.5
    out = ev.make_figure(ev.assemble(metrics), tmp_path / "pd2_convergence")
    assert all(p.exists() and p.stat().st_size > 0 for p in out)


def test_metric_windows_default_matches_protocol():
    from exp1_fjord_plumes.analysis.pd2_convergence_eval import metric_windows, WINDOW
    w = metric_windows()
    assert w["X"] == WINDOW == (12.0, 24.0)
    assert w["first"] == (12.0, 18.0) and w["second"] == (18.0, 24.0)
    assert w["samples"][0] == (12.0, 13.0) and w["samples"][-1] == (23.0, 24.0)
    assert len(w["samples"]) == 12


def test_metric_windows_shifted_for_extended_run():
    from exp1_fjord_plumes.analysis.pd2_convergence_eval import metric_windows
    w = metric_windows(24.0)
    assert w["X"] == (24.0, 36.0)
    assert w["first"] == (24.0, 30.0) and w["second"] == (30.0, 36.0)
    assert w["samples"][-1] == (35.0, 36.0)
