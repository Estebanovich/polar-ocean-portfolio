"""Logic of fig2_two_plume_melt: C4' edges, the both-stratification rule, fail-closed behaviour,
the pre-registered-window guard and the expected-case list (synthetic; no run data)."""
from __future__ import annotations

import json

import numpy as np
import pytest

from exp1_fjord_plumes.analysis import fig2_two_plume_melt as f2


@pytest.mark.parametrize("ratio,ok", [(0.4, True), (0.8, True), (0.3999, False), (0.8001, False), (0.52, True)])
def test_c4p_range_edges_inclusive(ratio, ok):
    # M0 = ratio * 2 * Ms exactly representable with Ms = 1
    assert f2.c4p_entry(2 * ratio, 1.0)["pass"] is ok


@pytest.mark.parametrize("m0,ms", [(float("nan"), 1.0), (1.0, float("nan")), (1.0, 0.0), (float("inf"), 1.0)])
def test_c4p_fail_closed(m0, ms):
    assert f2.c4p_entry(m0, ms)["pass"] is False


def test_both_stratifications_must_pass():
    ok = {"pass": True}
    bad = {"pass": False}
    v = f2.reduce_verdicts({"twolayer": ok, "sermilik": ok}, {"twolayer": ok, "sermilik": ok})
    assert v["C4"] == v["C4p"] == "PASS"
    v = f2.reduce_verdicts({"twolayer": ok, "sermilik": bad}, {"twolayer": bad, "sermilik": ok})
    assert v["C4"] == v["C4p"] == "FAIL"


@pytest.mark.parametrize("entry", [{}, {"error": "missing"}, {"pass": None}, {"pass": "yes"}])
def test_missing_or_malformed_stratification_fails_closed(entry):
    ok = {"pass": True}
    assert f2.reduce_verdicts({"twolayer": ok}, {"twolayer": ok, "sermilik": ok})["C4"] == "FAIL"
    assert f2.reduce_verdicts({"twolayer": ok, "sermilik": entry}, {"twolayer": ok, "sermilik": ok})["C4"] == "FAIL"


def test_expected_cases_from_registry():
    names = f2.expected_cases()
    assert len(names) == 14 and "b_3d_sermilik_single" in names and "b_3d_twolayer_x0_320" in names
    assert not any(("orlanski" in n or "short" in n or "stv0" in n) for n in names)


def test_missing_run_dir_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="b_3d_"):
        f2.collect(tmp_path)


def test_non_preregistered_window_gives_no_verdict(tmp_path, monkeypatch):
    cases = {}
    for s in f2.STRATS:
        for tag in ("single",) + tuple(f"x0_{x}" for x in f2.X0S):
            cases[(s, tag)] = {"run": f"b_3d_{s}_{tag}", "x0_actual": float(tag[3:]) if tag != "single" else 0.0,
                               "x0c": 57.14, "Q0_total": 20.0, "M_main": 1.0, "M_info": 1.0,
                               "rel_change_main_vs_info": 0.0}
    c4 = {"pass": True, "z_detach_model": 0.0, "post_hoc": {"pass_initial": True}}
    monkeypatch.setattr(f2, "collect", lambda root, window: cases)
    monkeypatch.setattr(f2, "c4_for", lambda run, window: dict(c4))
    monkeypatch.setattr(f2, "make_figure", lambda cases, window, stem, excess=None: None)
    # canonical window: verdict, canonical file name
    r = f2.main(tmp_path, f2.WIN_MAIN, tmp_path)
    assert r["verdict"]["C4"] == "PASS" and (tmp_path / "fig2_metrics.json").exists()
    (tmp_path / "fig2_metrics.json").unlink()
    # any other window: no verdict, no canonical file, labelled post_hoc_window, no pass flags
    r = f2.main(tmp_path, (1.0, 2.0), tmp_path)
    assert "verdict" not in r and r["post_hoc_window"] is True
    assert not (tmp_path / "fig2_metrics.json").exists()
    out = tmp_path / "fig2_post_hoc_window_1-2h_metrics.json"
    txt = out.read_text()
    assert "post_hoc_window" in out.name and '"pass' not in txt and "NOT the pre-registered" in txt
    assert json.loads(txt)["C4p"]["twolayer"]["ratio"] == pytest.approx(0.5)


def test_load_excess_needs_both_resolutions(tmp_path):
    assert f2.load_excess(tmp_path) is None
    (tmp_path / "fig2_excess_clean.json").write_text("{}")
    assert f2.load_excess(tmp_path) is None
    (tmp_path / "fig2_excess_dx5.json").write_text('{"R_E": {}}')
    e10, e5 = f2.load_excess(tmp_path)
    assert e10 == {} and e5 == {"R_E": {}}
