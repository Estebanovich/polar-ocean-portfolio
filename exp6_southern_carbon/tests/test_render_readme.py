"""scripts/render_readme.py (plan Task 10): every number in the README comes from a committed
metrics file, cases.yaml or the analysis constants; unknown placeholders and missing metrics
files fail unless --allow-pending."""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("exp6_render_readme", PKG / "scripts" / "render_readme.py")
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)

Y = ("1900", "1980", "2020", "2060")


def _synthetic(tmp):
    terms = {y: {"mld": 1e10, "moc": 3e10, "bio": -1e10, "sum_circ_bio": 3e10, "temp": -2e9,
                 "conductance": 5e4, "parent": {"g200": 0.3}, "branch": {}} for y in Y}
    hyp = {"h6a": {"dS": 0.08, "drift_gc": 0.004, "pass": True, "dS_pct_of_S1900": 80.0,
                   "s_nat": {y: -0.1 for y in Y}, "s_ant": {y: 0.02 for y in Y}, "zero_crossing_year": None},
           "h6b": {"within": {y: True for y in Y}, "n_within": 4, "pass": True, "ratio": {y: 1.07 for y in Y}},
           "s_gt": {y: -0.08 for y in Y}, "s_gc": {y: -0.1 for y in Y}, "pred": {y: 3e10 for y in Y},
           "meas": {y: 2.8e10 for y in Y}, "terms": terms}
    rep = {"years": 1000.0, "v6a": {"d_flux_century": 2.4e-5, "d_deep_dic_century": 1.4e-7, "pass": True},
           "moc_R_sv": 0.375, "moc_max_sv": 0.375, "upwelling_fraction_R": 0.813, "mld_winter_R_m": 752.4,
           "mld_month_max": 9, "flux_last_mol_m2_yr": 0.448}
    files = {"spinup": {"G0": rep, "G0s": rep}, "validation": {"V6b": {"max_frac": 3e-5, "pass": True}},
             "hypotheses": hyp,
             "figdata": {"sections": {"g200_R": {"1900": 0.3, "2060": 0.2}}},
             "crosscheck": {"H6c": {"dS": 0.07, "sigma": 0.01, "detectable": True, "sign_ok": True,
                                    "ratio": 0.9, "pass": True, "indicative": False, "n": 3}}}
    for n, d in files.items():
        (tmp / f"{n}.json").write_text(json.dumps(d))


def test_every_placeholder_has_a_value(tmp_path, monkeypatch):
    _synthetic(tmp_path)
    monkeypatch.setattr(rr, "FIG", tmp_path)
    text = rr.TEMPLATE.read_text()
    assert text.count("{{") == len(re.findall(r"\{\{(\w+)\}\}", text)), "malformed placeholder"
    missing = sorted(set(re.findall(r"\{\{(\w+)\}\}", text)) - set(rr.values()))
    assert not missing, missing


def test_synthetic_values(tmp_path, monkeypatch):
    _synthetic(tmp_path)
    monkeypatch.setattr(rr, "FIG", tmp_path)
    v = rr.values()
    assert v["dx_coarse_km"] == "50" and v["dx_fine_km"] == "10" and v["Lx_km"] == "1000" and v["Ly_km"] == "2000"
    assert v["R_rows"] == "19" and v["R_edge_km"] == "1000" and v["win_lo"] == "6" and v["win_hi"] == "15"
    assert v["v6a_dflux"] == "2.4 × 10⁻⁵" and v["v6a_verdict"] == "passes" and v["mld_month"] == "September"
    assert v["v6a_ddic"] == "1.4 × 10⁻⁴" and v["moc_edge_frac"] == "100" and v["upw_frac"] == "81"
    assert v["dS"] == "+0.080" and v["h6a_verdict"] == "passes" and v["zero_cross"] == "keeps its sign in every branch"
    assert v["g200_drop"] == "33" and v["mld_1"] == "+10.0" and v["temp_4"] == "-2.0" and v["ratio_2"] == "1.07"
    assert v["h6c_verdict"] == "passes" and v["n_members"] == "three" and v["accel"] == "4"


def test_unknown_placeholder_raises(tmp_path, monkeypatch):
    _synthetic(tmp_path)
    monkeypatch.setattr(rr, "FIG", tmp_path)
    t = tmp_path / "README.template.md"
    t.write_text("value {{no_such_value}}\n")
    monkeypatch.setattr(rr, "TEMPLATE", t)
    with pytest.raises(KeyError, match="no_such_value"):
        rr.render(allow_pending=False)


def test_missing_metrics_file_raises_unless_pending(tmp_path, monkeypatch):
    monkeypatch.setattr(rr, "FIG", tmp_path)
    with pytest.raises(FileNotFoundError):
        rr.values()
    v = rr.values(allow_pending=True)
    assert "dS" not in v and v["dx_coarse_km"] == "50"
    assert rr.PENDING in rr.render(allow_pending=True)


def test_template_has_no_hand_typed_numbers():
    left = rr.strip_allowed(rr.TEMPLATE.read_text())
    digits = re.findall(r".{0,30}\d.{0,30}", left)
    assert not digits, digits


def test_committed_readme_is_the_rendering():
    if not rr.OUT.exists():
        pytest.skip("README.md not rendered yet")
    pending = not all((rr.FIG / f"{n}.json").exists() for n in rr.SOURCES)
    assert rr.OUT.read_text() == rr.render(allow_pending=pending)
