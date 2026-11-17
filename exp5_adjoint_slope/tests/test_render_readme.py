"""scripts/render_readme.py (plan Task 7): every number in the README comes from a committed
metrics file or cases.yaml; unknown placeholders and missing metrics files fail."""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("exp5_render_readme", PKG / "scripts" / "render_readme.py")
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)


def test_every_placeholder_has_a_value():
    text = rr.TEMPLATE.read_text()
    assert text.count("{{") == len(re.findall(r"\{\{(\w+)\}\}", text)), "malformed placeholder"
    missing = sorted(set(re.findall(r"\{\{(\w+)\}\}", text)) - set(rr.values()))
    assert not missing, missing


def test_unknown_placeholder_raises(tmp_path, monkeypatch):
    t = tmp_path / "README.template.md"
    t.write_text("value {{no_such_value}}\n")
    monkeypatch.setattr(rr, "TEMPLATE", t)
    with pytest.raises(KeyError, match="no_such_value"):
        rr.render(allow_pending=False)


def test_missing_metrics_file_raises_unless_pending(tmp_path, monkeypatch):
    monkeypatch.setattr(rr, "FIG", tmp_path)
    with pytest.raises(FileNotFoundError):
        rr.values()
    assert "v1_n" not in rr.values(allow_pending=True)


def test_template_has_no_hand_typed_numbers():
    left = rr.strip_allowed(rr.TEMPLATE.read_text())
    digits = re.findall(r".{0,30}\d.{0,30}", left)
    assert not digits, digits


def test_committed_readme_is_the_rendering():
    assert rr.OUT.read_text() == rr.render(allow_pending=False)


def test_crosscheck_values_fill_h5b_and_v5(tmp_path, monkeypatch):
    import json
    import shutil
    for f in rr.FIG.glob("*.json"):
        shutil.copy(f, tmp_path / f.name)
    (tmp_path / "crosscheck.json").write_text(json.dumps({
        "J_control_TW": {"FC1": 4e-3, "FC2": 3e-3, "FC3": 1e-3},
        "J_pert_TW": {"FP1": 3.7e-3, "FP2": 4.1e-3, "FP3": 9.6e-3},
        "H5b": {"dj": 3.1e-3, "sigma_delta": 2.1e-3, "dj_pred": 1.4e-3, "ratio": 2.2,
                "detectable": False, "sign_ok": True, "within_factor": False, "pass": False,
                "indicative": False, "n": 3},
        "dj_pred_se_TW": 1.6e-4, "V5_eddy_fraction": {"FC1": 1.6, "FC2": 2.1, "FC3": 1.8},
        "EKE_spinup_m2s2": {"F1": [1e-2, 3.5e-2, 7e-4], "F1b": [6e-4, 5e-4, 4e-4]},
        "F1_equilibrium": {"last": 6.6e-4, "prev": 7.6e-4, "rel_change": 0.134, "pass": False},
        "F1b_equilibrium": {"last": 4.6e-4, "prev": 5.9e-4, "rel_change": 0.216, "pass": False}}))
    monkeypatch.setattr(rr, "FIG", tmp_path)
    v = rr.values()
    assert v["h5b_verdict"] == "not detectable with N = 3"
    assert v["h5b_dj"] == "3.10" and v["h5b_2sd"] == "4.20" and v["h5b_ratio"] == "2.20"
    assert v["h5b_jc_min"] == "1.00" and v["h5b_jp_max"] == "9.60"
    assert v["h5b_n"] == "three"
    assert v["v5_min"] == "160" and v["v5_max"] == "210"
    assert v["eke_peak"] == "3.5 × 10⁻²" and v["eke_final"] == "4.6 × 10⁻⁴"
    assert v["eke_drift_f1"] == "13" and v["eke_drift_f1b"] == "22"
    assert v["v5_verdict"] == "passes"
