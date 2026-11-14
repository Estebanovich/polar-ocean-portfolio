"""exp4_arctic_bayiw/scripts/render_readme.py (Task 8): every number in the README comes from a
committed metrics file (figures/*.json, or the case registry cases.yaml for the configuration),
never typed by hand; the committed README is the rendering of the template."""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parents[1]


def _rr():
    spec = importlib.util.spec_from_file_location("exp4_render_readme", PKG / "scripts" / "render_readme.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


rr = _rr()


def test_every_placeholder_has_a_value():
    text = rr.TEMPLATE.read_text()
    assert text.count("{{") == len(re.findall(r"\{\{(\w+)\}\}", text)), "malformed placeholder"
    names = set(re.findall(r"\{\{(\w+)\}\}", text))
    assert names, "template has no placeholders"
    missing = sorted(names - set(rr.values()))
    assert not missing, missing


def test_unknown_placeholder_raises(tmp_path, monkeypatch):
    t = tmp_path / "README.template.md"
    t.write_text("value {{no_such_value}}\n")
    monkeypatch.setattr(rr, "TEMPLATE", t)
    with pytest.raises(KeyError, match="no_such_value"):
        rr.render(allow_pending=False)


def test_every_number_from_json():
    """Outside placeholders the template has no digit except the whitelisted labels (run and
    criterion names, equation, section and figure labels, citation years, amendment numbers).
    Code spans, URLs, comments and the References section are skipped."""
    left = rr.strip_allowed(rr.TEMPLATE.read_text())
    offending = sorted({m.group(0) for m in re.finditer(r"\S*[0-9]\S*", left)})
    assert not offending, offending


def test_strip_allowed_catches_a_typed_result():
    assert re.search(r"[0-9]", rr.strip_allowed("H4b(i) is met with a ratio of 0.228"))
    assert re.search(r"[0-9]", rr.strip_allowed("the fjord is 119 m deep"))
    assert not re.search(r"[0-9]", rr.strip_allowed(
        "H4a, H4b(ii), F1c, R1k, V2a and Figure 3 (Gill, 1982; Eqs. 13–17; Amendment 2) {{x}}"))


def test_typographic_minus_signs_only():
    assert rr.typographic_minus("-1.17") == "−1.17"
    assert rr.typographic_minus("(-0.75") == "(−0.75"
    assert rr.typographic_minus("1.7-6.6") == "1.7-6.6"
    assert rr.typographic_minus("F1_out-1") == "F1_out-1"


def test_formatters():
    assert rr.gj(-1171008778.97) == "-1.17"
    assert rr.gj(6744107.9, 4) == "0.0067"
    assert rr.pct(-0.0606703, 2) == "-6.07"
    assert rr.rng([0.35, 0.52]) == "0.35–0.52"


def test_values_follow_the_json(monkeypatch):
    """A changed metrics value changes the README: nothing is cached or typed."""
    v = rr.values()
    assert v["H4a_low_share"] == "0.930"
    assert v["H4d_F1_eta"] == "-6.07"
    real = rr._load

    def fake(name):
        d = real(name)
        if name == "hypotheses.json":
            d["H4a"]["low_band_share"] = 0.5
        return d

    monkeypatch.setattr(rr, "_load", fake)
    assert rr.values()["H4a_low_share"] == "0.500"


def test_committed_readme_is_up_to_date():
    assert rr.OUT.read_text() == rr.render(allow_pending=False)


# ---- Task 8 fix round 1 -------------------------------------------------------------------------

def test_missing_bridge_station_raises(monkeypatch):
    real = rr._load

    def fake(name):
        d = real(name)
        if name == "hypotheses.json":
            d["post_hoc"]["bridge"]["R2"]["station_fraction"].pop("Seiche")
        return d

    monkeypatch.setattr(rr, "_load", fake)
    with pytest.raises(KeyError, match="Seiche"):
        rr.values()


def test_counts_in_words_follow_the_json(monkeypatch):
    v = rr.values()
    assert (v["n_failed"], v["n_failed_Word"], v["n_hyp_Word"], v["n_resolved"]) == ("three", "Three", "Four", "four")
    real = rr._load

    def fake(name):
        d = real(name)
        if name == "validation.json":
            d["summary"]["failed"] = ["V2a", "V4"]
        return d

    monkeypatch.setattr(rr, "_load", fake)
    assert rr.values()["n_failed_Word"] == "Two"


def test_factor_bounds_not_typed():
    assert "1/1.5" not in rr.TEMPLATE.read_text()
    assert re.search(r"[0-9]", rr.strip_allowed("between 1/1.5 and 1.5"))
