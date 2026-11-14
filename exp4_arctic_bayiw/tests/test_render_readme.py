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


# ---- final review: fixed prose words are checked against the JSON (Minor 11) ------------------

@pytest.mark.parametrize("name, mutate, what", [
    ("hypotheses.json", lambda d: d["H4a"]["bands_out_J"].update(low=1.0), "H4a"),
    ("hypotheses.json", lambda d: d["H4d"]["R1k_vs_R1"].update(outcome="refuted"), "H4d"),
    ("hypotheses.json", lambda d: d["H4c"]["i"].update(ratio=2.0), "neither tail"),
    ("hypotheses.json", lambda d: d["post_hoc"]["bridge"]["R2"]["flux_share_signed"].update(low=-0.1), "R2"),
    ("validation.json", lambda d: d["V3"].update({"pass": False}), "V3 passes"),
    ("posthoc_moving_reference.json", lambda d: d["runs"]["F1c"].update(low_band_out_J=-1.0), "F1c"),
])
def test_prose_mismatch_refuses_to_render(monkeypatch, name, mutate, what):
    real = rr._load

    def fake(n):
        d = real(n)
        if n == name:
            mutate(d)
        return d

    monkeypatch.setattr(rr, "_load", fake)
    with pytest.raises(ValueError, match=what):
        rr.values()


def test_amendment_1_ratio_is_rendered_from_the_json():
    v = rr.values()
    assert v["a1_R_ratio"] == "-0.17" and v["a1_R_outcome"] == "refuted"
    assert (v["r1_F_out"], v["r1_A_in"]) == ("0.28", "5.68")
    assert "`fjord_vmf_400`" in v["mod_configs"] and "`fjord_vmf_column`" in v["mod_configs"]


# ---- JPO-style rewrite (2026-11-14): guards behind the new fixed words and cross-references ----

@pytest.mark.parametrize("name, mutate, what", [
    ("hypotheses.json", lambda d: d["H4b"]["i"]["post_hoc_ke_bc_only"].update(ratio=0.9), "either number"),
    ("hypotheses.json", lambda d: d["H4d"]["F1_vs_F1c"].update(eta_kpp=-0.5), "small against the wind work"),
    ("hypotheses.json", lambda d: d["H4b"]["i"].update(ratio=0.4), "no simple fraction"),
    ("hypotheses.json", lambda d: d["H4d"]["F1_vs_F1c"]["numerator_J"].update(kpp=1.0), "H4d"),
    ("validation.json", lambda d: d["V4"]["V4a"].update(F1_out_J=1.0), "V4"),
    ("pd2_v4_diagnostics.json", lambda d: [x.update(gross_J=1.0) for x in d["F1"].values()], "gross exchange"),
])
def test_new_prose_guards_refuse_to_render(monkeypatch, name, mutate, what):
    real = rr._load

    def fake(n):
        d = real(n)
        if n == name:
            mutate(d)
        return d

    monkeypatch.setattr(rr, "_load", fake)
    with pytest.raises(ValueError, match=what):
        rr.values()


def test_fraction_words():
    assert rr.frac_words(0.228) == "about a quarter"
    assert rr.values()["H4b_frac_words"] == "about a quarter"
    with pytest.raises(ValueError):
        rr.frac_words(0.4)


def test_cross_references_resolve():
    rr.check_cross_refs(rr.TEMPLATE.read_text())
    head = "## 2. Methods\n\n### a. Model\n\n## Appendix A: Validation\n\n"
    rr.check_cross_refs(head + "See section 2a, appendix A and Gill 1982, §10.3.\n")
    for bad in ("see section 2b.", "see section 3a.", "see appendix B.", "as in §2.6.", "sections 2a and 2c"):
        with pytest.raises(ValueError, match="cross-references"):
            rr.check_cross_refs(head + bad + "\n")
