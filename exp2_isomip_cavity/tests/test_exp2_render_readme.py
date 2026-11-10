"""Tests for exp2_isomip_cavity/scripts/render_readme.py (Task 7): README numbers come only from
the figures/*.json metrics files, never from the void pre-R21 records."""
from __future__ import annotations

import re

import pytest

from exp2_isomip_cavity.scripts import render_readme as rr


def test_typographic_minus_signs_only():
    assert rr.typographic_minus("-1.5") == "−1.5"
    assert rr.typographic_minus("(-2.5") == "(−2.5"
    assert rr.typographic_minus("0.50-0.61") == "0.50-0.61"
    assert rr.typographic_minus("ens_T-0.5_zc250") == "ens_T-0.5_zc250"


def test_every_placeholder_has_a_value_from_the_metrics_files():
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


def test_void_records_are_never_read(monkeypatch):
    seen = []
    real = rr._load

    def spy(name):
        seen.append(name)
        return real(name)

    monkeypatch.setattr(rr, "_load", spy)
    rr.values()
    assert seen and not [n for n in seen if "void" in n or "/" in n]
    with pytest.raises(ValueError, match="void"):
        real("void_pre_r21/verify_4km.json")


def test_formatters():
    assert rr.pct(-0.00049271785512044, 2) == "-0.05"
    assert rr.pct(1.2458158225697127, 0, sign=True) == "+125"
    assert rr.rng([2.40, 2.48]) == "2.40–2.48"


def test_committed_readme_is_up_to_date():
    assert rr.OUT.read_text() == rr.render(allow_pending=False)
