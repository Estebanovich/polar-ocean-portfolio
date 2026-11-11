"""Tests for exp3_melt_emulator/scripts/render_readme.py (Task 10): every number of the README
comes from the figures/*.json metrics files, never typed by hand. Pure standard library: runs
under both pyocean and anemoi."""
from __future__ import annotations

import re
import subprocess
import sys

import pytest

from exp3_melt_emulator.scripts import render_readme as rr


def test_render_check_passes():
    """The committed README equals the rendered template, and the CLI --check agrees."""
    assert rr.OUT.read_text() == rr.render(allow_pending=False)
    p = subprocess.run([sys.executable, str(rr.EXP / "scripts" / "render_readme.py"), "--check"],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr


def test_unknown_placeholder_raises(tmp_path, monkeypatch):
    t = tmp_path / "README.template.md"
    t.write_text("value {{no_such_value}}\n")
    monkeypatch.setattr(rr, "TEMPLATE", t)
    with pytest.raises(KeyError, match="no_such_value"):
        rr.render(allow_pending=False)


def test_every_placeholder_has_a_value():
    text = rr.TEMPLATE.read_text()
    assert text.count("{{") == len(re.findall(r"\{\{(\w+)\}\}", text)), "malformed placeholder"
    names = set(re.findall(r"\{\{(\w+)\}\}", text))
    assert names, "template has no placeholders"
    assert not sorted(names - set(rr.values()))


def test_every_number_from_json():
    """Outside placeholders the template has no digits except the whitelisted tokens
    (rr.NUMBER_WHITELIST: equation and figure labels, method/hypothesis names, citation years,
    the two Exp. 2 design constants). Code spans, URLs, comments and the References section
    (bibliographic data) are not results and are skipped."""
    left = rr.strip_allowed(rr.TEMPLATE.read_text())
    offending = sorted({m.group(0) for m in re.finditer(r"\S*[0-9]\S*", left)})
    assert not offending, offending


def test_strip_allowed_catches_a_typed_result():
    assert re.search(r"[0-9]", rr.strip_allowed("M4 wins with rmse_map 3.83 m a-1"))
    assert not re.search(r"[0-9]", rr.strip_allowed(
        "M4 (Burgard et al., 2022, Eq. 14) in Figure 3c, §2.1; {{x_1}} `fig1.json`"))


def test_typographic_minus_signs_only():
    assert rr.typographic_minus("-4.8") == "−4.8"
    assert rr.typographic_minus("ens_T-0.5_zc250") == "ens_T-0.5_zc250"
    assert rr.typographic_minus("[-250, -150)") == "[−250, −150)"


def test_formatters():
    assert rr.pct(-0.04819419273906036, 1, sign=True) == "-4.8"
    assert rr.pct(4.129321558054767, 0, sign=True) == "+413"
    assert rr.f(3.8295206720983903) == "3.83"
