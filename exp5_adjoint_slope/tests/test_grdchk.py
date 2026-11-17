"""analysis.grdchk.parse on real grdchk output (plan Task 3). The fragments are copied from
MITgcm verification/tutorial_global_oce_biogeo/results/output_tap_{adj,tlm}.txt."""
from pathlib import Path

import pytest

from exp5_adjoint_slope.analysis import grdchk

DATA = Path(__file__).parent / "data"


def test_parse_adjoint_fragment():
    pts = grdchk.parse((DATA / "grdchk_tap_adj.txt").read_text())
    assert [p["ind"] for p in pts] == [1, 2]
    p = pts[0]
    assert p["kind"] == "adj"
    assert (p["i"], p["j"], p["k"]) == (61, 5, 1)
    assert p["fc_ref"] == pytest.approx(-6.05596159367e11)
    assert p["fc_plus"] == pytest.approx(-6.05596159457e11)
    assert p["grad"] == pytest.approx(-8.96642030128e06)
    assert p["fd_grad"] == pytest.approx(-8.96647949219e06)
    assert p["rel_err"] == pytest.approx(6.60139790920e-06)
    assert pts[1]["rel_err"] == pytest.approx(1.81285748702e-06)


def test_parse_tlm_fragment():
    pts = grdchk.parse((DATA / "grdchk_tap_tlm.txt").read_text())
    assert len(pts) == 1 and pts[0]["kind"] == "tlm"
    assert pts[0]["rel_err"] == pytest.approx(6.60139793096e-06)


def test_parse_empty():
    assert grdchk.parse("no gradient check here\n") == []


def test_parse_merged_fixed_width_fields():
    # GS, xx_theta at (20, 18, 15), 2026-11-16: the 3-D control index 24260 runs into the
    # point number in the second line ("124260"), which the column regex missed.
    pts = grdchk.parse((DATA / "grdchk_exp5_theta_adj.txt").read_text())
    assert len(pts) == 1
    p = pts[0]
    assert (p["i"], p["j"], p["k"]) == (20, 18, 15)
    assert p["grad"] == pytest.approx(3.85026515048e-04)
    assert p["fd_grad"] == pytest.approx(3.87559974166e-04)
    assert p["rel_err"] == pytest.approx(6.57996012930e-03)
