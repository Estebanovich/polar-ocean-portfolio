import json

import pytest

from exp6_southern_carbon.analysis import hypotheses

Y = (1900, 1980, 2020, 2060)


def test_assemble_is_json_ready_and_carries_verdicts():
    s_gt = {1900: -1.0, 1980: -0.5, 2020: 0.1, 2060: 0.5}
    s_gc = {1900: 0.0, 1980: 0.05, 2020: 0.0, 2060: 0.1}
    pred = {y: 1.0 for y in Y}
    meas = {1900: 1.1, 1980: 0.9, 2020: 1.2, 2060: 1.0}
    terms = {y: {"mld": 0.5, "moc": 0.7, "bio": -0.2, "sum_circ_bio": 1.0} for y in Y}
    out = hypotheses.assemble(s_gt, s_gc, pred, meas, terms)
    json.dumps(out)                                     # keys are strings, values plain floats/bools
    assert out["h6a"]["pass"] and out["h6a"]["dS"] == pytest.approx(1.5)
    assert out["h6a"]["dS_pct_of_S1900"] == pytest.approx(150.0)
    assert out["h6b"]["pass"] and out["s_gt"]["2060"] == 0.5 and out["terms"]["1980"]["moc"] == 0.7


def test_assemble_h6b_fails_with_one_branch_off():
    s = {y: 0.0 for y in Y}
    s[2060] = 1.0
    meas = {y: 1.0 for y in Y}
    pred = {**meas, 2020: 1.5}
    out = hypotheses.assemble(s, {y: 0.0 for y in Y}, pred, meas, {y: {} for y in Y})
    assert not out["h6b"]["pass"] and out["h6b"]["n_within"] == 3
