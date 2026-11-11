"""Mechanical evaluation of H3a-H3c (pd0_preregistro.md) on synthetic result entries."""
import json

import numpy as np
import pytest

from exp3_melt_emulator.analysis import hypotheses as hy

METHODS = tuple(f"M{i}" for i in range(8))


def _entry(kind, method, case, ire, rmse=1.0, gl=0.5, spread=None, variant="main",
           n1_pass=True):
    """One result entry in the cv_results.json / final_test.json layout."""
    true = 1000.0
    pred = true * (1.0 + ire)
    e = {"case": case, "method": method, "n1_pass": n1_pass,
         "metrics": {"rmse_map": rmse, "integrated_rel_error": ire, "rmse_gl": gl,
                     "integrated_pred": pred, "integrated_true": true}}
    if kind == "final":
        e["fold"] = "final"                       # final_test.json entries carry no 'kind'
    else:
        e.update(kind=kind, fold=kind, variant=variant, label=variant)
    if method in ("M6", "M7"):
        s = 0.1 if spread is None else spread
        # 5 members with population std = s * |pred| and mean = pred
        dev = np.array([-1.0, -1.0, 0.0, 1.0, 1.0]) * s * abs(pred) / np.sqrt(0.8)
        e["members"] = {"integrated": list(pred + dev), "integrated_mixture": pred,
                        "spread_rel": s}
    return e


def _results(final_ire, loro_rmse=None, h3c_points=None, n_zc=4, n_warm=2):
    """Synthetic cv + final dicts. ``final_ire`` maps method -> 3 integrated rel. errors;
    ``loro_rmse`` maps method -> rmse_map of all loro runs; ``h3c_points`` maps
    (method, case) -> (spread, ire) for the extrapolation cases."""
    loro_rmse = loro_rmse or {}
    h3c_points = h3c_points or {}
    cv = []
    for variant in ("main", "n1"):
        for m in METHODS:
            for i in range(5):
                cv.append(_entry("loro", m, f"loro_{i}", 0.05, rmse=loro_rmse.get(m, 2.0),
                                 variant=variant, n1_pass=i != 0))
            for kind, n in (("zc_out", n_zc), ("warm_out", n_warm)):
                for i in range(n):
                    case = f"{kind}_{i}"
                    sp, ire = h3c_points.get((m, case), (0.2, 0.05))
                    cv.append(_entry(kind, m, case, ire, spread=sp, variant=variant))
    final = []
    for m in METHODS:
        for i, ire in enumerate(final_ire.get(m, (0.1, 0.1, 0.1))):
            case = f"final_{i}"
            sp, ire = h3c_points.get((m, case), (0.2, ire))
            final.append(_entry("final", m, case, ire, spread=sp))
    return {"entries": cv, "amendments": [{"text": "a"}]}, {"entries": final}


def test_h3a_pass():
    cv, fin = _results({"M6": (0.10, 0.20, 0.30), "M7": (0.05, 0.10, 0.40)})
    r = hy.evaluate(cv, fin, expected=None)["H3a"]
    assert r["mean_abs_ire"]["M7"] == pytest.approx(0.55 / 3)
    assert r["mean_abs_ire"]["M6"] == pytest.approx(0.20)
    assert r["n_runs_M7_better"] == 2
    assert r["verdict"] == "PASS"


def test_h3a_fail_by_run_count():
    # M7 has the lower mean but is better in only 1 of 3 runs
    cv, fin = _results({"M6": (0.10, 0.10, 0.50), "M7": (0.12, 0.12, 0.05)})
    r = hy.evaluate(cv, fin, expected=None)["H3a"]
    assert r["mean_abs_ire"]["M7"] < r["mean_abs_ire"]["M6"]
    assert r["n_runs_M7_better"] == 1
    assert r["verdict"] == "FAIL"


def test_h3a_fail_by_mean():
    cv, fin = _results({"M6": (0.05, 0.05, 0.05), "M7": (-0.20, 0.10, 0.30)})
    assert hy.evaluate(cv, fin, expected=None)["H3a"]["verdict"] == "FAIL"


def test_h3a_inconclusive_below_floor():
    # mean difference 1.0 % < 1.3 %: inconclusive even though M7 is better in 3 of 3
    cv, fin = _results({"M6": (0.10, 0.10, 0.10), "M7": (0.09, 0.09, -0.09)})
    r = hy.evaluate(cv, fin, expected=None)["H3a"]
    assert r["n_runs_M7_better"] == 3
    assert r["abs_mean_difference"] == pytest.approx(0.01)
    assert r["verdict"] == "INCONCLUSIVE"


def test_h3b_threshold():
    cv, fin = _results({}, loro_rmse={"M6": 4.0, "M7": 4.4})          # ratio 1.10 exactly
    assert hy.evaluate(cv, fin, expected=None)["H3b"]["verdict"] == "PASS"
    cv, fin = _results({}, loro_rmse={"M6": 4.0, "M7": 4.41})
    r = hy.evaluate(cv, fin, expected=None)["H3b"]
    assert r["ratio_M7_over_M6"] == pytest.approx(4.41 / 4.0)
    assert r["verdict"] == "FAIL"


def test_h3b_ignores_n1_variant():
    cv, fin = _results({}, loro_rmse={"M6": 4.0, "M7": 4.0})
    for e in cv["entries"]:
        if e["variant"] == "n1" and e["method"] == "M7":
            e["metrics"]["rmse_map"] = 100.0
    assert hy.evaluate(cv, fin, expected=None)["H3b"]["verdict"] == "PASS"


def test_h3c_violation_and_confusion():
    pts = {("M7", "zc_out_1"): (0.03, 0.25),          # violation: low spread, large error
           ("M7", "final_2"): (0.08, 0.30),           # not H3c (spread >= 5 %), O26 miss
           ("M6", "warm_out_0"): (0.04, 0.20),        # error not > 20 %: no violation
           ("M6", "zc_out_0"): (0.30, 0.50)}          # flagged by spread: O26 hit
    cv, fin = _results({}, h3c_points=pts)
    res = hy.evaluate(cv, fin, expected=None)["H3c"]
    m7, m6 = res["M7"], res["M6"]
    assert m7["n_cases"] == m6["n_cases"] == 4 + 2 + 3
    assert m7["verdict"] == "FAIL"
    assert [v["case"] for v in m7["violations"]] == ["zc_out_1"]
    assert m6["verdict"] == "PASS" and m6["violations"] == []
    c7 = m7["ockenden_10pct_confusion"]
    assert c7["spread_le_10pct_error_gt_20pct"] == 2
    assert sorted(c7["cases_low_spread_high_error"]) == ["final_2", "zc_out_1"]
    assert m6["ockenden_10pct_confusion"]["spread_gt_10pct_error_gt_20pct"] == 1
    tot = sum(v for k, v in c7.items() if k.startswith("spread_"))
    assert tot == 9


def test_h3c_spread_is_recomputed_and_checked():
    cv, fin = _results({})
    e = next(x for x in fin["entries"] if x["method"] == "M7")
    e["members"]["spread_rel"] = 0.5                 # inconsistent with the member values
    with pytest.raises(ValueError, match="spread"):
        hy.evaluate(cv, fin, expected=None)


def test_expected_counts_enforced():
    cv, fin = _results({})
    with pytest.raises(ValueError, match="H3c"):
        hy.evaluate(cv, fin, expected={"final": 3, "loro": 5, "h3c": 33})
    hy.evaluate(cv, fin, expected={"final": 3, "loro": 5, "h3c": 9})


def test_ranking_tables():
    cv, fin = _results({"M2": (0.1, -0.2, 0.3)}, loro_rmse={"M3": 1.0})
    for e in cv["entries"]:
        if e["kind"] == "loro" and e["method"] == "M3":
            e["metrics"]["rmse_gl"] = 3.0 if e["case"] == "loro_0" else 4.0
    r = hy.evaluate(cv, fin, expected=None)
    t = r["ranking"]["main"]["loro"]
    assert t["order_by_mean_rmse_map"][0] == "M3"
    assert t["methods"]["M3"]["rms_rmse_gl"] == pytest.approx(np.sqrt((9 + 4 * 16) / 5))
    assert r["ranking"]["main"]["final"]["methods"]["M2"]["mean_abs_ire"] == pytest.approx(0.2)
    assert r["ranking"]["n1"]["label"] == "sensitivity"
    assert "final" not in r["ranking"]["n1"]["kinds"]
    json.dumps(r)                                     # serialisable
