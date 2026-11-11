"""Coverage, CRPS and skill summaries (pre-registered metrics, PD0 "Métricas") on synthetic
entries; the per-run values come from cv_results.json / final_test.json, nothing is refitted."""
import json

import pytest

from exp3_melt_emulator.analysis import probabilistic as pr
from exp3_melt_emulator.tests.test_hypotheses import _results


def _data():
    cv, fin = _results({"M6": (0.10, 0.20, 0.30)})
    for i, e in enumerate(cv["entries"] + fin["entries"]):
        e["skill"] = {"rmse_map": 0.1 * (i % 3), "rmse_gl": 0.2, "abs_integrated_rel_error": -0.5}
        if e["method"] in pr.PROB_METHODS:
            e["prob"] = {"coverage_68": 0.5 + 0.1 * (i % 2), "coverage_95": 0.9, "crps": 1.0 + i % 4,
                         "mean_sigma": 2.0}
    return cv, fin


def test_means_per_variant_kind_method():
    cv, fin = _data()
    r = pr.summarise(cv, fin)
    ents = [e for e in cv["entries"] if e["variant"] == "main" and e["kind"] == "loro"
            and e["method"] == "M6"]
    x = r["variants"]["main"]["loro"]["M6"]
    assert x["n_runs"] == len(ents)
    assert x["coverage_68"] == pytest.approx(sum(e["prob"]["coverage_68"] for e in ents) / len(ents))
    assert x["crps"] == pytest.approx(sum(e["prob"]["crps"] for e in ents) / len(ents))
    assert x["skill_rmse_map"] == pytest.approx(sum(e["skill"]["rmse_map"] for e in ents) / len(ents))
    assert "coverage_68" not in r["variants"]["main"]["loro"]["M2"]
    assert r["variants"]["main"]["final"]["M6"]["n_runs"] == 3
    assert "M0" not in r["variants"]["main"]["loro"]          # skill reference


def test_write(tmp_path):
    cv, fin = _data()
    (tmp_path / "cv.json").write_text(json.dumps(cv))
    (tmp_path / "f.json").write_text(json.dumps(fin))
    rec = pr.write(tmp_path / "cv.json", tmp_path / "f.json", tmp_path / "p.json")
    assert json.loads((tmp_path / "p.json").read_text()) == rec and "definition" in rec
