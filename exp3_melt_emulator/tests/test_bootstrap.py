"""Bootstrap intervals over runs (ruling R17) on synthetic result entries: resampling only,
no refits. numpy only (runs under pyocean and anemoi)."""
import json

import numpy as np
import pytest

from exp3_melt_emulator.analysis import bootstrap as bs
from exp3_melt_emulator.tests.test_hypotheses import METHODS, _entry, _results


def _cv_final():
    cv, fin = _results({"M6": (0.10, 0.20, 0.30), "M7": (0.05, 0.10, 0.40)})
    rng = np.random.default_rng(3)
    for e in cv["entries"] + fin["entries"]:          # make the metrics vary between runs
        e["metrics"]["rmse_map"] = float(rng.uniform(1, 5))
        e["metrics"]["rmse_gl"] = float(rng.uniform(0, 2))
    return cv, fin


def test_point_estimates_are_the_plain_means():
    cv, fin = _cv_final()
    r = bs.bootstrap(cv, fin, n_boot=200, seed=0)
    x = r["variants"]["main"]["loro"]["M4"]
    runs = [e for e in cv["entries"] if e["variant"] == "main" and e["kind"] == "loro"
            and e["method"] == "M4"]
    assert x["n_runs"] == len(runs)
    assert x["rmse_map"]["mean"] == pytest.approx(np.mean([e["metrics"]["rmse_map"] for e in runs]))
    gl = np.array([e["metrics"]["rmse_gl"] for e in runs])
    assert x["rmse_gl"]["mean"] == pytest.approx(gl.mean())
    assert x["rmse_gl_rms"]["mean"] == pytest.approx(np.sqrt((gl ** 2).mean()))
    for key in bs.STATS:
        lo, hi = x[key]["ci95"]
        assert lo <= x[key]["mean"] <= hi


def test_final_comes_from_final_test_and_only_in_main():
    cv, fin = _cv_final()
    r = bs.bootstrap(cv, fin, n_boot=200, seed=0)
    x = r["variants"]["main"]["final"]["M6"]
    assert x["n_runs"] == 3
    assert x["abs_integrated_rel_error"]["mean"] == pytest.approx(0.2)
    assert "final" not in r["variants"]["n1"]


def test_constant_metric_has_zero_width_and_seed_is_reproducible():
    cv, fin = _cv_final()
    for e in cv["entries"]:
        if e["kind"] == "zc_out" and e["method"] == "M2":
            e["metrics"]["rmse_map"] = 2.5
    a = bs.bootstrap(cv, fin, n_boot=300, seed=0)
    assert a["variants"]["main"]["zc_out"]["M2"]["rmse_map"]["ci95"] == [2.5, 2.5]
    assert a == bs.bootstrap(cv, fin, n_boot=300, seed=0)
    assert a != bs.bootstrap(cv, fin, n_boot=300, seed=1)


def test_resampling_is_paired_across_methods():
    """The same resampled runs are used for every method: identical per-run values give
    identical intervals, whatever the method."""
    cv, fin = _cv_final()
    for e in cv["entries"]:
        if e["kind"] == "loro" and e["variant"] == "main" and e["method"] == "M1":
            twin = next(t for t in cv["entries"] if t["kind"] == "loro" and t["variant"] == "main"
                        and t["method"] == "M0" and t["case"] == e["case"])
            e["metrics"] = dict(twin["metrics"])
    r = bs.bootstrap(cv, fin, n_boot=500, seed=0)["variants"]["main"]["loro"]
    assert r["M1"]["rmse_map"]["ci95"] == r["M0"]["rmse_map"]["ci95"]


def test_kinds_with_one_run_are_skipped():
    cv, fin = _cv_final()
    for m in METHODS:
        cv["entries"].append({**_entry("loro", m, "ref", 0.1), "kind": "profile_test",
                              "fold": "profile_test", "variant": "main"})
    r = bs.bootstrap(cv, fin, n_boot=50, seed=0)
    assert "profile_test" not in r["variants"]["main"]


def test_write(tmp_path):
    cv, fin = _cv_final()
    (tmp_path / "cv.json").write_text(json.dumps(cv))
    (tmp_path / "fin.json").write_text(json.dumps(fin))
    rec = bs.write(tmp_path / "cv.json", tmp_path / "fin.json", tmp_path / "b.json", n_boot=100)
    assert json.loads((tmp_path / "b.json").read_text()) == rec
    assert rec["n_boot"] == 100 and rec["seed"] == bs.SEED and "analyst" in rec["choice"]
