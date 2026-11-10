"""Ruling of the final review (Task 7): the analysis defaults of the Exp. 2 run root follow
pfcommon.paths.EXP2_RUNS (environment variable EXP2_RUNS), never a hard-coded disk path."""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parents[1] / "analysis"
REPO = ANALYSIS.parents[1]

PROBE = r"""
import json, sys
sys.path.insert(0, %r)
from exp2_isomip_cavity.analysis import ensemble, fig1_ocean0_validation as f1, jpl_compare as jc
from exp2_isomip_cavity.analysis import ensemble_n_check as nc, pd1_calibration_eval as ce
from exp2_isomip_cavity.analysis import export_exp3 as ex
print(json.dumps({"ensemble": str(ensemble.DEFAULT_RUNS), "fig1": str(f1.DEFAULT_RUN),
                  "jpl": str(jc.RUNS), "ncheck": str(nc.DEFAULT_RUNS),
                  "caleval": str(ce.DEFAULT_RUNS), "export": str(ex.DEFAULT_OUT)}))
"""


def test_defaults_follow_exp2_runs(tmp_path):
    root = tmp_path / "elsewhere"
    env = {**os.environ, "EXP2_RUNS": str(root)}
    out = subprocess.run([sys.executable, "-c", PROBE % str(REPO)], env=env, cwd=tmp_path,
                         capture_output=True, text=True, check=True).stdout
    d = json.loads(out.strip().splitlines()[-1])
    assert d["fig1"] == str(root / "ocean1_2km_r21")
    assert d["export"] == str(root / "exp3_dataset.npz")
    for k in ("ensemble", "jpl", "ncheck", "caleval"):
        assert d[k] == str(root), k


def test_no_hardcoded_run_disk_in_analysis_code():
    bad = []
    for f in ANALYSIS.glob("*.py"):
        for i, line in enumerate(f.read_text().splitlines(), 1):
            code = line.split("#", 1)[0]
            if re.search(r"Path\(\s*[\"']/Volumes/LaCie", code):
                bad.append(f"{f.name}:{i}")
    assert not bad, bad


def test_jpl_compare_takes_runs_root(tmp_path):
    from exp2_isomip_cavity.analysis import jpl_compare as jc
    assert jc.case_run("ocean1", tmp_path) == tmp_path / "ocean1_2km_r21"
    assert jc.case_run("ocean0", tmp_path) == tmp_path / "ocean0_2km_r21"
