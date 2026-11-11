"""reproduce.sh full mode never leaves a committed record rewritten (review Important 1):
analysis/compare_records.py and scripts/rerun_compare_restore.sh on synthetic records."""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from exp3_melt_emulator.analysis import compare_records as cr

SCRIPT = Path(cr.__file__).resolve().parents[1] / "scripts" / "rerun_compare_restore.sh"


def _cv(rmse=1.0, wall=10.0, created="2026-11-11T09:14:27"):
    return {"created": created, "session": {"wall_s": wall, "units_resumed": 9},
            "entries": [{"variant": "main", "fold": "loro_a", "method": "M4", "case": "a",
                         "metrics": {"rmse_map": rmse}}],
            "units": [{"unit": "main__loro_a__M4", "started": "t0", "finished": "t1",
                       "fit_s": 1.0, "wall_s": wall, "info": {"degree": 2}}]}


def test_compare_cv_ignores_session_metadata_only():
    assert cr.compare_cv(_cv(), _cv(wall=0.01, created="2030-01-01")) == []
    assert cr.compare_cv(_cv(), _cv(rmse=1.5)) == ["entry main/loro_a/M4/a differs"]
    other = _cv()
    other["units"][0]["info"]["degree"] = 1
    assert cr.compare_cv(_cv(), other) == ["unit main__loro_a__M4 differs"]


def test_compare_profile_ignores_fit_s():
    a = {"x": 1, "M7_refit_decomposition": {"fit_s": 78.0, "n": [2.4]}}
    b = {"x": 1, "M7_refit_decomposition": {"fit_s": 81.0, "n": [2.4]}}
    assert cr.compare_profile(a, b) == []
    b["M7_refit_decomposition"]["n"] = [2.5]
    assert cr.compare_profile(a, b) == ["M7_refit_decomposition differs"]


def _run(tmp_path, new_record: dict | None, fail: bool = False):
    rec = tmp_path / "cv_results.json"
    rec.write_text(json.dumps(_cv()))
    committed = rec.read_bytes()
    code = ("import json,sys; open(sys.argv[1],'w').write(sys.argv[2]); sys.exit(sys.argv[3] == '1')")
    cmd = [sys.executable, "-c", code, str(rec), json.dumps(new_record or {}), "1" if fail else "0"]
    env = {**os.environ, "PY": sys.executable, "RERUN_DIR": str(tmp_path)}
    p = subprocess.run(["bash", str(SCRIPT), "cv", str(rec), "--", *cmd], capture_output=True,
                       text=True, env=env)
    assert rec.read_bytes() == committed, "the committed record was not restored"
    return p


@pytest.mark.parametrize("new, code", [(_cv(wall=0.02, created="2030-01-01"), 0), (_cv(rmse=9.0), 1)])
def test_rerun_is_compared_and_committed_record_restored(tmp_path, new, code):
    p = _run(tmp_path, new)
    assert p.returncode == code, p.stdout + p.stderr
    assert json.loads((tmp_path / "cv_results.json.rerun").read_text()) == new


def test_failed_rerun_restores_and_propagates(tmp_path):
    p = _run(tmp_path, _cv(rmse=3.0), fail=True)
    assert p.returncode == 1 and "restored" in p.stderr
