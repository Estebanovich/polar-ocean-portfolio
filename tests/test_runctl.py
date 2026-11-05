import json
import subprocess
from pathlib import Path

import pytest

from pfcommon.runctl import require_ok, run_status

ROOT = Path(__file__).resolve().parents[1]
END = "PROGRAM MAIN: Execution ended Normally\n"


def _fake_run(tmp_path, state):
    (tmp_path / "STDOUT.0000").write_text("(PID.TID 0000.0001) time_tsnumber = 20\n" + END)
    (tmp_path / "run_status.json").write_text(
        json.dumps({"state": state, "last_iter": 20, "n_timesteps": 20, "wall_s": 1.5})
    )
    return tmp_path


def test_status_ok(tmp_path):
    d = _fake_run(tmp_path, "OK")
    assert run_status(d)["state"] == "OK"
    require_ok(d)


def test_status_failed_refused(tmp_path):
    d = _fake_run(tmp_path, "FAILED")
    with pytest.raises(RuntimeError, match="FAILED"):
        require_ok(d)


def test_status_missing_refused(tmp_path):
    with pytest.raises(RuntimeError):
        require_ok(tmp_path)


def test_runner_refuses_unmounted_lacie(tmp_path):
    target = Path("/Volumes/LaCie_NOPE/x")
    r = subprocess.run(
        ["bash", str(ROOT / "scripts" / "run_mitgcm.sh"), str(tmp_path), str(tmp_path), str(target), "1"],
        capture_output=True, text=True,
    )
    assert r.returncode == 2
    assert "not mounted" in r.stderr
    assert not target.exists()
