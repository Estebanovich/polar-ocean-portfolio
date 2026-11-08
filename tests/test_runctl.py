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


def test_short_run_marked_ok_refused(tmp_path):
    (tmp_path / "run_status.json").write_text(
        json.dumps({"state": "OK", "last_iter": 5, "n_timesteps": 20, "wall_s": 1.0}))
    with pytest.raises(RuntimeError, match="stopped"):
        require_ok(tmp_path)


@pytest.mark.parametrize("text", ["{not json", json.dumps({"last_iter": 1})])
def test_corrupt_status_refused(tmp_path, text):
    (tmp_path / "run_status.json").write_text(text)
    with pytest.raises(RuntimeError):
        require_ok(tmp_path)


def _run_stub(tmp_path, stub_body, data_text):
    build = tmp_path / "build"; build.mkdir()
    exe = build / "mitgcmuv"
    exe.write_text("#!/bin/bash\n" + stub_body)
    exe.chmod(0o755)
    inp = tmp_path / "input"; inp.mkdir()
    (inp / "data").write_text(data_text)
    run = tmp_path / "run"
    r = subprocess.run(["bash", str(ROOT / "scripts" / "run_mitgcm.sh"), str(build), str(inp), str(run), "1"],
                       capture_output=True, text=True)
    return r, json.loads((run / "run_status.json").read_text())


DATA = "&PARM03\n#nTimeSteps=999,\n nIter0=0,\n nTimeSteps=10,\n&\n"
STEPS = lambda n: "".join(f'echo " time_tsnumber                     =                    {i}"\n' for i in range(1, n + 1))
NORMAL = 'echo "PROGRAM MAIN: Execution ended Normally"\n'


def test_runner_short_run_with_normal_end_is_failed(tmp_path):
    r, st = _run_stub(tmp_path, STEPS(4) + NORMAL, DATA)
    assert r.returncode != 0
    assert st["state"] == "FAILED" and st["last_iter"] == 4 and st["n_timesteps"] == 10


def test_runner_abnormal_exit_is_failed(tmp_path):
    r, st = _run_stub(tmp_path, STEPS(10) + "exit 3\n", DATA)
    assert st["state"] == "FAILED"


def test_runner_complete_run_is_ok(tmp_path):
    r, st = _run_stub(tmp_path, STEPS(10) + NORMAL, DATA)
    assert r.returncode == 0
    assert st["state"] == "OK" and st["last_iter"] == 10
    require_ok(tmp_path / "run")


def test_runner_endtime_derives_steps(tmp_path):
    data = "&PARM03\n startTime=0.,\n endTime=200.,\n deltaT=20.,\n&\n"
    _, st = _run_stub(tmp_path, STEPS(10) + NORMAL, data)
    assert st["state"] == "OK" and st["n_timesteps"] == 10


def test_runner_unparsable_length_is_failed(tmp_path):
    _, st = _run_stub(tmp_path, STEPS(3) + NORMAL, "&PARM03\n&\n")
    assert st["state"] == "FAILED" and "reason" in st


NAN_MON = "(PID.TID 0000.0001) %MON dynstat_uvel_max             =                    NaN\n"


def test_status_ok_but_nan_in_log_refused(tmp_path):
    run = _fake_run(tmp_path, "OK")
    (run / "STDOUT.0000").write_text("(PID.TID 0000.0001) time_tsnumber = 20\n" + NAN_MON + END)
    with pytest.raises(RuntimeError, match="non-finite"):
        require_ok(run)


def test_runner_nan_in_log_is_failed(tmp_path):
    r, st = _run_stub(tmp_path, STEPS(10) + 'echo "%MON dynstat_eta_max = NaN"\n' + NORMAL, DATA)
    assert r.returncode != 0
    assert st["state"] == "FAILED" and "non-finite" in st["reason"]


# Ruling (2026-11-07, b_dx5_twolayer_x0_20): a run whose only defect is a missing normal-end line
# (MPI_Finalize failed after the last write) may be accepted by an explicit run_accepted.json.
NOEND = "no normal-end line in model log"


def _accept(d, **kw):
    rec = {"accepted": True, "ruling": "MPI_Finalize failed after the last write", "evidence": "sizes match"}
    rec.update(kw)
    (d / "run_accepted.json").write_text(json.dumps(rec))


def _noend_run(tmp_path, last_iter=20):
    (tmp_path / "STDOUT.0000").write_text("(PID.TID 0000.0001) time_tsnumber = 20\n")
    (tmp_path / "run_status.json").write_text(json.dumps(
        {"state": "FAILED", "last_iter": last_iter, "n_timesteps": 20, "wall_s": 1.0, "reason": NOEND}))
    return tmp_path


def test_noend_without_acceptance_refused(tmp_path):
    with pytest.raises(RuntimeError, match="FAILED"):
        require_ok(_noend_run(tmp_path))


def test_noend_with_acceptance_ok(tmp_path):
    d = _noend_run(tmp_path)
    _accept(d)
    require_ok(d)


def test_acceptance_ignored_for_short_run(tmp_path):
    d = _noend_run(tmp_path, last_iter=10)
    _accept(d)
    with pytest.raises(RuntimeError):
        require_ok(d)


def test_acceptance_ignored_for_other_failure(tmp_path):
    d = _fake_run(tmp_path, "FAILED")
    _accept(d)
    with pytest.raises(RuntimeError, match="FAILED"):
        require_ok(d)


def test_acceptance_needs_ruling_and_evidence(tmp_path):
    d = _noend_run(tmp_path)
    _accept(d, evidence="")
    with pytest.raises(RuntimeError, match="FAILED"):
        require_ok(d)


def test_acceptance_does_not_skip_nan_check(tmp_path):
    d = _noend_run(tmp_path)
    (d / "STDOUT.0000").write_text("(PID.TID 0000.0001) %MON theta_max = NaN\n")
    _accept(d)
    with pytest.raises(RuntimeError, match="non-finite"):
        require_ok(d)
