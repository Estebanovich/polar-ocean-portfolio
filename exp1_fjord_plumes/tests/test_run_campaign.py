"""Tests for exp1_fjord_plumes/scripts/run_campaign.py (non-MITgcm logic, Task 4).

MITgcm itself is never run here: the build and run shell commands go through an
injected ``sh`` recorder and ``pgrep`` is injected as a callable.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("run_campaign", PKG / "scripts" / "run_campaign.py")
rc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rc)

A_TWOLAYER = ["a_line2d_twolayer_q0.01", "a_line2d_twolayer_q0.03", "a_line2d_twolayer_q0.1"]


class Recorder:
    """Fake shell: records commands; a run command writes an OK run_status.json."""

    def __init__(self, state="OK"):
        self.calls = []
        self.state = state

    def __call__(self, cmd):
        self.calls.append([str(c) for c in cmd])
        if Path(cmd[0]).name == "run_mitgcm.sh":
            run_dir = Path(cmd[3])
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "run_status.json").write_text(json.dumps(
                {"state": self.state, "last_iter": 10, "n_timesteps": 10, "wall_s": 1}))
        return 0


def _kinds(calls):
    return [Path(c[0]).name for c in calls]


# ---- case globbing ------------------------------------------------------------

def test_select_cases_glob_in_cases_yaml_order():
    assert rc.select_cases("a_line2d", "a_line2d_twolayer_*") == A_TWOLAYER
    assert len(rc.select_cases("a_line2d", "*")) == 6


def test_select_cases_only_from_requested_config():
    # b_ cases match "*" but belong to another configuration
    assert all(n.startswith("a_line2d") for n in rc.select_cases("a_line2d", "*"))
    with pytest.raises(ValueError, match="b_3d_twolayer_x0_0"):
        rc.select_cases("a_line2d", "b_3d_twolayer_x0_0")


def test_select_cases_empty_glob_raises():
    with pytest.raises(ValueError, match="no case"):
        rc.select_cases("a_line2d", "a_line2d_nonexistent_*")


# ---- refuse to start if another MITgcm is running ------------------------------

def test_other_processes_ignores_empty_and_self():
    assert rc.other_mitgcm_processes(pgrep=lambda: "") == []
    lines = rc.other_mitgcm_processes(pgrep=lambda: "123 ./mitgcmuv\n")
    assert lines == ["123 ./mitgcmuv"]


def test_refuses_to_start_when_mitgcmuv_running(tmp_path):
    sh = Recorder()
    with pytest.raises(RuntimeError, match="mitgcmuv"):
        rc.run_campaign("a_line2d", "a_line2d_twolayer_*", 1,
                        runs_root=tmp_path / "runs", build_root=tmp_path / "builds",
                        sh=sh, pgrep=lambda: "4242 ./mitgcmuv")
    assert sh.calls == []                      # nothing built, nothing run
    assert not (tmp_path / "runs").exists()    # nothing written


# ---- build once, then run in series -------------------------------------------

def test_build_once_then_runs_in_series(tmp_path):
    sh = Recorder()
    out = rc.run_campaign("a_line2d", "a_line2d_twolayer_*", 1,
                          runs_root=tmp_path / "runs", build_root=tmp_path / "builds",
                          sh=sh, pgrep=lambda: "")
    kinds = _kinds(sh.calls)
    assert kinds == ["build_mitgcm.sh"] + ["run_mitgcm.sh"] * 3
    build = sh.calls[0]
    assert build[1].endswith("config/a_line2d") and build[2] == str(tmp_path / "builds" / "a_line2d")
    assert "--mpi" not in build                # np = 1: serial build
    for name, call in zip(A_TWOLAYER, sh.calls[1:]):
        # run_mitgcm.sh <build_dir> <input_dir> <run_dir> <np>
        assert call[1] == str(tmp_path / "builds" / "a_line2d")
        assert call[2] == str(tmp_path / "runs" / name / "input")
        assert call[3] == str(tmp_path / "runs" / name)
        assert call[4] == "1"
        assert (tmp_path / "runs" / name / "input" / "data").exists()  # generated before the run
    assert out == {n: "OK" for n in A_TWOLAYER}


def test_mpi_build_flag_when_np_gt_1(tmp_path):
    sh = Recorder()
    rc.run_campaign("a_line2d", "a_line2d_twolayer_q0.01", 4,
                    runs_root=tmp_path / "runs", build_root=tmp_path / "builds",
                    sh=sh, pgrep=lambda: "")
    assert sh.calls[0][-1] == "--mpi" and sh.calls[1][-1] == "4"


def test_np_above_five_refused(tmp_path):
    with pytest.raises(ValueError, match="5"):
        rc.run_campaign("a_line2d", "*", 6, runs_root=tmp_path, build_root=tmp_path,
                        sh=Recorder(), pgrep=lambda: "")


def test_failed_case_is_reported_and_series_continues(tmp_path):
    sh = Recorder(state="FAILED")
    out = rc.run_campaign("a_line2d", "a_line2d_twolayer_*", 1,
                          runs_root=tmp_path / "runs", build_root=tmp_path / "builds",
                          sh=sh, pgrep=lambda: "")
    assert out == {n: "FAILED" for n in A_TWOLAYER}
    assert _kinds(sh.calls).count("run_mitgcm.sh") == 3


def test_overrides_and_suffix_for_short_run(tmp_path):
    sh = Recorder()
    rc.run_campaign("a_line2d", "a_line2d_twolayer_q0.03", 1,
                    runs_root=tmp_path / "runs", build_root=tmp_path / "builds",
                    sh=sh, pgrep=lambda: "", overrides={"endTime": 7200.0},
                    run_suffix="_short2h")
    inp = tmp_path / "runs" / "a_line2d_twolayer_q0.03_short2h" / "input"
    meta = json.loads((inp / "case.json").read_text())
    assert meta["params"]["nTimeSteps"] == round(7200.0 / meta["params"]["deltaT"])
    assert f"nTimeSteps={meta['params']['nTimeSteps']}," in (inp / "data").read_text()


# ---- SIZE.h must match the generated grid -------------------------------------

def test_check_size_matches_real_config(tmp_path):
    from exp1_fjord_plumes.gendata.make_case import make_case
    d = make_case("a_line2d_twolayer_q0.03", tmp_path / "in")
    rc.check_size(d, PKG / "config" / "a_line2d" / "code" / "SIZE.h")  # no error


def test_check_size_mismatch_raises(tmp_path):
    from exp1_fjord_plumes.gendata.make_case import make_case
    d = make_case("a_line2d_twolayer_q0.03", tmp_path / "in")
    real = (PKG / "config" / "a_line2d" / "code" / "SIZE.h").read_text()
    bad = tmp_path / "SIZE.h"
    bad.write_text(real.replace("sNx = 193", "sNx = 100"))
    with pytest.raises(RuntimeError, match="Nx"):
        rc.check_size(d, bad)
