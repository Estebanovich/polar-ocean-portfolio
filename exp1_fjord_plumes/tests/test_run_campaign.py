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


# ---- --skip-ok: resume a campaign after the first-case gate (R30) --------------

def _fake_build(build_dir: Path, mpi: bool) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    (build_dir / "Makefile").write_text("FC = mpif77\nDEFINES = -DALLOW_USE_MPI\n" if mpi
                                        else "FC = gfortran\nDEFINES = \n")
    exe = build_dir / "mitgcmuv"
    exe.write_text("")
    import os
    import time
    t = time.time() + 10  # newer than every file in config/<cfg>/code
    os.utime(exe, (t, t))


def test_skip_ok_skips_finished_case_and_up_to_date_build(tmp_path):
    runs, builds = tmp_path / "runs", tmp_path / "builds"
    sh = Recorder()
    rc.run_campaign("a_line2d", "a_line2d_twolayer_q0.01", 1, runs_root=runs,
                    build_root=builds, sh=sh, pgrep=lambda: "")
    _fake_build(builds / "a_line2d", mpi=False)
    sh2 = Recorder()
    out = rc.run_campaign("a_line2d", "a_line2d_twolayer_*", 1, runs_root=runs,
                          build_root=builds, sh=sh2, pgrep=lambda: "", skip_ok=True)
    kinds = _kinds(sh2.calls)
    assert kinds == ["run_mitgcm.sh"] * 2          # no rebuild, first case skipped
    assert [c[3] for c in sh2.calls] == [str(runs / n) for n in A_TWOLAYER[1:]]
    assert out == {n: "OK" for n in A_TWOLAYER}


def test_skip_ok_reruns_failed_case(tmp_path):
    runs, builds = tmp_path / "runs", tmp_path / "builds"
    rc.run_campaign("a_line2d", "a_line2d_twolayer_q0.01", 1, runs_root=runs,
                    build_root=builds, sh=Recorder(state="FAILED"), pgrep=lambda: "")
    _fake_build(builds / "a_line2d", mpi=False)
    sh2 = Recorder()
    rc.run_campaign("a_line2d", "a_line2d_twolayer_q0.01", 1, runs_root=runs,
                    build_root=builds, sh=sh2, pgrep=lambda: "", skip_ok=True)
    assert _kinds(sh2.calls) == ["run_mitgcm.sh"]


def test_skip_ok_rebuilds_when_binary_stale_or_wrong_mpi(tmp_path):
    builds = tmp_path / "builds"
    _fake_build(builds / "a_line2d", mpi=False)
    sh = Recorder()
    rc.run_campaign("a_line2d", "a_line2d_twolayer_q0.01", 2, runs_root=tmp_path / "runs",
                    build_root=builds, sh=sh, pgrep=lambda: "", skip_ok=True)
    assert _kinds(sh.calls)[0] == "build_mitgcm.sh"   # serial binary, np = 2 wants MPI
    import os
    exe = builds / "a_line2d" / "mitgcmuv"
    _fake_build(builds / "a_line2d", mpi=True)
    os.utime(exe, (0, 0))                              # older than the code/ files
    sh2 = Recorder()
    rc.run_campaign("a_line2d", "a_line2d_twolayer_q0.03", 2, runs_root=tmp_path / "runs",
                    build_root=builds, sh=sh2, pgrep=lambda: "", skip_ok=True)
    assert _kinds(sh2.calls)[0] == "build_mitgcm.sh"


def test_skip_ok_refuses_ok_run_with_other_parameters(tmp_path):
    runs, builds = tmp_path / "runs", tmp_path / "builds"
    rc.run_campaign("a_line2d", "a_line2d_twolayer_q0.01", 1, runs_root=runs,
                    build_root=builds, sh=Recorder(), pgrep=lambda: "",
                    overrides={"endTime": 7200.0})
    sh2 = Recorder()
    with pytest.raises(RuntimeError, match="endTime"):
        rc.run_campaign("a_line2d", "a_line2d_twolayer_q0.01", 1, runs_root=runs,
                        build_root=builds, sh=sh2, pgrep=lambda: "", skip_ok=True)
    assert sh2.calls == []


def test_cli_has_skip_ok():
    assert "--skip-ok" in Path(rc.__file__).read_text()


# ---- two-stage cases of (c): stage 2 restarts from the stage-1 pickups --------
# The real scripts/run_mitgcm.sh runs a stub mitgcmuv that reads nIter0 and
# nTimeSteps from `data`, refuses a restart without pickup.<nIter0> and the
# Orlanski pickup in its cwd, and writes the pickups of its last iteration.

STUB = r'''#!/bin/bash
get() { grep -v '^ *#' data | grep -iE "^ *$1 *=" | head -1 | sed -E 's/^[^=]*= *//; s/[ ,].*$//'; }
N0=$(get nIter0); NT=$(get nTimeSteps)
if [ "$N0" -gt 0 ]; then
  P=$(printf %010d "$N0")
  for f in pickup.$P.data pickup.$P.meta pickup_orlanskiE.$P.001.001.data; do
    [ -f "$f" ] && [ ! -L "$f" ] || { echo "missing or linked $f"; exit 4; }
  done
fi
END=$((N0 + NT)); E=$(printf %010d "$END")
echo " time_tsnumber = $END"
echo "pickup" > pickup.$E.data; echo "meta" > pickup.$E.meta
[ -n "$STUB_NO_ORLANSKI" ] || echo "orl" > pickup_orlanskiE.$E.001.001.data
[ -n "$STUB_FAIL_STAGE1" ] && [ "$N0" -eq 0 ] && exit 3
echo "PROGRAM MAIN: Execution ended Normally"
'''


def _stub_build(build_dir: Path) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    exe = build_dir / "mitgcmuv"
    exe.write_text(STUB)
    exe.chmod(0o755)


class StubShell:
    """Build commands are no-ops (stub mitgcmuv already in place); run commands
    go through the REAL scripts/run_mitgcm.sh."""

    def __init__(self, env=None):
        self.calls, self.env = [], env or {}

    def __call__(self, cmd):
        import os
        import subprocess
        self.calls.append([str(c) for c in cmd])
        if Path(cmd[0]).name == "build_mitgcm.sh":
            return 0
        return subprocess.run(["bash"] + [str(c) for c in cmd], capture_output=True,
                              env={**os.environ, **self.env}).returncode


SHORT2 = {"endTime": 21600.0, "wind_stage1_s": 10800.0}   # 180 + 180 steps of 60 s


def _two_stage(tmp_path, env=None, case="c_fjord_wind"):
    builds = tmp_path / "builds"
    _stub_build(builds / "c_fjord3d")
    sh = StubShell(env)
    out = rc.run_campaign("c_fjord3d", case, 1, runs_root=tmp_path / "runs",
                          build_root=builds, sh=sh, pgrep=lambda: "", overrides=SHORT2,
                          run_suffix="_t")
    return out, tmp_path / "runs" / f"{case}_t", sh


def test_two_stage_case_restarts_from_copied_pickups(tmp_path):
    out, run, sh = _two_stage(tmp_path)
    assert out == {"c_fjord_wind": "OK"}
    runs = [c for c in sh.calls if Path(c[0]).name == "run_mitgcm.sh"]
    assert [c[3] for c in runs] == [str(run), str(run / "stage2")]
    inp2 = run / "input_stage2"
    assert "nIter0=180," in (inp2 / "data").read_text() and not (inp2 / "data.stage2").exists()
    for f in ("pickup.0000000180.data", "pickup.0000000180.meta",
              "pickup_orlanskiE.0000000180.001.001.data"):
        for d in (inp2, run / "stage2"):
            assert (d / f).is_file() and not (d / f).is_symlink()
    s1 = json.loads((run / "run_status.json").read_text())
    s2 = json.loads((run / "stage2" / "run_status.json").read_text())
    fin = json.loads((run / "run_status_final.json").read_text())
    assert s1["state"] == s2["state"] == fin["state"] == "OK"
    assert s1["last_iter"] == 180 and s2["last_iter"] == 360
    assert fin["last_iter"] == fin["expected_last_iter"] == 360
    assert fin["stages"]["1"] == s1 and fin["stages"]["2"] == s2


def test_stage2_not_run_when_stage1_failed(tmp_path):
    out, run, sh = _two_stage(tmp_path, env={"STUB_FAIL_STAGE1": "1"})
    assert out == {"c_fjord_wind": "FAILED"}
    assert not (run / "stage2").exists()
    fin = json.loads((run / "run_status_final.json").read_text())
    assert fin["state"] == "FAILED" and "stage 1" in fin["reason"]


def test_stage2_refused_without_orlanski_pickup(tmp_path):
    """C3 lesson from (a): the restart is bit-exact only with the Orlanski pickup."""
    out, run, sh = _two_stage(tmp_path, env={"STUB_NO_ORLANSKI": "1"})
    assert out == {"c_fjord_wind": "FAILED"}
    fin = json.loads((run / "run_status_final.json").read_text())
    assert "orlanski" in fin["reason"].lower() and not (run / "stage2").exists()


def test_final_state_needs_stage2_last_iter(tmp_path):
    """Stage-2 run_status OK but last_iter != nIter0 + nTimeSteps -> FAILED."""
    runs, builds = tmp_path / "runs", tmp_path / "builds"

    class Rec(Recorder):
        def __call__(self, cmd):
            r = super().__call__(cmd)
            if Path(cmd[0]).name == "run_mitgcm.sh":
                d = Path(cmd[3])
                for f in ("pickup.0000000180.data", "pickup.0000000180.meta",
                          "pickup_orlanskiE.0000000180.001.001.data"):
                    (d / f).write_text("x")
                (d / "run_status.json").write_text(json.dumps(
                    {"state": "OK", "last_iter": 180 if d.name != "stage2" else 300,
                     "n_timesteps": 180, "wall_s": 1}))
            return r

    out = rc.run_campaign("c_fjord3d", "c_fjord_wind", 1, runs_root=runs, build_root=builds,
                          sh=Rec(), pgrep=lambda: "", overrides=SHORT2)
    assert out == {"c_fjord_wind": "FAILED"}
    fin = json.loads((runs / "c_fjord_wind" / "run_status_final.json").read_text())
    assert fin["expected_last_iter"] == 360 and fin["last_iter"] == 300


def test_single_stage_c_case_has_no_stage2(tmp_path):
    out, run, sh = _two_stage(tmp_path, case="c_fjord_plume")
    assert out == {"c_fjord_plume": "OK"}
    assert not (run / "stage2").exists() and not (run / "input_stage2").exists()
    assert json.loads((run / "run_status_final.json").read_text())["state"] == "OK"


def test_skip_ok_two_stage_uses_final_status(tmp_path):
    out, run, sh = _two_stage(tmp_path)
    assert rc.finished_ok(run)
    fin = json.loads((run / "run_status_final.json").read_text())
    fin["state"] = "FAILED"
    (run / "run_status_final.json").write_text(json.dumps(fin))
    assert not rc.finished_ok(run)          # stage 1 OK alone is not enough


def test_cli_set_overrides():
    assert rc.parse_sets(["wind_stage1_s=43200", "endTime=86400.5"]) == {
        "wind_stage1_s": 43200, "endTime": 86400.5}
    with pytest.raises(ValueError):
        rc.parse_sets(["novalue"])


def test_set_rejects_unknown_key_and_keeps_int_types(tmp_path):
    names = ["c_fjord_all"]
    with pytest.raises(ValueError, match="not a parameter"):
        rc.check_overrides(names, {"wind_stage_1": 43200})
    ov = rc.check_overrides(names, rc.parse_sets(["n_plume_cols=6", "wind_stage1_s=43200",
                                                  "endTime=86400"]))
    assert ov == {"n_plume_cols": 6, "wind_stage1_s": 43200.0, "endTime": 86400.0}
    assert isinstance(ov["n_plume_cols"], int) and isinstance(ov["wind_stage1_s"], float)
    with pytest.raises(ValueError, match="integer"):
        rc.check_overrides(names, {"n_plume_cols": 6.5})
    with pytest.raises(ValueError, match="not a parameter"):
        rc.run_campaign("c_fjord3d", "c_fjord_all", 1, runs_root=tmp_path, build_root=tmp_path,
                        sh=Recorder(), pgrep=lambda: "", overrides={"bogus": 1})
