"""run_campaign.py --dry-run plan: run directories, stages, skip decisions and whether a build
would be triggered, without generating inputs, building or running anything."""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

from pfcommon import campaign as cp

PKG = Path(__file__).resolve().parents[1]


def _rc():
    spec = importlib.util.spec_from_file_location("exp4_rc_plan", PKG / "scripts" / "run_campaign.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _fake_build(build_root: Path, cfg: str, mpi=True, newer=True):
    b = build_root / cfg
    b.mkdir(parents=True)
    (b / "Makefile").write_text("DEFINES = -DALLOW_USE_MPI\n" if mpi else "DEFINES =\n")
    exe = b / cp.EXE_NAME
    exe.write_bytes(b"x")
    newest = max(f.stat().st_mtime for f in (PKG / "config" / cfg / "code").iterdir() if f.is_file())
    t = newest + 10 if newer else newest - 10
    os.utime(exe, (t, t))


def test_plan_reports_build_and_stages_for_a_smoke_run(tmp_path):
    rc = _rc()
    p = rc.plan_campaign("bay_thesis_kpp", "R[23]", 4, overrides={"t_end_h": 2.0},
                         run_suffix="_smoke2h", runs_root=tmp_path / "runs",
                         build_root=tmp_path / "b", skip_ok=True)
    assert p["build"]["action"] == "build"                      # no executable yet
    assert [c["case"] for c in p["cases"]] == ["R2", "R3"]
    c = p["cases"][0]
    assert c["run_dir"] == str(tmp_path / "runs" / "R2_smoke2h")
    assert c["t_end_h"] == 2.0 and c["stages"] == [[0, 240]] and c["action"] == "run"
    assert not (tmp_path / "runs").exists()                      # nothing written


def test_plan_no_rebuild_when_up_to_date_with_skip_ok_only(tmp_path):
    rc = _rc()
    _fake_build(tmp_path / "b", "fjord_vmf")
    kw = dict(runs_root=tmp_path / "runs", build_root=tmp_path / "b")
    assert rc.plan_campaign("fjord_vmf", "F5", 5, skip_ok=True, **kw)["build"]["action"] == "none (up to date)"
    assert rc.plan_campaign("fjord_vmf", "F5", 5, skip_ok=False, **kw)["build"]["action"] == "build"
    old = tmp_path / "old"
    _fake_build(old, "fjord_vmf", newer=False)
    assert rc.plan_campaign("fjord_vmf", "F5", 5, skip_ok=True, runs_root=tmp_path / "runs",
                            build_root=old)["build"]["action"] == "build"


def test_plan_skips_a_finished_run_with_the_same_parameters(tmp_path):
    rc = _rc()
    _fake_build(tmp_path / "b", "fjord_vmf")
    r = tmp_path / "runs" / "F5"
    (r / "input").mkdir(parents=True)
    spec = rc.load_cases()
    params = {**spec["configs"]["fjord_vmf"], **spec["cases"]["F5"]}
    (r / "input" / "case.json").write_text(json.dumps({"params": params}))
    (r / "run_status_final.json").write_text(json.dumps({"state": "OK", "last_iter": 8400,
                                                          "expected_last_iter": 8400}))
    p = rc.plan_campaign("fjord_vmf", "F5", 5, skip_ok=True, runs_root=tmp_path / "runs",
                         build_root=tmp_path / "b")
    assert p["cases"][0]["action"] == "skip (already OK)"
    assert p["build"]["action"] == "none (nothing to run)"


def test_dry_run_prints_names_on_stdout_and_the_plan_on_stderr(capsys, tmp_path):
    rc = _rc()
    assert rc.main(["--config", "bay_thesis_kpp", "--cases", "R[23]", "--np", "4", "--end-time", "7200",
                    "--run-suffix", "_smoke2h", "--skip-ok", "--dry-run",
                    "--runs-root", str(tmp_path), "--build-root", str(tmp_path)]) == 0
    out, err = capsys.readouterr()
    assert out.split() == ["R2", "R3"]
    assert "R2_smoke2h" in err and "build:" in err and "stages [[0, 240]]" in err
