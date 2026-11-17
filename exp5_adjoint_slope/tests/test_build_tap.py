"""scripts/build_mitgcm_tap.sh (Exp. 5 Task 0): builds the forward, Tapenade adjoint or
Tapenade tangent-linear executable with <config>/code_tap as -mods (genmake2 -tap)."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "build_mitgcm_tap.sh"


def _run(args, tmp_path, dry=True):
    env = dict(os.environ)
    if dry:
        env["DRYRUN"] = "1"
    return subprocess.run(["bash", str(SCRIPT), *args], capture_output=True, text=True, env=env,
                          cwd=tmp_path)


def _cfg(tmp_path):
    cfg = tmp_path / "cfg"
    (cfg / "code_tap").mkdir(parents=True)
    return cfg


def test_no_args_prints_usage(tmp_path):
    r = _run([], tmp_path)
    assert r.returncode != 0 and "usage" in r.stderr


def test_dryrun_adj_builds_tap_adj(tmp_path):
    cfg = _cfg(tmp_path)
    r = _run([str(cfg), str(tmp_path / "b"), "adj"], tmp_path)
    assert r.returncode == 0, r.stderr
    out = r.stdout
    assert "genmake2" in out and "-tap" in out
    assert f"-mods {cfg.resolve()}/code_tap" in out
    assert "make tap_adj" in out and "mitgcmuv_tap_adj" in out
    assert "TAPENADE_HOME=" in out


def test_dryrun_tlm_and_fwd(tmp_path):
    cfg = _cfg(tmp_path)
    tlm = _run([str(cfg), str(tmp_path / "b1"), "tlm"], tmp_path).stdout
    assert "make tap_tlm" in tlm
    fwd = _run([str(cfg), str(tmp_path / "b2"), "fwd", "--mpi"], tmp_path).stdout
    assert "make tap_" not in fwd and "-mpi" in fwd and "mitgcmuv" in fwd


def test_unknown_target_fails(tmp_path):
    r = _run([str(_cfg(tmp_path)), str(tmp_path / "b"), "foo"], tmp_path)
    assert r.returncode != 0 and "target" in r.stderr
