"""exp4_arctic_bayiw/reproduce.sh (Task 8): the quick mode recomputes every cache whose key misses
sources, both modes refuse while a model run is active, and the full mode can be dry-run without
touching anything (backup and comparison of the committed JSONs included in its plan)."""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parents[1]
SH = PKG / "reproduce.sh"


def _run(args, tmp_path, fake_pgrep: str = "", env_extra=None):
    """Run the script with RUNS in tmp_path and a fake pgrep on PATH printing ``fake_pgrep``."""
    bin_ = tmp_path / "bin"
    bin_.mkdir(exist_ok=True)
    pg = bin_ / "pgrep"
    pg.write_text(f"#!/bin/sh\nprintf '%s' '{fake_pgrep}'\n")
    pg.chmod(0o755)
    runs = tmp_path / "runs"
    runs.mkdir(exist_ok=True)
    env = {**os.environ, "PATH": f"{bin_}:{os.environ['PATH']}", "RUNS": str(runs), **(env_extra or {})}
    return subprocess.run(["bash", str(SH), *args], capture_output=True, text=True, env=env, timeout=120)


def _figures_state():
    return {p.name: p.stat().st_mtime_ns for p in (PKG / "figures").iterdir() if p.is_file()}


def test_script_parses():
    assert subprocess.run(["bash", "-n", str(SH)]).returncode == 0


def test_quick_recomputes_caches_whose_keys_miss_sources():
    s = SH.read_text()
    # partition caches: keys omit helpers outside the hashed functions (and validate.py for kf/faces)
    assert re.search(r'partition\.py"? compute .*--recompute', s)
    # thesis_contrast library_kf: CACHE_SOURCES omits thesis_contrast.py / flux.py (Task 5 parked item)
    assert re.search(r'thesis_contrast\.py"? contrast .*--recompute', s)


def test_no_bare_model_binary_name_in_the_script():
    """The process check goes through pfcommon.campaign (process names), never a pattern that could
    match the checking command itself."""
    s = SH.read_text()
    assert "other_mitgcm_processes" in s
    assert not re.search(r"(?<!\[m\])itgcmuv", s.replace("[m]itgcmuv", ""))


@pytest.mark.parametrize("args", [["--quick"], []])
def test_refuses_while_a_model_run_is_active(tmp_path, args):
    before = _figures_state()
    r = _run(args, tmp_path, fake_pgrep="4242 mpirun\n")
    assert r.returncode == 3, r.stdout + r.stderr
    assert "active" in r.stderr
    assert _figures_state() == before


def test_full_dry_run_prints_the_plan_and_writes_nothing(tmp_path):
    before = _figures_state()
    r = _run(["--dry-run"], tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    out = r.stdout
    for needle in ("backup", "run_campaign.py --config bay_thesis ", "--config fjord_vmf ",
                   "--config fjord_vmf_400 ", "column_test.py", "restart_test.py",
                   "partition.py compute", "--recompute", "render_readme.py", "compare"):
        assert needle in out, needle
    # campaigns before the analysis, the comparison last
    assert out.index("run_campaign.py") < out.index("partition.py compute") < out.rindex("compare")
    assert _figures_state() == before
    assert not any((tmp_path / "runs").iterdir())


def test_quick_dry_run_has_no_campaign(tmp_path):
    r = _run(["--quick", "--dry-run"], tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "run_campaign.py" not in r.stdout and "partition.py compute" in r.stdout
    assert "backup" not in r.stdout


def test_unknown_option(tmp_path):
    assert _run(["--bogus"], tmp_path).returncode == 2
