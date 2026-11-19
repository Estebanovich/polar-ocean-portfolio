"""reproduce.sh (plan Task 10): --check lists every step in order and runs nothing."""
from __future__ import annotations

import subprocess
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
SH = PKG / "reproduce.sh"

ORDER = ["build", "spin-up G0", "spin-up G0s", "coarse campaign", "fine campaign", "evaluation",
         "figures", "README"]


def test_script_parses():
    assert subprocess.run(["bash", "-n", str(SH)]).returncode == 0


def test_check_lists_steps_in_order_and_runs_nothing(tmp_path):
    fig = PKG / "figures"
    before = {p.name: p.stat().st_mtime_ns for p in fig.iterdir() if p.is_file()} if fig.exists() else {}
    r = subprocess.run(["bash", str(SH), "--check"], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr
    steps = [ln[4:].split(":")[0] for ln in r.stdout.splitlines() if ln.startswith("=== ")]
    assert steps == ORDER
    assert all(ln.startswith(("=== ", "+ ", "  ")) for ln in r.stdout.splitlines()), r.stdout
    after = {p.name: p.stat().st_mtime_ns for p in fig.iterdir() if p.is_file()} if fig.exists() else {}
    assert before == after
