"""Run-status bookkeeping written by scripts/run_mitgcm.sh."""
from __future__ import annotations

import json
from pathlib import Path


def run_status(run_dir: Path) -> dict:
    """Read run_status.json: {"state": OK|FAILED|RUNNING, last_iter, n_timesteps, wall_s}."""
    f = Path(run_dir) / "run_status.json"
    if not f.exists():
        raise RuntimeError(f"{run_dir}: no run_status.json (run did not finish through run_mitgcm.sh)")
    return json.loads(f.read_text())


def require_ok(run_dir: Path) -> None:
    """Raise RuntimeError unless the run finished with state OK."""
    state = run_status(run_dir)["state"]
    if state != "OK":
        raise RuntimeError(f"{run_dir}: run state is {state}, refusing to use it")
