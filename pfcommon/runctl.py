"""Run-status bookkeeping written by scripts/run_mitgcm.sh."""
from __future__ import annotations

import json
import re
from pathlib import Path

_KEYS = ("state", "last_iter", "n_timesteps")
# Same pattern as scripts/run_mitgcm.sh: monitor / solver lines with NaN or Infinity (gfortran).
_NONFINITE = re.compile(r"(%MON|cg2d|cg3d).*(NaN|Infinity)")


def run_status(run_dir: Path) -> dict:
    """Read run_status.json: {"state": OK|FAILED|RUNNING, last_iter, n_timesteps, wall_s}."""
    f = Path(run_dir) / "run_status.json"
    if not f.exists():
        raise RuntimeError(f"{run_dir}: no run_status.json (run did not finish through run_mitgcm.sh)")
    try:
        st = json.loads(f.read_text())
    except json.JSONDecodeError as e:
        raise RuntimeError(f"{run_dir}: corrupt run_status.json ({e})") from e
    if not isinstance(st, dict) or any(k not in st for k in _KEYS):
        raise RuntimeError(f"{run_dir}: run_status.json lacks one of {_KEYS}")
    return st


def require_ok(run_dir: Path) -> None:
    """Raise RuntimeError unless the run finished OK and reached its last time step."""
    st = run_status(run_dir)
    if st["state"] != "OK":
        raise RuntimeError(f"{run_dir}: run state is {st['state']}, refusing to use it")
    if st["last_iter"] < st["n_timesteps"]:
        raise RuntimeError(
            f"{run_dir}: run stopped at iteration {st['last_iter']} < nTimeSteps {st['n_timesteps']}"
        )
    for name in ("STDOUT.0000", "output.txt"):
        log = Path(run_dir) / name
        if log.exists():
            with log.open(errors="replace") as fh:
                for n, line in enumerate(fh, 1):
                    if _NONFINITE.search(line):
                        raise RuntimeError(
                            f"{run_dir}: non-finite values in {name} line {n} (run diverged): {line.strip()[:120]}"
                        )
            break
