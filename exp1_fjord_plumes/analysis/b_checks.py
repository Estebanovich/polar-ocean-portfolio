"""Validity checks for the config (b) campaign (Task 6; PD1 note §6).

Eta retention (volume budget of the discharge)
----------------------------------------------
addMass adds the discharge volume Q_total [m3 s-1] inside the box; the
meltwater enters as a virtual salt flux without volume (pkg/icefront), so with
a linear free surface the wet-area mean surface elevation obeys

    A dEta/dt = Q_total - Q_OB                       (volume budget; A = wet area)

If nothing left through the open boundary, Eta would follow the
"trapped-volume" line Eta_trapped(t) = Q_total t / A [m]. Two fractions are
reported (definitions made here, not from a reference):

    retained_cum       = Eta(t_end) / Eta_trapped(t_end)
    retained_last_hour = [Eta(t_end) - Eta(t_end - 1 h)] / (Q_total * 1 h / A)

``retained_cum`` includes the spin-up, while the outflow has not yet reached
the open boundary, so it is informative only. The run is flagged when
``retained_last_hour > RETENTION_MAX``: in the last hour the open boundary
passes less than half of the discharge (justification in
notes/pd1_b_twoplume3d.md §6).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import read_bin, read_mds
from pfcommon.runctl import require_ok

RETENTION_MAX = 0.5  # flag threshold on retained_last_hour (PD1 note §6)


def _eta_series(run_dir: Path, dt: float) -> dict[float, np.ndarray]:
    """{t [s]: Eta (ny, nx)} from ETAN snapshots and the Eta state dumps."""
    out: dict[float, np.ndarray] = {}
    for name in ("ETAN", "Eta"):
        for f in run_dir.glob(f"{name}.[0-9]*.meta"):
            m = re.fullmatch(rf"{name}\.(\d+)\.meta", f.name)
            if m:  # iteration 0 (Eta = 0) is the start of the budget
                it = int(m.group(1))
                out.setdefault(it * dt, read_mds(run_dir / name, it)[0][0])
    return dict(sorted(out.items()))


def eta_retention(run_dir: Path) -> dict:
    """Retained fraction of the discharge volume (module docstring). Refuses runs not OK."""
    run_dir = Path(run_dir)
    require_ok(run_dir)
    meta = json.loads((run_dir / "case.json").read_text())
    p = meta["params"]
    ny, nx, dt = p["ny"], p["nx"], p["deltaT"]
    Q = float(meta["Q0_total"])
    wet = read_bin(run_dir / "topog.bin", (ny, nx)) < 0
    cell = p["dx"] * p["dy"]  # uniform grid in (b)
    A = wet.sum() * cell  # m2
    series = _eta_series(run_dir, dt)
    if not series:
        raise FileNotFoundError(f"no ETAN/Eta output in {run_dir}")
    mean = {t: float(e[wet].mean()) for t, e in series.items()}
    t_end = max(mean)
    t_prev = min(mean, key=lambda t: abs(t - (t_end - 3600.0)))
    trapped = Q * t_end / A
    last = ((mean[t_end] - mean[t_prev]) / (Q * (t_end - t_prev) / A)
            if t_end > t_prev else float("nan"))
    return {"t_end_h": t_end / 3600.0, "t_prev_h": t_prev / 3600.0, "wet_area_m2": A,
            "eta_mean_m": mean[t_end], "eta_trapped_m": trapped,
            "retained_cum": mean[t_end] / trapped, "retained_last_hour": last,
            "threshold": RETENTION_MAX, "flag": bool(last > RETENTION_MAX)}
