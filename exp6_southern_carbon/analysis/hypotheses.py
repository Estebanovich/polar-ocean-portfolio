"""H6a and H6b (spec Sec. 7; frozen in notes/pd0_preregistro.md).

H6a (Esteban: signed, sign and detection only): dS = S_GT(2060) - S_GT(1900) > 0 (the wind
response moves toward uptake) and dS > 2 |S_GC(2060) - S_GC(1900)| (control drift, V6c).
S_nat = S_GC, S_ant = S_GT - S_GC (Lovenduski et al. 2008).
H6b: |pred - meas| <= 0.25 |meas| in every GT branch."""
from __future__ import annotations

import numpy as np


def h6a(s_gt: dict[int, float], s_gc: dict[int, float]) -> dict:
    ds = s_gt[2060] - s_gt[1900]
    drift = abs(s_gc[2060] - s_gc[1900])
    years = sorted(s_gt)
    s0 = np.sign(s_gt[1900])
    cross = next((y for y in years if y > 1900 and np.sign(s_gt[y]) != s0), None)
    return {"dS": float(ds), "drift_gc": float(drift), "pass": bool(ds > 0 and ds > 2 * drift),
            "s_nat": {y: float(s_gc[y]) for y in sorted(s_gc)},
            "s_ant": {y: float(s_gt[y] - s_gc[y]) for y in years if y in s_gc},
            "zero_crossing_year": cross}


def h6b(pred: dict[int, float], meas: dict[int, float], tol: float = 0.25) -> dict:
    within = {y: bool(abs(pred[y] - meas[y]) <= tol * abs(meas[y])) for y in sorted(meas)}
    return {"within": within, "n_within": int(sum(within.values())), "pass": bool(all(within.values())),
            "ratio": {y: float(pred[y] / meas[y]) if meas[y] else float("nan") for y in sorted(meas)}}
