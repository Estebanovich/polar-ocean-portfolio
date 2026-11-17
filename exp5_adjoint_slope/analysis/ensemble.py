"""Post hoc ensemble criteria (PD2 amendment, notes/pd0_preregistro.md; written before EA and the
ensemble forwards were seen). The coarse model is chaotic, so the adjoint is averaged over the
EA21-EA30 windows (Lea et al. 2000) and checked against paired forward ensembles.

dp, dm are per-member responses J(+du) - J0 and J(-du) - J0 (TW), J0 the member's own trajectory.
"""
from __future__ import annotations

import numpy as np

H5C_RANGE = (30.0, 180.0)        # days, frozen H5c range (PD0)


def mean_se(x) -> tuple[float, float]:
    """Member mean and its standard error (ddof = 1)."""
    x = np.asarray(x, dtype=float)
    return float(x.mean()), float(x.std(ddof=1) / np.sqrt(x.size))


def v2e(dp, dm, tol: float = 0.1) -> dict:
    """Linearity on member means: |<d+ + d->| / |<d+ - d->| < tol."""
    dp, dm = np.asarray(dp, float), np.asarray(dm, float)
    val = abs((dp + dm).mean()) / abs((dp - dm).mean())
    return {"pass": bool(val < tol), "value": float(val), "tol": tol}


def v3e(pred: float, dp, dm, tol: float = 0.30) -> dict:
    """Ensemble-adjoint prediction against the forward response <(d+ - d-)/2>; passes only if
    within tol and the response exceeds two standard errors."""
    half = (np.asarray(dp, float) - np.asarray(dm, float)) / 2.0
    fwd, se = mean_se(half)
    rel = abs(pred - fwd) / abs(fwd)
    det = bool(abs(fwd) > 2.0 * se)
    return {"pass": bool(rel < tol and det), "value": float(rel), "tol": tol, "pred": float(pred),
            "fwd": fwd, "se": se, "detectable": det}


def t_half(leads, s) -> float | None:
    """Lead (days) at which |s| first falls to half its maximum after the maximum, linearly
    interpolated between records; None if it never does within the leads given."""
    leads, a = np.asarray(leads, float), np.abs(np.asarray(s, float))
    k = int(a.argmax())
    h = 0.5 * a[k]
    for n in range(k + 1, a.size):
        if a[n] <= h:
            return float(leads[n - 1] + (leads[n] - leads[n - 1]) * (a[n - 1] - h) / (a[n - 1] - a[n]))
    return None


def h5c(leads, s) -> dict:
    """H5c on the ensemble-mean sensitivity: half-decay between 30 and 180 days."""
    t = t_half(leads, s)
    if t is None:
        v = f"not testable beyond {max(leads):.0f} days"
    else:
        v = "pass" if H5C_RANGE[0] <= t <= H5C_RANGE[1] else "fail"
    return {"verdict": v, "t_half_days": t, "leads_days": [float(x) for x in leads]}
