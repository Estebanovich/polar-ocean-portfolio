"""V1-V3 (notes/pd0_preregistro.md, frozen): each returns {"pass", "value", ...} for the JSON."""
from __future__ import annotations

N_POINTS_MIN = 10


def v1(points: list[dict], tol: float = 0.01) -> dict:
    """grdchk: relative error < tol at every point, and at least 10 points."""
    errs = [float(p["rel_err"]) for p in points]
    worst = max(errs) if errs else float("nan")
    ok = len(errs) >= N_POINTS_MIN and worst < tol
    return {"pass": bool(ok), "value": worst, "n_points": len(errs), "tol": tol}


def v2(dj_plus: float, dj_minus: float, tol: float = 0.1) -> dict:
    """Linearity: |dJ+ + dJ-| / |dJ+ - dJ-| < tol."""
    val = abs(dj_plus + dj_minus) / abs(dj_plus - dj_minus)
    return {"pass": bool(val < tol), "value": val, "dj_plus": dj_plus, "dj_minus": dj_minus,
            "tol": tol}


def v3(dj_pred: float, dj_fwd: float, tol: float = 0.30) -> dict:
    """Adjoint-predicted vs forward dJ: |pred - fwd| / |fwd| < tol."""
    val = abs(dj_pred - dj_fwd) / abs(dj_fwd)
    return {"pass": bool(val < tol), "value": val, "dj_pred": dj_pred, "dj_fwd": dj_fwd,
            "tol": tol}


def top_points(sens, n: int, min_sep: int = 2, floor: float = 1e-10) -> list[tuple[int, ...]]:
    """Indices of the n largest |sens| (V1 points, PD0: "where the GA sensitivity is large"),
    at least min_sep cells apart in every horizontal direction from the ones already chosen, and
    above floor * max|sens| (plan Review Focus 4: a point without signal is replaced)."""
    import numpy as np
    a = np.abs(np.asarray(sens, dtype=float))
    thr = floor * a.max()
    chosen: list[tuple[int, ...]] = []
    for flat in np.argsort(a, axis=None)[::-1]:
        if len(chosen) == n:
            break
        idx = tuple(int(v) for v in np.unravel_index(flat, a.shape))
        if a[idx] <= thr:
            break
        if all(max(abs(idx[-1] - c[-1]), abs(idx[-2] - c[-2])) >= min_sep for c in chosen):
            chosen.append(idx)
    return chosen
