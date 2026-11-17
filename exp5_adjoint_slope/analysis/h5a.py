"""H5a (notes/pd0_preregistro.md, PD0 + PD1b scales + PD2): curl against magnitude.

Two time-constant, x-uniform tau_x(y) patterns with the same RMS over the shelf and slope
(y < 275 km): (i) a uniform increase of the easterlies, (ii) a dipole at the shelf break
(y0 = 125 km), the derivative of a Gaussian of half-width w (37.5 km; robustness 18.75 and
75 km). dJ_pred = sum_y S(y) dtau(y), with S the ensemble-adjoint sensitivity summed over x and
over the monthly records (TW per N m-2). H5a holds if |dJ(ii)| > |dJ(i)| (frozen); "robust" if
it holds at the three widths. PD2 adds, without changing the verdict, whether the paired member
difference |dJ(ii)| - |dJ(i)| exceeds two standard errors."""
from __future__ import annotations

import numpy as np

Y_MAX = 275e3            # m, RMS region and extent of (i) (PD1b)
Y0 = 125e3               # m, shelf break (PD1b)
WIDTHS = (37.5e3, 18.75e3, 75e3)
RMS = 0.01               # N m-2, common RMS over y < Y_MAX


def _norm(p, yc):
    m = np.asarray(yc) < Y_MAX
    return RMS * p / np.sqrt((p[m] ** 2).mean())


def uniform(yc) -> np.ndarray:
    """(i): stronger easterlies (tau_x more negative) over y < 275 km."""
    yc = np.asarray(yc, float)
    return _norm(np.where(yc < Y_MAX, -1.0, 0.0), yc)


def dipole(yc, w: float) -> np.ndarray:
    """(ii): x (y - y0)/w exp(-((y - y0)/w)^2), derivative of a Gaussian centred at the break."""
    yc = np.asarray(yc, float)
    x = (yc - Y0) / w
    return _norm(x * np.exp(-x ** 2), yc)


def dj(sens_y, pattern) -> float:
    return float(np.dot(np.asarray(sens_y, float), np.asarray(pattern, float)))


def evaluate(sens_members, yc) -> dict:
    """sens_members: (n_members, ny), each summed over x and records."""
    s = np.atleast_2d(np.asarray(sens_members, float))
    u = uniform(yc)
    out = {"dj_uniform_mean": float(np.mean([dj(m, u) for m in s])), "widths": {}}
    for w in WIDTHS:
        d = dipole(yc, w)
        ji, jii = np.array([dj(m, u) for m in s]), np.array([dj(m, d) for m in s])
        diff = np.abs(jii) - np.abs(ji)
        se = float(diff.std(ddof=1) / np.sqrt(diff.size)) if diff.size > 1 else float("nan")
        mi, mii = float(ji.mean()), float(jii.mean())
        out["widths"][f"{w / 1e3:g}"] = {
            "dj_uniform": mi, "dj_dipole": mii, "pass": bool(abs(mii) > abs(mi)),
            "paired_diff": float(diff.mean()), "paired_se": se,
            "significant": bool(diff.mean() > 2 * se)}
    out["robust"] = all(v["pass"] for v in out["widths"].values())
    return out
