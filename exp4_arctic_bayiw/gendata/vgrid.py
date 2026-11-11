"""Vertical grid of the Van Mijenfjorden configuration (fjord_vmf).

The level thicknesses grow smoothly from ``dz_top`` at the surface as a power of the level index,
dz_k = dz_top + (dz_bot - dz_top) (k/(nz-1))^p, so that the top ~60 m (the summer pycnocline of
Van Mijenfjorden) hold about 24 nearly uniform 1.5-5 m cells, as wanted for KPP, and the cells
below grow (ratio of neighbours <= 1.12) to reach the maximum depth of the model box.
All depths in metres, positive downward.
"""
from __future__ import annotations

import numpy as np


def faces(delR) -> np.ndarray:
    """Depth of the cell interfaces (rF in MITgcm, positive down): [0, dz1, dz1+dz2, ...] (m)."""
    d = np.asarray(delR, dtype=float)
    if d.ndim != 1 or np.any(d <= 0):
        raise ValueError("delR must be a 1-D array of positive thicknesses")
    return np.concatenate([[0.0], np.cumsum(d)])


def stretched_delR(nz: int, dz_top: float, H: float, p: float) -> np.ndarray:
    """Level thicknesses dz_k = dz_top + (dz_bot - dz_top) (k/(nz-1))^p, k = 0..nz-1 (m).

    The top cell is exactly ``dz_top``; dz_bot follows from sum(dz) = H, which is linear in dz_bot.
    The thicknesses are rounded to 1 mm and the rounding residual is added to the last (thickest)
    level so that sum(dz) = H exactly.
    """
    k = np.arange(nz, dtype=float)
    s = (k / (nz - 1)) ** p
    if H <= nz * dz_top:
        raise ValueError("H must exceed nz * dz_top")
    dz_bot = dz_top + (H - nz * dz_top) / s.sum()
    dz = np.round(dz_top + (dz_bot - dz_top) * s, 3)
    dz[-1] = round(dz[-1] + (H - dz.sum()), 3)
    if np.any(np.diff(dz) < 0):
        raise ValueError("rounding broke monotonicity")
    return dz
