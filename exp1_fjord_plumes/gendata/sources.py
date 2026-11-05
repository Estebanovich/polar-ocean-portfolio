"""Subglacial-discharge sources (pkg/addmass) and their placement."""
from __future__ import annotations

import numpy as np


def discharge_addmass(shape, cells, Q0_total, rho0, *, hFacC, front_i,
                      face_offset=0, require_bottom=True):
    """Mass source array [kg s^-1 per cell], shape (nz, ny, nx).

    cells : list of (k, j, i); duplicates accumulate. The total is
    rho0 * Q0_total (Q0_total in m^3 s^-1; negative gives a sink). Each cell is
    given an equal share. Raises ValueError, naming the offending cell, if a
    cell (Review Focus 1)
      * is land (hFacC = 0) or outside the grid,
      * is not adjacent to the ice face: i != front_i + face_offset (the front
        column holds R_icefront, i.e. the wet cell touching the ice, so the
        default offset is 0),
      * is above the bottom wet cell of its column (if require_bottom).
    """
    nz, ny, nx = shape
    hFacC = np.asarray(hFacC)
    if hFacC.shape != tuple(shape):
        raise ValueError(f"hFacC shape {hFacC.shape} != {tuple(shape)}")
    if not cells:
        raise ValueError("no source cells given")
    out = np.zeros(shape)
    share = rho0 * Q0_total / len(cells)
    for (k, j, i) in cells:
        c = (k, j, i)
        if not (0 <= k < nz and 0 <= j < ny and 0 <= i < nx):
            raise ValueError(f"source cell (k,j,i)={c} is outside the grid {tuple(shape)}")
        if hFacC[k, j, i] <= 0:
            raise ValueError(f"source cell (k,j,i)={c} is land (hFacC=0)")
        if i != front_i + face_offset:
            raise ValueError(f"source cell (k,j,i)={c} is not adjacent to the ice "
                             f"face (needs i={front_i + face_offset})")
        if require_bottom:
            wet = np.nonzero(hFacC[:, j, i] > 0)[0]
            if k != wet.max():
                raise ValueError(f"source cell (k,j,i)={c} is above the bottom "
                                 f"layer of its column (k={wet.max()})")
        out[k, j, i] += share
    return out


def place_sources(y_centers, x0_nominal, dy):
    """Snap two symmetric sources to the grid.

    Separation = n dy with n = floor(x0_nominal/dy + 0.5) (round half up),
    sources at j1 = jmid - n//2 and j2 = j1 + n, jmid = ny//2. Returns
    ([j1, j2], x0_actual) where x0_actual = y_centers[j2] - y_centers[j1] [m],
    the separation the analysis must use (Review Focus 2). x0 = 0 gives
    j1 = j2 (both discharges in the same column).
    """
    y = np.asarray(y_centers, dtype=float)
    n = int(np.floor(x0_nominal / dy + 0.5))
    j1 = len(y) // 2 - n // 2
    j2 = j1 + n
    if j1 < 0 or j2 >= len(y):
        raise ValueError(f"x0={x0_nominal} m does not fit in {len(y)} columns")
    return [j1, j2], float(y[j2] - y[j1])
