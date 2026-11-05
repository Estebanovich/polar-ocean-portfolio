"""Grid and ice-front fields for the Exp. 1 MITgcm configurations."""
from __future__ import annotations

import numpy as np


def stretched_x(dx0: float, dxmax: float, r: float, L: float) -> np.ndarray:
    """Cell widths [m] growing geometrically from the ice face.

    dx_n = min(dx0 r^n, dxmax); cells are added until the cumulative length
    reaches L, so sum(dx) >= L (overshoot < dxmax). The first cell (index 0)
    is next to the face.
    """
    if dx0 <= 0 or dxmax < dx0 or r < 1.0 or L <= 0:
        raise ValueError("need 0 < dx0 <= dxmax, r >= 1, L > 0")
    out, total, dx = [], 0.0, dx0
    while total < L - 1e-9:
        out.append(dx)
        total += dx
        dx = min(dx * r, dxmax)
    return np.array(out)


def icefront_fields(nx, ny, nz, dx, dy, front_i, H):
    """Inputs of pkg/icefront (2-D, shape (ny, nx)).

    R_icefront : depth of the submerged ice wall [m, positive] in the wet
    column ``front_i`` (the cell whose western face is the ice), 0 elsewhere
    (icefront_description.tex:12-14; icefront_init_fixed.F:73-75 reads it and
    :78-86 takes ABS). nz is only used to check that H is resolved.
    icefrontlength : "ratio of the horizontal length of vertical ice front to
    horizontal grid area" (icefront_description.tex:15-16) = dy/(dx dy) =
    1/dx [m^-1] in the front column; 0 elsewhere (it is the switch tested at
    icefront_thermodynamics.F:89). dx may be an array (per column) or scalar.
    """
    dx = np.broadcast_to(np.asarray(dx, dtype=float), (nx,))
    if not 0 <= front_i < nx:
        raise ValueError(f"front_i={front_i} outside 0..{nx - 1}")
    R = np.zeros((ny, nx))
    length = np.zeros((ny, nx))
    R[:, front_i] = H
    length[:, front_i] = 1.0 / dx[front_i]
    return R, length
