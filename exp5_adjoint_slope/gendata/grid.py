"""Model grid of the Exp. 5 channel (spec Sec. 3.1): zonally reentrant, walls at y = 0 and y = Ly.

Conventions follow MITgcm: tracer cells are centred at yc; the v point of row j sits on the
southern face of cell j, at yg[j] = j*dy (pkg-independent, model/src/ini_cartesian_grid.F).
Depths z are positive downward (m)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def vertical_spacing(nr: int = 35, top_m: float = 10.0, total_m: float = 3000.0) -> np.ndarray:
    """Layer thicknesses (m), growing linearly from `top_m` so that they sum to `total_m`.

    With nr = 35, 10 m at the top and 3000 m total, the bottom layer is 2*3000/35 - 10 = 161 m.
    """
    bottom = 2.0 * total_m / nr - top_m
    if bottom <= top_m:
        raise ValueError("total_m too small for a downward-growing spacing")
    return np.linspace(top_m, bottom, nr)


@dataclass
class Grid:
    nx: int
    ny: int
    dx_m: float
    dz_m: np.ndarray
    yc: np.ndarray = field(init=False)
    yg: np.ndarray = field(init=False)
    xc: np.ndarray = field(init=False)
    zf: np.ndarray = field(init=False)
    zc: np.ndarray = field(init=False)

    def __post_init__(self):
        self.dz_m = np.asarray(self.dz_m, dtype=float)
        j = np.arange(self.ny)
        self.yg = j * self.dx_m
        self.yc = (j + 0.5) * self.dx_m
        self.xc = (np.arange(self.nx) + 0.5) * self.dx_m
        self.zf = np.concatenate([[0.0], np.cumsum(self.dz_m)])
        self.zc = 0.5 * (self.zf[:-1] + self.zf[1:])

    @property
    def nr(self) -> int:
        return self.dz_m.size
