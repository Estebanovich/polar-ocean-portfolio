"""Grid of the MITgcm tutorial_reentrant_channel (Abernathey et al. 2011 channel): Cartesian,
Lx = 1000 km, Ly = 2000 km, 49 levels from the hyperbolic-tangent design of Stewart et al. (2017)
(verification/tutorial_reentrant_channel/input/gendata.50km.py:13-25)."""
from __future__ import annotations

import numpy as np

DR = np.array([5.48716549, 6.19462098, 6.99291201, 7.89353689,
               8.90937723, 10.05483267, 11.34595414, 12.80056778,
               14.43837763, 16.28102917, 18.35210877, 20.67704362,
               23.28285446, 26.1976981, 29.45012046, 33.06792588,
               37.07656002, 41.496912, 46.34247864, 51.61592052,
               57.30518684, 63.37960847, 69.78661289, 76.44996107,
               83.27047568, 90.13003112, 96.89898027, 103.44631852,
               109.65099217, 115.4122275, 120.65692923, 125.34295968,
               129.45821977, 133.01641219, 136.05088105, 138.60793752,
               140.74074276, 142.50436556, 143.95220912, 145.133724,
               146.09317287, 146.86917206, 147.49475454, 147.99774783,
               148.40131516, 148.72455653, 148.98310489, 149.18968055,
               149.35458582])


class Grid:
    def __init__(self, nx: int, ny: int, dx_m: float):
        self.nx, self.ny, self.dx_m = int(nx), int(ny), float(dx_m)
        self.dr = DR.copy()
        self.nr = DR.size
        self.rf = np.insert(np.cumsum(DR), 0, 0.0)      # face depths, m (positive down)
        self.z = 0.5 * (self.rf[1:] + self.rf[:-1])      # centre depths, m
        self.H = -float(DR.sum())                        # bottom, m (negative, as in the tutorial)
        self.xc = (np.arange(self.nx) + 0.5) * self.dx_m
        self.yc = (np.arange(self.ny) + 0.5) * self.dx_m
