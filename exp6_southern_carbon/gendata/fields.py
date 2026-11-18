"""Physical input fields of the tutorial channel.

g50 (50 km): exact port of verification/tutorial_reentrant_channel/input/gendata.50km.py:27-75
(bit-identical >f4 files, tests/test_fields.py). f10 (10 km): 2x2 block mean of the exact port of
gendata.5km.py:27-79, so the 10-km channel is the tutorial's eddy-permitting set-up coarsened;
the tutorial's index-based ridge and sponge have no general-resolution formula.
Seasonal SST relaxation and the wind factor are Exp. 6 additions (spec Sec. 3.1)."""
from __future__ import annotations

import numpy as np
from numpy import exp, pi, sin

from .grid import Grid

TMAX, TMIN, H_T, TAUX_MAX, BUMP_MAX = 10.0, -2.0, 500.0, 0.2, 2000.0


def _temp3d(sst, z, H):
    zscale = ((exp(-z / H_T) - exp(H / H_T)) / (1 - exp(H / H_T))).reshape(-1, 1, 1)
    return (sst - TMIN)[None] * zscale + TMIN


def coarse50_fields() -> dict[str, np.ndarray]:
    g = Grid(20, 40, 50e3)
    nx, ny, nr = g.nx, g.ny, g.nr
    bathy = g.H * np.ones([ny, nx])
    bump = np.zeros([ny, nx])
    r1 = BUMP_MAX * sin(np.linspace(0, pi, 9))
    r2 = np.reshape(np.linspace(0, 1, 6), [6, 1])
    bump[:, 5:14] = r1
    bump[13:19, :] = bump[13:19, :] * np.flip(r2[:])
    bump[20:26, :] = bump[20:26, :] * r2[:]
    bump[18:21, :] = 0.0
    bathy = bathy + bump
    bathy[0, :] = 0.0
    taux = np.zeros([ny, nx])
    taux[:, :] = np.reshape(TAUX_MAX * sin(np.linspace(0, pi, ny)), [ny, 1])
    mask = np.zeros([nr, ny, nx])
    mask[1:, -1, :] = 1.0
    mask[1:, -2, :] = 0.25
    sst = np.zeros([ny, nx]) + np.reshape(TMIN + (TMAX - TMIN) / 40 * np.linspace(0.5, 39.5, ny), [ny, 1])
    return {"bathy": bathy, "taux": taux, "mask": mask, "sst": sst, "temp": _temp3d(sst, g.z, g.H)}


def fine5_fields() -> dict[str, np.ndarray]:
    g = Grid(200, 400, 5e3)
    nx, ny, nr = g.nx, g.ny, g.nr
    bathy = g.H * np.ones([ny, nx])
    bump = np.zeros([ny, nx])
    r1 = BUMP_MAX * sin(np.linspace(0, pi, 80))
    r2 = np.reshape(np.linspace(0, 1, 52), [52, 1])
    bump[:, 55:135] = r1
    bump[134:186, :] = bump[134:186, :] * np.flip(r2[:])
    bump[204:256, :] = bump[204:256, :] * r2[:]
    bump[185:205, :] = 0.0
    bathy = bathy + bump
    bathy[0:10, :] = 0.0
    taux = np.zeros([ny, nx])
    taux[5:395, :] = np.reshape(TAUX_MAX * sin(np.linspace(0, pi, 390)), [390, 1])
    mask = np.zeros([nr, ny, nx])
    mask[1:, -10:, :] = 1.0
    mask[1:, -20:-10, :] = 0.25
    sst = np.zeros([ny, nx]) + np.reshape(TMIN + (TMAX - TMIN) / 40 * np.linspace(0.05, 39.95, ny), [ny, 1])
    return {"bathy": bathy, "taux": taux, "mask": mask, "sst": sst, "temp": _temp3d(sst, g.z, g.H)}


def _block2(a: np.ndarray) -> np.ndarray:
    *lead, ny, nx = a.shape
    return a.reshape(*lead, ny // 2, 2, nx // 2, 2).mean(axis=(-3, -1))


def tutorial_fields(g: Grid) -> dict[str, np.ndarray]:
    if (g.nx, g.ny, g.dx_m) == (20, 40, 50e3):
        return coarse50_fields()
    if (g.nx, g.ny, g.dx_m) == (100, 200, 10e3):
        return {k: _block2(v) for k, v in fine5_fields().items()}
    raise ValueError(f"no tutorial fields for {g.nx}x{g.ny} at {g.dx_m} m")


def bathymetry(g: Grid) -> np.ndarray:
    return tutorial_fields(g)["bathy"]


def taux(g: Grid, factor: float = 1.0) -> np.ndarray:
    return factor * tutorial_fields(g)["taux"]


def sst_relax(g: Grid) -> np.ndarray:
    return tutorial_fields(g)["sst"]


def t_relax_mask(g: Grid) -> np.ndarray:
    return tutorial_fields(g)["mask"]


def temperature(g: Grid) -> np.ndarray:
    return tutorial_fields(g)["temp"]


def sst_relax_monthly(g: Grid, amp: float = 1.5, month_min: int = 8) -> np.ndarray:
    """(12, ny, nx): SST_relax + amp cos(2 pi (m - month_min + 6)/12), m = 1..12; minimum in month_min."""
    m = np.arange(1, 13)
    s = amp * np.cos(2 * np.pi * (m - month_min + 6) / 12.0)
    return sst_relax(g)[None] + s[:, None, None]
