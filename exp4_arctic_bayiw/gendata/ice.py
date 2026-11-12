"""Initial sea-ice fields of the Exp. 4 fjord runs (Task 4; spec section 3.5).

Ice of thickness h = 1 m and concentration A = 1 (Hoyland 2009 [verificar]) covers the wet
cells whose centre lies inside the inner-fjord polygon (east of Akseloya; the profile polygon
of cases.yaml, ``configs.fjord_vmf.profiles.polygon_utm``). Outside it there is open water.
pkg/seaice reads them as ``HeffFile`` (effective thickness, m) and ``AreaFile`` (fraction)
(pkg/seaice/seaice_init_varia.F:291-292 and below; defaults in pkg/seaice/seaice_readparms.F:512-515).
Thermodynamics is off, so h and A only change by advection (free drift, F4).
"""
from __future__ import annotations

import numpy as np
import shapely


def mask(xc, yc, wet, polygon_utm) -> np.ndarray:
    """(ny, nx) bool: wet cells whose centre (xc[i], yc[j], EPSG:32633 m) is strictly inside
    ``polygon_utm`` (list of [x, y] vertices)."""
    XC, YC = np.meshgrid(np.asarray(xc, float), np.asarray(yc, float))
    inside = shapely.contains_xy(shapely.Polygon(polygon_utm), XC, YC)
    return np.asarray(inside, bool) & np.asarray(wet, bool)


def fields(mask, h: float = 1.0, A: float = 1.0) -> dict:
    """{'HEFF.bin': h on the mask, 0 elsewhere [m]; 'AREA.bin': A on the mask [0-1]}."""
    if not 0.0 <= A <= 1.0:
        raise ValueError(f"ice concentration A = {A} outside [0, 1]")
    if h < 0:
        raise ValueError(f"ice thickness h = {h} < 0")
    m = np.asarray(mask, bool)
    return {"HEFF.bin": np.where(m, float(h), 0.0), "AREA.bin": np.where(m, float(A), 0.0)}
