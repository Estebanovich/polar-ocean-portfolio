"""Test data helpers for Exp. 2: synthetic 1 km ISOMIP+ geometry files."""
from __future__ import annotations

from pathlib import Path

import numpy as np

REAL_GEOM = Path(__file__).resolve().parents[1] / "data_external" / "Ocean1_input_geom_v1.01.nc"


def analytic_bed(X, Y):
    """Independent copy of AD16 Eq. 1-3 (Table 1 values) for the synthetic files."""
    xt = X / 300e3
    bx = -150.0 - 728.8 * xt**2 + 343.91 * xt**4 - 50.57 * xt**6
    by = (500.0 / (1 + np.exp(-2 * (Y - 40e3 - 24e3) / 4e3))
          + 500.0 / (1 + np.exp(2 * (Y - 40e3 + 24e3) / 4e3)))
    return np.maximum(bx + by, -720.0)


def synthetic_fields():
    """1 km fields on the file's grid: an ice shelf thinning to the front at 640 km,
    grounded where the draft would lie below the bed, and a thin (80 m) strip at
    630-640 km that must calve (AD16 sect. 3.1.2)."""
    x = 320.5e3 + 1e3 * np.arange(480)
    y = 0.5e3 + 1e3 * np.arange(80)
    X, Y = np.meshgrid(x, y)
    bed = analytic_bed(X, Y)
    shelf = -650.0 + 500.0 * (X - 320e3) / 320e3          # -650 m at 320 km, -150 m at 640 km
    lower = np.where(X < 640e3, np.maximum(shelf, bed), 0.0)
    thick = np.where(X < 640e3, -lower * 1028.0 / 918.0, 0.0)
    thick = np.where((X >= 630e3) & (X < 640e3), 80.0, thick)
    upper = lower + thick
    grounded = (X < 640e3) & (lower <= bed)
    floating = (X < 640e3) & ~grounded
    return dict(x=x, y=y, bedrockTopography=bed, lowerSurface=lower, upperSurface=upper,
                groundedMask=grounded.astype(float), floatingMask=floating.astype(float),
                openOceanMask=(X >= 640e3).astype(float))


def write_geom(path: Path, f: dict, dims=("y", "x"), xscale=1.0) -> Path:
    import netCDF4
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("y", len(f["y"]))
        nc.createDimension("x", len(f["x"]))
        nc.createVariable("x", "f8", ("x",))[:] = f["x"] * xscale
        nc.createVariable("y", "f8", ("y",))[:] = f["y"] * xscale
        for k in ("bedrockTopography", "lowerSurface", "upperSurface",
                  "groundedMask", "floatingMask", "openOceanMask"):
            arr = f[k] if dims == ("y", "x") else f[k].T
            nc.createVariable(k, "f8", dims)[:] = arr
    return path


