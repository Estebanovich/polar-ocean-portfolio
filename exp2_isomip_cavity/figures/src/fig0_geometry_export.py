"""Figure 0d, step 1: export the ISOMIP+ geometry for the Blender render (fig0_geometry_blender.py).

Reads the official ISOMIP+ geometry (Cornford and Asay-Davis, 2016; data_external/Ocean1_input_geom_v1.01.nc,
1 km grid, dims (y, x), x = 320-800 km, y = 0-80 km) and writes, into WORK (default: a fig0_geometry folder in
the system temporary directory, not the repository):

  geom.npz       float32 arrays (y, x) in metres: bed (bedrockTopography), lower and upper ice surfaces
                 (lowerSurface, upperSurface), and the masks floating, grounded, open (0/1), plus x, y.

and, into this directory (committed, a few kB, used by fig0_schematic.py for panel a):

  isomip_outline.json   the grounding line x(y) and the calving-front x of the geometry file, in km.

Only the input geometry is read: no model output enters Figure 0.
Usage: fig0_geometry_export.py [--nc PATH] [--work DIR]
"""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
NC = EXP / "data_external" / "Ocean1_input_geom_v1.01.nc"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nc", default=str(NC))
    ap.add_argument("--work", default=str(Path(tempfile.gettempdir()) / "fig0_geometry"))
    a = ap.parse_args()
    import netCDF4
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(a.nc) as d:
        g = {k: np.asarray(d[v][:], dtype=np.float32) for k, v in (
            ("bed", "bedrockTopography"), ("lower", "lowerSurface"), ("upper", "upperSurface"),
            ("floating", "floatingMask"), ("grounded", "groundedMask"), ("open", "openOceanMask"))}
        g["x"] = np.asarray(d["x"][:], dtype=np.float32)
        g["y"] = np.asarray(d["y"][:], dtype=np.float32)
    np.savez_compressed(work / "geom.npz", **g)
    # grounding line: the first floating cell along x in each row; front: the last floating cell (1 km cells)
    x_km = g["x"] / 1e3
    gl = [float(x_km[np.argmax(row > 0)] - 0.5) for row in g["floating"]]
    front = float(max(x_km[np.nonzero(row > 0)[0][-1]] for row in g["floating"]) + 0.5)
    out = {"source": "Ocean1_input_geom_v1.01.nc (Cornford and Asay-Davis, 2016)",
           "y_km": [float(v) for v in g["y"] / 1e3], "grounding_line_x_km": gl, "front_x_km": front,
           "x_range_km": [320.0, 800.0], "y_range_km": [0.0, 80.0]}
    (HERE / "isomip_outline.json").write_text(json.dumps(out, indent=0))
    print("wrote", work / "geom.npz", "and", HERE / "isomip_outline.json")


if __name__ == "__main__":
    main()
