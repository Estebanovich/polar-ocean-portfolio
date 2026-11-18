"""Readers for Exp. 6 runs (MDS diagnostics, pickups, grid)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import read_bin, read_mds


def case(run: Path) -> dict:
    return json.loads((Path(run) / "input" / "case.json").read_text()) if (Path(run) / "input" / "case.json").exists() \
        else json.loads((Path(run) / "case.json").read_text())


def grid(run: Path) -> dict:
    c = case(run)
    nr, ny, nx = c["nr"], c["ny"], c["nx"]
    r = Path(run)
    return {"hfac": read_bin(r / "hFacC.data", (nr, ny, nx)), "rac": read_bin(r / "RAC.data", (ny, nx)),
            "drf": read_bin(r / "DRF.data", (nr,)), "rc": -read_bin(r / "RC.data", (nr,)),
            "yc": read_bin(r / "YC.data", (ny, nx))[:, 0], "nr": nr, "ny": ny, "nx": nx}


def iters(run: Path, prefix: str) -> list[int]:
    return sorted(int(f.name.split(".")[1]) for f in Path(run).glob(f"{prefix}.*.meta") if not f.name.startswith("._"))


def fields(run: Path, prefix: str, itr: int) -> dict[str, np.ndarray]:
    a, meta = read_mds(Path(run) / prefix, itr)
    return {name: a[i] for i, name in enumerate(meta["fldList"])}


def year_iter(run: Path, year: int) -> int:
    """Iteration at the end of model year `year` (1-based) of the run."""
    c = case(run)
    return int(c["nIter0"]) + int(year) * int(c["steps_per_year"])
