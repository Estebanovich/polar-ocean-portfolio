"""Readers for Exp. 6 runs (MDS diagnostics, pickups, grid)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import read_bin, read_mds


def case(run: Path) -> dict:
    return json.loads((Path(run) / "input" / "case.json").read_text()) if (Path(run) / "input" / "case.json").exists() \
        else json.loads((Path(run) / "case.json").read_text())


def grid_dir(run: Path) -> Path:
    """Directory holding the grid files of `run`: the run itself, or, for a run restarted from a
    pickup (MITgcm writes hFacC/RAC/... only on a cold start), the nearest ancestor along
    pickup_from (sibling directories under the same runs root)."""
    r, seen = Path(run), set()
    while not (r / "hFacC.data").exists():
        parent = case(r).get("pickup_from")
        if not parent or parent in seen:
            raise FileNotFoundError(f"no grid files in {Path(run)} or its pickup_from ancestors")
        seen.add(parent)
        r = r.parent / parent
    return r


def grid(run: Path) -> dict:
    c = case(run)
    nr, ny, nx = c["nr"], c["ny"], c["nx"]
    r = grid_dir(run)
    return {"hfac": read_bin(r / "hFacC.data", (nr, ny, nx)), "rac": read_bin(r / "RAC.data", (ny, nx)),
            "drf": read_bin(r / "DRF.data", (nr,)), "rc": -read_bin(r / "RC.data", (nr,)),
            "yc": read_bin(r / "YC.data", (ny, nx))[:, 0], "dxg": read_bin(r / "DXG.data", (ny, nx)), "nr": nr, "ny": ny, "nx": nx}


def iters(run: Path, prefix: str) -> list[int]:
    return sorted(int(f.name.split(".")[1]) for f in Path(run).glob(f"{prefix}.*.meta") if not f.name.startswith("._"))


_2D = ("EtaN", "dEtaHdt", "EtaH", "DIC_pH2d")


def fields(run: Path, prefix: str, itr: int) -> dict[str, np.ndarray]:
    """{name: array}. Diagnostics: one record per field. Pickups written as 2-D records
    (nrecords > nFlds): Nr records per 3-D field, one per 2-D field (_2D)."""
    a, meta = read_mds(Path(run) / prefix, itr)
    names = meta["fldList"]
    if a.shape[0] == len(names):
        return {n: a[i] for i, n in enumerate(names)}
    nr = int(case(run)["nr"])
    out, k = {}, 0
    for n in names:
        m = 1 if n in _2D else nr
        out[n] = a[k] if m == 1 else a[k:k + m]
        k += m
    if k != a.shape[0]:
        raise ValueError(f"{prefix}.{itr}: {a.shape[0]} records, fields need {k}")
    return out


def year_iter(run: Path, year: int) -> int:
    """Iteration at the end of model year `year` (1-based) of the run."""
    c = case(run)
    return int(c["nIter0"]) + int(year) * int(c["steps_per_year"])
