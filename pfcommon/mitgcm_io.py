"""Binary and MDS I/O for MITgcm (big-endian float64, as in the thesis setup)."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

_DTYPES = {"float64": ">f8", "float32": ">f4"}


def write_bin(path: Path, arr: np.ndarray) -> None:
    """Write `arr` as big-endian float64 in C order; shape (nz, ny, nx) or (ny, nx)."""
    np.ascontiguousarray(arr, dtype=">f8").tofile(Path(path))


def read_bin(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    """Read a big-endian float64 file written by `write_bin` (returns native float64)."""
    data = np.fromfile(Path(path), dtype=">f8")
    if data.size != int(np.prod(shape)):
        raise ValueError(f"{path}: {data.size} values, expected shape {shape}")
    return data.reshape(shape).astype(np.float64)


def _parse_meta(text: str) -> dict:
    def block(name):
        m = re.search(rf"{name}\s*=\s*[\[{{](.*?)[\]}}]\s*;", text, re.S)
        return m.group(1) if m else None

    meta = {}
    dims = block("dimList")
    if dims is None:
        raise ValueError("meta file has no dimList")
    nums = [int(v) for v in re.findall(r"-?\d+", dims)]
    # dimList = (global size, first index, last index) per dimension
    meta["dimList"] = nums
    meta["shape"] = tuple(nums[i] for i in range(len(nums) - 3, -1, -3))  # (nz?, ny, nx)
    prec = block("dataprec")
    meta["dataprec"] = re.findall(r"'([^']+)'", prec)[0].strip() if prec else "float64"
    nrec = block("nrecords")
    meta["nrecords"] = int(re.findall(r"\d+", nrec)[0]) if nrec else 1
    fld = block("fldList")
    meta["fldList"] = [s.strip() for s in re.findall(r"'([^']*)'", fld)] if fld else []
    its = block("timeStepNumber")
    meta["timeStepNumber"] = int(re.findall(r"\d+", its)[0]) if its else None
    return meta


def read_mds(prefix: Path, itr: int | None = None) -> tuple[np.ndarray, dict]:
    """Read MDS output `<prefix>.<itr:010d>.{meta,data}`.

    Returns (array of shape (nrec, nz, ny, nx) with size-1 dimensions squeezed
    (the record axis is always kept), meta dict). `itr=None` takes the last iteration.
    """
    prefix = Path(prefix)
    if itr is None:
        files = sorted(prefix.parent.glob(prefix.name + ".[0-9]*.meta"))
        if not files:
            raise FileNotFoundError(f"no MDS files for {prefix}")
        stem = files[-1].parent / files[-1].name[: -len(".meta")]
    else:
        stem = prefix.parent / f"{prefix.name}.{itr:010d}"
    meta = _parse_meta(Path(f"{stem}.meta").read_text())
    dt = _DTYPES[meta["dataprec"]]
    data = np.fromfile(Path(f"{stem}.data"), dtype=dt)
    full = (meta["nrecords"], *meta["shape"])
    if data.size != int(np.prod(full)):
        raise ValueError(f"{stem}.data: {data.size} values, expected shape {full}")
    arr = data.reshape(full).astype(np.float64)
    keep = [0] + [i for i in range(1, arr.ndim) if arr.shape[i] != 1]
    return arr.reshape([arr.shape[i] for i in keep]), meta
