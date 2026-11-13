"""Mandatory ice-open-boundary distance check of the Exp. 4 ice runs (notes/pd1_seaice.md section 8).

Why. The fjord_vmf build defines DISABLE_SEAICE_OBCS (code/SEAICE_OPTIONS.h), which removes the
open-boundary calls of pkg/seaice/seaice_model.F:190-195 and :309-314 (OBCS_ADJUST_UVICE,
OBCS_APPLY_SEAICE); in free drift (F4) pkg/seaice/seaice_dynsolver.F:324-328 still sets the ice
velocity on the boundary to OB*uice/OB*vice = 0 (pkg/obcs/obcs_init_variables.F:221, 241, 283).
Ice reaching the boundaries would therefore go unnoticed. The check: over every seaice record of
every stage of the run, the ice (SIarea > 0) stays at >= n_sponge cells (10, the sponge width) of
every wet open-boundary cell (W column, N row, S row of data.obcs). A run that fails is declared
affected.

Enforcement. ``check_and_stamp`` writes the result to ``figures/ice_ob_check.json`` keyed by the
run directory, with the sha256 of its run status file; ``require_ice_ob_check`` raises unless a
current (same status) passing stamp exists. ``load.run`` calls it for every run whose case.json
has ice, so no analysis (validate.py, Task 7) can read an ice run that was not checked.

Record coverage. A seaice snapshot labelled n holds the state at the start of step n + 1
(model/src/do_statevars_diags.F:96), so the state after the last step (one deltaT = 120 s later
than the last record) is not checked; with ice speeds <= 0.22 m s^-1 that is < 27 m, 0.13 cell.

CLI: ``validate.py ice-check RUN [RUN ...]`` (exit 1 if any run fails).
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import _parse_meta

PKG = Path(__file__).resolve().parents[1]
STAMPS = PKG / "figures" / "ice_ob_check.json"


def _case(run_dir: Path) -> dict | None:
    f = Path(run_dir) / "input" / "case.json"
    return json.loads(f.read_text()) if f.exists() else None


def is_ice_run(run_dir: Path) -> bool:
    """True if the run's case.json (run_campaign layout: <run>/input/case.json) has ice."""
    cj = _case(run_dir)
    return bool(cj and (cj.get("actual") or {}).get("ice"))


def stage_dirs(run_dir: Path) -> list[Path]:
    """The run directory and its stage<k> subdirectories (exactly ``stage`` + digits, so a
    renamed ``stage2_old`` is ignored), in stage order."""
    run_dir = Path(run_dir)
    return [run_dir] + sorted((p for p in run_dir.iterdir() if p.is_dir() and re.fullmatch(r"stage\d+", p.name)),
                              key=lambda p: int(p.name[5:]))


def ob_distance_field(wet, ob_sides) -> np.ndarray:
    """Euclidean distance (cells) of every cell to the nearest wet open-boundary cell (W column,
    S row, N row; j = 0 is the south row), as analysis/pd1_timing.ob_distance_field."""
    from scipy import ndimage
    wet = np.asarray(wet, bool)
    ob = np.zeros_like(wet)
    if "W" in ob_sides:
        ob[:, 0] = wet[:, 0]
    if "E" in ob_sides:
        ob[:, -1] = wet[:, -1]
    if "S" in ob_sides:
        ob[0, :] = wet[0, :]
    if "N" in ob_sides:
        ob[-1, :] = wet[-1, :]
    return ndimage.distance_transform_edt(~ob)


def _area(data: Path, meta: dict) -> np.ndarray:
    """SIarea of one seaice record, memory-mapped (only that field is read)."""
    if meta.get("dataprec") != "float64":
        raise ValueError(f"{data}: dataprec {meta.get('dataprec')}, expected float64")
    k = meta["fldList"].index("SIarea")
    ny, nx = meta["shape"][-2:]
    m = np.memmap(data, dtype=">f8", mode="r", shape=(meta["nrecords"], ny, nx))
    return np.array(m[k])


def ice_ob_distance(run_dir: Path, min_cells: int | None = None) -> dict:
    """Smallest distance (cells, between cell centres) from ice (SIarea > 0) to a wet open-boundary
    cell over every seaice record of every stage. passes: distance >= min_cells (default the
    case's n_sponge) at every record, AND the records cover the run: every stage holds at least
    one record and the last record read is within one output interval (the record spacing) of
    expected_last_iter (run_status_final.json, else last_iter of the last stage). The final
    state itself is never a diagnostic (module docstring), hence the one-interval allowance."""
    run_dir = Path(run_dir)
    cj = _case(run_dir)
    if cj is None:
        raise RuntimeError(f"{run_dir}: no input/case.json")
    p = cj["params"]
    n = int(p["n_sponge"]) if min_cells is None else int(min_cells)
    ny, nx = int(p["ny"]), int(p["nx"])
    wet = np.fromfile(run_dir / "bathy.bin", ">f8").reshape(ny, nx) < 0
    dist = ob_distance_field(wet, p["ob_sides"])
    dmin, at, nrec, seen, empty = None, None, 0, set(), []
    for d in stage_dirs(run_dir):
        metas = sorted(d.glob("seaice.*.meta"))
        if not metas:
            empty.append(d.name)
        for m in metas:
            it = int(re.match(r"seaice\.(\d{10})\.meta", m.name).group(1))
            if it in seen:                                  # a stage boundary written twice
                continue
            seen.add(it)
            a = _area(m.with_suffix(".data"), _parse_meta(m.read_text()))
            nrec += 1
            if not (a > 0).any():
                continue
            dd = float(dist[a > 0].min())
            if dmin is None or dd < dmin:
                dmin, at = dd, it
    its = sorted(seen)
    fin = run_dir / "run_status_final.json"
    if fin.exists():
        expected = int(json.loads(fin.read_text())["expected_last_iter"])
    else:
        expected = int(json.loads((stage_dirs(run_dir)[-1] / "run_status.json").read_text())["last_iter"])
    step = int(np.median(np.diff(its))) if len(its) > 1 else None
    reason = []
    if empty:
        reason.append("no seaice record in " + ", ".join(empty))
    if not its or step is None or its[-1] < expected - step:
        reason.append(f"records do not cover the run (last {its[-1] if its else None}, "
                      f"expected {expected}, interval {step})")
    dist_ok = dmin is None or dmin >= n
    if not dist_ok:
        reason.append(f"ice within {dmin} cells of an open boundary at iteration {at}")
    return {"min_cells": dmin, "at_iter": at, "required_cells": n, "records": nrec,
            "first_iter": its[0] if its else None, "last_iter_read": its[-1] if its else None,
            "expected_last_iter": expected, "record_interval_iter": step,
            "passes": bool(not reason), "reason": "; ".join(reason)}


def status_sha(run_dir: Path) -> str:
    """sha256 of the run's final status (run_status_final.json, else run_status.json) and of the
    sorted (name, size, mtime) of every seaice.*.data file of every stage: a stamp goes stale if
    the run or its sea-ice output changes."""
    run_dir = Path(run_dir)
    f = run_dir / "run_status_final.json"
    if not f.exists():
        f = run_dir / "run_status.json"
    h = hashlib.sha256(f.read_bytes())
    files = sorted((str(p.relative_to(run_dir)), p.stat().st_size, p.stat().st_mtime_ns)
                   for d in stage_dirs(run_dir) for p in d.glob("seaice.*.data"))
    h.update(json.dumps(files).encode())
    return h.hexdigest()


def check_and_stamp(run_dir: Path, stamps: Path | None = None) -> dict:
    """Run ice_ob_distance and record it in ``stamps`` (default figures/ice_ob_check.json)."""
    stamps = Path(stamps or STAMPS)
    res = ice_ob_distance(run_dir)
    res.update(status_sha256=status_sha(run_dir), source="notes/pd1_seaice.md section 8")
    allr = json.loads(stamps.read_text()) if stamps.exists() else {}
    allr[str(Path(run_dir).resolve())] = res
    stamps.parent.mkdir(parents=True, exist_ok=True)
    stamps.write_text(json.dumps(allr, indent=2))
    return res


def require_ice_ob_check(run_dir: Path, stamps: Path | None = None) -> None:
    """RuntimeError unless ``run_dir`` has a passing, current ice-OB stamp."""
    stamps = Path(stamps or STAMPS)
    key = str(Path(run_dir).resolve())
    allr = json.loads(stamps.read_text()) if stamps.exists() else {}
    rec = allr.get(key)
    if rec is None:
        raise RuntimeError(f"{run_dir}: ice run without the mandatory ice-open-boundary check; "
                           f"run `validate.py ice-check {run_dir}` first (notes/pd1_seaice.md 8)")
    if rec.get("status_sha256") != status_sha(run_dir):
        raise RuntimeError(f"{run_dir}: the ice-open-boundary stamp is stale (the run changed); "
                           "rerun validate.py ice-check")
    if not rec.get("passes"):
        raise RuntimeError(f"{run_dir}: ice-open-boundary check failed: {rec.get('reason')} "
                           f"(ice within {rec.get('min_cells')} cells, required >= "
                           f"{rec.get('required_cells')}); the run is affected (notes/pd1_seaice.md 8)")
