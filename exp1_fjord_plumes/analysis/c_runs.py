"""Config (c) run access shared by fig3_exchange_partition.py and animate_c.py.

* ``case_parts``: the directories of one case, in time order. Two-stage cases
  (wind, all; ``data.stage2`` copied into the run dir) are ``[run, run/stage2]``
  and are accepted only if ``run_status_final.json`` is OK (both stages OK and
  the last iteration reached, scripts/run_campaign.py); single-stage cases need
  an OK run_status.json (``pfcommon.runctl.require_ok``).
* ``section_area`` / ``exchange_series``: Q_ex and heat toward the front on the
  cross-fjord section at the U points of column ``i`` (western face of cell i),
  from the 3-h means ``mean3h`` (UVEL, THETA; data.diagnostics). The section
  must lie outside the OBCS sponge: OB at column nx-2 (OB_Ieast=-2), sponge in
  the 10 columns before it, i.e. i >= nx - 12 (= 188) is refused
  (notes/pd1_c_fjord3d.md). T on the face = mean of the two adjacent cells.
  The time of a 3-h mean is the centre of its window: iteration*deltaT - 1.5 h.
* ``c5p_check``: criterion C5' (T, S stay inside the range of the initial
  field, the plume outflow T_out/S_out and the prescribed OB series
  OBEt/OBEs; the heaved OB profile reaches slightly beyond the initial range),
  from the ``%MON dynstat_{theta,salt}_{min,max}`` lines of every stage.
* ``w_over_rd``: W/Rd with the two-layer first-mode speed
  c = sqrt(g' h1 h2 / H), g' from the top and bottom initial T/S (linear EOS),
  h1 = depth of the interface of max |d rho/dz| of the initial profile
  (200 m, the pycnocline centre of continuous_sermilik, in (c)), Rd = c/f0.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from exp1_fjord_plumes.analysis.exchange import exchange_flux, heat_toward_front
from exp1_fjord_plumes.theory.plume import LinearEOS
from pfcommon.mitgcm_io import read_bin, read_mds
from pfcommon.runctl import require_ok

SPONGE = 10            # data.obcs spongeThickness
MEAN_S = 10800.0       # averaging period of mean3h (data.diagnostics)
C5P_TOL = 1e-6         # round-off tolerance of the "exact" C5'


def params(run_dir: Path) -> dict:
    return json.loads((Path(run_dir) / "case.json").read_text())


def case_parts(run_dir: Path) -> list[Path]:
    run_dir = Path(run_dir)
    if (run_dir / "data.stage2").exists():
        f = run_dir / "run_status_final.json"
        if not f.exists():
            raise RuntimeError(f"{run_dir}: two-stage case without run_status_final.json")
        fin = json.loads(f.read_text())
        if fin.get("state") != "OK" or fin.get("last_iter") != fin.get("expected_last_iter"):
            raise RuntimeError(f"{run_dir}: final state {fin.get('state')} "
                               f"(last_iter {fin.get('last_iter')}), refusing to use it")
        return [run_dir, run_dir / "stage2"]
    require_ok(run_dir)
    return [run_dir]


def mds_iters(dirs: list[Path], name: str) -> list[tuple[int, Path]]:
    """[(iteration, directory)] of ``name.<it>.meta`` over all parts, sorted."""
    out: dict[int, Path] = {}
    for d in dirs:
        for f in Path(d).glob(f"{name}.[0-9]*.meta"):
            m = re.fullmatch(rf"{re.escape(name)}\.(\d+)\.meta", f.name)
            if m:
                out.setdefault(int(m.group(1)), Path(d))
    return sorted(out.items())


def depth(run_dir: Path) -> np.ndarray:
    p = params(run_dir)["params"]
    return -read_bin(Path(run_dir) / "topog.bin", (p["ny"], p["nx"]))


def section_area(run_dir: Path, i: int) -> np.ndarray:
    """Face areas dA (nz, ny) [m2] of the U points of column i; 0 on land."""
    p = params(run_dir)["params"]
    nx, nz = p["nx"], p["nz"]
    first_sponge = nx - 2 - SPONGE
    if not 1 <= i < first_sponge:
        raise ValueError(f"section i={i} must be in 1..{first_sponge - 1} "
                         f"(OBCS sponge covers i >= {first_sponge})")
    dep = depth(run_dir)
    zb = (np.arange(nz) + 1) * p["dz"]                       # cell bottoms, depth
    wet = (zb[:, None] <= dep[None, :, i - 1] + 1e-9) & (zb[:, None] <= dep[None, :, i] + 1e-9)
    return wet * p["dy"] * p["dz"]


def exchange_series(run_dir: Path, i: int, T_ref: float) -> dict:
    """Q_ex(t) [m3/s] and heat toward the front(t) [W] on section i."""
    dirs = case_parts(run_dir)
    p = params(run_dir)["params"]
    dA = section_area(run_dir, i)
    wet = dA > 0
    t, q, h = [], [], []
    for it, d in mds_iters(dirs, "mean3h"):
        a = read_mds(d / "mean3h", it)[0]                    # (2, nz, ny, nx)
        u = np.where(wet, a[0][:, :, i], 0.0)
        T = np.where(wet, 0.5 * (a[1][:, :, i - 1] + a[1][:, :, i]), T_ref)
        t.append(it * p["deltaT"] - MEAN_S / 2)
        q.append(exchange_flux(u, dA))
        h.append(heat_toward_front(u, T, dA, T_ref=T_ref))
    return {"t_s": np.array(t), "Q_ex": np.array(q), "heat_W": np.array(h)}


_MON = re.compile(r"%MON dynstat_(theta|salt)_(min|max)\s*=\s*(\S+)")


def monitor_range(dirs: list[Path]) -> dict:
    """{'T': [min, max], 'S': [min, max]} over the monitor lines of all parts."""
    vals = {("theta", "min"): [], ("theta", "max"): [], ("salt", "min"): [], ("salt", "max"): []}
    for d in dirs:
        for name in ("STDOUT.0000", "output.txt"):
            f = Path(d) / name
            if f.exists() and "%MON" in f.read_text():
                for m in _MON.finditer(f.read_text()):
                    vals[(m.group(1), m.group(2))].append(float(m.group(3)))
                break
    if not all(vals.values()):
        raise FileNotFoundError(f"no %MON dynstat T/S lines in {dirs}")
    return {"T": [min(vals[("theta", "min")]), max(vals[("theta", "max")])],
            "S": [min(vals[("salt", "min")]), max(vals[("salt", "max")])]}


def c5p_check(run_dir: Path) -> dict:
    run_dir = Path(run_dir)
    meta = params(run_dir)
    p = meta["params"]
    shape = (p["nz"], p["ny"], p["nx"])
    wet3 = np.broadcast_to(depth(run_dir)[None] > 0, shape)
    allowed = {}
    for key, ini, ob, src in (("T", "T.init.bin", "OBEt.bin", "T_out"),
                              ("S", "S.init.bin", "OBEs.bin", "S_out")):
        a = read_bin(run_dir / ini, shape)[wet3]
        lo, hi = float(a.min()), float(a.max())
        if (run_dir / ob).exists():
            o = np.fromfile(run_dir / ob, dtype=">f8")
            lo, hi = min(lo, float(o.min())), max(hi, float(o.max()))
        if meta.get(src) is not None:
            lo, hi = min(lo, float(meta[src])), max(hi, float(meta[src]))
        allowed[key] = [lo, hi]
    mon = monitor_range(case_parts(run_dir))
    ok = all(mon[k][0] >= allowed[k][0] - C5P_TOL and mon[k][1] <= allowed[k][1] + C5P_TOL
             for k in ("T", "S"))
    return {"pass": bool(ok), "monitor": mon, "allowed": allowed, "tol": C5P_TOL}


def w_over_rd(run_dir: Path, h1: float | None = None) -> dict:
    p = params(run_dir)["params"]
    shape = (p["nz"], p["ny"], p["nx"])
    T = read_bin(Path(run_dir) / "T.init.bin", shape)[:, :, 0]
    S = read_bin(Path(run_dir) / "S.init.bin", shape)[:, :, 0]
    jw = int(np.argmax(depth(run_dir)[:, 0] > 0))
    e = LinearEOS()
    if h1 is None:   # pycnocline: interface of max |d rho/dz| (200 m in (c))
        rho = e.rho(T[:, jw], S[:, jw])
        h1 = float((np.argmax(np.abs(np.diff(rho))) + 1) * p["dz"])
    drho = e.rho(T[-1, jw], S[-1, jw]) - e.rho(T[0, jw], S[0, jw])
    gp = e.g * drho / e.rho0
    H = p["nz"] * p["dz"]
    c = float(np.sqrt(gp * h1 * (H - h1) / H))
    W = float((depth(run_dir)[:, 0] > 0).sum() * p["dy"])
    Rd = c / p["f0"]
    return {"W_m": W, "Rd_m": Rd, "W_over_Rd": W / Rd, "c_m_s": c, "gprime": gp, "h1_m": h1}
