#!/opt/anaconda3/envs/pyocean/bin/python
r"""PD1a evaluation of the Gamma_T calibration runs (Exp. 2, Task 4 phase A).

Usage::

    pd1_calibration_eval.py --runs-root /Volumes/LaCie/portfolio_runs/exp2 \
        [--out exp2_isomip_cavity/figures/pd1_calibration_eval.json]

Mechanical against ``notes/pd0_preregistro.md`` (K1, V1, N1-N3, resolution test) and the PD2
sampling amendment (``notes/pd1_calibracion.md`` sect. 4.2). Also the controller rulings:

* R15 (supercooling, not a criterion): cells of the end state with T < T_f(S, p), with the
  AD16 liquidus T_f = lam1 S + lam2 + lam3 p (AD16 Eq. 25, Table 4) and p = rho0 g |z_c| at
  the cell centre (as ``metrics.gade_line_check``; the wet-part centre of a partial cell is
  lower, so p is slightly underestimated there: T_f a few mK warmer). Located as
  (i) wall rows: the first and last wet rows of the ISOMIP+ domain (y = 0 and y = 80 km
  walls), (ii) any lateral dry neighbour at the same level (side walls of the trough, ice
  front, grounded ice), (iii) thin top cells: top wet cell of the column with hFacC < 0.5.
* R17 (CFL, informative): maxima over the whole run of the %MON advcfl_* lines of the model
  log; > 0.5 is flagged.

Resolution test (pd0, descriptive, no threshold): Ocean0 2 km vs 4 km, each with its own
Gamma_T,cal: relative differences (4 km - 2 km) / 2 km of the whole-shelf mean melt, max
overturning and cavity max |psi_bt| (last 6 months); and the spatial correlation of the melt
maps, the 2 km map (melt in non-shelf cells = 0) averaged area-weighted over 2 x 2 blocks onto
the 4 km grid, over the 4 km cells where either map has shelf. |shelf melt difference| > 25 %
is flagged (pd0: highlighted in the README as an ensemble limitation).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp2_isomip_cavity.analysis import calibrate as cal  # noqa: E402
from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402
from exp2_isomip_cavity.theory import melt as ml  # noqa: E402

V1_RANGE = (0.011, 0.2)          # pd0 V1: Y26 Table 2 COM range
N1_MAX = 0.03                    # pd0 N1
N2_MAX = 0.01                    # pd0 N2
CFL_MAX = 0.5                    # Task 4 brief: advection / internal waves
SHELF_DIFF_FLAG = 0.25           # pd0 resolution test
THIN_HFAC = 0.5                  # R15: thin top cells
CFL_KEYS = ("advcfl_uvel_max", "advcfl_vvel_max", "advcfl_wvel_max", "advcfl_W_hf_max")
GAMMAS = ("0.0075", "0.015", "0.03", "0.06")


# ---- R17: CFL maxima ---------------------------------------------------------------------

def cfl_maxima(log_text: str) -> dict:
    """Max over all %MON lines of each advective CFL monitor (model log). Raises if a key
    is missing (a log without monitor output cannot be evaluated)."""
    out = {}
    for k in CFL_KEYS:
        vals = [abs(float(v.replace("D", "E")))
                for v in re.findall(rf"%MON {k}\s*=\s*(\S+)", log_text)]
        if not vals:
            raise ValueError(f"no %MON {k} lines in the log")
        out[k] = max(vals)
    out["max"] = max(out[k] for k in CFL_KEYS)
    out["n_monitor_lines"] = len(re.findall(r"%MON advcfl_uvel_max", log_text))
    out["above_0.5"] = bool(out["max"] > CFL_MAX)
    return out


def run_log(run_dir: Path) -> Path:
    p = Path(run_dir) / "STDOUT.0000"
    return p if p.is_file() and p.stat().st_size > 0 else Path(run_dir) / "output.txt"


# ---- R15: supercooling -------------------------------------------------------------------

def freezing_point(S, z_c, rho0=1027.51, g=9.81):
    """AD16 Eq. 25 liquidus [degC] at p = rho0 g |z| [Pa] (AD16 Table 4 constants)."""
    p = rho0 * g * np.abs(z_c)
    return ml.LAMBDA1 * np.asarray(S) + ml.LAMBDA2 + ml.LAMBDA3 * p


def supercooling_stats(T, S, hfac, rc, wall_rows=(0, None), rho0=1027.51, g=9.81) -> dict:
    """R15 statistics of T - T_f(S, p) over wet cells (hfac > 0). ``T``, ``S``, ``hfac``:
    (nz, ny, nx); ``rc``: (nz,) cell-centre z [m]; ``wall_rows``: (first, last) wet row
    index of the ISOMIP+ domain (last None = last row of the array)."""
    T, S, hfac = (np.asarray(a, float) for a in (T, S, hfac))
    nz, ny, nx = T.shape
    wet = hfac > 0
    dT = T - freezing_point(S, np.asarray(rc, float)[:, None, None], rho0, g)
    sc = wet & (dT < 0)
    j0, j1 = wall_rows[0], (ny - 1 if wall_rows[1] is None else wall_rows[1])
    rows = np.zeros_like(wet)
    rows[:, [j0, j1], :] = True
    pad = np.pad(wet, ((0, 0), (1, 1), (1, 1)), constant_values=False)
    lateral_dry = ~(pad[:, :-2, 1:-1] & pad[:, 2:, 1:-1] & pad[:, 1:-1, :-2] & pad[:, 1:-1, 2:])
    ktop = np.argmax(wet, axis=0)
    top = np.zeros_like(wet)
    kk, jj, ii = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    top[...] = (kk == ktop[None]) & wet.any(axis=0)[None]
    thin_top = top & (hfac < THIN_HFAC)

    def part(mask):
        m = sc & mask
        return {"count": int(m.sum()),
                "min_dT": float(dT[m].min()) if m.any() else None,
                "n_below_minus_0.1": int((m & (dT < -0.1)).sum())}

    out = {"n_wet": int(wet.sum()), "count": int(sc.sum()),
           "min_dT": float(dT[sc].min()) if sc.any() else None,
           "n_below_minus_0.1": int((sc & (dT < -0.1)).sum()),
           "wall_rows": part(rows), "interior_rows": part(~rows),
           "lateral_dry_neighbour": part(lateral_dry & ~rows),
           "no_dry_neighbour": part(~lateral_dry & ~rows),
           "thin_top_cells": part(thin_top), "other_cells": part(~thin_top),
           "top_cells_any_hfac": part(top)}
    if sc.any():
        k, j, i = np.unravel_index(np.where(sc, dT, np.inf).argmin(), dT.shape)
        out["coldest"] = {"k": int(k), "j": int(j), "i": int(i), "z_c": float(rc[k]),
                          "hFacC": float(hfac[k, j, i]), "T": float(T[k, j, i]),
                          "S": float(S[k, j, i]), "dT": float(dT[k, j, i])}
    return out


# ---- resolution test ---------------------------------------------------------------------

def coarsen_2x2(field2, rA2, ny_wet4: int, nx4: int) -> np.ndarray:
    """Area-weighted 2 x 2 block mean of a 2 km field onto the 4 km grid (rows 0..ny_wet4-1;
    the land row is not used)."""
    f = np.asarray(field2, float)[: 2 * ny_wet4, : 2 * nx4]
    a = np.asarray(rA2, float)[: 2 * ny_wet4, : 2 * nx4]
    fa = (f * a).reshape(ny_wet4, 2, nx4, 2).sum(axis=(1, 3))
    return fa / a.reshape(ny_wet4, 2, nx4, 2).sum(axis=(1, 3))


def melt_map(rd, win) -> np.ndarray:
    m = mt.melt_rate_m_per_a(rd.mean_field("shiDiag", "SHIfwFlx", win))
    return np.where(rd.shelf, m, 0.0)


def resolution_test(rd2, rd4) -> dict:
    w2, w4 = mt.calendar_window(rd2, "last6"), mt.calendar_window(rd4, "last6")
    q = {}
    for name, f in (("shelf_mean_melt", lambda r, w: mt.mean_melt(r, w, "shelf")),
                    ("overturning_max", mt.overturning_max),
                    ("psi_bt_max_cavity", lambda r, w: mt.psi_bt_max(r, w, "cavity"))):
        a, b = f(rd2, w2), f(rd4, w4)
        q[name] = {"2km": a, "4km": b, "rel_diff_4_vs_2": (b - a) / a}
    ny4 = int(rd4.params["ny_wet_domain"])       # rows 0..ny4-1; land row(s) after them
    m2 = coarsen_2x2(melt_map(rd2, w2), rd2.rA, ny4, rd4.nx)
    m4 = melt_map(rd4, w4)[:ny4]
    sh2 = coarsen_2x2(rd2.shelf.astype(float), rd2.rA, ny4, rd4.nx) > 0
    sel = sh2 | rd4.shelf[:ny4]
    r = float(np.corrcoef(m2[sel], m4[sel])[0, 1])
    q["melt_map_correlation"] = {"r": r, "n_cells": int(sel.sum()),
                                 "method": "2 km map (0 off-shelf) area-averaged 2x2 onto "
                                           "4 km; cells where either has shelf"}
    q["flag_shelf_melt_diff_gt_25pct"] = bool(
        abs(q["shelf_mean_melt"]["rel_diff_4_vs_2"]) > SHELF_DIFF_FLAG)
    return q


# ---- per-run evaluation ------------------------------------------------------------------

def evaluate_run(run_dir: Path, supercool: bool = False) -> dict:
    rd = mt.load_run(run_dir)                     # N3: require_ok first
    win = mt.calendar_window(rd, "last6")
    n1 = mt.equilibrium_drift(rd, "2yr")
    n2 = mt.tracer_budget_residual(rd)
    out = {"run_dir": str(run_dir), "gammaT": float(rd.params["gammaT"]),
           "deltaT": rd.params.get("deltaT"),
           "deep_melt_last6": mt.mean_melt(rd, win, "deep300"),
           "shelf_melt_last6": mt.mean_melt(rd, win, "shelf"),
           "N1": {"value": n1, "pass": bool(n1 < N1_MAX)},
           "N2": {"value": n2, "pass": bool(abs(n2) < N2_MAX)},
           "N3": {"pass": True, "source": "metrics.load_run -> pfcommon.runctl.require_ok"},
           "R17_cfl": cfl_maxima(run_log(run_dir).read_text(errors="replace"))}
    if supercool:
        last = rd.stage_dirs[-1]
        itr = int(rd.params["nIter0"]) + int(rd.params["nTimeSteps"])
        T = mt._read_state(last, "T", itr)
        S = mt._read_state(last, "S", itr)
        nyw = int(rd.params["ny_wet_domain"])    # ISOMIP+ rows 0..nyw-1 (land row after)
        out["R15_supercooling_end_state"] = supercooling_stats(
            T, S, rd.hfac, rd.rc, wall_rows=(0, nyw - 1),
            rho0=float(rd.params.get("rho_ref", 1027.51)),
            g=float(rd.params.get("gravity", 9.81)))
        Tm = rd.mean_field("dynDiag", "THETA", win)
        Sm = rd.mean_field("dynDiag", "SALT", win)
        out["R15_supercooling_last6_mean"] = supercooling_stats(
            Tm, Sm, rd.hfac, rd.rc, wall_rows=(0, nyw - 1))
    return out, rd


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="PD1a evaluation of the calibration runs")
    ap.add_argument("--runs-root", type=Path, default=Path("/Volumes/LaCie/portfolio_runs/exp2"))
    ap.add_argument("--out", type=Path,
                    default=REPO / "exp2_isomip_cavity" / "figures" / "pd1_calibration_eval.json")
    a = ap.parse_args(argv)
    res = {"runs": {}, "K1": {}, "V1": None}
    rds = {}
    for res_lbl, pre, ver in (("4km", "cal4", "ocean0_4km"), ("2km", "cal2", "ocean0_2km")):
        for g in GAMMAS:
            name = f"{pre}_gT{g}"
            res["runs"][name], _ = evaluate_run(a.runs_root / name)
        res["runs"][ver], rds[res_lbl] = evaluate_run(a.runs_root / ver, supercool=True)
        res["K1"][res_lbl] = {**cal.k1_verdict(res["runs"][ver]["deep_melt_last6"]),
                              "gammaT_cal": res["runs"][ver]["gammaT"],
                              "verification_run": ver}
    g2 = res["K1"]["2km"]["gammaT_cal"]
    v1 = V1_RANGE[0] <= g2 <= V1_RANGE[1]
    res["V1"] = {"gammaT_cal_2km": g2, "range": list(V1_RANGE), "pass": v1,
                 "class": None if v1 else "(c) numerical choice (PD2 amendment, "
                 "notes/pd1_calibracion.md 4.2: z-level, 20 m top cells; Gwyther 2020 L30-31)"}
    res["resolution_test"] = resolution_test(rds["2km"], rds["4km"])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=2, default=float) + "\n")
    print(json.dumps({"K1": res["K1"], "V1": res["V1"]}, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
