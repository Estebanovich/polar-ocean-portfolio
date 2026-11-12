"""Timing, CFL and disk numbers of the Exp. 4 PD1 short runs (Task 4) -> figures/pd1_timing.json.

For each run directory (scripts/run_campaign.py layout): wall time (run_status.json), steps,
s/step, the MITgcm main-loop timer (STDOUT.0000, "MAIN LOOP"), the maximum over all monitor
records of the advective CFL numbers MITgcm prints (pkg/monitor/monitor.F:148-163: trAdv_CFL_*,
advcfl_uvel/vvel/wvel/W_hf, computed with deltaT), the maximum speeds, the first time the
maximum vertical velocity exceeds W_NOISE (onset of the grid-scale noise, notes/pd1_config.md
section 2; threshold post hoc), the size of the
output written per diagnostics record, and, for runs with sea ice, the stress the ocean
receives (diagnostic oceTAUX/oceTAUY, after pkg/seaice) inside and outside the ice.
For an aborted run (no run_status OK) the numbers cover the steps it made.

Usage: pd1_timing.py RUN_DIR [RUN_DIR ...] [--out figures/pd1_timing.json]
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

from pfcommon.mitgcm_io import read_mds  # noqa: E402

PKG = Path(__file__).resolve().parents[1]
MON_KEYS = ("trAdv_CFL_u_max", "trAdv_CFL_v_max", "trAdv_CFL_w_max", "advcfl_uvel_max",
            "advcfl_vvel_max", "advcfl_wvel_max", "advcfl_W_hf_max", "dynstat_uvel_max",
            "dynstat_uvel_min", "dynstat_vvel_max", "dynstat_vvel_min", "dynstat_wvel_max",
            "dynstat_eta_max", "dynstat_eta_min")
NUM = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"


def monitor_series(log: str) -> dict:
    """{key: [(iteration, value), ...]} of the %MON lines (the iteration of each block is the
    preceding time_tsnumber)."""
    out: dict = {k: [] for k in MON_KEYS}
    it = None
    for line in log.splitlines():
        m = re.search(rf"%MON (\w+)\s*=\s*({NUM})", line)
        if not m:
            continue
        k, v = m.group(1), float(m.group(2).replace("D", "E"))
        if k == "time_tsnumber":
            it = int(v)
        elif k in out:
            out[k].append((it, v))
    return out


def main_loop_wall(log: str) -> float | None:
    m = re.search(r'Seconds in section "MAIN LOOP\s*\[THE_MAIN_LOOP\]":\s*\n.*User time:.*\n.*System time:.*\n'
                  rf".*Wall clock time:\s*({NUM})", log)
    return float(m.group(1)) if m else None


W_NOISE = 5e-3     # m/s; post hoc threshold (set after seeing the runs): the smooth response of
                   # the 6 h runs stays below 1.5e-3 m/s, the grid-scale noise grows past 1e-2


def stage_info(run: Path, dt: float | None = None) -> dict:
    st = json.loads((run / "run_status.json").read_text())
    log = (run / "STDOUT.0000").read_text(errors="replace")
    mon = monitor_series(log)
    n = int(st["n_timesteps"])
    out = {"state": st["state"], "n_timesteps": n, "last_iter": st["last_iter"], "wall_s": st["wall_s"],
           # steps actually made (an aborted run stops before n_timesteps; nIter0 = 0 here)
           "s_per_step": st["wall_s"] / st["last_iter"] if st.get("last_iter") else None, "main_loop_wall_s": main_loop_wall(log)}
    out["monitor_max"] = {k: (max(abs(v) for _, v in s) if s else None) for k, s in mon.items()}
    cfl = {k: max(v for _, v in mon[k]) if mon[k] else None
           for k in ("advcfl_uvel_max", "advcfl_vvel_max", "advcfl_wvel_max", "advcfl_W_hf_max")}
    out["cfl_max"] = max(v for v in cfl.values() if v is not None) if any(cfl.values()) else None
    on = [it for it, v in mon["dynstat_wvel_max"] if abs(v) > W_NOISE]
    out["noise_onset_iter"] = on[0] if on else None
    out["noise_onset_h"] = on[0] * dt / 3600 if (on and dt) else None
    out["wvel_max_series"] = mon["dynstat_wvel_max"]
    out["cfl_max_iter"] = None
    if out["cfl_max"] is not None:
        for k in cfl:
            for it, v in mon[k]:
                if v == out["cfl_max"]:
                    out["cfl_max_iter"], out["cfl_max_key"] = it, k
    return out


def output_bytes(run: Path) -> dict:
    """Bytes per record of each diagnostics stream (largest .data of the stream)."""
    out = {}
    for f in run.glob("*.data"):
        m = re.match(r"(.+)\.(\d{10})\.data$", f.name)
        if not m or m.group(1).startswith("pickup"):
            continue
        out[m.group(1)] = max(out.get(m.group(1), 0), f.stat().st_size)
    return out


def ice_stress(run: Path) -> dict | None:
    """oceTAUX (field 0 of the surf stream) at every record, at U points under full ice (A = 1 on
    both sides, from the seaice stream of the same record if written, else AREA.bin) and in
    open water (A = 0)."""
    cj = json.loads((run / "input" / "case.json").read_text()) if (run / "input" / "case.json").exists() \
        else json.loads((run / "case.json").read_text())
    if not cj["actual"].get("ice"):
        return None
    ny, nx = cj["params"]["ny"], cj["params"]["nx"]
    area = np.fromfile(run / "AREA.bin", ">f8").reshape(ny, nx)
    wet = np.fromfile(run / "bathy.bin", ">f8").reshape(ny, nx) < 0
    au = np.zeros_like(area)
    au[:, 1:] = 0.5 * (area[:, 1:] + area[:, :-1])
    inside = (au == 1.0) & wet
    outside = (au == 0.0) & wet
    recs = sorted(int(re.match(r"surf\.(\d{10})\.data", f.name).group(1)) for f in run.glob("surf.*.data"))
    res = []
    for it in recs:
        d, _ = read_mds(run / "surf", it)
        tx = d[0]
        ins, outs = inside, outside
        if (run / f"seaice.{it:010d}.data").exists():            # ice that moved (free drift)
            a_now = read_mds(run / "seaice", it)[0][0]
            au_now = np.zeros_like(a_now)
            au_now[:, 1:] = 0.5 * (a_now[:, 1:] + a_now[:, :-1])
            ins, outs = (au_now >= 1.0 - 1e-12) & wet, (au_now == 0.0) & wet
        res.append({"iter": it, "cells_inside": int(ins.sum()),
                    "max_abs_taux_inside_ice": float(np.abs(tx[ins]).max()) if ins.any() else None,
                    "max_abs_taux_open_water": float(np.abs(tx[outs]).max())})
    sea = None
    sif = sorted(run.glob("seaice.*.data"))
    if sif:
        it = int(re.match(r"seaice\.(\d{10})\.data", sif[-1].name).group(1))
        d, _ = read_mds(run / "seaice", it)
        sea = {"iter": it, "max_abs_uice": float(np.abs(d[2]).max()), "max_abs_vice": float(np.abs(d[3]).max()),
               "area_sum": float(d[0].sum())}
    edge = ice_distance_to_open_boundaries(run)
    drag = []
    for it in recs:
        if (run / f"seaice.{it:010d}.data").exists() and (run / f"state3d.{it:010d}.data").exists() and it > 0:
            drag.append(drag_check(run, it, wet))
    return {"inside_cells_u": int(inside.sum()), "records": res, "seaice_last": sea,
            "drag_check": drag, "open_boundary_distance": edge}


def ob_distance_field(wet, ob_sides) -> np.ndarray:
    """Distance (cells) of every cell to the nearest wet open-boundary cell (W column, S row,
    N row; j = 0 is the south row)."""
    from scipy import ndimage
    wet = np.asarray(wet, bool)
    ob = np.zeros_like(wet)
    if "W" in ob_sides:
        ob[:, 0] = wet[:, 0]
    if "S" in ob_sides:
        ob[0, :] = wet[0, :]
    if "N" in ob_sides:
        ob[-1, :] = wet[-1, :]
    return ndimage.distance_transform_edt(~ob)


def ice_distance_to_open_boundaries(run: Path, min_cells: int | None = None) -> dict:
    """Smallest distance (cells, Euclidean between cell centres) between ice (SIarea > 0) and
    an open-boundary cell (a wet cell of the W column, N row or S row, as OB_Iwest/OB_Jnorth/
    OB_Jsouth of data.obcs) over every seaice record of ``run``. With DISABLE_SEAICE_OBCS the
    ice gets no OB condition (code/obcs_check.F header), so ice reaching the boundaries would
    go unnoticed: the check passes if the distance stays >= ``min_cells`` (default n_sponge,
    the sponge width). Required for F2/F4 after every run (Task 6)."""
    cj = json.loads((run / "input" / "case.json").read_text())
    p = cj["params"]
    n = int(p["n_sponge"]) if min_cells is None else int(min_cells)
    wet = np.fromfile(run / "bathy.bin", ">f8").reshape(int(p["ny"]), int(p["nx"])) < 0
    dist = ob_distance_field(wet, p["ob_sides"])
    recs = sorted(int(re.match(r"seaice\.(\d{10})\.data", f.name).group(1)) for f in run.glob("seaice.*.data"))
    dmin, at = None, None
    for it in recs:
        a = read_mds(run / "seaice", it)[0][0]
        if not (a > 0).any():
            continue
        dd = float(dist[a > 0].min())
        if dmin is None or dd < dmin:
            dmin, at = dd, it
    return {"min_cells": dmin, "at_iter": at, "required_cells": n, "records": len(recs),
            "passes": bool(dmin is not None and dmin >= n)}


RHO_CONST = 1026.96     # kg m^-3, rhoConst = rhoNil of the fitted EOS (cases.yaml eos)
C_W = 5.5e-3            # SEAICE_waterDrag (PD0 section 8)
DWAT_MIN = 0.25         # kg m^-2 s^-1, SEAICEdWatMin default (pkg/seaice/seaice_readparms.F:403)


RHO_ICE = 910.0         # kg m^-3, SEAICE_rhoIce default (pkg/seaice/seaice_readparms.F:396)
GRAV = 9.81


def free_drift_velocity(run: Path, it: int, wet) -> tuple:
    """Free-drift ice velocity of step ``it`` at U and V points, recomputed as
    pkg/seaice/seaice_freedrift.F:65-170 does, from the diagnosed state at the start of the
    step: surface velocity (state3d), free surface (state2d, tilt term), ice h and A of the
    seaice stream (fields at the start of step it) and the wind stress at t = it deltaT
    (linear interpolation of the records, model/src/external_fields_load.F). Forcing on the
    ice: FORCEX0 = TAUX A_u - m_u g d(eta)/dx (seaice_dynsolver.F:266-296, scaleSurfStress and
    tilt on), TAUX = fu (SEAICE_drag = OCEAN_drag). The SIuice/SIvice diagnostics are written at
    the start of the step (model/src/do_statevars_diags.F:96), before this update."""
    import sys as _s
    if str(PKG.parent) not in _s.path:
        _s.path.insert(0, str(PKG.parent))
    from exp4_arctic_bayiw.gendata import forcing
    cj = json.loads((run / "input" / "case.json").read_text())
    w, dt = cj["actual"]["wind"], cj["actual"]["deltaT_s"]
    f0, dx = cj["actual"]["f0_s-1"], cj["actual"]["grid"]["dx_m"]
    tr = forcing.record_times(w["n_records"], w["record_period_s"])
    mag = np.interp(it * dt, tr, forcing.pulse(tr, tau_max=w["tau_max"], sigma_h=w["sigma_h"],
                                                  center_h=w["center_h"]))
    tx0, ty0 = mag * w["unit_vector_xy"][0], mag * w["unit_vector_xy"][1]
    s3 = read_mds(run / "state3d", it)[0]
    eta = read_mds(run / "state2d", it)[0][0]
    si = read_mds(run / "seaice", it)[0]
    a, h = si[0], si[1]
    u1, v1 = s3[0][0], s3[1][0]
    mw = np.zeros_like(wet, float)
    mw[:, 1:] = (wet[:, 1:] & wet[:, :-1])
    ms = np.zeros_like(wet, float)
    ms[1:, :] = (wet[1:, :] & wet[:-1, :])

    def xm(q):                     # q(i-1) at i
        r = np.zeros_like(q)
        r[:, 1:] = q[:, :-1]
        return r

    def ym(q):
        r = np.zeros_like(q)
        r[1:, :] = q[:-1, :]
        return r

    def xp(q):
        r = np.zeros_like(q)
        r[:, :-1] = q[:, 1:]
        return r

    def yp(q):
        r = np.zeros_like(q)
        r[:-1, :] = q[1:, :]
        return r

    fx = (tx0 * mw) * 0.5 * (a + xm(a)) - RHO_ICE * 0.5 * (h + xm(h)) * GRAV * (eta - xm(eta)) / dx * mw
    fy = (ty0 * ms) * 0.5 * (a + ym(a)) - RHO_ICE * 0.5 * (h + ym(h)) * GRAV * (eta - ym(eta)) / dx * ms
    tauc_x, tauc_y = 0.5 * (fx + xp(fx)), 0.5 * (fy + yp(fy))
    mcor = RHO_ICE * h * f0
    uc, vc = 0.5 * (u1 + xp(u1)), 0.5 * (v1 + yp(v1))
    rx, ry = -tauc_x - mcor * vc, -tauc_y + mcor * uc
    rn, ra = np.hypot(rx, ry), np.arctan2(ry, rx)
    t1 = 1.0 / RHO_CONST / C_W
    t2, t3 = (t1 * mcor) ** 2, (t1 * rn) ** 2
    soln = np.where(t3 > 0, np.sqrt(0.5 * (np.sqrt(t2 * t2 + 4 * t3) - t2)), 0.0)
    sola = np.where(soln > 0, ra - np.arctan2(mcor * soln, C_W * RHO_CONST * soln ** 2), 0.0)
    uic, vic = uc - soln * np.cos(sola), vc - soln * np.sin(sola)
    return 0.5 * (xm(uic) + uic) * mw, 0.5 * (ym(vic) + vic) * ms


def drag_check(run: Path, it: int, wet) -> dict:
    """oceTAUX under full ice against the drag MITgcm should apply after the PD1 change of
    seaice_dynsolver.F: DWATN = max(rho Cw |du_c|, SEAICEdWatMin) at cell centres
    (pkg/seaice/seaice_oceandrag_coeffs.F:84-107; |du_c| from the two neighbouring U and V
    points), tau_x = (DWATN_i + DWATN_i-1)/2 (u_ice - u_1) at U points
    (pkg/seaice/seaice_ocean_stress.F:109-119, no turning angle)."""
    s3 = read_mds(run / "state3d", it)[0]
    si = read_mds(run / "seaice", it)[0]
    tx = read_mds(run / "surf", it)[0][0]
    cj = json.loads((run / "input" / "case.json").read_text())
    if cj["actual"]["ice"]["mode"] == "free_drift":
        ui, vi = free_drift_velocity(run, it, wet)      # ice velocity of this step (see above)
    else:
        ui, vi = si[2], si[3]                            # landfast: 0
    du, dv = ui - s3[0][0], vi - s3[1][0]
    dup = np.zeros_like(du)
    dup[:, :-1] = du[:, 1:]
    dvp = np.zeros_like(dv)
    dvp[:-1, :] = dv[1:, :]
    mag = np.sqrt(0.25 * ((du + dup) ** 2 + (dv + dvp) ** 2))
    dw = np.maximum(RHO_CONST * C_W * mag, DWAT_MIN) * wet
    dwm = np.zeros_like(dw)
    dwm[:, 1:] = dw[:, :-1]
    fu = 0.5 * (dw + dwm) * du
    a = si[0]
    au = np.zeros_like(a)
    au[:, 1:] = 0.5 * (a[:, 1:] + a[:, :-1])
    m = (au >= 1.0 - 1e-12) & wet
    m[:, 0] = False
    err = np.abs(tx[m] - fu[m])
    # interior: >= 2 cells from the ice edge and the coast, where the reconstruction of the
    # masks (SIMaskU, A_u, h at the edge) does not matter
    from scipy import ndimage
    inner = ndimage.binary_erosion(m, iterations=2)
    ei = np.abs(tx[inner] - fu[inner]) if inner.any() else np.array([np.nan])
    return {"iter": it, "cells": int(m.sum()), "max_abs_tau": float(np.abs(tx[m]).max()),
            "interior_cells": int(inner.sum()), "interior_max_abs_err": float(ei.max()),
            "interior_p99_abs_err": float(np.percentile(ei, 99)),
            "interior_median_abs_err": float(np.median(ei)),
            "max_abs_expected": float(np.abs(fu[m]).max()), "max_abs_err": float(err.max()),
            "max_rel_err_to_max_tau": float(err.max() / max(np.abs(tx[m]).max(), 1e-30)),
            "ice_velocity": "recomputed free drift of this step" if cj["actual"]["ice"]["mode"] == "free_drift"
            else "SIuice/SIvice (landfast)",
            "p99_abs_err": float(np.percentile(err, 99)),
            "floor_fraction": float((RHO_CONST * C_W * mag < DWAT_MIN)[m].mean()),
            "max_abs_du": float(np.abs(du[m]).max())}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, default=PKG / "figures" / "pd1_timing.json")
    a = ap.parse_args(argv)
    res = json.loads(a.out.read_text()) if a.out.exists() else {}
    for run in a.runs:
        cj = json.loads((run / "input" / "case.json").read_text())
        e = {"run_dir": str(run), "case": cj["params"]["case"], "config": cj["config"],
             "deltaT_s": cj["actual"]["deltaT_s"], "grid": cj["actual"]["grid"],
             "edits": cj["actual"]["edits"], "binary_sha256": cj.get("binary_sha256"),
             **stage_info(run, cj["actual"]["deltaT_s"]), "bytes_per_record": output_bytes(run), "ice_stress": ice_stress(run)}
        res[run.name] = e
        print(run.name, {k: e[k] for k in ("state", "n_timesteps", "wall_s", "s_per_step", "cfl_max")})
    a.out.write_text(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
