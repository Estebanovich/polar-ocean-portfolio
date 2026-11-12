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
    drag = []
    for it in recs:
        if (run / f"seaice.{it:010d}.data").exists() and (run / f"state3d.{it:010d}.data").exists() and it > 0:
            drag.append(drag_check(run, it, wet))
    return {"inside_cells_u": int(inside.sum()), "records": res, "seaice_last": sea,
            "drag_check": drag}


RHO_CONST = 1026.96     # kg m^-3, rhoConst = rhoNil of the fitted EOS (cases.yaml eos)
C_W = 5.5e-3            # SEAICE_waterDrag (PD0 section 8)
DWAT_MIN = 0.25         # kg m^-2 s^-1, SEAICEdWatMin default (pkg/seaice/seaice_readparms.F:403)


def drag_check(run: Path, it: int, wet) -> dict:
    """oceTAUX under full ice against the drag MITgcm should apply after the PD1 change of
    seaice_dynsolver.F: DWATN = max(rho Cw |du_c|, SEAICEdWatMin) at cell centres
    (pkg/seaice/seaice_oceandrag_coeffs.F:84-107; |du_c| from the two neighbouring U and V
    points), tau_x = (DWATN_i + DWATN_i-1)/2 (u_ice - u_1) at U points
    (pkg/seaice/seaice_ocean_stress.F:109-119, no turning angle)."""
    s3 = read_mds(run / "state3d", it)[0]
    si = read_mds(run / "seaice", it)[0]
    tx = read_mds(run / "surf", it)[0][0]
    du, dv = si[2] - s3[0][0], si[3] - s3[1][0]
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
    return {"iter": it, "cells": int(m.sum()), "max_abs_tau": float(np.abs(tx[m]).max()),
            "max_abs_expected": float(np.abs(fu[m]).max()), "max_abs_err": float(err.max()),
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
