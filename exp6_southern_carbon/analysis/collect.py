"""Read the coarse campaign into the inputs of H6a and H6b (plan Task 8; notes/pd0_preregistro.md).

Branch and parent share the absolute iteration count (the branch starts from the parent's pickup
at its nIter0), so annual means are keyed by the absolute iteration that closes each model year
and a branch year is compared with the parent year that ends at the same iteration.

Per branch (window = branch years 6-15, PD0):
  S         mean of F_branch - F_parent, F = channel-mean DICCFLX in mol C m-2 yr-1, + into ocean
  mld_inv   maximum of the monthly climatology of MXLDEPTH averaged over R (m)
  moc_sv    upper-cell residual MOC leaving R, max_l psi(l, j_edge) (Sv; analysis/spinup.py)
  measured  branch - parent of adv + mix + bio of the DIC budget of the mixed layer of R (mol C yr-1)
Per parent (decade centred on the branch year for the gradients, PD0; window years for e):
  g200, ddic1000, dpo4dz200, dpo4_1000 (analysis/gradients.py, umol kg-1 [m-1]), area_R (m2),
  r_cp = 117 (pkg/dic/dic_readparms.F:198, not overridden in data.dic),
  e = export / supply of PO4 in the mixed layer of R = -gchem / (adv + mix) of the PO4 budget."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import read_bin

from . import budget, flux, gradients, io, spinup

R_CP = 117.0
J_EDGE = 20
WINDOW = (6, 15)


def annual_means(iters, values, start: int, steps_per_year: int) -> dict[int, float]:
    """Monthly records -> {absolute end iteration of the model year: annual mean}; years with
    fewer than 12 records are dropped."""
    it = np.asarray(iters, dtype=np.int64)
    v = np.asarray(values, float)
    yr = np.ceil((it - int(start)) / int(steps_per_year)).astype(int)
    out = {}
    for y in np.unique(yr):
        sel = yr == y
        if sel.sum() == 12:
            out[int(start) + int(y) * int(steps_per_year)] = float(v[sel].mean())
    return out


def branch_sensitivity(f_branch: dict, f_parent: dict, n0: int, steps_per_year: int, window=WINDOW) -> float:
    keys = [int(n0) + k * int(steps_per_year) for k in range(window[0], window[1] + 1)]
    return float(np.mean([f_branch[i] - f_parent[i] for i in keys]))


def winter_mld(monthly_R) -> float:
    m = np.asarray(monthly_R, float)
    return float(m.reshape(-1, 12).mean(axis=0).max())


def export_fraction(po4_budget: dict) -> float:
    return float(-po4_budget["gchem"] / (po4_budget["adv"] + po4_budget["mix"]))


# ---- readers (run directories) ----

def region(run) -> np.ndarray:
    return spinup._region(io.grid(run), J_EDGE)


def _start(run) -> tuple[int, int]:
    c = io.case(run)
    return int(c["nIter0"]), int(c["steps_per_year"])


def flux_series(run) -> dict[int, float]:
    g = io.grid(run)
    n0, spy = _start(run)
    its = io.iters(run, "surf_mon")
    wet = g["hfac"][0] > 0
    vals = [flux.channel_mean(flux.flux_into_ocean(io.fields(run, "surf_mon", i)["DICCFLX"]), g["rac"], wet)
            for i in its]
    return annual_means(its, vals, n0, spy)


def _window_iters(run, n0_branch: int, spy: int, prefix: str, window=WINDOW) -> list[int]:
    lo, hi = n0_branch + (window[0] - 1) * spy, n0_branch + window[1] * spy
    return [i for i in io.iters(run, prefix) if lo < i <= hi]


def circulation(run, n0_branch: int, window=WINDOW) -> dict:
    """mld_inv and moc_sv of `run` over the branch window (branch years 6-15 after n0_branch)."""
    g = io.grid(run)
    _, spy = _start(run)
    R = region(run)
    mon = _window_iters(run, n0_branch, spy, "surf_mon", window)
    mld = [io.fields(run, "surf_mon", i)["MXLDEPTH"][R].mean() for i in mon]
    dxg = read_bin(Path(run) / "DXG.data", (g["ny"], g["nx"]))
    ann = _window_iters(run, n0_branch, spy, "layers_ann", window)
    moc = [spinup.moc_out_of_R(spinup.residual_moc(io.fields(run, "layers_ann", i)["LaVH1TH"], dxg), J_EDGE)
           for i in ann]
    return {"mld_inv": winter_mld(mld), "moc_sv": float(np.mean(moc))}


def _year_of(run, n0_branch: int, k: int) -> int:
    n0, spy = _start(run)
    return (n0_branch - n0) // spy + k


def ml_budget(run, n0_branch: int, tracer: int, window=WINDOW) -> dict:
    """Mean over the window of the mixed-layer budget of R (budget.dic_budget, layer='ml')."""
    R = region(run)
    bs = [budget.dic_budget(run, _year_of(run, n0_branch, k), R, "ml", tracer)
          for k in range(window[0], window[1] + 1)]
    return {key: float(np.mean([b[key] for b in bs])) for key in bs[0]}


def parent_state(run, n0_branch: int) -> dict:
    """Gradients over R in the parent decade centred on the branch year, area_R, r_cp, e."""
    g = io.grid(run)
    _, spy = _start(run)
    R = region(run)
    dec = [i for i in io.iters(run, "carbon_ann") if n0_branch - 5 * spy < i <= n0_branch + 5 * spy]
    if len(dec) != 10:
        raise ValueError(f"{run}: {len(dec)} annual records in the decade centred on iteration {n0_branch}")
    f = [io.fields(run, "carbon_ann", i) for i in dec]
    dic = np.mean([x["TRAC01"] for x in f], axis=0)
    po4 = np.mean([x["TRAC03"] for x in f], axis=0)
    z = g["rc"]
    p4 = ml_budget(run, n0_branch, tracer=3)
    return {"g200": gradients.g200(dic, z, R), "ddic1000": gradients.ddic1000(dic, z, R),
            "dpo4dz200": gradients.g200(po4, z, R), "dpo4_1000": gradients.ddic1000(po4, z, R),
            "area_R": float(g["rac"][R].sum()), "r_cp": R_CP, "e": export_fraction(p4)}


def measured(branch_run, parent_run, n0_branch: int, window=WINDOW) -> float:
    """Branch - parent of adv + mix + bio of the DIC budget of the mixed layer of R (mol C yr-1)."""
    b, p = ml_budget(branch_run, n0_branch, 1, window), ml_budget(parent_run, n0_branch, 1, window)
    s = lambda d: d["adv"] + d["mix"] + d["bio"]
    return float(s(b) - s(p))
