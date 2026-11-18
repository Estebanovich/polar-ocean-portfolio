"""DIC budget of a region over one model year (spec Sec. 3.5, 6 V6b; plan Task 5).

Terms (mol C yr-1, positive = gain of the region):
  adv   convergence of ADVxTr01/ADVyTr01/ADVrTr01 (resolved + GM bolus advection; pkg/ptracers
        diagnostics, fluxes in mol s-1 through the west, south and upper faces, ADVr positive up)
  mix   convergence of DFxETr01/DFyETr01/DFrETr01/DFrITr01 (explicit, GM-Redi and implicit vertical)
  forc  ForcTr01 x volume (pkg/ptracers forcing tendency: the rbcs sponge; no E-P here)
  gchem the tendency pkg/dic applies in GCHEM_FORCING_SEP: GDIC = car + R_CP GPO4 (+ SURC at
        the top; pkg/dic/dic_biotic_forcing.F:311-374), diagnosed as DICGDIC (and DICGPO4 for
        PO4) by the EXP6 blocks of code/dic_biotic_forcing.F and code/dic_diagnostics_init.F.
        (TRAC01 - GC_Tr01 is not usable: TRAC01 is filled at the start of the step, so the
        difference is minus the transport and forcing increment, and the budget would close
        by construction.)
  surf  the air-sea part of gchem, from DICTFLX (top-level tendency, mol m-3 s-1)
  bio   gchem - surf
  tend  inventory change between the year's first and last pickup_ptracers
  residual = tend - (adv + mix + forc + gchem)."""
from __future__ import annotations

import numpy as np

from . import io

YEAR_S = 31536000.0
_SUF = {1: "Tr01", 3: "Tr03"}


def _conv(fx, fy, fr):
    """Per-cell convergence (nr, ny, nx) of face fluxes; x periodic, y walls, fr positive up."""
    cx = fx - np.roll(fx, -1, axis=-1)
    fyn = np.concatenate([fy[:, 1:, :], np.zeros_like(fy[:, :1, :])], axis=1)
    cy = fy - fyn
    frb = np.concatenate([fr[1:], np.zeros_like(fr[:1])], axis=0)    # flux through the bottom face
    cr = frb - fr
    return cx + cy + cr


def region_cells(g: dict, region_mask2d, layer: str = "col", mld2d=None) -> np.ndarray:
    m = np.asarray(region_mask2d, bool)[None] & (g["hfac"] > 0)
    if layer == "ml":
        m &= g["rc"][:, None, None] <= np.asarray(mld2d, float)[None]
        m[0] |= np.asarray(region_mask2d, bool) & (g["hfac"][0] > 0)       # always the top level
    return m


def dic_budget(run, year: int, region_mask2d=None, layer: str = "col", tracer: int = 1) -> dict:
    g = io.grid(run)
    c = io.case(run)
    spy = int(c["steps_per_year"])
    i1, i0 = io.year_iter(run, year), io.year_iter(run, year - 1)
    if region_mask2d is None:
        region_mask2d = np.ones((g["ny"], g["nx"]), bool)
    vol = g["rac"][None] * g["drf"][:, None, None] * g["hfac"]
    d = io.fields(run, "carbon_ann", i1)
    s = _SUF[tracer]
    adv = _conv(d["ADVx" + s], d["ADVy" + s], d["ADVr" + s])
    mix = _conv(d["DFxE" + s], d["DFyE" + s], d["DFrE" + s] + d["DFrI" + s])
    forc = d["Forc" + s] * vol
    gchem = d["DICGDIC" if tracer == 1 else "DICGPO4"] * vol
    mld = None
    if layer == "ml":
        mon = [io.fields(run, "surf_mon", it)["MXLDEPTH"] for it in io.iters(run, "surf_mon") if i0 < it <= i1]
        mld = np.mean(mon, axis=0)
    m = region_cells(g, region_mask2d, layer, mld)
    p0, p1 = (io.fields(run, "pickup_ptracers", i) for i in (i0, i1))
    key = f"pTr0{tracer}"
    tend = float(((p1[key] - p0[key]) * vol)[m].sum())
    out = {k: float((v * YEAR_S)[m].sum()) for k, v in (("adv", adv), ("mix", mix), ("forc", forc), ("gchem", gchem))}
    if tracer == 1:
        surf = io.fields(run, "surf_mon", i1)["DICTFLX"] if "DICTFLX" in io.fields(run, "surf_mon", i1) else None
        mon = [io.fields(run, "surf_mon", it)["DICTFLX"] for it in io.iters(run, "surf_mon") if i0 < it <= i1]
        tf = np.mean(mon, axis=0) * vol[0]
        top = m[0]
        out["surf"] = float((tf * YEAR_S)[top].sum())
        out["bio"] = out["gchem"] - out["surf"]
    out["tend"] = tend
    out["residual"] = tend - (out["adv"] + out["mix"] + out["forc"] + out["gchem"])
    return out
