"""Spin-up report for G0 / G0s (plan Task 6): V6a, the residual MOC and the analysis region R.

Residual MOC in temperature layers from pkg/layers LaVH1TH (v h per layer, m2 s-1; includes the
GM bolus velocity because layers_bolus defaults to useGMRedi, pkg/layers/layers_readparms.F:68):
psi[l, j] = sum_x sum_{l' >= l} vh dxG (Sv), the northward transport warmer than bound l.
The upper cell upwells inside R (south wall): moc_out_of_R = max_l psi[l, j_edge], the light
water leaving R across its northern edge (the decomposition's MOC, analysis/decomp.py).
upwelling_fraction: share of the upper-cell growth d(max_l psi)/dy that lies south of j_edge."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import read_bin

from . import drift, io
from .flux import YEAR_S

SV = 1e6


def residual_moc(vh, dxg) -> np.ndarray:
    t = (np.asarray(vh, float) * np.asarray(dxg, float)[None]).sum(axis=-1)      # (nl, ny)
    return np.flip(np.cumsum(np.flip(t, 0), 0), 0) / SV


def moc_out_of_R(psi, j_edge: int) -> float:
    return float(np.max(np.asarray(psi)[:, j_edge]))


def upwelling_fraction(psi_upper, j_edge: int) -> float:
    d = np.clip(np.diff(np.asarray(psi_upper, float)), 0.0, None)
    return float(d[:j_edge - 1].sum() / d.sum()) if d.sum() > 0 else float("nan")


def run_years(ann_iters, steps_per_year: int) -> np.ndarray:
    """Model years 1..N of the annual records, counted from the start of the run directory (the
    first record closes year 1). Independent of case.json nIter0, which a --continue stage resets."""
    it = np.asarray(ann_iters, dtype=np.int64)
    return (it - (it[0] - int(steps_per_year))) / int(steps_per_year)


def _region(g, j_edge):
    m = np.zeros((g["ny"], g["nx"]), bool)
    m[1:j_edge] = True
    return m & (g["hfac"][0] > 0)


def report(run, j_edge: int = 20, last_years: int = 10) -> dict:
    run = Path(run)
    g = io.grid(run)
    c = io.case(run)
    dxg = read_bin(run / "DXG.data", (g["ny"], g["nx"]))
    R = _region(g, j_edge)
    area = g["rac"]
    vol = g["rac"][None] * g["drf"][:, None, None] * g["hfac"]
    deep = (g["rc"] >= 1000.0)[:, None, None] & (g["hfac"] > 0)
    ann = io.iters(run, "carbon_ann")
    years = list(run_years(ann, int(c["steps_per_year"])))
    deep_dic = [float((io.fields(run, "carbon_ann", it)["TRAC01"] * vol)[deep].sum() / vol[deep].sum()) for it in ann]
    mon = io.iters(run, "surf_mon")
    key = "DICCFLX" if "DICCFLX" in io.fields(run, "surf_mon", mon[0]) else "DICTFLX"
    wet0 = g["hfac"][0] > 0
    fl = []
    for it in mon:
        f = io.fields(run, "surf_mon", it)[key]
        if key == "DICTFLX":
            f = f * g["drf"][0] * g["hfac"][0]
        fl.append(float((f * area)[wet0].sum() / area[wet0].sum() * YEAR_S))
    fl_ann = np.array(fl[: 12 * len(ann)]).reshape(len(ann), 12).mean(1) if len(fl) >= 12 * len(ann) else np.array(fl)
    dec = max(1, min(10, len(ann)))
    ny_ = len(ann) // dec
    yd = np.array(years[: ny_ * dec]).reshape(ny_, dec).mean(1)
    fd = fl_ann[: ny_ * dec].reshape(ny_, dec).mean(1)
    dd = np.array(deep_dic[: ny_ * dec]).reshape(ny_, dec).mean(1)
    lay = io.fields(run, "layers_ann", ann[-1])
    psi = residual_moc(lay["LaVH1TH"], dxg)
    psi_upper = psi.max(axis=0)
    mld = np.array([io.fields(run, "surf_mon", it)["MXLDEPTH"] for it in mon[-12 * last_years:]])
    clim = mld.reshape(-1, 12, *mld.shape[1:]).mean(0) if mld.shape[0] % 12 == 0 else mld
    mld_R = clim[:, R].mean(axis=1)
    return {"years": float(years[-1]), "v6a": drift.v6a(yd, fd, dd, run_years=years[-1]), "moc_R_sv": moc_out_of_R(psi, j_edge),
            "moc_max_sv": float(psi.max()), "upwelling_fraction_R": upwelling_fraction(psi_upper, j_edge),
            "mld_winter_R_m": float(mld_R.max()), "mld_month_max": int(np.argmax(mld_R)) + 1,
            "flux_last_mol_m2_yr": float(fl_ann[-1]), "deep_dic_last_mol_m3": float(deep_dic[-1])}
