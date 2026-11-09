#!/opt/anaconda3/envs/pyocean/bin/python
"""Evaluate the pre-registered smoke criteria S1-S5 of notes/pd1_config.md sect. 6 on one run.

Usage::

    pd1_smoke_check.py <run_dir> [--json OUT]

S1 (N3) require_ok; S2 (N2) meltwater-tracer budget < 1 %; S3 T, S inside the initial
range +- 0.1; S4 shelf-mean melt > 0; S5 MITgcm melt vs theory.melt.three_equation < 1 %
where m > 1 m/yr, on each snapshot. Plus the informative hFacC and timing checks.

Melt sign (notes/pd1_config.md sect. 1.3): SHIfwFlx [kg m^-2 s^-1] is negative for
melting, m [m/yr] = -SHIfwFlx / rho_fw * 365 * 86400.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pfcommon.mitgcm_io import read_bin, read_mds  # noqa: E402
from pfcommon.runctl import require_ok, run_status  # noqa: E402
from exp2_isomip_cavity.gendata import geometry as geo  # noqa: E402
from exp2_isomip_cavity.gendata import make_case as mc  # noqa: E402
from exp2_isomip_cavity.theory import melt as ml  # noqa: E402

RHO_FW = ml.RHO_FW
TOL_T_S = 0.1          # S3 margin


def _grid(rd: Path, name: str, shape) -> np.ndarray:
    """Model grid file <name>.data (big-endian float64, written with useSingleCpuIO)."""
    return np.fromfile(rd / f"{name}.data", dtype=">f8").reshape(shape).astype(float)


def _fields(prefix: Path, itr=None) -> tuple[dict, dict]:
    arr, meta = read_mds(prefix, itr)
    names = meta["fldList"]
    if arr.shape[0] != len(names):
        raise ValueError(f"{prefix}: {arr.shape[0]} records for fields {names}")
    return {n: arr[k] for k, n in enumerate(names)}, meta


def boundary_layer_mean(fld: np.ndarray, hfac: np.ndarray, ktop: np.ndarray, dz: float,
                        clip0: bool = False) -> np.ndarray:
    """Losch (2008) boundary-layer average below the ice, as shelfice_thermodynamics.F:341-389
    (original): top wet cell K plus drKp1 = min(dz (1 - hFacC(K)), dz hFacC(K+1)) of K+1.
    ``clip0`` applies MAX(S, 0) as for salinity. Columns with ktop < 0 give NaN."""
    nz, ny, nx = fld.shape
    out = np.full((ny, nx), np.nan)
    f = np.maximum(fld, 0.0) if clip0 else fld
    for j, i in zip(*np.nonzero(ktop >= 0)):
        k = ktop[j, i]
        top = f[k, j, i]
        if k < nz - 1:
            dr = max(min(dz * (1.0 - hfac[k, j, i]), dz * hfac[k + 1, j, i]), 0.0)
            out[j, i] = (top * dz * hfac[k, j, i] + f[k + 1, j, i] * dr) / (dz * hfac[k, j, i] + dr)
        else:
            out[j, i] = top
    return out


def evaluate(rd: Path) -> dict:
    rd = Path(rd)
    res: dict = {"run_dir": str(rd)}
    # ---- S1 (N3)
    try:
        require_ok(rd)
        res["S1"] = {"pass": True}
    except RuntimeError as e:
        res["S1"] = {"pass": False, "reason": str(e)}
        return res
    st = run_status(rd)
    meta = json.loads((rd / "case.json").read_text())
    p = meta["params"]
    nx, ny, nz, dz = p["nx"], p["ny"], p["nz"], float(p["dz"])
    nsteps, dt = int(p["nTimeSteps"]), float(p["deltaT"])
    T_run = nsteps * dt
    res["timing"] = {"wall_s": st["wall_s"], "n_steps": nsteps,
                     "s_per_step": st["wall_s"] / nsteps,
                     "two_year_h": st["wall_s"] / nsteps * (63072000.0 / dt) / 3600.0}

    hfac = _grid(rd, "hFacC", (nz, ny, nx))
    rA = _grid(rd, "RAC", (ny, nx))
    drF = _grid(rd, "DRF", (nz,))
    bathy = read_bin(rd / "bathy.bin", (ny, nx))
    draft = read_bin(rd / "icetopo.bin", (ny, nx))
    wet3 = hfac > 0
    wet = wet3.any(axis=0)
    shelf = wet & (draft < 0)
    ktop = np.where(shelf, np.argmax(wet3, axis=0), -1)

    # ---- informative: model hFacC vs gendata's reproduction of ini_masks_etc.F
    h_py = geo.mitgcm_hfac(bathy, draft, np.full(nz, dz), p["hFacMin"])[0]
    res["hfac_check"] = {"max_abs_diff": float(np.abs(h_py - hfac).max())}

    # ---- S2 (N2): tracer budget
    shi_mean, _ = _fields(rd / "shiDiag")
    dyn_mean, _ = _fields(rd / "dynDiag")
    fw_mean = shi_mean["SHIfwFlx"][-1] if shi_mean["SHIfwFlx"].ndim == 3 else shi_mean["SHIfwFlx"]
    trac_mean = dyn_mean["TRAC01"][-1] if dyn_mean["TRAC01"].ndim == 4 else dyn_mean["TRAC01"]
    it_end = int(p["nIter0"]) + nsteps
    C_end = read_mds(rd / "PTRACER01", it_end)[0][0]
    eta_end = read_mds(rd / "Eta", it_end)[0][0]
    vol0 = hfac * drF[:, None, None] * rA[None]
    hcol = (hfac * drF[:, None, None]).sum(axis=0)
    rstar = np.where(hcol > 0, (eta_end + hcol) / np.where(hcol > 0, hcol, 1.0), 1.0)
    inv = float((C_end * vol0 * rstar[None]).sum())
    inv_norstar = float((C_end * vol0).sum())
    src = float((-fw_mean * rA).sum() * T_run / RHO_FW)
    mask = read_bin(rd / "rbcs_mask.bin", (nz, ny, nx))
    sink = float((mask / p["tauRelax"] * trac_mean * vol0).sum() * T_run)
    err = abs(inv - (src - sink)) / src
    res["S2"] = {"pass": bool(err < 0.01), "rel_err": err, "inventory_m3": inv,
                 "source_m3": src, "rbcs_sink_m3": sink,
                 "rel_err_without_rstar (informative)": abs(inv_norstar - (src - sink)) / src}

    # ---- S3: T, S range
    # range of the initial profile (analytic, 0-720 m: for the linear / ramp profiles the
    # extremes are the end points; notes/pd1_config.md sect. 6 gives T in [-2.0, 1.1],
    # S in [33.7, 34.8] for Ocean0)
    zz = np.linspace(0.0, -nz * dz, 7201)
    Tp, Sp, *_ = mc._profile(p["init_profile"], zz, p)
    lo_T, hi_T = Tp.min() - TOL_T_S, Tp.max() + TOL_T_S
    lo_S, hi_S = Sp.min() - TOL_T_S, Sp.max() + TOL_T_S
    sets = {"final": (read_mds(rd / "T", it_end)[0][0], read_mds(rd / "S", it_end)[0][0]),
            "mean30d": (dyn_mean["THETA"].reshape(-1, nz, ny, nx)[-1],
                        dyn_mean["SALT"].reshape(-1, nz, ny, nx)[-1])}
    snaps = sorted((rd).glob("snapDiag.*.meta"))
    snap_fields = []
    for f in snaps:
        itr = int(f.name.split(".")[1])
        fl, _ = _fields(rd / "snapDiag", itr)          # 2-D list
        fl3, _ = _fields(rd / "snap3D", itr)           # 3-D list, same time step
        fl.update(fl3)
        snap_fields.append((itr, fl))
        sets[f"snap_{itr}"] = (fl["THETA"], fl["SALT"])
    s3 = {}
    ok3 = True
    for k, (Tf, Sf) in sets.items():
        tmin, tmax = float(Tf[wet3].min()), float(Tf[wet3].max())
        smin, smax = float(Sf[wet3].min()), float(Sf[wet3].max())
        good = lo_T <= tmin and tmax <= hi_T and lo_S <= smin and smax <= hi_S
        ok3 &= good
        s3[k] = {"T": [tmin, tmax], "S": [smin, smax], "pass": bool(good)}
    res["S3"] = {"pass": bool(ok3), "T_bounds": [lo_T, hi_T], "S_bounds": [lo_S, hi_S], **s3}

    # ---- S4: shelf-mean melt > 0
    m_mean = ml.melt_m_per_yr(fw_mean)
    mbar = float((m_mean * rA)[shelf].sum() / rA[shelf].sum())
    res["S4"] = {"pass": bool(mbar > 0), "shelf_mean_melt_m_per_yr": mbar,
                 "n_shelf_columns": int(shelf.sum()),
                 "deep_mean_melt_m_per_yr (z_draft < -300 m, informative)":
                     float((m_mean * rA)[shelf & (draft < -300)].sum() / rA[shelf & (draft < -300)].sum())}

    # ---- S5: MITgcm melt vs three_equation on snapshots
    mass_file = sorted(rd.glob("shelficemassinit*.data"))
    mass = np.fromfile(mass_file[0], dtype=">f8").reshape(ny, nx).astype(float)
    p_pa = mass * p["gravity"]                     # pLoc [dbar] = mass g SItodBar, line 289
    s5 = {}
    ok5 = bool(snap_fields)
    for itr, fl in snap_fields:
        Tbl = boundary_layer_mean(fl["THETA"], hfac, ktop, dz)
        Sbl = boundary_layer_mean(fl["SALT"], hfac, ktop, dz, clip0=True)
        us = fl["SHIuStar"]
        m_model = ml.melt_m_per_yr(fl["SHIfwFlx"])
        sel = shelf & (m_model > 1.0)
        mt = np.full((ny, nx), np.nan)
        mt[sel] = ml.three_equation(Tbl[sel], Sbl[sel], p_pa[sel], us[sel], p["gammaT"])[0] \
            / RHO_FW * ml.SECONDS_PER_YEAR
        rel = np.abs(mt[sel] - m_model[sel]) / np.abs(m_model[sel])
        good = bool(sel.any() and rel.max() < 0.01)
        ok5 &= good
        s5[f"snap_{itr}"] = {"n_cells": int(sel.sum()), "max_rel_err": float(rel.max()) if sel.any() else None,
                             "median_rel_err": float(np.median(rel)) if sel.any() else None,
                             "model_melt_range_m_per_yr": [float(m_model[shelf].min()), float(m_model[shelf].max())],
                             "pass": good}
    res["S5"] = {"pass": ok5, "n_snapshots": len(snap_fields), **s5}
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)
    res = evaluate(a.run_dir)
    txt = json.dumps(res, indent=2, default=float)
    print(txt)
    if a.json:
        a.json.write_text(txt)
    return 0 if all(res.get(k, {}).get("pass") for k in ("S1", "S2", "S3", "S4", "S5")) else 1


if __name__ == "__main__":
    sys.exit(main())
