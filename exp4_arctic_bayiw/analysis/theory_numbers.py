"""Linear internal-wave theory numbers for the PD0 of Exp. 4 (Task 2).

For the representative Van Mijenfjorden profiles (summer, late winter), the thesis bay with its
"realistic" profile and its f, and the thesis bay at 78 N with the summer VMF profile (bridge R2),
this script computes c_1-c_3, the internal Rossby radius, the Burger number, the internal Kelvin
circuit period, the Merian seiche and the inertial period, and the spectral bands proposed for the
PD0. It writes figures/theory.json, figures/profiles_vmf.json (the profiles on the model grid) and
the control figure figures/profiles_vmf.{png,pdf}.

Geometry of the inner fjord: model grid of Task 1b (gendata/bathy.py, ``inner_fjord_mask``).
Geometry of the thesis bay: the thesis bathyFile and delX/delYFile (read-only thesis tree).

Usage: theory_numbers.py [--out figures/theory.json] [--data DIR] [--thesis-input DIR] [--no-figure]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import contourpy  # noqa: E402
import gsw  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import dijkstra  # noqa: E402

from exp4_arctic_bayiw.gendata import profiles as prof  # noqa: E402
from exp4_arctic_bayiw.gendata import vgrid  # noqa: E402
from exp4_arctic_bayiw.theory import modes  # noqa: E402

PKG = Path(__file__).resolve().parents[1]
CASES = PKG / "cases.yaml"
FIG = PKG / "figures"
THESIS_INPUT = Path("/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/input")
THESIS_DATA = PKG / "thesis_copy" / "BayIW_Rectan" / "run_expand" / "data.stage1"
G = modes.G
LAT_FJORD = 78.0          # plan: f-plane at 78 N, f = 2 Omega sin 78 = 1.4266e-4 s^-1
RECORD_H = 240.0          # spec 3.4: free decay to ~240 h (PD0 decision 2: 280 h if PD1 allows)
RECORDS_H = (240.0, 280.0)
DT_OUT_S = 1200.0         # plan: records every 1200 s
N_MODES = 4


# ----------------------------------------------------------------------------- vertical structure
def modes_column(delR, T, S, H: float, eos=prof.EOS_THESIS):
    """Nodes (z, m, positive up, 0 to -H) and N^2 (s^-2) for theory.modes.vertical_modes.

    The nodes are the model interfaces shallower than H plus z = -H. N^2 at an interior interface
    between cell centres k-1 and k uses the linear EOS of the model (model/src/find_rho.F:104-107):
    N^2 = g (tAlpha dT/dz - sBeta dS/dz), z positive up (Gill 1982 eq. 3.7.4 for a linear EOS
    [verificar]); at z = 0 and z = -H N^2 is copied from the nearest interior interface.
    """
    f = vgrid.faces(delR)
    zc = -0.5 * (f[1:] + f[:-1])
    T, S = np.asarray(T, float), np.asarray(S, float)
    n2f = G * (eos["tAlpha"] * np.diff(T) - eos["sBeta"] * np.diff(S)) / np.diff(zc)   # at f[1:-1]
    zf = -f[1:-1]
    keep = zf > -H + 1e-6
    z = np.concatenate([[0.0], zf[keep], [-H]])
    n2_in = n2f[keep]
    if n2_in.size == 0:
        raise ValueError("H is shallower than the first interface")
    n2 = np.concatenate([[n2_in[0]], n2_in, [n2_in[-1]]])
    return z, n2


def modes_column_teos10(delR, T, S, H: float, lat: float, lon: float):
    """As ``modes_column`` but with TEOS-10 N^2 (gsw.Nsquared; IOC et al. 2010), for the EOS check.

    ``lon`` is the mean longitude of the casts (absolute salinity anomaly). Negative N^2 values
    (statically unstable under TEOS-10 although stable under the model EOS) are clipped to 0 and
    counted: returns (z, N2, n_negative).
    """
    f = vgrid.faces(delR)
    zc = -0.5 * (f[1:] + f[:-1])
    p = gsw.p_from_z(zc, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_t(SA, T, p)
    n2f, _ = gsw.Nsquared(SA, CT, p, lat)
    zf = -f[1:-1]
    keep = zf > -H + 1e-6
    n2_in = n2f[keep]
    nneg = int(np.sum(n2_in < 0))
    n2_in = np.maximum(n2_in, 0.0)
    return (np.concatenate([[0.0], zf[keep], [-H]]), np.concatenate([[n2_in[0]], n2_in, [n2_in[-1]]]),
            nneg)


# ----------------------------------------------------------------------------- basin geometry
def mask_geometry(basin, wet, mouth, dx: float) -> dict:
    """Plan-view scales of a basin given as a mask on a uniform grid (cell size dx, m).

    - ``staircase_coast_m``: cell edges between basin and land (not wet) cells (an upper bound of
      the coast length on a staircase grid).
    - ``contour_*``: the 0.5 contour of the basin mask (marching squares, contourpy) on the cell
      centres, which cuts the staircase corners. Each segment is "coast" if the non-basin cells
      next to it are all land and "mouth" if one of them is wet (the sound openings).
      ``contour_outer_m`` is the longest closed loop (the basin circuit; islands are inner loops).
    - ``along_axis_length_m``: longest water path (8-connected cells, Dijkstra) from the basin
      cells next to the mouth cells to any basin cell, plus one cell (mouth edge to far wall).
    - ``axis_angle_deg``: principal axis of the basin cell centres, degrees counter-clockwise from
      east, in (-90, 90].
    """
    basin, wet, mouth = (np.asarray(a, bool) for a in (basin, wet, mouth))
    ny, nx = basin.shape
    pad = np.pad(basin, 1)
    padwet = np.pad(wet, 1)
    stair = 0
    for ax_, sh in ((0, 1), (0, -1), (1, 1), (1, -1)):
        nb = np.roll(pad, sh, axis=ax_)
        nbw = np.roll(padwet, sh, axis=ax_)
        stair += int(np.sum(pad & ~nb & ~nbw))
    # contour of the basin mask (index coordinates, padded so loops close)
    gen = contourpy.contour_generator(z=pad.astype(float), name="serial")
    lines = gen.lines(0.5)
    total = coast = 0.0
    loops = []
    for ln in lines:
        seg = np.diff(ln, axis=0)
        lens = np.hypot(seg[:, 0], seg[:, 1]) * dx
        mid = 0.5 * (ln[1:] + ln[:-1])
        is_mouth = np.zeros(len(mid), bool)
        for k, (x, y) in enumerate(mid):
            i0, j0 = int(np.floor(x)), int(np.floor(y))
            for jj in (j0, j0 + 1):
                for ii in (i0, i0 + 1):
                    if 0 <= jj < ny + 2 and 0 <= ii < nx + 2 and not pad[jj, ii] and padwet[jj, ii]:
                        is_mouth[k] = True
        total += lens.sum()
        coast += lens[~is_mouth].sum()
        loops.append((lens.sum(), lens[~is_mouth].sum()))
    outer = max(loops) if loops else (0.0, 0.0)
    # along-axis water path
    idx = -np.ones(basin.shape, int)
    J, I = np.nonzero(basin)
    idx[J, I] = np.arange(J.size)
    rows, cols, w = [], [], []
    for dj, di, d in ((0, 1, dx), (1, 0, dx), (1, 1, dx * np.sqrt(2)), (1, -1, dx * np.sqrt(2))):
        j2, i2 = J + dj, I + di
        ok = (j2 >= 0) & (j2 < ny) & (i2 >= 0) & (i2 < nx)
        ok[ok] &= basin[j2[ok], i2[ok]]
        rows += list(idx[J[ok], I[ok]])
        cols += list(idx[j2[ok], i2[ok]])
        w += [d] * int(ok.sum())
    graph = coo_matrix((w, (rows, cols)), shape=(J.size, J.size)).tocsr()
    near_mouth = np.zeros_like(basin)
    for sh, ax_ in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
        near_mouth |= np.roll(mouth, sh, axis=ax_)
    src = idx[basin & near_mouth]
    if src.size == 0:
        raise ValueError("no basin cell touches the mouth")
    dist = dijkstra(graph, directed=False, indices=src, min_only=True)
    L = float(np.max(dist[np.isfinite(dist)]) + dx)
    # principal axis
    X = np.c_[I * dx, J * dx].astype(float)
    X -= X.mean(axis=0)
    w_, v_ = np.linalg.eigh(X.T @ X)
    ax_vec = v_[:, np.argmax(w_)]
    ang = float(np.degrees(np.arctan2(ax_vec[1], ax_vec[0])))
    ang = ang - 180.0 if ang > 90 else (ang + 180.0 if ang <= -90 else ang)
    area = float(basin.sum() * dx * dx)
    return {"area_m2": area, "staircase_coast_m": float(stair * dx), "contour_total_m": float(total),
            "contour_coast_m": float(coast), "contour_outer_m": float(outer[0]),
            "contour_outer_coast_m": float(outer[1]), "along_axis_length_m": L,
            "mean_width_m": area / L, "axis_angle_deg": ang}


def fjord_geometry(cfg: dict, data_dir: Path) -> dict:
    """Inner Van Mijenfjorden on the Task 1b model grid (bathy.build + bathy.inner_fjord_mask)."""
    from exp4_arctic_bayiw.gendata import bathy

    _, grid = bathy.build(cfg, data_dir)
    inner, outer, disks = bathy.inner_fjord_mask(grid)
    both = np.logical_or.reduce(list(disks.values()))
    dx = grid["dx_m"]
    geo = mask_geometry(inner, grid["wet"], both & grid["wet"], dx)
    d = grid["depth"][inner]
    geo.update({"max_depth_m": float(d.max()), "mean_depth_m": float(d.mean()),
                "volume_m3": float(d.sum() * dx * dx), "cells": int(inner.sum()), "dx_m": dx,
                "method": "inner_fjord_mask of gendata/bathy.py (Task 1b); perimeters from "
                          "mask_geometry (0.5 contour of the mask, marching squares)"})
    # control of the profile polygon against the inner-fjord mask
    import shapely
    from shapely.geometry import Polygon

    poly = Polygon(cfg["profiles"]["polygon_utm"])
    XC, YC = np.meshgrid(grid["xc"], grid["yc"])
    ins = shapely.contains_xy(poly, XC, YC)
    geo["polygon_control"] = {"inner_cells_inside": int((inner & ins).sum()),
                              "inner_cells": int(inner.sum()),
                              "other_sea_cells_inside_outside_disks": int((grid["wet"] & ~inner & ~both & ins).sum())}
    geo["_grid"] = grid
    geo["_inner"] = inner
    return geo


def parse_namelist_list(text: str, name: str) -> list[float]:
    m = re.search(rf"^\s*{name}\s*=(.*?)(?=^\s*[A-Za-z_]\w*\s*=|^\s*&)", text, re.S | re.M)
    if not m:
        raise KeyError(name)
    return [float(v) for v in re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?", m.group(1).replace("D", "E"))]


def thesis_bay(input_dir: Path, data_file: Path = THESIS_DATA) -> dict:
    """Thesis bay geometry, f0, delR and the 'realistic' (Feb) T/S column, read-only.

    The bay is the set of cells with 0 < depth < 999 m (the shelf is 1000 m, the land 0). The grid
    is non-uniform (delX/delYFile), so L, W are sums of the cell sizes over the bay rows/columns;
    the bay is a rectangle closed on three sides: coast = 2L + W, closed loop = 2L + 2W.
    """
    txt = data_file.read_text()
    f0 = parse_namelist_list(txt, "f0")[0]
    delR = np.array(parse_namelist_list(txt, "delR"))
    tAlpha = parse_namelist_list(txt, "tAlpha")[0]
    sBeta = parse_namelist_list(txt, "sBeta")[0]
    dx = np.fromfile(input_dir / "bahia_01_expand_dx.bin", ">f8")
    dy = np.fromfile(input_dir / "bahia_01_expand_dy.bin", ">f8")
    b = -np.fromfile(input_dir / "bahia_01_expand_bat.bin", ">f8").reshape(dy.size, dx.size)
    bay = (b > 0) & (b < 999.0)
    J, I = np.nonzero(bay)
    L = float(dy[np.unique(J)].sum())
    W = float(dx[np.unique(I)].sum())
    area = float((np.outer(dy, dx) * bay).sum())
    vol = float((np.outer(dy, dx) * b * bay).sum())
    nz, ny, nx = delR.size, dy.size, dx.size
    n = nz * ny * nx
    T = np.fromfile(input_dir / "feb_temp_50zlev_560x352.bin", ">f8", count=n).reshape(nz, ny, nx)[:, 0, 0]
    S = np.fromfile(input_dir / "feb_salt_50zlev_560x352.bin", ">f8", count=n).reshape(nz, ny, nx)[:, 0, 0]
    return {"f0_s-1": f0, "delR": delR.tolist(), "eos": {"tAlpha": tAlpha, "sBeta": sBeta},
            "length_m": L, "width_m": W, "max_depth_m": float(b[bay].max()), "mean_depth_m": vol / area,
            "area_m2": area, "coast_m": 2 * L + W, "closed_loop_m": 2 * L + 2 * W,
            "dx_bay_m": [float(dx[I].min()), float(dx[I].max())], "dy_bay_m": [float(dy[J].min()), float(dy[J].max())],
            "T": T.tolist(), "S": S.tolist(),
            "source": f"{input_dir}/bahia_01_expand_*.bin, feb_*_50zlev_560x352.bin (column 0,0); {data_file.name}"}


# ----------------------------------------------------------------------------- profiles
def _median_ignore_warnings(fn, *a, **k):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return fn(*a, **k)


def vmf_profiles(cfg: dict, data_dir: Path, eos: dict | None = None) -> dict:
    """Selection, QC and representative profiles per season. ``eos`` (default: the fitted EOS of
    cases.yaml, ``configs.fjord_vmf.eos``) is the one used for the static-stability adjustment."""
    pc = cfg["profiles"]
    f = vgrid.faces(cfg["delR"])
    zc = -0.5 * (f[1:] + f[:-1])
    allp = prof.load_unis_hd(data_dir / pc["source"])
    eos = cfg["eos"] if eos is None else eos
    out = {"n_files": len(allp), "z_c_m": zc.tolist(), "seasons": {}}
    for season, s in pc["seasons"].items():
        sel = prof.select(allp, pc["polygon_lonlat"], set(s["months"]), tuple(s["years"]))
        sel_all_years = prof.select(allp, pc["polygon_lonlat"], set(s["months"]), None)
        reasons = defaultdict(int)
        ok = []
        for p in sel:
            q, why = prof.qc_check(p, **pc["qc"])
            if q is not None:
                ok.append(q)
            else:
                reasons[why] += 1
        if len(ok) < pc["min_profiles"]:
            raise RuntimeError(f"{season}: {len(ok)} profiles after QC < {pc['min_profiles']}")
        info = {}
        T, S = _median_ignore_warnings(prof.representative, ok, zc, eos=eos,
                                       min_count=pc["min_count_per_level"], info=info)
        # sensitivity 1: each year weighs the same (median of the annual representative profiles)
        by_year = defaultdict(list)
        for q in ok:
            by_year[q["time"].year].append(q)
        Ty, Sy = [], []
        for yr, ps in sorted(by_year.items()):
            t_, s_ = _median_ignore_warnings(prof.representative, ps, zc, eos=eos, min_count=1)
            Ty.append(t_)
            Sy.append(s_)
        Tyb = np.median(Ty, axis=0)
        Syb, _ = prof.stabilize(zc, Tyb, np.median(Sy, axis=0), eos)
        # sensitivity 2: QC inversion tolerance
        tol_sens = {}
        for tol in (1e-3, 2e-3, 4e-3, 1e-2):
            kw = dict(pc["qc"], inversion_tol_kg_m3=tol)
            tol_sens[f"{tol:g}"] = int(sum(prof.qc(p, **kw) is not None for p in sel))
        lon_mean = float(np.mean([q["lon"] for q in ok]))
        lat_mean = float(np.mean([q["lat"] for q in ok]))
        years = sorted({q["time"].year for q in ok})
        out["seasons"][season] = {
            "months": s["months"], "years_window": s["years"], "n_selected": len(sel),
            "n_selected_all_years": len(sel_all_years), "n_qc_pass": len(ok),
            "qc_rejections": dict(reasons), "qc_pass_vs_inversion_tol": tol_sens,
            "years_present": years, "n_per_year": {str(k): len(v) for k, v in sorted(by_year.items())},
            "n_per_month": {str(m): sum(q["time"].month == m for q in ok) for m in s["months"]},
            "instrument_types": sorted({str(q.get("instrument_type")) for q in ok}),
            "lon_mean": lon_mean, "lat_mean": lat_mean,
            "T": T.tolist(), "S": S.tolist(), "info": info,
            "T_year_balanced": Tyb.tolist(), "S_year_balanced": Syb.tolist(),
            "_casts": ok,
        }
    return out


def fit_eos(cfg: dict, data_dir: Path, zmax_m: float = 120.0) -> dict:
    """Linear EOS for the fjord and bridge runs (PD0 decision 3): ``prof.fit_linear_eos`` on the
    summer and winter representative profiles (levels above ``zmax_m``, both seasons pooled).
    The representative profiles used for the fit are stabilised with the thesis EOS (the
    adjustment changes S by < 0.004); test_real_fitted_eos_in_cases_yaml checks that re-fitting
    on the final profiles changes the coefficients by < 1 %."""
    pr = vmf_profiles(cfg, data_dir, eos=prof.EOS_THESIS)
    zc = np.array(pr["z_c_m"])
    k = zc > -zmax_m
    T = np.concatenate([np.array(pr["seasons"][s]["T"])[k] for s in ("summer", "winter")])
    S = np.concatenate([np.array(pr["seasons"][s]["S"])[k] for s in ("summer", "winter")])
    lon = float(np.mean([pr["seasons"][s]["lon_mean"] for s in ("summer", "winter")]))
    lat = float(np.mean([pr["seasons"][s]["lat_mean"] for s in ("summer", "winter")]))
    fit = prof.fit_linear_eos(T, S, lon=lon, lat=lat)
    fit.update({"fit_lon": lon, "fit_lat": lat, "zmax_m": zmax_m})
    return fit


# ----------------------------------------------------------------------------- numbers
def case_numbers(z, n2, H_modes, f, L, perimeter, H_merian, W=None, n_modes=N_MODES) -> dict:
    """Scales of one configuration (theory/modes.py). Units: SI, periods in hours.

    With the basin width W (m) it adds the width Burger number R_n / W and the lowest
    cross-basin Poincare frequency of each vertical mode, omega_n = sqrt(f^2 + c_n^2 (pi/W)^2)
    (channel of width W, one cross-channel half wavelength; Gill 1982 sec. 10.3, eq. 10.3.5, with c_n
    for (gH)^(1/2); thesis library Gill_2016 .md l. 14593-14600),
    given as omega_n / f.
    """
    c, _ = modes.vertical_modes(z, n2, H_modes, n_modes=n_modes)
    c = [float(x) for x in c]
    out = {
        "f_s-1": float(f),
        "inertial_period_h": float(2 * np.pi / abs(f) / 3600),                     # T_i = 2 pi / f
        "H_modes_m": float(H_modes), "L_m": float(L), "perimeter_m": float(perimeter),
        "H_merian_m": float(H_merian),
        "c_n": c,                                                                   # m s^-1
        "rossby_radius_m": [modes.rossby_radius(x, f) for x in c],                  # R_n = c_n / f
        "burger": [modes.burger(x, f, L) for x in c],                               # Bu = R / L
        "kelvin_circuit_period_h": [modes.kelvin_circuit_period(perimeter, x) / 3600 for x in c],
        "c_barotropic_m_s-1": float(np.sqrt(G * H_merian)),
        "merian_closed_open_h": [modes.merian_period(L, H_merian, n, True) / 3600 for n in (1, 2, 3)],
        "merian_closed_closed_h": [modes.merian_period(L, H_merian, n, False) / 3600 for n in (1, 2, 3)],
        "N_max_s-1": float(np.sqrt(np.max(n2))),
    }
    if W is not None:
        out["W_m"] = float(W)
        out["burger_width"] = [modes.rossby_radius(x, f) / W for x in c]
        out["poincare_cross_1_over_f"] = [float(np.sqrt(f**2 + (x * np.pi / W) ** 2) / abs(f)) for x in c]
    return out


def spectral_bands(f_cph: float, record_h: float, dt_out_s: float, merian_lowest_cph: float,
                   ni=(0.85, 1.20), seiche_factor: float = 0.5) -> dict:
    """Non-overlapping frequency bands (cph), relative to f and to the lowest Merian frequency.

    low [0, 2 df): unresolved (df = 1/record); kelvin_subinertial [2 df, ni0 f); near_inertial
    [ni0 f, ni1 f); superinertial [ni1 f, seiche_factor * merian_lowest); seiche [.., Nyquist].
    Raises ValueError if the edges are not increasing or a resolved band has < 2 frequency bins.
    """
    df = 1.0 / record_h
    nyq = 1.0 / (2.0 * dt_out_s / 3600.0)
    edges = [0.0, 2 * df, ni[0] * f_cph, ni[1] * f_cph, seiche_factor * merian_lowest_cph, nyq]
    names = ["low", "kelvin_subinertial", "near_inertial", "superinertial", "seiche"]
    if np.any(np.diff(edges) <= 0):
        raise ValueError(f"band edges not increasing: {edges}")
    bands = []
    for nm, lo, hi in zip(names, edges[:-1], edges[1:]):
        nb = (hi - lo) / df
        if nm != "low" and nb < 2:
            raise ValueError(f"band {nm} has {nb:.2f} < 2 bins of df = {df:.4g} cph")
        bands.append({"name": nm, "lo_cph": lo, "hi_cph": hi, "n_bins": nb,
                      "period_h": [1 / hi if hi > 0 else None, 1 / lo if lo > 0 else None]})
    return {"df_cph": df, "nyquist_cph": nyq, "record_h": record_h, "dt_out_s": dt_out_s,
            "f_cph": f_cph, "ni_factors": list(ni), "seiche_factor": seiche_factor,
            "merian_lowest_cph": merian_lowest_cph, "bands": bands}


def band_of(freq_cph: float, bands: dict) -> str:
    for b in bands["bands"]:
        if b["lo_cph"] <= freq_cph < b["hi_cph"]:
            return b["name"]
    return "above_nyquist"


# ----------------------------------------------------------------------------- figure
def figure(prof_out: dict, cases: dict, thesis: dict | None, stem: Path):
    import matplotlib.pyplot as plt
    from pfcommon import style

    style.apply()
    zc = np.array(prof_out["z_c_m"])
    fig, axs = plt.subplots(1, 4, figsize=(7.0, 4.6), sharey=False)
    cols = {"summer": style.OKABE_ITO[1], "winter": style.OKABE_ITO[0]}
    for season, s in prof_out["seasons"].items():
        c = cols[season]
        for q in s["_casts"]:
            axs[0].plot(q["T"], q["z"], color=c, lw=0.2, alpha=0.08, rasterized=True)
            axs[1].plot(q["S"], q["z"], color=c, lw=0.2, alpha=0.08, rasterized=True)
        T, S = np.array(s["T"]), np.array(s["S"])
        lab = f"{season} (n = {s['n_qc_pass']})"
        axs[0].plot(T, zc, color=c, lw=1.6, label=lab)
        axs[1].plot(S, zc, color=c, lw=1.6)
        axs[0].plot(s["T_year_balanced"], zc, color=c, lw=0.9, ls="--")
        axs[1].plot(s["S_year_balanced"], zc, color=c, lw=0.9, ls="--")
        zl = s["info"]["last_data_z_m"]
        for a in axs[:3]:
            a.axhline(zl, color=c, lw=0.6, ls=":")
        cs = cases[f"vmf_{season}"]
        z, n2 = cs["_z"], cs["_n2"]
        axs[2].plot(np.sqrt(np.maximum(n2, 0)) * 3600 / (2 * np.pi), z, color=c, lw=1.2)
        z2, n22 = cs["_z_teos"], cs["_n2_teos"]
        axs[2].plot(np.sqrt(np.maximum(n22, 0)) * 3600 / (2 * np.pi), z2, color=c, lw=0.8, ls="--")
        for k in range(3):
            phi = cs["_phi"][k]
            axs[3].plot(phi / np.max(np.abs(phi)), z, color=c, lw=1.0, ls=["-", "--", ":"][k])
    axs[0].set_xlabel("T (°C)")
    axs[1].set_xlabel("S (PSU)")
    axs[2].set_xlabel("N (cph)")
    axs[0].set_ylabel("z (m)")
    for a in axs[:2]:
        a.set_ylim(-125, 0)
    axs[2].set_ylim(-125, 0)
    axs[3].set_ylim(-125, 0)
    axs[3].set_xlabel("w-mode (norm.)")
    axs[1].set_xlim(24, 35.5)
    fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.03))
    axs[2].text(0.97, 0.03, "solid: fitted linear EOS\ndashed: TEOS-10", transform=axs[2].transAxes,
                ha="right", va="bottom", fontsize=7)
    axs[3].text(0.03, 0.97, "solid n = 1\ndashed n = 2\ndotted n = 3", transform=axs[3].transAxes,
                ha="left", va="top", fontsize=6.5, bbox=dict(fc="w", ec="none", alpha=0.8))
    axs[1].text(0.03, 0.97, "thin: QC casts\nthick: median\ndashed: year-balanced\ndotted: last datum",
                transform=axs[1].transAxes, ha="left", va="top", fontsize=6.5)
    for a, t in zip(axs, "abcd"):
        a.set_title(f"({t})", loc="left")
    fig.suptitle("Inner Van Mijenfjorden, UNIS HD (1990–2019): representative profiles", fontsize=9)
    fig.tight_layout()
    return style.save(fig, stem)


# ----------------------------------------------------------------------------- main
def _strip(d):
    if isinstance(d, dict):
        return {k: _strip(v) for k, v in d.items() if not k.startswith("_")}
    if isinstance(d, list):
        return [_strip(v) for v in d]
    if isinstance(d, (np.floating, np.integer)):
        return d.item()
    return d


def run(data_dir: Path, thesis_input: Path) -> tuple[dict, dict, dict | None]:
    cfg = yaml.safe_load(CASES.read_text())["configs"]["fjord_vmf"]
    eos = cfg["eos"]                       # linear EOS fitted to VMF (PD0 decision 3)
    delR = np.array(cfg["delR"])
    geo = fjord_geometry(cfg, data_dir)
    P = geo["contour_outer_m"]
    L = geo["along_axis_length_m"]
    pr_out = vmf_profiles(cfg, data_dir)
    f78 = modes.coriolis(LAT_FJORD)
    cases = {}
    for season in ("summer", "winter"):
        s = pr_out["seasons"][season]
        T, S = np.array(s["T"]), np.array(s["S"])
        for tag, H in (("", geo["max_depth_m"]), ("_Hmean", geo["mean_depth_m"])):
            z, n2 = modes_column(delR, T, S, H, eos)
            c = case_numbers(z, n2, H, f78, L, P, geo["mean_depth_m"], W=geo["mean_width_m"])
            c["eos"] = "fitted_vmf"
            zt, n2t, nneg = modes_column_teos10(delR, T, S, H, LAT_FJORD, s["lon_mean"])
            ct, _ = modes.vertical_modes(zt, n2t, H, n_modes=N_MODES)
            c["c_n_teos10_eos"] = [float(x) for x in ct]
            c["teos10_negative_N2_interfaces_clipped"] = nneg
            # sensitivity: the thesis EOS on the same T/S (the PD0 draft numbers)
            zth, n2th = modes_column(delR, T, S, H, prof.EOS_THESIS)
            cth, _ = modes.vertical_modes(zth, n2th, H, n_modes=N_MODES)
            c["c_n_thesis_eos"] = [float(x) for x in cth]
            zy, n2y = modes_column(delR, s["T_year_balanced"], s["S_year_balanced"], H, eos)
            cy, _ = modes.vertical_modes(zy, n2y, H, n_modes=N_MODES)
            c["c_n_year_balanced"] = [float(x) for x in cy]
            c["profile"] = f"vmf_{season}"
            if tag == "":
                _, phi = modes.vertical_modes(z, n2, H, n_modes=N_MODES)
                c.update({"_z": z, "_n2": n2, "_z_teos": zt, "_n2_teos": n2t, "_phi": phi})
            cases[f"vmf_{season}{tag}"] = c
    thesis = None
    if thesis_input.is_dir():
        thesis = thesis_bay(thesis_input)
        tdel = np.array(thesis["delR"])
        Hb = thesis["max_depth_m"]
        z, n2 = modes_column(tdel, thesis["T"], thesis["S"], Hb, thesis["eos"])
        c = case_numbers(z, n2, Hb, thesis["f0_s-1"], thesis["length_m"], thesis["coast_m"], Hb,
                         W=thesis["width_m"])
        c["kelvin_circuit_period_closed_loop_h"] = [thesis["closed_loop_m"] / x / 3600 for x in c["c_n"]]
        c["profile"] = "thesis realistic (feb_temp/feb_salt, column 0,0)"
        c["eos"] = "thesis"
        # ruling R11: the thesis f0 is the constant of its beta plane (beta = 1e-11, MITgcm default)
        c["inertial_period_note"] = ("f0 of the thesis beta plane (f = f0 + beta yC, "
                                     "model/src/ini_cori.F:77, beta = 1e-11 default): 1-3 % below "
                                     "the local f of the bay (7.05e-5 to 7.16e-5 s^-1 at rows "
                                     "j 324-347, R1 case.json), so this period is 1-3 % longer "
                                     "than the local inertial period (24.4-24.8 h)")
        cases["thesis_bay_thesis_f"] = c
        # bridge R2: VMF summer profile on the thesis levels (below 240 m the last value)
        tf = vgrid.faces(tdel)
        tzc = -0.5 * (tf[1:] + tf[:-1])
        s = pr_out["seasons"]["summer"]
        zc = np.array(pr_out["z_c_m"])
        Tb = np.interp(-tzc, -zc, s["T"])
        Sb = np.interp(-tzc, -zc, s["S"])
        Sb, nadj = prof.stabilize(tzc, Tb, Sb, eos)
        z, n2 = modes_column(tdel, Tb, Sb, Hb, eos)
        c = case_numbers(z, n2, Hb, f78, thesis["length_m"], thesis["coast_m"], Hb, W=thesis["width_m"])
        c["kelvin_circuit_period_closed_loop_h"] = [thesis["closed_loop_m"] / x / 3600 for x in c["c_n"]]
        c["profile"] = "vmf_summer interpolated to the thesis levels"
        c["eos"] = "fitted_vmf"
        Sb_th, _ = prof.stabilize(tzc, Tb, np.interp(-tzc, -zc, s["S"]), prof.EOS_THESIS)
        zth, n2th = modes_column(tdel, Tb, Sb_th, Hb, prof.EOS_THESIS)
        cth, _ = modes.vertical_modes(zth, n2th, Hb, n_modes=N_MODES)
        c["c_n_thesis_eos"] = [float(x) for x in cth]
        c["levels_adjusted_after_interp"] = nadj
        cases["thesis_bay_78N_vmf_summer"] = c
        pr_out["thesis_levels_vmf_summer"] = {"T": Tb.tolist(), "S": Sb.tolist(), "delR": tdel.tolist()}
    # spectral bands: 240 h and 280 h (PD0 decision 2; T_end fixed at the PD1)
    bands = {}
    f78_cph = f78 / (2 * np.pi) * 3600
    fj_merian_low = 1.0 / cases["vmf_summer"]["merian_closed_open_h"][0]
    for rec in RECORDS_H:
        bands[f"fjord_78N_{rec:.0f}h"] = spectral_bands(f78_cph, rec, DT_OUT_S, fj_merian_low)
    bands["fjord_78N_free_decay_180h"] = spectral_bands(f78_cph, 180.0, DT_OUT_S, fj_merian_low)
    if thesis is not None:
        bay_merian_low = 1.0 / cases["thesis_bay_thesis_f"]["merian_closed_open_h"][0]
        f0_cph = thesis["f0_s-1"] / (2 * np.pi) * 3600
        # R1k runs 240 h (spec 3.4); with the 120 h of R1 the near-inertial band is under-resolved
        bands["bay_thesis_f_240h"] = spectral_bands(f0_cph, RECORD_H, 900.0, bay_merian_low)
        try:
            spectral_bands(f0_cph, 120.0, 900.0, bay_merian_low)
            bands["bay_thesis_f_120h_check"] = "resolved"
        except ValueError as e:
            bands["bay_thesis_f_120h_check"] = f"not resolved: {e}"
        for rec in RECORDS_H:
            bands[f"bay_78N_{rec:.0f}h"] = spectral_bands(f78_cph, rec, DT_OUT_S, bay_merian_low)
    # where the theoretical periods fall
    placement = {}
    for key, bkey in (("vmf_summer", "fjord_78N_240h"), ("vmf_winter", "fjord_78N_240h"),
                      ("vmf_summer_Hmean", "fjord_78N_240h"), ("vmf_winter_Hmean", "fjord_78N_240h"),
                      ("thesis_bay_thesis_f", "bay_thesis_f_240h"),
                      ("thesis_bay_78N_vmf_summer", "bay_78N_240h")):
        if key not in cases or bkey not in bands:
            continue
        c, b = cases[key], bands[bkey]
        TK1 = c["kelvin_circuit_period_h"][0]
        placement[key] = {
            "bands": bkey,
            "inertial": band_of(1 / c["inertial_period_h"], b),
            "kelvin_mode1_harmonics": {str(n): {"period_h": TK1 / n, "band": band_of(n / TK1, b)} for n in range(1, 11)},
            "kelvin_circuits_in_record": b["record_h"] / TK1,
            "kelvin_circuits_in_280h": 280.0 / TK1,
            "kelvin_mode1_band_280h": band_of(1 / TK1, bands[bkey.replace("240h", "280h")])
            if bkey.replace("240h", "280h") in bands else None,
            "merian_closed_open_1": band_of(1 / c["merian_closed_open_h"][0], b),
            "merian_closed_closed_1": band_of(1 / c["merian_closed_closed_h"][0], b),
        }
    resolution = {k: {"R1_over_dx_200m": c["rossby_radius_m"][0] / geo["dx_m"],
                      "R1_over_dx_400m": c["rossby_radius_m"][0] / (2 * geo["dx_m"]),
                      "R1_over_dx_100m": c["rossby_radius_m"][0] / (geo["dx_m"] / 2)}
                  for k, c in cases.items() if k.startswith("vmf_")}
    # proposed wind direction: along the principal axis of the inner fjord, down-fjord (blowing
    # towards the mouth, i.e. towards the west end of the axis)
    a = geo["axis_angle_deg"]                    # axis direction, CCW from east, (-90, 90]
    toward = a + 180.0                           # the western end of the axis
    u = np.cos(np.radians(toward))
    v = np.sin(np.radians(toward))
    wind = {"axis_angle_ccw_from_east_deg": a,
            "down_fjord_stress_unit_vector_xy": [float(u), float(v)],
            "down_fjord_from_direction_deg_true_utm": float((90.0 - a) % 360.0),
            "note": "direction 'from' measured clockwise from UTM grid north; meteorological convention"}
    theory = {
        "description": "Linear internal-wave scales for the Exp. 4 PD0 (analysis/theory_numbers.py)",
        "lat_fjord_deg": LAT_FJORD,
        "fjord_geometry": {k: v for k, v in geo.items() if not k.startswith("_")},
        "thesis_bay_geometry": {k: v for k, v in (thesis or {}).items() if k not in ("T", "S", "delR")},
        "eos": {"fitted_vmf": {k: v for k, v in eos.items()}, "thesis": prof.EOS_THESIS,
                "check_on_final_profiles": _eos_check(pr_out, eos)},
        "v2b_range_m_s-1": {
            "formula": "[0.8 c1(H mean), 1.2 c1(H max)], summer (PD0 4)",
            "fitted_vmf": [0.8 * cases["vmf_summer_Hmean"]["c_n"][0], 1.2 * cases["vmf_summer"]["c_n"][0]],
            "thesis_eos": [0.8 * cases["vmf_summer_Hmean"]["c_n_thesis_eos"][0],
                           1.2 * cases["vmf_summer"]["c_n_thesis_eos"][0]]},
        "resolution": resolution,
        "cases": cases, "bands": bands, "placement": placement, "wind": wind,
        "profiles_summary": {k: {kk: vv for kk, vv in v.items()
                                 if kk not in ("T", "S", "T_year_balanced", "S_year_balanced", "_casts", "info")}
                             | {"info": {kk: vv for kk, vv in v["info"].items() if kk != "count_per_level"}}
                             for k, v in pr_out["seasons"].items()},
    }
    return theory, pr_out, thesis


def _eos_check(pr_out: dict, eos: dict, zmax_m: float = 120.0) -> dict:
    """Fitted linear sigma against TEOS-10 sigma0 on the final representative profiles."""
    zc = np.array(pr_out["z_c_m"])
    k = zc > -zmax_m
    T = np.concatenate([np.array(pr_out["seasons"][s]["T"])[k] for s in ("summer", "winter")])
    S = np.concatenate([np.array(pr_out["seasons"][s]["S"])[k] for s in ("summer", "winter")])
    SA = gsw.SA_from_SP(S, 0.0, eos["fit_lon"], eos["fit_lat"])
    sig = gsw.sigma0(SA, gsw.CT_from_t(SA, T, 0.0))
    r = prof.linear_sigma(T, S, eos) - sig
    refit = prof.fit_linear_eos(T, S, lon=eos["fit_lon"], lat=eos["fit_lat"])
    return {"max_abs_err_kg_m3": float(np.max(np.abs(r))), "rms_err_kg_m3": float(np.sqrt(np.mean(r**2))),
            "refit_tAlpha": refit["tAlpha"], "refit_sBeta": refit["sBeta"],
            "refit_rel_change": [refit["tAlpha"] / eos["tAlpha"] - 1, refit["sBeta"] / eos["sBeta"] - 1]}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=FIG / "theory.json")
    ap.add_argument("--data", type=Path, default=PKG / "data_external")
    ap.add_argument("--thesis-input", type=Path, default=THESIS_INPUT)
    ap.add_argument("--no-figure", action="store_true")
    ap.add_argument("--fit-eos", action="store_true",
                    help="only fit the linear EOS (PD0 decision 3) and print it for cases.yaml")
    a = ap.parse_args(argv)
    if a.fit_eos:
        cfg = yaml.safe_load(CASES.read_text())["configs"]["fjord_vmf"]
        print(json.dumps(fit_eos(cfg, a.data), indent=1))
        return
    theory, pr_out, thesis = run(a.data, a.thesis_input)
    if not a.no_figure:
        for p in figure(pr_out, theory["cases"], thesis, FIG / "profiles_vmf"):
            print(p)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(_strip(theory), indent=1) + "\n")
    pj = {"z_c_m": pr_out["z_c_m"], "delR": yaml.safe_load(CASES.read_text())["configs"]["fjord_vmf"]["delR"]}
    for k, s in pr_out["seasons"].items():
        pj[k] = {"T": s["T"], "S": s["S"], "n": s["n_qc_pass"], "years": s["years_present"],
                 "months": s["months"], "info": s["info"]}
    if "thesis_levels_vmf_summer" in pr_out:
        pj["summer_on_thesis_levels"] = pr_out["thesis_levels_vmf_summer"]
    (a.out.parent / "profiles_vmf.json").write_text(json.dumps(_strip(pj), indent=1) + "\n")
    print(a.out)
    for k, c in theory["cases"].items():
        print(f"{k:28s} c = {[round(x, 3) for x in c['c_n'][:3]]} m/s  R1 = {c['rossby_radius_m'][0]/1e3:.2f} km  "
              f"Bu1 = {c['burger'][0]:.3f}  TK1 = {c['kelvin_circuit_period_h'][0]:.1f} h  "
              f"Ti = {c['inertial_period_h']:.2f} h  Merian co = {c['merian_closed_open_h'][0]:.2f} h")


if __name__ == "__main__":
    main()
