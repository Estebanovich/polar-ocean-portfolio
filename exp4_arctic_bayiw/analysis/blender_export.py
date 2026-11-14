#!/opt/anaconda3/envs/pyocean/bin/python
r"""Scene export for the DESCRIPTIVE Blender animation of Exp. 4 (Task 10): one isopycnal surface
over the model bathymetry of Van Mijenfjorden, written as plain .npz/.json files that
analysis/blender_render.py reads inside Blender (Blender's own Python has numpy only).

Usage (pyocean env; LaCie read-only)::

    blender_export.py A  --out DIR [--t-end-h 96]    # shot A: F1, summer profile
    blender_export.py B  --out DIR [--t-end-h 96]    # shot B: F3 open water | F2 landfast ice
    blender_export.py probe --runs F3 F2 [--every-h 3] [--t-end-h 96]   # isopycnal candidates
    blender_export.py restyle --scene DIR [--limit L] [--mask-sponge]  # after the style gate

No number of this animation enters the README claims; no hypothesis or criterion uses it.

Field. The isopycnal sigma is extracted from every 3-D density snapshot (one every 1200 s) as the
depth z_s(x, y, t) of the first crossing rho = sigma from the surface down, linear between the
level centres, the wet centre for a partial bottom cell (isopycnal_depth; NaN where the column is denser than sigma at the top cell, i.e.
outcropped, or lighter at its last wet cell, i.e. the isopycnal meets the bed). Its displacement
from t = 0 is zeta = z_s(t) - z_s(0) [m, + up] (displacement). The GEOMETRY of the surface is the
raw z_s(t). Its COLOUR is the inertially filtered zeta' = zeta minus its centred running mean over
exactly one inertial period T_f = 2 pi/|f| (animate.inertial_filter, the Task 9 filter; the t = 0
reference cancels in zeta'), on one fixed symmetric scale for every frame of the shot. Frames
within T_f/2 of the record ends have no full window and are not exported (no padding); the record
read extends to t_end + T_f/2 so the filtered frames reach t_end, and they start at T_f/2
(6.1 h in the fjord). A requested still inside the dropped edge is moved to the first kept frame
and the shift is written to meta.json (select_frames).

Isopycnal choice (rule per shot; the numbers are written to meta.json):

* Shot A (F1, summer): sigma = the density half-way between the two levels of the largest N_b^2
  of the fixed t = 0 reference rho_b (budget.background_profile, Enmienda 1) below the mixed layer
  (max_n2_sigma with min_depth_m = the depth of the uniform top layer). ``probe`` checks that the
  pulse does not mix it away (coverage of the wet columns over the record).
* Shot B (F3 | F2, winter): the max-N_b^2 rule picks the top cell, whose 0.03 kg m^-3 surface step
  the pulse mixes away within 3 h (Task 9). The winter profile is a staircase. sigma is the
  SHALLOWEST stratified interface (N_b^2 >= WINTER_N2_MIN) no deeper than the Akselsundet sill
  (figures/bathy_report.json, so the surface spans the inner fjord under the ice) whose isopycnal
  keeps >= WINTER_MIN_SHARE of its t = 0 columns at every probed time (every 3 h, 0-96 h) in both
  runs (probe + choose_winter_sigma; Esteban, fix round 1). One sigma for both runs (same profile);
  the probe table is written to meta.json.

Scene. Domain cropped to Bellsund + Van Mijenfjorden: the open-boundary rows and their sponges of
data.obcs (sponge_mask: the OB row plus spongeThickness cells, W, N and S) removed from the frames
(crop_outside); shot B shows them as a grey strip (restyle --mask-sponge). Seafloor = - Depth (model grid). Land =
IBCAO v5.2 100 m elevation (Jakobsson et al. 2024) area-averaged onto the model land cells (the
same source and averaging as gendata/bathy.to_model_grid), 0 at the coast. One vertical
exaggeration VE for seafloor, land, isopycnal and ice (stated in the overlay). Landfast ice (F2):
the AREA.bin cells (>= 0.5) with the thickness h_m of case.json.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

# Blender's Python imports this module for quad_faces/lut_colour: numpy-only at the top level.

VE = 20.0                       # vertical exaggeration below sea level (seafloor, isopycnal, ice)
VE_LAND = 2.0                   # vertical exaggeration of the land relief (IBCAO up to ~1000 m)
PCTL = 99.5                     # colour limit: nice ceiling of this percentile of |zeta'|
WINTER_N2_MIN = 5e-6            # [s^-2] interfaces weaker than this are not candidates (winter)
WINTER_MIN_SHARE = 0.98         # the winter isopycnal keeps >= this share of its columns (all runs)
G = 9.81
SHOTS = {"A": ["F1"], "B": ["F3", "F2"]}


# ---- pure pieces --------------------------------------------------------------------------------

def isopycnal_depth(rho, RC, hfac, sigma: float, drF=None) -> np.ndarray:
    """Depth z [m, negative] of the first crossing rho = sigma from the top cell down, linear between
    the level centres; rho, hfac (nz, ny, nx). Only the wet levels (hfac > 0) count. NaN where the
    top wet cell is already >= sigma (outcropped), where no wet level reaches sigma (the isopycnal
    meets the bed) and on land. With drF, the centre of a partial bottom cell is its wet centre
    RF_k - hFac drF_k / 2 (RF_k = RC_k + drF_k / 2, the top face; partial cells of
    model/src/ini_masks_etc.F:105-121 fill the cell from its top); without, the nominal RC."""
    rho = np.asarray(rho, float)
    wet = np.asarray(hfac) > 0
    r = np.where(wet, rho, -np.inf)                         # dry levels never reach sigma
    ge = r >= sigma
    k = np.argmax(ge, axis=0)                               # first level >= sigma (0 if none)
    found = ge.any(axis=0) & (k > 0)
    kk = np.clip(k, 1, rho.shape[0] - 1)
    r1 = np.take_along_axis(rho, kk[None], 0)[0]
    r0 = np.take_along_axis(rho, (kk - 1)[None], 0)[0]
    RC = np.asarray(RC, float)
    if drF is None:
        zc = np.broadcast_to(RC[:, None, None], rho.shape)
    else:
        dr = np.asarray(drF, float)[:, None, None]
        zc = RC[:, None, None] + dr / 2 - np.asarray(hfac, float) * dr / 2
    z1 = np.take_along_axis(zc, kk[None], 0)[0]
    z0 = np.take_along_axis(zc, (kk - 1)[None], 0)[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        w = (sigma - r0) / (r1 - r0)
    z = z0 + w * (z1 - z0)
    return np.where(found, z, np.nan)


def displacement(z, z0) -> np.ndarray:
    """Isopycnal displacement zeta = z - z0 [m], + up (z negative downward)."""
    return np.asarray(z, float) - np.asarray(z0, float)


def max_n2_sigma(rho_b, RC, min_depth_m: float, rho0: float, max_depth_m: float | None = None,
                 g: float = G) -> dict:
    """The interface of largest N_b^2 = (g/rho0) (rho_{k+1} - rho_k)/(RC_k - RC_{k+1}) [s^-2]
    (Gill 1982 eq. 3.7.4 [verificar], finite difference between the level centres) whose depth
    -(RC_k + RC_{k+1})/2 lies in [min_depth_m, max_depth_m]; sigma = (rho_k + rho_{k+1})/2."""
    rb = np.asarray(rho_b, float)
    RC = np.asarray(RC, float)
    n2 = g / rho0 * (rb[1:] - rb[:-1]) / (RC[:-1] - RC[1:])
    zi = 0.5 * (RC[:-1] + RC[1:])
    ok = -zi >= min_depth_m
    if max_depth_m is not None:
        ok &= -zi <= max_depth_m
    if not ok.any():
        raise ValueError("no interface in the depth range")
    k = int(np.flatnonzero(ok)[np.argmax(n2[ok])])
    return {"sigma": float(0.5 * (rb[k] + rb[k + 1])), "k_upper": k, "k_lower": k + 1,
            "z_m": float(zi[k]), "n2_s-2": float(n2[k])}


def tf_s(f: float) -> float:
    """Inertial period 2 pi/|f| [s]."""
    return 2.0 * np.pi / abs(f)


def filtered_colour(x, t_s, f: float) -> tuple:
    """(kept times, x') with the Task 9 one-inertial-period running-mean filter
    (animate.inertial_filter): frames within T_f/2 of the ends are dropped."""
    from exp4_arctic_bayiw.analysis import animate
    return animate.inertial_filter(x, t_s, tf_s(f))


def select_frames(t_s, onset_s: float, hours, tol_s: float = 1.0) -> list:
    """For each requested time [h after the onset], the nearest available frame time; ``shifted``
    marks a request that is not an available frame (e.g. inside the dropped filter edge)."""
    t = np.asarray(t_s, float)
    out = []
    for h in hours:
        want = onset_s + 3600.0 * float(h)
        k = int(np.argmin(np.abs(t - want)))
        out.append({"requested_h": float(h), "index": k, "t_s": float(t[k]),
                    "t_h": float((t[k] - onset_s) / 3600.0), "shifted": bool(abs(t[k] - want) > tol_s)})
    return out


def time_label(t_s: float, onset_s: float, f: float) -> str:
    """'t = H h = N inertial periods' after the pulse onset (animate.clock)."""
    from exp4_arctic_bayiw.analysis import animate
    h, n = animate.clock(t_s, onset_s, f)
    return f"t = {h:.1f} h = {n:.2f} inertial periods"


def symmetric_limit(frames: dict, pctl: float = PCTL) -> float:
    """Fixed symmetric colour limit L: animate.nice_ceil of the pctl-th percentile of |x| pooled over
    every finite value of every frame of every run of the shot."""
    from exp4_arctic_bayiw.analysis import animate
    a = np.concatenate([np.abs(np.ravel(np.asarray(v, float))) for v in frames.values()])
    a = a[np.isfinite(a)]
    return animate.nice_ceil(float(np.percentile(a, pctl)))


def lut_colour(v, L: float, lut, nan_rgb=(0.6, 0.6, 0.6)) -> np.ndarray:
    """RGB (..., 3) of v on the fixed scale [-L, L] through the colour table lut (n, 3), linear
    between the entries; values beyond +-L clip to the end colours; NaN -> nan_rgb."""
    v = np.asarray(v, float)
    lut = np.asarray(lut, float)
    n = lut.shape[0]
    s = (np.clip(np.nan_to_num(v, nan=0.0), -L, L) + L) / (2.0 * L) * (n - 1)
    out = np.stack([np.interp(s, np.arange(n), lut[:, c]) for c in range(3)], axis=-1)
    out[np.isnan(v)] = nan_rgb
    return out


def quad_faces(valid) -> np.ndarray:
    """(nf, 4) vertex indices (j nx + i, counter-clockwise seen from above) of the quads between
    four neighbouring cell centres that are all valid; valid (ny, nx) bool."""
    v = np.asarray(valid, bool)
    ny, nx = v.shape
    q = v[:-1, :-1] & v[:-1, 1:] & v[1:, 1:] & v[1:, :-1]
    j, i = np.nonzero(q)
    a = j * nx + i
    return np.stack([a, a + 1, a + nx + 1, a + nx], axis=1).astype(np.int64)


def _lut(name: str, n: int = 256) -> list:
    from matplotlib import colormaps
    from pfcommon import style
    cm = style.CMAPS["diverging"] if name == "diverging" else colormaps[name]
    return [list(map(float, cm(x)[:3])) for x in np.linspace(0, 1, n)]


@dataclass
class Meta:
    shot: str
    runs: list
    sigma: dict
    ve: float
    ve_land: float
    colour_limit_m: float
    colour_field: str
    onset_s: float
    f: dict
    frames: list
    crop: dict
    extra: dict = field(default_factory=dict)


def write_meta(meta: Meta, path: Path) -> None:
    """meta.json with the colour tables (the diverging map of pfcommon.style, cmocean balance)."""
    d = asdict(meta)
    d["lut_diverging"] = _lut("diverging")
    Path(path).write_text(json.dumps(d, indent=1))


def obcs_params(data_obcs: str) -> dict:
    """{'sponge': useOBCSsponge, 'thickness': spongeThickness (0 without the sponge), 'ob': {key:
    first non-zero index}} of a data.obcs namelist (OB_Iwest, OB_Ieast, OB_Jsouth, OB_Jnorth;
    1-based, negative counts from the far end)."""
    import re
    txt = re.sub(r"(?m)^\s*#.*$", "", data_obcs)

    def val(key):
        m = re.search(rf"\b{key}\s*=\s*([^\n]*)", txt, flags=re.I)
        return m.group(1).strip().rstrip(",") if m else None

    sponge = (val("useOBCSsponge") or ".FALSE.").upper().strip(".") == "TRUE"
    ob = {}
    for key in ("OB_Iwest", "OB_Ieast", "OB_Jsouth", "OB_Jnorth"):
        v = val(key)
        idx = [int(x) for x in re.findall(r"(?:\d+\*)?(-?\d+)", v)] if v else []
        idx = [i for i in idx if i != 0]
        if idx:
            ob[key] = idx[0]
    return {"sponge": sponge, "thickness": int(float(val("spongeThickness") or 0)) if sponge else 0,
            "ob": ob}


def sponge_mask(data_obcs: str, shape) -> np.ndarray:
    """(ny, nx) bool: the open-boundary rows/columns of data.obcs (obcs_params) plus, with
    useOBCSsponge, the spongeThickness cells inward of each (pkg/obcs/obcs_sponge.F:68-69,
    j = OB - jsl for jsl = 1..spongeThickness; the weight of the last one is zero, :74-77)."""
    ny, nx = shape
    p = obcs_params(data_obcs)
    T = p["thickness"]
    m = np.zeros((ny, nx), bool)
    for key, i in p["ob"].items():
        n = nx if key.startswith("OB_I") else ny
        k = i - 1 if i > 0 else n + i                        # 0-based index of the OB row
        inward = 1 if key in ("OB_Iwest", "OB_Jsouth") else -1
        ks = [k + inward * j for j in range(T + 1) if 0 <= k + inward * j < n]
        if key.startswith("OB_I"):
            m[:, ks] = True
        else:
            m[ks, :] = True
    return m


def crop_outside(mask) -> dict:
    """The largest index box of whole rows and columns with no masked cell: the masked edge strips
    (sponge_mask) removed, so the frames never overlap the grey mask."""
    m = np.asarray(mask, bool)
    rows = np.flatnonzero(~m.all(axis=1))
    cols = np.flatnonzero(~m.all(axis=0))
    j0, j1, i0, i1 = rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
    while m[j0, i0:i1].any() and j0 < j1:
        j0 += 1
    while m[j1 - 1, i0:i1].any() and j1 > j0:
        j1 -= 1
    while m[j0:j1, i0].any() and i0 < i1:
        i0 += 1
    while m[j0:j1, i1 - 1].any() and i1 > i0:
        i1 -= 1
    return {"j0": int(j0), "j1": int(j1), "i0": int(i0), "i1": int(i1)}


def ice_cells(area_path: Path, shape) -> np.ndarray:
    """Landfast-ice cells: AREA.bin (big-endian float64, the actual model input) >= 0.5."""
    return np.fromfile(area_path, ">f8").reshape(shape) >= 0.5


def top_layer_depth(rho_b, drF) -> float:
    """Bottom face [m] of the uniform top layer of rho_b (levels equal to the top one within
    1e-6 kg m^-3); the top cell alone if the profile is stratified from the surface."""
    rb = np.asarray(rho_b, float)
    k = int(np.flatnonzero(np.abs(rb - rb[0]) > 1e-6)[0]) - 1
    return float(np.sum(np.asarray(drF, float)[:k + 1]))


def choose_winter_sigma(rows: list, sill_m: float, min_share: float = 0.98,
                        n2_min: float | None = None) -> dict:
    """The winter isopycnal (rule of Esteban's fix round 1): among the probe rows (stratified
    interfaces, N_b^2 >= n2_min) no deeper than the sill, the SHALLOWEST whose isopycnal keeps at
    least min_share of its t = 0 columns at every probed time in every run. ValueError if none."""
    n2_min = WINTER_N2_MIN if n2_min is None else n2_min
    ok = [r for r in rows if r["n2_s-2"] >= n2_min and -r["z_m"] <= sill_m
          and all(v["share_min"] >= min_share for v in r["runs"].values())]
    if not ok:
        raise ValueError("no interface above the sill keeps the share in every run")
    return min(ok, key=lambda r: -r["z_m"])


def saturation_stats(frames: dict, L: float) -> dict:
    """{run: {mean, worst, worst_index}}: per frame the fraction of the finite pixels with
    |x| >= L, its mean over the frames and the worst frame."""
    out = {}
    for n, a in frames.items():
        a = np.asarray(a, float)
        fr = []
        for x in a:
            x = x[np.isfinite(x)]
            fr.append(float(np.mean(np.abs(x) >= L)) if x.size else 0.0)
        out[n] = {"mean": float(np.mean(fr)), "worst": float(np.max(fr)), "worst_index": int(np.argmax(fr))}
    return out


# ---- model runs ---------------------------------------------------------------------------------

def _sponge(run_dir: Path, shape) -> np.ndarray:
    return sponge_mask((Path(run_dir) / "data.obcs").read_text(), shape)


def _crop(run_dir: Path, shape) -> dict:
    """Index box of the frames: the grid without the open boundaries and sponges of data.obcs."""
    return crop_outside(_sponge(run_dir, shape))


def sill_depth_m() -> float:
    """Akselsundet sill depth on the model grid [m] (figures/bathy_report.json)."""
    rep = json.loads((REPO / "exp4_arctic_bayiw" / "figures" / "bathy_report.json").read_text())
    return float(rep["sounds"]["Akselsundet"]["sill_depth_grid_m"])


def _sl(c: dict):
    return slice(c["j0"], c["j1"]), slice(c["i0"], c["i1"])


def land_elevation(case: dict, land) -> np.ndarray:
    """IBCAO v5.2 elevation [m, >= 0] area-averaged onto the model cells (the averaging of
    gendata/bathy.to_model_grid), set to 0 on the model's water cells."""
    from exp4_arctic_bayiw.gendata import bathy
    p = case["params"]
    src = bathy.load_ibcao(REPO / "exp4_arctic_bayiw" / "data_external" / p["bathy"]["ibcao_dir"],
                           bathy.bbox_utm_to_lonlat(p["bbox_utm"]))
    g = bathy.to_model_grid(src, p["bbox_utm"], p["dx_m"])
    c = g["src_cell"].ravel()
    z = src["z"].ravel()
    ins = c >= 0
    ny, nx = np.asarray(land).shape
    num = np.bincount(c[ins], weights=np.maximum(z[ins], 0.0), minlength=ny * nx)
    den = np.bincount(c[ins], minlength=ny * nx)
    elev = (num / np.maximum(den, 1)).reshape(ny, nx)
    return np.where(land, elev, 0.0)


def open_run(run_dir: Path) -> dict:
    from exp4_arctic_bayiw.analysis import animate, budget, load
    ds = load.run(run_dir)
    hfac = ds["hFacC"].values
    wet3 = hfac > 0
    rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, wet3)     # Enmienda 1
    r = {"ds": ds, "hfac": hfac, "rho_b": rho_b, "land": ~wet3.any(axis=0),
         "f": animate._f0(run_dir), "onset_s": animate._onset(run_dir),
         "case": animate._case(run_dir), "run_dir": Path(run_dir)}
    ice = r["case"].get("actual", {}).get("ice")
    r["ice"] = (ice_cells(r["run_dir"] / "AREA.bin", r["land"].shape)
                if ice and ice.get("mode") not in (None, "none") else None)
    return r


def iso_frames(r: dict, t_s, sigma: float, crop: dict) -> np.ndarray:
    """(nt, ny_c, nx_c) isopycnal depth at the record times t_s (one 3-D snapshot each), cropped."""
    ds = r["ds"]
    tm = ds["time"].values
    RC, drF = ds["RC"].values, ds["drF"].values
    js, is_ = _sl(crop)
    hf = r["hfac"][:, js, is_]
    out = []
    for t in t_s:
        n = int(np.argmin(np.abs(tm - t)))
        if abs(tm[n] - t) > 1e-6:
            raise ValueError(f"{r['run_dir']}: no record at t = {t} s")
        rho = ds["RHO"].isel(time=n).values[:, js, is_]
        out.append(isopycnal_depth(rho, RC, hf, sigma, drF=drF).astype(np.float32))
    return np.stack(out)


def record_times(r: dict, t_end_h: float) -> np.ndarray:
    """Record times from the onset to onset + t_end_h + T_f/2 (or the end of the record)."""
    t = r["ds"]["time"].values
    t1 = r["onset_s"] + t_end_h * 3600.0 + tf_s(r["f"]) / 2 + 1e-6
    return t[(t >= r["onset_s"] - 1e-6) & (t <= t1)]


def _interface(rb, RC, rho0: float, k: int) -> dict:
    return {"sigma": float(0.5 * (rb[k] + rb[k + 1])), "k_upper": k, "k_lower": k + 1,
            "z_m": float(0.5 * (RC[k] + RC[k + 1])),
            "n2_s-2": float(G / rho0 * (rb[k + 1] - rb[k]) / (RC[k] - RC[k + 1])),
            "step_kg_m3": float(rb[k + 1] - rb[k])}


def sigma_for_shot(shot: str, runs: dict, every_h: float = 3.0, t_end_h: float = 96.0) -> dict:
    """The sigma of a shot by the rule of the module docstring, applied in code. Shot B runs the
    probe (every_h, 0..t_end_h) and keeps its table under ``probe``."""
    r = next(iter(runs.values()))
    rb, RC = r["rho_b"], r["ds"]["RC"].values
    rho0 = float(r["ds"].attrs["rho0"])
    if shot == "A":
        ml = top_layer_depth(rb, r["ds"]["drF"].values)
        s = max_n2_sigma(rb, RC, min_depth_m=ml, rho0=rho0)
        s = _interface(rb, RC, rho0, s["k_upper"])
        s["rule"] = ("half-way density of the max-N_b^2 interface of the t = 0 reference below the "
                     f"uniform top layer ({ml:.2f} m)")
        return s
    sill = sill_depth_m()
    rows = probe(runs, every_h, t_end_h, max_depth_m=sill)
    c = choose_winter_sigma(rows, sill)
    s = _interface(rb, RC, rho0, c["k_upper"])
    s["rule"] = (f"winter staircase: the shallowest interface (N_b^2 >= {WINTER_N2_MIN:g} s^-2) no deeper "
                 f"than the Akselsundet sill ({sill:.2f} m) whose isopycnal keeps >= "
                 f"{100 * WINTER_MIN_SHARE:g} % of its t = 0 columns at every probed time "
                 f"(every {every_h:g} h, 0-{t_end_h:g} h) in every run")
    s["probe"] = rows
    return s


def probe(runs: dict, every_h: float, t_end_h: float, max_depth_m: float = 60.0) -> list:
    """For each stratified interface (N_b^2 >= WINTER_N2_MIN, depth <= max_depth_m) of the t = 0
    reference: per run, the share of the t = 0 columns that still hold the isopycnal (all wet
    columns and, if the run has ice, the ice cells), its minimum over the times, and the median
    raw displacement [m] at the last time."""
    r0 = next(iter(runs.values()))
    rb, RC = r0["rho_b"], r0["ds"]["RC"].values
    rho0 = float(r0["ds"].attrs["rho0"])
    n2 = G / rho0 * (rb[1:] - rb[:-1]) / (RC[:-1] - RC[1:])
    zi = 0.5 * (RC[:-1] + RC[1:])
    cand = [k for k in range(len(n2)) if n2[k] >= WINTER_N2_MIN and -zi[k] <= max_depth_m]
    crop = _crop(r0["run_dir"], r0["land"].shape)
    sigmas = [0.5 * (rb[k] + rb[k + 1]) for k in cand]
    js, is_ = _sl(crop)
    rows = [{"k_upper": k, "z_m": float(zi[k]), "n2_s-2": float(n2[k]), "sigma": float(s), "runs": {}}
            for k, s in zip(cand, sigmas)]
    for n, r in runs.items():
        ds = r["ds"]
        t = ds["time"].values
        step = max(1, int(round(every_h * 3600 / float(np.median(np.diff(t))))))
        idx = np.flatnonzero((t >= r["onset_s"]) & (t <= r["onset_s"] + t_end_h * 3600 + 1e-6))[::step]
        hf = r["hfac"][:, js, is_]
        drF = ds["drF"].values
        z = np.empty((len(cand), len(idx)) + hf.shape[1:], np.float32)
        for q, n_t in enumerate(idx):                       # one 3-D read per time, all candidates
            rho = ds["RHO"].isel(time=int(n_t)).values[:, js, is_]
            for c, s in enumerate(sigmas):
                z[c, q] = isopycnal_depth(rho, RC, hf, s, drF=drF)
        ice = r["ice"][js, is_] if r["ice"] is not None else None
        for c, row in enumerate(rows):
            have0 = np.isfinite(z[c, 0])
            cov = [float(np.mean(np.isfinite(zz[have0]))) if have0.any() else 0.0 for zz in z[c]]
            d = {"columns_t0": int(have0.sum()), "share_min": min(cov), "share_last": cov[-1],
                 "median_zeta_last_m": float(np.nanmedian(z[c, -1] - z[c, 0]))}
            if ice is not None:
                ic = ice & have0
                d["ice_columns_t0"] = int(ic.sum())
                d["ice_share_min"] = min(float(np.mean(np.isfinite(zz[ic]))) if ic.any() else 0.0
                                         for zz in z[c])
            row["runs"][n] = d
    for row in rows:
        print(json.dumps(row), flush=True)
    return rows


def write_terrain(out: Path, runs: dict, crop: dict, mask_sponge: bool) -> dict:
    """terrain.npz: the data crop, or with ``mask_sponge`` the whole grid plus the open-boundary
    and sponge mask of data.obcs (identical in every run, asserted). Returns the meta keys to set."""
    names = list(runs)
    r0 = runs[names[0]]
    ds, land = r0["ds"], r0["land"]
    for r in runs.values():
        if not np.array_equal(r["land"], land):
            raise ValueError("the runs of a shot must share the grid")
    elev = land_elevation(r0["case"], land)
    js, is_ = _sl(crop) if not mask_sponge else (slice(None), slice(None))
    ter = {"xc_m": ds["XC"].values[0, is_], "yc_m": ds["YC"].values[js, 0],
           "depth_m": ds["Depth"].values[js, is_], "land": land[js, is_], "elev_m": elev[js, is_]}
    for n, r in runs.items():
        if r["ice"] is not None:
            ter[f"ice_{n}"] = r["ice"][js, is_]
    keys = {"elev_max_m": float(ter["elev_m"].max()), "depth_max_m": float(ter["depth_m"].max())}
    if mask_sponge:
        masks = [_sponge(r["run_dir"], land.shape) for r in runs.values()]
        if any(not np.array_equal(m, masks[0]) for m in masks):
            raise ValueError("the sponge differs between the runs of the shot")
        ter["sponge"] = masks[0]
        keys["view_offset"] = {"j": crop["j0"], "i": crop["i0"]}
        keys["sponge_cells"] = obcs_params((r0["run_dir"] / "data.obcs").read_text())["thickness"]
    np.savez_compressed(Path(out) / "terrain.npz", **ter)
    return keys


def export(shot: str, root: Path, out: Path, t_end_h: float, still_h, ve: float = VE,
           mask_sponge: bool = False) -> dict:
    """Write terrain.npz, frames/<run>/<NNNN>.npz (isopycnal depth and colour per kept frame),
    raw_<run>.npz (the raw isopycnal depth at every record read, for re-filtering without the run
    disk) and meta.json of one shot into ``out``."""
    names = SHOTS[shot]
    runs = {n: open_run(root / n) for n in names}
    onsets = {r["onset_s"] for r in runs.values()}
    fs = {r["f"] for r in runs.values()}
    if len(onsets) != 1 or len(fs) != 1:
        raise ValueError("the runs of a shot must share the onset and f")
    onset, f = onsets.pop(), fs.pop()
    r0 = runs[names[0]]
    crop = _crop(r0["run_dir"], r0["land"].shape)
    sig = sigma_for_shot(shot, runs, t_end_h=t_end_h)
    out = Path(out)
    (out / "frames").mkdir(parents=True, exist_ok=True)
    tkeys = write_terrain(out, runs, crop, mask_sponge)
    zs, zf, tk, cover = {}, {}, None, {}
    for n, r in runs.items():
        t = record_times(r, t_end_h)
        z = iso_frames(r, t, sig["sigma"], crop)
        np.savez_compressed(out / f"raw_{n}.npz", t_s=t, z_m=z)
        zeta = displacement(z, z[0][None])
        tkn, zfn = filtered_colour(zeta, t, f)
        keep = np.isin(t, tkn)
        zs[n], zf[n] = z[keep], zfn
        tk = tkn if tk is None else tk
        if not np.array_equal(tk, tkn):
            raise ValueError("kept frame times differ between the runs")
        have0 = np.isfinite(z[0])
        cover[n] = {"columns_t0": int(have0.sum()),
                    "share_min": float(min(np.mean(np.isfinite(zz[have0])) for zz in z)),
                    "raw_zeta_p99_5_m": float(np.nanpercentile(np.abs(zeta), PCTL))}
    tk = tk[tk <= onset + t_end_h * 3600.0 + 1e-6]
    zs = {n: v[:tk.size] for n, v in zs.items()}
    zf = {n: v[:tk.size] for n, v in zf.items()}
    L = symmetric_limit(zf)
    sat = {n: float(np.mean(np.abs(v[np.isfinite(v)]) >= L)) for n, v in zf.items()}
    frames = []
    for k, t in enumerate(tk):
        for n in names:
            (out / "frames" / n).mkdir(exist_ok=True)
            np.savez_compressed(out / "frames" / n / f"{k:04d}.npz", z_m=zs[n][k], colour_m=zf[n][k])
        frames.append({"index": k, "t_s": float(t), "label": time_label(t, onset, f)})
    stills = select_frames(tk, onset, still_h)
    from exp4_arctic_bayiw.analysis import animate
    t_tau = np.arange(onset, onset + t_end_h * 3600.0 + 1e-6, animate.TAU_STEP_S)
    taus = {n: animate._tau_peak(r["ds"], t_tau) for n, r in runs.items()}
    tau_same = all(np.allclose(v, taus[names[0]]) for v in taus.values())
    probe_rows = sig.pop("probe", None)
    meta = Meta(shot=shot, runs=names, sigma={n: sig for n in names}, ve=ve, ve_land=VE_LAND, colour_limit_m=L,
                colour_field="filt", onset_s=float(onset), f={n: float(f) for n in names},
                frames=frames, crop=crop,
                extra={"stills": stills, "coverage": cover, "saturated_fraction": sat,
                       "tau_t_s": t_tau.tolist(), "tau_peak_Nm2": taus[names[0]].tolist(),
                       "tau_same_in_all_runs": bool(tau_same),
                       "T_f_h": tf_s(f) / 3600.0, "t_end_h": t_end_h,
                       "ice_h_m": {n: (r["case"]["actual"]["ice"].get("h_m") if r["ice"] is not None
                                       else None) for n, r in runs.items()},
                       "probe": probe_rows,
                       **{k: v for k, v in tkeys.items() if k != "view_offset"}})
    write_meta(meta, out / "meta.json")
    if "view_offset" in tkeys:
        d = json.loads((out / "meta.json").read_text())
        d["view_offset"] = tkeys["view_offset"]
        (out / "meta.json").write_text(json.dumps(d, indent=1))
    return asdict(meta)


def restyle(scene: Path, root: Path, limit: float | None, mask_sponge: bool) -> dict:
    """Apply the style-gate decisions (Esteban, 2026-11-14) to an exported scene without re-reading
    the frames: ``limit``, a fixed colour limit set by decision (the percentile value is kept as
    ``colour_limit_p99_5_m``); ``mask_sponge``, the terrain on the whole grid with the open
    boundaries and sponges of data.obcs shown in grey (write_terrain)."""
    scene = Path(scene)
    meta = json.loads((scene / "meta.json").read_text())
    if limit is not None:
        meta["extra"].setdefault("colour_limit_p99_5_m", meta["colour_limit_m"])
        meta["colour_limit_m"] = float(limit)
        meta["extra"]["colour_limit_rule"] = "fixed by decision at the style gate (Esteban, 2026-11-14)"
    if mask_sponge:
        runs = {n: open_run(root / n) for n in meta["runs"]}
        keys = write_terrain(scene, runs, meta["crop"], True)
        meta["view_offset"] = keys.pop("view_offset")
        meta["extra"].update(keys)
    (scene / "meta.json").write_text(json.dumps(meta, indent=1))
    return meta


def main(argv=None) -> int:
    from pfcommon.paths import EXP4_RUNS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("what", choices=["A", "B", "probe", "restyle"])
    ap.add_argument("--scene", type=Path)
    ap.add_argument("--limit", type=float, default=None)
    ap.add_argument("--mask-sponge", action="store_true")
    ap.add_argument("--runs-root", type=Path, default=EXP4_RUNS)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--runs", nargs="+", default=["F3", "F2"])
    ap.add_argument("--t-end-h", type=float, default=96.0)
    ap.add_argument("--every-h", type=float, default=3.0)
    ap.add_argument("--still-h", type=float, nargs="+", default=[6.0, 24.0, 48.0])
    ap.add_argument("--ve", type=float, default=VE)
    a = ap.parse_args(argv)
    if a.what == "probe":
        probe({n: open_run(a.runs_root / n) for n in a.runs}, a.every_h, a.t_end_h)
        return 0
    if a.what == "restyle":
        m = restyle(a.scene, a.runs_root, a.limit, a.mask_sponge)
        print(json.dumps({k: m[k] for k in ("sigma", "colour_limit_m")} | {"view_offset": m.get("view_offset")},
                         indent=1))
        return 0
    m = export(a.what, a.runs_root, a.out, a.t_end_h, a.still_h, a.ve)
    m = {**m, "sigma": {n: {k: v for k, v in s_.items() if k != "probe"} for n, s_ in m["sigma"].items()}}
    m["extra"] = {k: v for k, v in m["extra"].items() if k not in ("tau_t_s", "tau_peak_Nm2", "probe")}
    print(json.dumps({k: m[k] for k in ("shot", "sigma", "colour_limit_m", "crop")} | {"extra": m["extra"]},
                     indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
