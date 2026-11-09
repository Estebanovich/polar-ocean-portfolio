"""Generate all MITgcm input for one Exp. 2 (ISOMIP+ cavity) case.

``make_case(name, out_dir, template_dir=None, ...)`` reads ``cases.yaml``, writes the
binary inputs (big-endian float64 via ``pfcommon.mitgcm_io.write_bin``), renders the
namelist templates and writes ``case.json`` with the *actual* parameters. It returns the
case.json content.

Grid: x 320-800 km, y 0-80 km (AD16 sect. 3.1.1), cells of dx_km (2 or 4 km), 36 levels of
20 m. Binary files, written under FIXED names that templates reference literally:

  bathy.bin                (ny, nx)      bed [m], 0 on land
                                         (draft' = s draft on floating columns only, R9)
  icetopo.bin              (ny, nx)      ice draft [m], 0 on land and open ocean
                                         (both snapped to MITgcm's partial cells when
                                         the case has hFacMin, so the model keeps them)
  shelficeLoadAnomaly.bin  (ny, nx)      pressure-load anomaly [Pa] (profiles.load_anomaly)
  T.init.bin S.init.bin    (nz, ny, nx)  initial state, horizontally uniform
  T.restore.bin S.restore.bin (nz, ny, nx) restoring targets (AD16 Eq. 18-19)
  rbcs_mask.bin            (nz, ny, nx)  restoring mask, (x - 790 km) / 10 km >= 0 (AD16 Eq. 20)

Templates: every regular, non-hidden file in ``template_dir`` (default
``exp2_isomip_cavity/config/<cfg>/input_templates/``, ``<cfg>`` = the case's config, or
its ``template_config``) is rendered with ``string.Template.substitute``; an unknown
placeholder raises KeyError. Variables: every scalar of the merged case parameters plus
``case``, ``config``, ``nx``, ``ny``, ``nz``, ``dx``, ``dy`` [m], ``f0`` [s^-1],
``gammaS``, ``nIter0`` (0) and ``nTimeSteps`` (= endTime / deltaT, which must be integral).
Every quoted ``*.bin`` name in a rendered file must be one of the files above
(RuntimeError, before anything is written).

Ensemble: ``load_cases`` expands the ``ensemble:`` block of cases.yaml into cases
``ens_T{Tbot}_zc{z_c}_s{s}``. The geometry group is a Latin hypercube in (Tbot, z_c) per s,
drawn with numpy's default_rng(seed); values are rounded to the precision of the name
(0.01 degC, 1 m), so name and actual parameters agree.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path
from string import Template

import numpy as np
import yaml

from pfcommon.mitgcm_io import write_bin

from . import geometry as geo
from . import profiles as pr

PKG = Path(__file__).resolve().parents[1]
CASES_FILE = PKG / "cases.yaml"
GEOM_FILE = PKG / "data_external" / "Ocean1_input_geom_v1.01.nc"
BED_SOURCE = "GFZ Ocean1_input_geom_v1.01.nc bedrockTopography"
DRAFT_SOURCE = "GFZ Ocean1_input_geom_v1.01.nc lowerSurface, calved where thickness < 100 m"
OMEGA = 2 * np.pi / 86164.0          # s^-1; MITgcm rotationPeriod default (model/src/set_defaults.F:115)
LX, LY = 480e3, 80e3                 # AD16 sect. 3.1.1
BINS = ("bathy.bin", "icetopo.bin", "shelficeLoadAnomaly.bin", "T.init.bin", "S.init.bin",
        "T.restore.bin", "S.restore.bin", "rbcs_mask.bin")


# ---- registry ---------------------------------------------------------------------------

def fmt_num(v: float) -> str:
    """Case-name number: at most 2 decimals, at least 1 (-0.5, 0.0, 0.37)."""
    s = f"{float(v) + 0.0:.2f}".rstrip("0")
    return s + "0" if s.endswith(".") else s


def ens_name(Tbot: float, z_c: float, s: float) -> str:
    return f"ens_T{fmt_num(Tbot)}_zc{z_c:.0f}_s{s:.1f}"


def _lhs(rng, n: int, lo: float, hi: float, decimals: int) -> list[float]:
    """One sample in each of n equal strata of [lo, hi], in random order, rounded to
    ``decimals`` but kept inside its stratum."""
    step = 10.0 ** -decimals
    edges = lo + (hi - lo) * np.arange(n + 1) / n
    strata = rng.permutation(n)
    out = []
    for k in strata:
        v = edges[k] + (edges[k + 1] - edges[k]) * rng.random()
        v = round(v, decimals)
        if k < n - 1 and v >= edges[k + 1]:
            v = round(v - step, decimals)
        out.append(float(min(max(v, lo), hi)))
    return out


def _expand_ensemble(spec: dict) -> dict:
    e = spec["ensemble"]
    cfg, d = e["config"], float(e["d"])
    new = {}

    def add(group, T, z, s, **extra):
        name = ens_name(T, z, s)
        if name in new or name in spec["cases"]:
            raise ValueError(f"duplicate ensemble case {name}")
        new[name] = {"config": cfg, "ensemble_group": group, "init_profile": "ramp",
                     "restore_profile": "ramp", "Tbot": float(T), "z_c": float(z), "d": d,
                     "s": float(s), **extra}

    f = e["factorial"]
    for T in f["Tbot"]:
        for z in f["z_c"]:
            add("factorial", T, z, f["s"])
    g = e["geometry"]
    rng = np.random.default_rng(g["seed"])
    for s in g["s"]:
        Ts = _lhs(rng, g["n_per_s"], *g["Tbot_range"], decimals=2)
        zs = _lhs(rng, g["n_per_s"], *g["z_c_range"], decimals=0)
        for T, z in zip(Ts, zs):
            add("geometry", T, z, s, lhs_seed=g["seed"])
    x = e["extrapolation"]
    for T in x["Tbot"]:
        for z in x["z_c"]:
            add("extrapolation", T, z, x["s"])
    return new


def load_cases(path: Path | None = None) -> dict:
    """cases.yaml with the ensemble block expanded into ``cases`` (cases.yaml order first)."""
    spec = yaml.safe_load(Path(path or CASES_FILE).read_text())
    if "ensemble" in spec:
        spec["cases"].update(_expand_ensemble(spec))
    return spec


def ensemble_cases(spec: dict) -> list[str]:
    """The ensemble members: the reference case, then the expanded ens_ cases."""
    ref = spec["ensemble"]["reference"]
    return [ref] + [n for n, c in spec["cases"].items()
                    if n != ref and c.get("ensemble_group") is not None]


# ---- helpers ----------------------------------------------------------------------------

def grid_centres(dx_km: float) -> tuple[np.ndarray, np.ndarray]:
    """Cell centres (x, y) [m] of the dx_km grid over 320-800 x 0-80 km."""
    dx = float(dx_km) * 1e3
    nx, ny = LX / dx, LY / dx
    if nx != int(nx) or ny != int(ny):
        raise ValueError(f"dx_km = {dx_km} does not divide the 480 x 80 km domain")
    return geo.X_MIN + dx * (np.arange(int(nx)) + 0.5), dx * (np.arange(int(ny)) + 0.5)


def _profile(kind: str, z: np.ndarray, p: dict):
    """(T, S, Tbot, z_c, d) of a profile kind; ocean0 is the ramp with z_c 360 m, d 720 m."""
    if kind == "ocean0":
        T, S = pr.ocean0(z)
        return T, S, pr.TBOT_WARM, -pr.ZB_DEEP / 2, -pr.ZB_DEEP
    if kind == "cold":
        T, S = pr.cold(z)
        return T, S, pr.TBOT_COLD, None, None
    if kind == "ramp":
        T, S = pr.ramp(z, p["Tbot"], p["z_c"], d=p["d"])
        return T, S, float(p["Tbot"]), float(p["z_c"]), float(p["d"])
    raise ValueError(f"unknown profile {kind!r}")


def _git_state() -> dict:
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PKG, capture_output=True,
                              text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=PKG,
                                    capture_output=True, text=True, check=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        return {"git_commit": "unknown", "git_dirty": None}
    return {"git_commit": head, "git_dirty": dirty}


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---- make_case --------------------------------------------------------------------------

def make_case(name: str, out_dir: Path, template_dir: Path | None = None,
              cases_file: Path | None = None, overrides: dict | None = None,
              geom_path: Path | None = None) -> dict:
    """Write all inputs of case ``name`` into ``out_dir``; return the case.json content."""
    spec = load_cases(cases_file)
    if name not in spec["cases"]:
        raise KeyError(f"unknown case {name!r}")
    case = spec["cases"][name]
    cfg = case["config"]
    p = {**spec["configs"][cfg], **{k: v for k, v in case.items() if k != "config"}}
    p.update(overrides or {})
    p["case"], p["config"] = name, cfg
    out_dir = Path(out_dir)
    if p.get("gammaT") is None:
        raise ValueError(f"{name}: gammaT is not calibrated yet (null in cases.yaml); "
                         "pass it with --set gammaT=VALUE")

    # ---- grid and time
    xc, yc = grid_centres(p["dx_km"])
    nx, ny, nz, dz = xc.size, yc.size, int(p["nz"]), float(p["dz"])
    dzv = np.full(nz, dz)
    zc = -(np.arange(nz) + 0.5) * dz
    nsteps = p["endTime"] / p["deltaT"]
    if abs(nsteps - round(nsteps)) > 1e-9:
        raise ValueError(f"{name}: endTime {p['endTime']} / deltaT {p['deltaT']} is not an integer")
    p.update(nx=nx, ny=ny, nz=nz, dx=p["dx_km"] * 1e3, dy=p["dx_km"] * 1e3,
             nIter0=0, nTimeSteps=int(round(nsteps)),
             f0=float(2 * OMEGA * np.sin(np.deg2rad(p["lat_deg"]))),   # f plane, AD16 sect. 3.1.1
             gammaS=p["gammaT"] / 35.0)                                # AD16 Table 4

    # ---- geometry (ruling R7: bed and draft from the same file, same block average)
    gpath = Path(geom_path or GEOM_FILE)
    bed = geo.load_bed(str(gpath), xc, yc)
    draft = geo.load_draft(str(gpath), xc, yc)
    floating = geo.load_floating(str(gpath), xc, yc)
    # Task 2: minimum column and connectivity on MITgcm's own partial cells (hFacMin),
    # bed and draft snapped to them (geometry.mitgcm_hfac, model/src/ini_masks_etc.F:105-195)
    hfac_min = p.get("hFacMin")
    bathy, icetopo, wet = geo.apply_scale_and_mask(bed, draft, p["s"], dzv,
                                                   min_cells=p["min_cells"], xc_m=xc,
                                                   x_connect_m=p["x_restore0_m"],
                                                   floating=floating, hfac_min=hfac_min)
    if hfac_min is None:
        frac = geo.cell_fractions(bathy, icetopo, dzv) * wet[None]
    else:
        frac = geo.mitgcm_hfac(bathy, icetopo, dzv, hfac_min)[0] * wet[None]

    # ---- profiles, load anomaly, restoring
    Ti, Si, *_ = _profile(p["init_profile"], zc, p)
    Tr, Sr, Tbot, z_c, d = _profile(p["restore_profile"], zc, p)
    pf = pr.load_anomaly(Ti, Si, dzv, rho_const=p["rho_ref"], g=p["gravity"])
    pload = np.where(wet, pr.load_anomaly_at_draft(pf, dzv, icetopo), 0.0)
    shape3 = (nz, ny, nx)
    ramp = np.clip((xc - p["x_restore0_m"]) / (p["x_restore1_m"] - p["x_restore0_m"]), 0.0, None)
    bins = {
        "bathy.bin": bathy, "icetopo.bin": icetopo, "shelficeLoadAnomaly.bin": pload,
        "T.init.bin": np.broadcast_to(Ti[:, None, None], shape3),
        "S.init.bin": np.broadcast_to(Si[:, None, None], shape3),
        "T.restore.bin": np.broadcast_to(Tr[:, None, None], shape3),
        "S.restore.bin": np.broadcast_to(Sr[:, None, None], shape3),
        "rbcs_mask.bin": np.broadcast_to(ramp[None, None, :], shape3),   # AD16 Eq. 20
    }
    assert set(bins) == set(BINS)

    # ---- templates
    if template_dir is None:
        template_dir = PKG / "config" / p.get("template_config", cfg) / "input_templates"
    template_dir = Path(template_dir)
    if not template_dir.is_dir():
        raise FileNotFoundError(f"template dir {template_dir} does not exist")
    tfiles = [f for f in sorted(template_dir.iterdir()) if f.is_file() and not f.name.startswith(".")]
    scalars = {k: v for k, v in p.items() if not isinstance(v, (list, dict))}
    rendered = {f.name: Template(f.read_text()).substitute(scalars) for f in tfiles}
    missing = []
    for fname, text in rendered.items():
        for line in text.splitlines():
            if line.lstrip().startswith(("#", "!")):
                continue
            missing += [f"{fname}: {ref}" for ref in re.findall(r"'([^']+\.bin)'", line) if ref not in bins]
    if missing:
        raise RuntimeError("namelists name files that were not written: " + "; ".join(missing))

    gl = geo.grounding_line(icetopo, wet, xc, yc)
    actual = {
        "s": float(p["s"]), "gammaT": float(p["gammaT"]), "gammaS": float(p["gammaS"]),
        "init_profile": p["init_profile"], "restore_profile": p["restore_profile"],
        "Tbot": Tbot, "z_c": z_c, "d": d,
        "n_wet_columns": int(wet.sum()), "n_wet_cells": int((frac > 0).sum()),
        "wet_volume_m3": float(frac.sum() * dz * p["dx"] * p["dy"]),
        **gl,
        "draft_min_m": float(icetopo[wet].min()) if wet.any() else None,
        "deepest_wet_ice_base_depth_m": (float(-icetopo[wet & (icetopo < 0)].min())
                                         if (wet & (icetopo < 0)).any() else None),
        "n_floating_columns": int(floating.sum()),
        "hFacMin": hfac_min,
        "geometry_snapped_to_hfac": hfac_min is not None,
        "floating_rule": "strict majority of 1 km floatingMask sub-cells; s scales only these (R9)",
        "bed_source": BED_SOURCE, "draft_source": DRAFT_SOURCE,
        "geom_file": gpath.name, "geom_sha256": _sha256(gpath),
        "load_anomaly_profile": p["init_profile"],
    }
    meta = {"params": p, "actual": actual, **_git_state()}

    out_dir.mkdir(parents=True, exist_ok=True)
    for fname, arr in bins.items():
        write_bin(out_dir / fname, arr)
    for fname, text in rendered.items():
        (out_dir / fname).write_text(text)
    (out_dir / "case.json").write_text(json.dumps(meta, indent=2, default=float))
    return json.loads(json.dumps(meta, default=float))
