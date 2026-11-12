"""Generate all MITgcm input for one Exp. 4 case: the thesis bay (bay_thesis, Task 3) and the
Van Mijenfjorden cases (fjord_vmf, Task 4; ``make_fjord_case``, documented below it).

``make_case(name, out_dir, ...)`` reads ``cases.yaml``, writes the binary inputs, renders the
namelist templates of ``config/bay_thesis/input_templates``, applies the namelist edit sets of
the case and writes ``case.json`` with the *actual* parameters. It returns the case.json
content.

The thesis is the MSc thesis run BayIW_Rectan/run_expand (MITgcm, 560 x 352 x 50, stretched
horizontal grid, beta-plane). Its inputs are read from the read-only thesis tree
(``thesis_dir`` of the config):

  bathy.bin     (ny, nx)       bahia_01_expand_bat.bin (R3: nobahia_01_expand_bat.bin), byte copy
  delX.bin      (nx,)          bahia_01_expand_dx.bin, byte copy
  delY.bin      (ny,)          bahia_01_expand_dy.bin, byte copy
  T.init.bin    (nz, ny, nx)   thesis_realistic: feb_temp_50zlev_560x352.bin, byte copy;
  S.init.bin                   vmf_summer_thesis_levels: horizontally uniform, see below
  tauy.bin      (nt, ny, nx)   meridional wind stress [N m^-2], ``thesis_wind_patch``

Namelists: every template is rendered with ``string.Template.substitute`` (``nIter0``,
``n_stage1_steps``, ``n_free_steps``). ``data`` is the forced first stage (thesis
data.stage1), ``data.free`` the free stage (thesis data.stage2), written only if the run is
longer than ``forced_h``; ``data.kpp`` only for cases with useKPP. Then the edit sets named by
the case's ``edits`` (``namelist_sets`` of the config, plus ``eos_vmf`` built from
``configs.fjord_vmf.eos``) are applied in order with ``set_nml_key``; an edit of ``data``
applies to ``data`` and ``data.free`` (two stages of one run). Every quoted ``*.bin`` name in
a live namelist line must be a written file (RuntimeError).

R2/R3 profile (``vmf_summer_thesis_levels``): the Van Mijenfjorden summer representative
profile (figures/profiles_vmf.json, Task 2) interpolated linearly to the thesis cell centres,
the last value repeated below its deepest level, then made statically stable with the fitted
linear EOS (gendata.profiles.stabilize), as analysis/theory_numbers.py does for the theory
numbers (``summer_on_thesis_levels``).
"""
from __future__ import annotations

import copy
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path
from string import Template

import numpy as np
import yaml

from pfcommon.mitgcm_io import write_bin

from . import bathy as bth
from . import forcing
from . import ice as icemod
from . import profiles as prof
from . import vgrid

PKG = Path(__file__).resolve().parents[1]
CASES_FILE = PKG / "cases.yaml"
PROFILES_JSON = PKG / "figures" / "profiles_vmf.json"
BINS = ("bathy.bin", "delX.bin", "delY.bin", "T.init.bin", "S.init.bin", "tauy.bin")
THESIS_TS = {"T.init.bin": "feb_temp_50zlev_560x352.bin", "S.init.bin": "feb_salt_50zlev_560x352.bin"}
BETA_DEFAULT = 1.0e-11      # m^-1 s^-1, MITgcm default (model/src/set_defaults.F:112)
STAGE_FILES = {"data": ("data", "data.free")}


# ---- namelists ------------------------------------------------------------------------------

def _live(line: str) -> bool:
    s = line.strip()
    return bool(s) and not s.startswith(("#", "!"))


def _is_end(line: str) -> bool:
    s = line.strip().upper()
    return s in ("&", "/", "&END", "$END")


_KEY = re.compile(r"\s*([A-Za-z_]\w*(?:\([^)]*\))?)\s*=")


def _norm_key(k: str) -> str:
    return re.sub(r"\s+", "", k).lower()


def _norm_val(v: str) -> str:
    return re.sub(r"\s+", " ", v).strip().rstrip(",").strip()


def parse_namelist(text: str) -> dict[str, dict[str, str]]:
    """{GROUP: {key (lower case, no blanks): value text (blanks collapsed, no trailing
    comma)}} of a MITgcm namelist file. Lines starting with # or ! are comments (MITgcm
    strips them before reading). One key per line; a value may continue on the next lines
    (the thesis tRef/sRef/delR lists)."""
    out: dict[str, dict[str, str]] = {}
    group, key, buf = None, None, []

    def flush():
        if group is not None and key is not None:
            out[group][key] = _norm_val(" ".join(buf))

    for line in text.splitlines():
        if not _live(line):
            continue
        s = line.strip()
        if group is None:
            m = re.match(r"&(\w+)", s)
            if m:
                group, key, buf = m.group(1).upper(), None, []
                out.setdefault(group, {})
            continue
        if _is_end(line):
            flush()
            group, key, buf = None, None, []
            continue
        m = _KEY.match(line)
        if m:
            flush()
            key, buf = _norm_key(m.group(1)), [line[m.end():]]
        else:
            buf.append(line)
    return out


def parse_namelist_list(text: str, name: str) -> list[float]:
    """Numbers of the (possibly multi-line) value of ``name`` in any group."""
    for g in parse_namelist(text).values():
        v = g.get(name.lower())
        if v is not None:
            return [float(x) for x in re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?",
                                                 v.replace("D", "E").replace("d", "e"))]
    raise KeyError(name)


def set_nml_key(text: str, group: str, key: str, value: str | None) -> str:
    """Set ``key = value`` in namelist ``group`` of ``text`` (replace its line and its
    continuation lines, or insert it before the group terminator); ``value=None`` removes
    the key. Comment lines are kept. KeyError if the group is missing."""
    lines = text.splitlines(keepends=True)
    g0 = g1 = None
    for i, ln in enumerate(lines):
        if g0 is None:
            if _live(ln) and re.match(rf"\s*&{re.escape(group)}\s*$", ln, re.IGNORECASE):
                g0 = i
        elif _live(ln) and _is_end(ln):
            g1 = i
            break
    if g0 is None or g1 is None:
        raise KeyError(f"namelist group &{group} not found")
    k_norm = _norm_key(key)
    start = None
    for i in range(g0 + 1, g1):
        m = _KEY.match(lines[i]) if _live(lines[i]) else None
        if m and _norm_key(m.group(1)) == k_norm:
            start = i
            break
    new = [] if value is None else [f" {key}={value},\n"]
    if start is None:
        if value is None:
            return text
        return "".join(lines[:g1] + new + lines[g1:])
    end = start + 1
    while end < g1 and not lines[end].lstrip().startswith(("#", "!")) and not _KEY.match(lines[end]) \
            and not _is_end(lines[end]):
        end += 1
    return "".join(lines[:start] + new + lines[end:])


# ---- registry ---------------------------------------------------------------------------------

def load_cases(path: Path | None = None) -> dict:
    return yaml.safe_load(Path(path or CASES_FILE).read_text())


def case_params(spec: dict, name: str, overrides: dict | None = None) -> dict:
    """Merged parameters: config, then case, then overrides (KeyError for an unknown case)."""
    if name not in spec["cases"]:
        raise KeyError(f"unknown case {name!r}")
    case = spec["cases"][name]
    p = copy.deepcopy(spec["configs"][case["config"]])
    p.update(copy.deepcopy(case))
    p.update(overrides or {})
    p["case"] = name
    return p


def edit_sets(spec: dict, p: dict) -> list[tuple[str, dict]]:
    """[(set name, {file: {group: {key: value}}})] for the case's ``edits``, in order."""
    sets = dict(p.get("namelist_sets", {}))
    eos = spec["configs"]["fjord_vmf"]["eos"]
    nz = int(p["nz"])
    sets["eos_vmf"] = {"data": {"PARM01": {
        "tAlpha": f"{eos['tAlpha']}", "sBeta": f"{eos['sBeta']}",
        "rhoNil": f"{eos['rhoNil']}", "rhoConst": f"{eos['rhoNil']}",   # rhoConst = rhoNil
        "tRef": f"{nz}*{eos['tRef']}", "sRef": f"{nz}*{eos['sRef']}"}}}
    out = []
    for e in p.get("edits", []) or []:
        if e not in sets:
            raise KeyError(f"case {p['case']}: unknown namelist edit set {e!r}")
        out.append((e, sets[e]))
    return out


# ---- inputs -----------------------------------------------------------------------------------

def thesis_wind_patch(*, nt: int = 180, end_time_h: float = 60.0, max_tau: float = 0.10,
                      sigma_t_h: float = 1.0, center_t_h: float = 2.5, nx: int, ny: int,
                      cx: int, cy: int, radius: int) -> np.ndarray:
    """Meridional wind stress (nt, ny, nx) [N m^-2] of the thesis
    (input/wind_patch.py, make_meridional_wind_patch): a half Gaussian in time,
    tau = max_tau exp(-(t - t_c)^2 / (2 sigma_t^2)) for t <= t_c and 0 after, on
    t = linspace(0, end_time_h, nt) [h] (nt records written 1200 s apart in the model:
    externForcingPeriod of data); times a Gaussian patch exp(-r^2 / (2 (radius/2)^2)) in grid
    cells, cut at r > radius; times sin(90 deg)."""
    t = np.linspace(0, end_time_h, nt)
    g_t = max_tau * np.exp(-(t - center_t_h) ** 2 / (2 * sigma_t_h ** 2))
    g_t[t > center_t_h] = 0
    sig = radius / 2
    x = np.arange(nx)[:, None]
    y = np.arange(ny)[None, :]
    r2 = (x - cx) ** 2 + (y - cy) ** 2
    g_xy = np.exp(-r2 / (2 * sig ** 2))
    g_xy[r2 > radius ** 2] = 0
    tau = g_xy[:, :, None] * g_t[None, None, :] * np.sin(np.radians(90))   # (nx, ny, nt)
    return tau.transpose(2, 1, 0)


def thesis_cell_centres(delR) -> np.ndarray:
    f = vgrid.faces(delR)
    return -0.5 * (f[1:] + f[:-1])                 # m, positive up, surface first


def summer_on_thesis_levels(delR, eos: dict, profiles_json: Path = PROFILES_JSON):
    pj = json.loads(Path(profiles_json).read_text())
    zc = np.array(pj["z_c_m"])
    s = pj["summer"]
    tzc = thesis_cell_centres(delR)
    T = np.interp(-tzc, -zc, s["T"])               # constant beyond the last level
    S = np.interp(-tzc, -zc, s["S"])
    S, nadj = prof.stabilize(tzc, T, S, eos)
    return T, S, nadj


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_state() -> dict:
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PKG, capture_output=True,
                             text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=PKG,
                                    capture_output=True, text=True).stdout.strip())
        return {"git_commit": sha or None, "git_dirty": dirty}
    except OSError:
        return {"git_commit": None, "git_dirty": None}


def stage_plan(p: dict) -> list[list[int]]:
    """[[nIter0, nTimeSteps]]: stage 1 forced up to forced_h, then one free stage."""
    dt = float(p["deltaT"])
    n_tot = round(float(p["t_end_h"]) * 3600.0 / dt)
    n_f = round(float(p["forced_h"]) * 3600.0 / dt)
    if abs(n_tot * dt - float(p["t_end_h"]) * 3600.0) > 1e-6 or n_tot <= 0:
        raise ValueError(f"t_end_h = {p['t_end_h']} h is not a positive multiple of deltaT")
    n1 = min(n_f, n_tot)
    return [[0, n1]] + ([[n1, n_tot - n1]] if n_tot > n1 else [])


# ---- main -------------------------------------------------------------------------------------

def make_case(name: str, out_dir: Path, *, cases_file: Path | None = None,
              template_dir: Path | None = None, overrides: dict | None = None,
              write_bins=True) -> dict:
    """Write the inputs of case ``name`` to ``out_dir``; return the case.json content.

    ``write_bins``: True (all), or an iterable of the names in BINS to write (tests use a
    subset; the namelists and case.json are always written). ``overrides`` replaces
    top-level parameters (``run_campaign --set`` / ``--end-time``)."""
    spec = load_cases(cases_file)
    p = case_params(spec, name, overrides)
    cfg = spec["cases"][name]["config"]
    if cfg.startswith("fjord_vmf"):
        return make_fjord_case(name, out_dir, spec=spec, p=p, template_dir=template_dir,
                               write_bins=write_bins)
    sets = edit_sets(spec, p)
    nx, ny, nz = int(p["nx"]), int(p["ny"]), int(p["nz"])
    thesis_in = Path(p["thesis_dir"]) / "input"
    want = set(BINS) if write_bins is True else set(write_bins or ())
    out_dir = Path(out_dir)

    # ---- namelists
    plan = stage_plan(p)
    n1 = plan[0][1]
    n_free = plan[1][1] if len(plan) > 1 else 0
    if template_dir is None:
        template_dir = PKG / "config" / p.get("template_config", "bay_thesis") / "input_templates"
    tfiles = [f for f in sorted(Path(template_dir).iterdir()) if f.is_file() and not f.name.startswith(".")]
    subs = {"nIter0": 0, "n_stage1_steps": n1, "n_free_steps": n_free}
    rendered = {f.name: Template(f.read_text()).substitute(subs) for f in tfiles}
    if not n_free:
        rendered.pop("data.free", None)
    for _, es in sets:
        for fname, groups in es.items():
            for target in STAGE_FILES.get(fname, (fname,)):
                if target not in rendered:
                    continue
                for g, kv in groups.items():
                    for k, v in kv.items():
                        rendered[target] = set_nml_key(rendered[target], g, k, v)
    pkgs = parse_namelist(rendered["data.pkg"]).get("PACKAGES", {})
    use_kpp = pkgs.get("usekpp", ".FALSE.").upper() in (".TRUE.", "T", ".T.")
    if not use_kpp:
        rendered.pop("data.kpp", None)
    missing = []
    for fname, text in rendered.items():
        for line in text.splitlines():
            if _live(line):
                missing += [f"{fname}: {r}" for r in re.findall(r"'([^']+\.bin)'", line) if r not in BINS]
    if missing:
        raise RuntimeError("namelists name files that are not inputs: " + "; ".join(missing))

    nml = parse_namelist(rendered["data"])
    par1 = nml["PARM01"]
    delR = np.array(parse_namelist_list(rendered["data"], "delR"))
    if delR.size != nz:
        raise ValueError(f"delR has {delR.size} levels, nz = {nz}")
    f0 = float(par1["f0"])
    beta = float(par1["beta"]) if "beta" in par1 else BETA_DEFAULT
    eos = {"eosType": par1["eostype"].strip("'"), "tAlpha": float(par1["talpha"]),
           "sBeta": float(par1["sbeta"]), "rhoNil": float(par1["rhonil"]),
           "rhoConst": float(par1.get("rhoconst", par1["rhonil"])),   # model/src/ini_parms.F:483
           "tRef": par1["tref"], "sRef": par1["sref"]}

    # ---- binaries
    src = {"bathy.bin": thesis_in / p["bathy_file"], "delX.bin": thesis_in / p["delx_file"],
           "delY.bin": thesis_in / p["dely_file"]}
    if p["profile"] == "thesis_realistic":
        src.update({k: thesis_in / v for k, v in THESIS_TS.items()})
        T = np.fromfile(src["T.init.bin"], ">f8", count=nz * ny * nx)[:: ny * nx]
        S = np.fromfile(src["S.init.bin"], ">f8", count=nz * ny * nx)[:: ny * nx]
        prof_info = {"source": "thesis feb_{temp,salt}_50zlev_560x352.bin (byte copy)"}
    elif p["profile"] == "vmf_summer_thesis_levels":
        e = spec["configs"]["fjord_vmf"]["eos"]
        T, S, nadj = summer_on_thesis_levels(delR, e)
        prof_info = {"source": "figures/profiles_vmf.json summer, linear interpolation to the "
                               "thesis cell centres, last value below, stabilize (fitted EOS)",
                     "levels_adjusted": nadj}
    else:
        raise ValueError(f"unknown profile {p['profile']!r}")
    dx = np.fromfile(src["delX.bin"], ">f8")
    dy = np.fromfile(src["delY.bin"], ">f8")
    bathy = np.fromfile(src["bathy.bin"], ">f8")
    if dx.size != nx or dy.size != ny or bathy.size != nx * ny:
        raise ValueError("thesis grid files do not match nx, ny")
    bathy = bathy.reshape(ny, nx)
    w = p["wind"]
    tau = thesis_wind_patch(nt=w["nt"], end_time_h=w["end_time_h"], max_tau=w["max_tau"],
                            sigma_t_h=w["sigma_t_h"], center_t_h=w["center_t_h"], nx=nx, ny=ny,
                            cx=w["cx"], cy=w["cy"], radius=w["radius"]) if "tauy.bin" in want else None

    # ---- actual parameters
    yg0 = parse_namelist_list(rendered["data"], "ygOrigin")[0]
    yc = yg0 + np.cumsum(dy) - 0.5 * dy                        # Cartesian yC (m)
    fy = f0 + beta * yc                                        # model/src/ini_cori.F:77
    bay = (bathy < 0) & (bathy > -999.0)                       # bay: shallower than the 1000 m shelf
    jb = np.unique(np.nonzero(bay)[0]) if bay.any() else np.array([], int)
    actual = {
        "f0_s-1": f0, "beta_m-1s-1": beta,
        "f_domain_s-1": [float(fy.min()), float(fy.max())],
        "f_bay_s-1": [float(fy[jb].min()), float(fy[jb].max())] if jb.size else None,
        "inertial_period_bay_h": ([float(2 * np.pi / fy[jb].max() / 3600), float(2 * np.pi / fy[jb].min() / 3600)]
                                  if jb.size else None),
        "bay_rows_j": [int(jb.min()), int(jb.max())] if jb.size else None,
        "profile": p["profile"], "profile_info": prof_info,
        "T_column": [float(x) for x in T], "S_column": [float(x) for x in S],
        "eos": eos, "kpp": use_kpp, "surface": p["surface"],
        "wind": {**w, "max_tau_N_m-2": float(tau.max()) if tau is not None else None,
                 "record_period_s": float(parse_namelist(rendered["data"])["PARM03"]["externforcingperiod"])},
        "stages": plan, "edits": [e for e, _ in sets],
        "column_depth_m": float(delR.sum()), "bathy_max_depth_m": float(-bathy.min()),
        "effective_shelf_depth_m": float(min(delR.sum(), -bathy.min())),
        "sources": {k: {"path": str(v), "sha256": _sha256(v)} for k, v in src.items()},
    }
    meta = {"params": p, "actual": actual, "config": cfg, **_git_state()}

    out_dir.mkdir(parents=True, exist_ok=True)
    for k in ("bathy.bin", "delX.bin", "delY.bin"):
        if k in want:
            shutil.copyfile(src[k], out_dir / k)              # byte copy of the thesis file
    for k in ("T.init.bin", "S.init.bin"):
        if k not in want:
            continue
        if p["profile"] == "thesis_realistic":
            shutil.copyfile(src[k], out_dir / k)
        else:
            col = T if k == "T.init.bin" else S
            write_bin(out_dir / k, np.broadcast_to(col[:, None, None], (nz, ny, nx)))
    if tau is not None:
        write_bin(out_dir / "tauy.bin", tau)
    for fname, text in rendered.items():
        (out_dir / fname).write_text(text)
    (out_dir / "case.json").write_text(json.dumps(meta, indent=2, default=float))
    return json.loads(json.dumps(meta, default=float))


# ==== Van Mijenfjorden (fjord_vmf, Task 4) =====================================================
#
# make_fjord_case writes, for the cases of the fjord_vmf* configs (cases.yaml):
#
#   bathy.bin   (ny, nx)          Task 1b bathymetry, gendata.bathy.build at dx_m (fjord_vmf_100:
#                                 the 200 m grid split 2 x 2; fjord_vmf_100s: a sub-box of it;
#                                 timing only, bathy_source)
#   T.init.bin  (nz, ny, nx)      horizontally uniform summer / winter representative profile
#   S.init.bin                    (figures/profiles_vmf.json, Task 2, on the delR cell centres)
#   taux.bin    (nrec, ny, nx)    uniform pulse (gendata.forcing) towards wind_angle_deg,
#   tauy.bin                      records at (n - 1/2) * 1200 s, nrec = forced/1200 s + 1
#   HEFF.bin    (ny, nx)          ice cases: h and A inside the profile polygon (gendata.ice)
#   AREA.bin
#   OB?_?TidAm.bin, OB?_?TidPh.bin  V3: amplitude (m/s) and phase lag (s) of the inward M2 normal
#                                 velocity per OB point (pkg/obcs/OBCS_PARAMS.h:162)
#
# Namelists: the templates of config/fjord_vmf/input_templates rendered with the config values
# (string.Template.substitute), then the case's edit sets (namelist_sets) applied in order.
# data.free (free stage) = data with the wind forcing removed and nIter0/nTimeSteps of stage 2.
# data.kpp / data.ggl90 / data.seaice are kept only if data.pkg uses the package.

TIDE_FILES = ("OBW_uTidAm.bin", "OBW_uTidPh.bin", "OBN_vTidAm.bin", "OBN_vTidPh.bin",
              "OBS_vTidAm.bin", "OBS_vTidPh.bin")
FJORD_BINS = ("bathy.bin", "T.init.bin", "S.init.bin", "taux.bin", "tauy.bin", "HEFF.bin",
              "AREA.bin") + TIDE_FILES
DATA_EXTERNAL = PKG / "data_external"
WIND_KEYS = {"PARM03": ("periodicExternalForcing", "externForcingPeriod", "externForcingCycle"),
             "PARM05": ("zonalWindFile", "meridWindFile")}
_GRID_CACHE: dict = {}


def _true(v) -> bool:
    return str(v).strip().upper() in (".TRUE.", "T", ".T.", "TRUE")


def fjord_grid(p: dict, data_dir: Path = DATA_EXTERNAL) -> dict:
    """xc, yc (m), depth (m, > 0 at sea), wet of the case grid (cached per bathymetry config)."""
    keys = ("bbox_utm", "dx_m", "ob_sides", "n_sponge", "hFacMin", "hFacMinDr", "min_depth_m",
            "max_step_cells", "delR", "bathy", "bathy_source", "parent_bbox_utm")
    key = json.dumps({k: p.get(k) for k in keys}, sort_keys=True, default=str)
    if key in _GRID_CACHE:
        return _GRID_CACHE[key]
    if p.get("bathy_source") == "subbox_of_refine_200m":
        parent = {k: copy.deepcopy(p[k]) for k in keys if k in p}
        parent.update(bbox_utm=list(p["parent_bbox_utm"]), bathy_source="refine_200m",
                      ob_sides=["W", "N", "S"])
        gp = fjord_grid(parent, data_dir)
        px0, _, py0, _ = p["parent_bbox_utm"]
        x0, x1, y0, y1 = p["bbox_utm"]
        d = float(p["dx_m"])
        i0, i1 = round((x0 - px0) / d), round((x1 - px0) / d)
        j0, j1 = round((y0 - py0) / d), round((y1 - py0) / d)
        depth = gp["depth"][j0:j1, i0:i1].copy()
        xc, yc = bth.model_axes(p["bbox_utm"], d)
        g = {"xc": xc, "yc": yc, "depth": depth, "wet": depth > 0, "source": "subbox_of_refine_200m"}
    elif p.get("bathy_source") == "refine_200m":
        base = {k: copy.deepcopy(p[k]) for k in keys if k in p and k not in ("bathy_source", "parent_bbox_utm")}
        base["dx_m"] = 2 * float(p["dx_m"])
        g2 = fjord_grid(base, data_dir)
        depth = np.kron(g2["depth"], np.ones((2, 2)))
        xc, yc = bth.model_axes(p["bbox_utm"], float(p["dx_m"]))
        g = {"xc": xc, "yc": yc, "depth": depth, "wet": depth > 0, "source": "refine_200m"}
    else:
        _, gb = bth.build(p, data_dir)
        g = {"xc": gb["xc"], "yc": gb["yc"], "depth": gb["depth"], "wet": gb["wet"],
             "source": "gendata.bathy.build"}
    _GRID_CACHE[key] = g
    return g


def season_profile(profile: str, profiles_json: Path = PROFILES_JSON):
    """(T, S) columns of 'vmf_summer' / 'vmf_winter' on the fjord delR cell centres."""
    season = {"vmf_summer": "summer", "vmf_winter": "winter"}.get(profile)
    if season is None:
        raise ValueError(f"unknown fjord profile {profile!r}")
    pj = json.loads(Path(profiles_json).read_text())
    return np.array(pj[season]["T"], float), np.array(pj[season]["S"], float), pj["delR"]


def tide_v3(depth, wet, delR, dx: float, ob_sides, period_s: float, eta0_m: float,
            hfac_min: float, hfac_min_dr: float) -> dict:
    """Uniform inward barotropic M2 normal velocity U0 on the open boundaries, from the volume
    balance of a basin whose surface rises eta0 sin(omega t) everywhere:
    U0 * A_OB = omega * eta0 * A_wet  (A_OB: wet face area of the OB cells with partial cells,
    model/src/ini_masks_etc.F:105-121 via gendata.bathy.partial_cells; A_wet: wet surface area).
    The OB velocity is u_OB = U0 sin(omega t), starting at 0. pkg/obcs reads an amplitude Am
    (m/s) and a phase lag Ph (s) (pkg/obcs/OBCS_PARAMS.h:162), forms Cs = Am cos(omega Ph) and
    Sn = Am sin(omega Ph) (pkg/obcs/obcs_init_fixed.F:536-575) and adds Cs cos(omega t) +
    Sn sin(omega t) (pkg/obcs/obcs_add_tides.F), so Am = +-U0 and Ph = T/4. Sign: inward is +u
    at W, +v at S, -v at N. (Task 4 fix round 1: the first version wrote U0 into the phase
    file and no amplitude file, so the tide was zero.)"""
    faces = vgrid.faces(delR)
    col = bth.partial_cells(depth, faces, hfac_min, hfac_min_dr)       # (ny, nx) water column (m)
    col = np.where(wet, col, 0.0)
    a_ob = 0.0
    if "W" in ob_sides:
        a_ob += col[:, 0].sum() * dx
    if "N" in ob_sides:
        a_ob += col[-1, :].sum() * dx
    if "S" in ob_sides:
        a_ob += col[0, :].sum() * dx
    a_wet = float(np.asarray(wet).sum() * dx * dx)
    omega = 2 * np.pi / float(period_s)
    u0 = omega * float(eta0_m) * a_wet / a_ob
    w = np.asarray(wet)
    q = float(period_s) / 4.0
    files = {"OBW_uTidAm.bin": np.where(w[:, 0], u0, 0.0), "OBW_uTidPh.bin": np.where(w[:, 0], q, 0.0),
             "OBN_vTidAm.bin": np.where(w[-1, :], -u0, 0.0), "OBN_vTidPh.bin": np.where(w[-1, :], q, 0.0),
             "OBS_vTidAm.bin": np.where(w[0, :], u0, 0.0), "OBS_vTidPh.bin": np.where(w[0, :], q, 0.0)}
    return {"U0_m_s-1": float(u0), "ob_face_area_m2": float(a_ob), "wet_area_m2": a_wet,
            "eta0_m": float(eta0_m), "period_s": float(period_s), "files": files}


def _fmt_list(v, per_line: int = 8) -> str:
    items = [f"{x:g}" for x in v]
    return ",\n   ".join(", ".join(items[i:i + per_line]) for i in range(0, len(items), per_line))


def _apply_sets(rendered: dict, sets) -> dict:
    for _, es in sets:
        for fname, groups in es.items():
            if fname not in rendered:
                raise KeyError(f"edit set names {fname}, which is not a template")
            for g, kv in groups.items():
                for k, v in kv.items():
                    rendered[fname] = set_nml_key(rendered[fname], g, k, v)
    return rendered


def free_stage_data(data_text: str, n_iter0: int, n_steps: int) -> str:
    """data of the free stage: the forced data without the wind forcing (WIND_KEYS)."""
    t = data_text
    for g, ks in WIND_KEYS.items():
        for k in ks:
            t = set_nml_key(t, g, k, None)
    t = set_nml_key(t, "PARM03", "nIter0", str(int(n_iter0)))
    t = set_nml_key(t, "PARM03", "nTimeSteps", str(int(n_steps)))
    return t.replace("stage 1 (forced)", "free stage (data.free: data without the wind forcing)", 1)


def make_fjord_case(name: str, out_dir: Path, *, spec: dict, p: dict, template_dir=None,
                    write_bins=True) -> dict:
    """Write the inputs of fjord case ``name`` (see the block comment above)."""
    cfg = spec["cases"][name]["config"]
    sets = edit_sets(spec, p)
    names = [e for e, _ in sets]
    nx, ny, nz, dx = int(p["nx"]), int(p["ny"]), int(p["nz"]), float(p["dx_m"])
    x0, x1, y0, y1 = p["bbox_utm"]
    if round((x1 - x0) / dx) != nx or round((y1 - y0) / dx) != ny or len(p["delR"]) != nz:
        raise ValueError(f"{cfg}: nx, ny, nz do not match bbox_utm / dx_m / delR")
    want = set(FJORD_BINS) if write_bins is True else set(write_bins or ())
    out_dir = Path(out_dir)

    # ---- inputs
    g = fjord_grid(p)
    T, S, delR_prof = season_profile(p["profile"])
    if not np.allclose(delR_prof, p["delR"]):
        raise ValueError("profiles_vmf.json delR differs from the config delR")
    plan = stage_plan(p)
    n1 = plan[0][1]
    n_free = plan[1][1] if len(plan) > 1 else 0
    dt = float(p["deltaT"])
    w = p["wind"]
    period = float(w["record_period_s"])
    if abs(period / dt - round(period / dt)) > 1e-9:
        raise ValueError(f"deltaT = {dt} s does not divide the record period {period} s")
    angle = p.get("wind_angle_deg")
    use_wind = "no_wind" not in names
    if use_wind and angle is None:
        raise ValueError(f"case {name}: wind on but wind_angle_deg is null")
    nrec = forcing.cycle_records(p["forced_h"], period)
    bins: dict = {}
    wind_info = None
    if use_wind:
        tr = forcing.record_times(nrec, period)
        if want & {"taux.bin", "tauy.bin"}:
            tx, ty = forcing.uniform_stress(tr, angle, (ny, nx), tau_max=w["tau_max"],
                                            sigma_h=w["sigma_h"], center_h=w["center_h"])
            bins["taux.bin"], bins["tauy.bin"] = tx, ty
        mag = forcing.pulse(tr, tau_max=w["tau_max"], sigma_h=w["sigma_h"], center_h=w["center_h"])
        c, s_ = forcing._cos_sin(angle)
        wind_info = {**w, "angle_deg": float(angle), "unit_vector_xy": [c, s_], "n_records": nrec,
                     "record_times_s": [float(tr[0]), float(tr[-1])],
                     "max_tau_N_m-2": float(mag.max()),
                     "peak_record_time_s": float(tr[int(np.argmax(mag))]),
                     "cycle_s": nrec * period}
    ice_info = None
    if p.get("ice_mode"):
        poly = p["ice"]["polygon"]
        if isinstance(poly, str):
            poly = p[poly]["polygon_utm"]
        m = icemod.mask(g["xc"], g["yc"], g["wet"], poly)
        bins.update(icemod.fields(m, h=p["ice"]["h_m"], A=p["ice"]["A"]))
        ice_info = {"mode": p["ice_mode"], "cells": int(m.sum()), "area_km2": float(m.sum() * dx * dx / 1e6),
                    "h_m": float(p["ice"]["h_m"]), "A": float(p["ice"]["A"]),
                    "polygon": "configs.fjord_vmf.profiles.polygon_utm"}
    tide_info = None
    if "tide_m2" in names:
        td = tide_v3(g["depth"], g["wet"], p["delR"], dx, p["ob_sides"], p["tide"]["period_s"],
                     p["tide"]["eta0_m"], p["hFacMin"], p["hFacMinDr"])
        bins.update(td.pop("files"))
        tide_info = td

    # ---- namelists
    eos = spec["configs"]["fjord_vmf"]["eos"]
    zc = 0.5 * (vgrid.faces(p["delR"])[1:] + vgrid.faces(p["delR"])[:-1])
    lev = [int(k) for k in p["wlev_levels"]]
    subs = {"nIter0": 0, "n_stage1_steps": n1, "deltaT": f"{dt:g}",
            "staggerTimeStep": ".TRUE." if int(p.get("stagger_time_step", 0)) else ".FALSE.",
            "pChkptFreq": f"{float(p['forced_h']) * 3600:g}",
            "tAlpha": f"{eos['tAlpha']}", "sBeta": f"{eos['sBeta']}", "rhoNil": f"{eos['rhoNil']}",
            "tRef": f"{nz}*{eos['tRef']}", "sRef": f"{nz}*{eos['sRef']}", "f0": f"{p['f0']}",
            "hFacMin": f"{p['hFacMin']}", "hFacMinDr": f"{p['hFacMinDr']}",
            "forcing_period": f"{period:g}", "forcing_cycle": f"{nrec * period:g}",
            "delX": f"{nx}*{dx:.1f}", "delY": f"{ny}*{dx:.1f}", "xgOrigin": f"{float(x0):.1f}",
            "ygOrigin": f"{float(y0):.1f}", "delR": _fmt_list(p["delR"]),
            "nx": nx, "ny": ny, "n_sponge": int(p["n_sponge"]),
            "relax_bound": f"{max(60.0, 2 * dt):g}",
            "wlev_levels": ", ".join(f"{k}." for k in lev),
            "wlev_depths": ", ".join(f"{zc[k - 1]:.1f}" for k in lev)}
    if template_dir is None:
        template_dir = PKG / "config" / p.get("template_config", "fjord_vmf") / "input_templates"
    tfiles = [f for f in sorted(Path(template_dir).iterdir()) if f.is_file() and not f.name.startswith(".")]
    rendered = {f.name: Template(f.read_text()).substitute(subs) for f in tfiles}
    rendered = _apply_sets(rendered, sets)
    if n_free:
        rendered["data.free"] = free_stage_data(rendered["data"], n1, n_free)
    pk = parse_namelist(rendered["data.pkg"]).get("PACKAGES", {})
    use = {k: _true(pk.get(f"use{k}", ".FALSE.")) for k in ("kpp", "ggl90", "seaice")}
    for k, f in (("kpp", "data.kpp"), ("ggl90", "data.ggl90"), ("seaice", "data.seaice")):
        if not use[k]:
            rendered.pop(f, None)
    if use["seaice"] != bool(ice_info):
        raise RuntimeError(f"case {name}: useSEAICE = {use['seaice']} but ice_mode = {p.get('ice_mode')}")
    if use["kpp"] and use["ggl90"]:
        raise RuntimeError(f"case {name}: both KPP and GGL90 on")
    written = {"bathy.bin", "T.init.bin", "S.init.bin"} | set(bins) | ({"taux.bin", "tauy.bin"} if use_wind else set())
    missing = []
    for fname, text in rendered.items():
        for line in text.splitlines():
            if _live(line):
                missing += [f"{fname}: {r}" for r in re.findall(r"'([^']+\.bin)'", line) if r not in written]
    if missing:
        raise RuntimeError("namelists name files that are not inputs: " + "; ".join(missing))
    par1 = parse_namelist(rendered["data"])["PARM01"]

    # ---- actual parameters
    src_files = {"profiles_vmf.json": PROFILES_JSON, "cases.yaml": CASES_FILE}
    if g["source"] in ("gendata.bathy.build", "refine_200m", "subbox_of_refine_200m"):
        b = p["bathy"]
        src_files.update({k: DATA_EXTERNAL / b["ibcao_dir"] / v for k, v in bth.IBCAO_FILES.items()})
        src_files["coastline"] = DATA_EXTERNAL / b["coastline"]
    closure = "KPP" if use["kpp"] else ("GGL90" if use["ggl90"] else "constant (thesis)")
    actual = {
        "f0_s-1": float(par1["f0"]), "beta_m-1s-1": float(par1["beta"]),
        "profile": p["profile"], "T_column": [float(x) for x in T], "S_column": [float(x) for x in S],
        "eos": {k: par1[k.lower()] for k in ("eosType", "tAlpha", "sBeta", "rhoNil", "rhoConst", "tRef", "sRef")},
        "closure": closure, "packages": use, "surface": p["surface"],
        "wind": wind_info, "ice": ice_info, "tide": tide_info,
        "stages": plan, "edits": names, "deltaT_s": dt,
        "grid": {"nx": nx, "ny": ny, "nz": nz, "dx_m": dx, "bbox_utm": list(p["bbox_utm"]),
                 "wet_cells": int(g["wet"].sum()), "max_depth_m": float(g["depth"].max()),
                 "bathy_source": g["source"]},
        "sources": {k: {"path": str(v), "sha256": _sha256(v)} for k, v in src_files.items()},
    }
    meta = {"params": p, "actual": actual, "config": cfg, **_git_state()}

    out_dir.mkdir(parents=True, exist_ok=True)
    if "bathy.bin" in want:
        bth.write_bathy(out_dir / "bathy.bin", g["depth"], g["wet"])
        actual["bathy_sha256"] = _sha256(out_dir / "bathy.bin")
    for k, col in (("T.init.bin", T), ("S.init.bin", S)):
        if k in want:
            write_bin(out_dir / k, np.broadcast_to(col[:, None, None], (nz, ny, nx)))
    for k, arr in bins.items():
        if k in want:
            write_bin(out_dir / k, arr)
    for fname, text in rendered.items():
        (out_dir / fname).write_text(text)
    (out_dir / "case.json").write_text(json.dumps(meta, indent=2, default=float))
    return json.loads(json.dumps(meta, default=float))
