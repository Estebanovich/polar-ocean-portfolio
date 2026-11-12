"""Generate all MITgcm input for one Exp. 4 case of the thesis bay (bay_thesis, Task 3).

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
