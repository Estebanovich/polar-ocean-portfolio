"""Generate all MITgcm input for one Exp. 1 case.

``make_case(name, out_dir, template_dir=None)`` reads ``cases.yaml``, writes
the binary inputs (big-endian float64 via ``pfcommon.mitgcm_io.write_bin``),
renders the namelist templates and writes ``case.json`` (actual parameters).

Template mechanism
------------------
Every regular, non-hidden file in ``template_dir`` (default
``exp1_fjord_plumes/config/<cfg>/input_templates/``, ``<cfg>`` = a_line2d,
b_twoplume3d or c_fjord3d from the name prefix a_/b_/c_) is rendered with
``string.Template.substitute`` and written under the same name to ``out_dir``.
Placeholders are ``${var}``; a literal dollar sign is ``$$``. An unknown
placeholder raises KeyError (no silent blanks). Available variables: every
key of the merged case parameters (cases.yaml config + case entry; e.g.
``${deltaT}``, ``${f0}``, ``${rho0}``, ``${tAlpha}``, ``${sBeta}``) plus the
derived ``${case}``, ``${nx}``, ``${ny}``, ``${nz}``, ``${nTimeSteps}``
(= round(endTime/deltaT)). Binary files are written under FIXED names that
templates reference literally:

  delX.bin delY.bin delZ.bin topog.bin T.init.bin S.init.bin
  R_icefront.bin icefrontlength.bin addmass.bin

(all (ny, nx) or (nz, ny, nx), big-endian float64; delX/delY/delZ 1-D.)
Afterwards every quoted ``*.bin`` name found in a rendered file must exist in
``out_dir``; otherwise RuntimeError.

Conventions: front column i = front_i (0) holds the ice wall on its western
face; x increases away from it. Flat bottom at -H. Discharge enters the bottom
wet cell(s) of the front column (a, b); for (c) a source Q0 at the bottom
and a sink of the same size at the plume neutral level (point_plume, half
cone) are written to addmass.bin (net zero; the entrained flux is internal to
the parameterisation, see report). Wind and shelf forcing are parameters only
(case.json), not files.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from string import Template

import numpy as np
import yaml

from pfcommon.mitgcm_io import write_bin

from . import profiles
from .grids import icefront_fields, stretched_x
from .sources import discharge_addmass, place_sources

PKG = Path(__file__).resolve().parents[1]
CFG_DIRS = {"a_": "a_line2d", "b_": "b_twoplume3d", "c_": "c_fjord3d"}


def load_cases(path: Path | None = None) -> dict:
    return yaml.safe_load(Path(path or PKG / "cases.yaml").read_text())


def _profile(strat: str, zc, p):
    if strat == "twolayer":
        z_int = p["z_int"] if "z_int" in p else p["H"] - p["H2"]
        return profiles.two_layer(zc, z_int, profiles.T_PW, profiles.S_PW,
                                  profiles.T_AW, profiles.S_AW, p["dz_int"])
    if strat == "sermilik":
        return profiles.continuous_sermilik(zc)
    raise ValueError(f"unknown stratification {strat!r}")


def _plume_levels(p, zc, T1d, S1d):
    """Neutral/detachment level of the (c) point plume (half cone)."""
    from ..theory.melt import MeltParams
    from ..theory.plume import LinearEOS, point_plume
    eos = LinearEOS(rho0=p["rho0"], tAlpha=p["tAlpha"], sBeta=p["sBeta"])
    zg = np.concatenate([[-p["H"]], zc[::-1]])  # base first, increasing
    Tg = np.interp(zg, zc[::-1], T1d[::-1])
    Sg = np.interp(zg, zc[::-1], S1d[::-1])
    sol = point_plume(zg, Tg, Sg, p["Q0"], p["T0_src"], p["S0_src"], MeltParams(),
                      "mitgcm", alpha=p["alpha"], eos=eos)
    z = sol.z_neutral if np.isfinite(sol.z_neutral) else sol.z_detach
    return float(z), float(sol.z_detach), float(sol.z_neutral)


def make_case(name: str, out_dir: Path, template_dir: Path | None = None,
              cases_file: Path | None = None) -> Path:
    """Write all inputs for case ``name`` into ``out_dir``; return ``out_dir``."""
    spec = load_cases(cases_file)
    if name not in spec["cases"]:
        raise KeyError(f"unknown case {name!r}")
    case = spec["cases"][name]
    cfg = case["config"]
    p = {**spec["configs"][cfg], **{k: v for k, v in case.items() if k != "config"}}
    p["case"], p["config"] = name, cfg
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- grid
    dz = p["dz"]
    nz = int(round(p["H"] / dz))
    if cfg == "a_line2d":
        dxv = stretched_x(p["dx0"], p["dxmax"], p["r"], p["L"])
        nx, ny, dy = len(dxv), p["ny"], p["dy"]
    else:
        nx, ny, dy = p["nx"], p["ny"], p["dy"]
        dxv = np.full(nx, p["dx"])
    dyv, dzv = np.full(ny, dy), np.full(nz, dz)
    zc = -(np.arange(nz) + 0.5) * dz
    yc = (np.arange(ny) + 0.5) * dy
    front_i = p["front_i"]
    p.update(nx=nx, ny=ny, nz=nz, nTimeSteps=int(round(p["endTime"] / p["deltaT"])))

    write_bin(out_dir / "delX.bin", dxv)
    write_bin(out_dir / "delY.bin", dyv)
    write_bin(out_dir / "delZ.bin", dzv)
    write_bin(out_dir / "topog.bin", np.full((ny, nx), -p["H"]))
    T1d, S1d = _profile(p["strat"], zc, p)
    write_bin(out_dir / "T.init.bin", np.broadcast_to(T1d[:, None, None], (nz, ny, nx)))
    write_bin(out_dir / "S.init.bin", np.broadcast_to(S1d[:, None, None], (nz, ny, nx)))
    R, length = icefront_fields(nx, ny, nz, dxv, dy, front_i, p["H"])
    write_bin(out_dir / "R_icefront.bin", R)
    write_bin(out_dir / "icefrontlength.bin", length)

    # ---- discharge
    hFacC = np.ones((nz, ny, nx))
    rho0, kb = p["rho0"], nz - 1
    info: dict = {}
    if cfg == "a_line2d":
        Q0 = p["q0"] * dy  # per cell: rho0 q0 dy (ny = 1)
        cells = [(kb, j, front_i) for j in range(ny)]
        add = discharge_addmass((nz, ny, nx), cells, Q0, rho0, hFacC=hFacC, front_i=front_i)
        info.update(Q0_total=Q0, source_cells=cells)
    elif cfg == "b_twoplume3d":
        if p["nsrc"] == 2:
            js, x0_actual = place_sources(yc, p["x0"], dy)
            info["x0_nominal"], info["x0_actual"] = p["x0"], x0_actual
        else:
            js = [ny // 2]
            info["x0_nominal"], info["x0_actual"] = None, None
        cells = [(kb, j, front_i) for j in js]
        Q0 = p["Q0"] * len(cells)
        add = discharge_addmass((nz, ny, nx), cells, Q0, rho0, hFacC=hFacC, front_i=front_i)
        info.update(Q0_total=Q0, source_cells=cells)
    else:
        add = np.zeros((nz, ny, nx))
        info.update(Q0_total=p["Q0"], source_cells=[], sink_cells=[])
        if p["Q0"] > 0:
            j = ny // 2
            src = [(kb, j, front_i)]
            z_sink, z_det, z_neu = _plume_levels(p, zc, T1d, S1d)
            ksink = int(np.argmin(np.abs(zc - z_sink)))
            snk = [(ksink, j, front_i)]
            add = (discharge_addmass((nz, ny, nx), src, p["Q0"], rho0, hFacC=hFacC,
                                     front_i=front_i)
                   + discharge_addmass((nz, ny, nx), snk, -p["Q0"], rho0, hFacC=hFacC,
                                       front_i=front_i, require_bottom=False))
            info.update(source_cells=src, sink_cells=snk, z_sink=z_sink,
                        z_detach=z_det, z_neutral=z_neu)
    write_bin(out_dir / "addmass.bin", add)

    # ---- templates
    if template_dir is None:
        template_dir = PKG / "config" / cfg / "input_templates"
    template_dir = Path(template_dir)
    if not template_dir.is_dir():
        raise FileNotFoundError(f"template dir {template_dir} does not exist")
    for f in sorted(template_dir.iterdir()):
        if f.is_file() and not f.name.startswith("."):
            (out_dir / f.name).write_text(Template(f.read_text()).substitute(p))
    missing = []
    for f in sorted(template_dir.iterdir()):
        if f.is_file() and not f.name.startswith("."):
            for line in (out_dir / f.name).read_text().splitlines():
                if line.lstrip().startswith(("#", "!")):
                    continue
                for ref in re.findall(r"'([^']+\.bin)'", line):
                    if not (out_dir / ref).is_file():
                        missing.append(f"{f.name}: {ref}")
    if missing:
        raise RuntimeError("namelists name files that were not written: "
                           + "; ".join(missing))

    meta = {"params": p, "dx_min": float(dxv.min()), "dx_max": float(dxv.max()),
            "sum_addmass_kg_s": float(add.sum()), "sum_abs_addmass_kg_s": float(np.abs(add).sum()),
            **info}
    (out_dir / "case.json").write_text(json.dumps(meta, indent=2, default=float))
    return out_dir
