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
face; x increases away from it. Flat bottom at -H. In (a), (b) the discharge
enters the bottom wet cell(s) of the front column. In (c) the plume is
buoyant-plume theory as addMass entrainment sinks and outflow sources whose
net is rho0 Q0, plus RBCS for the outflow T/S (``plume_representation`` =
"bpt_sink_source_rbcs"; ``_c_forcing``), and (c) also writes
rbcs_mask/rbcs_T/rbcs_S.bin, OBEt/OBEs.bin and, for wind cases, windx.bin
plus a second-stage namelist ``data.stage2`` (``_c_stages``).
Optional case key ``land_columns: [[j, i], ...]`` makes columns land, and
``land_j_rows: [j, ...]`` / ``land_i_cols: [i, ...]`` make whole rows (all i)
or whole columns (all j) land, e.g. the (b) side walls; hFacC for
the discharge checks is derived from the generated topography. A case/config with
``dx0`` uses the stretched grid of (a) (``dz`` from the case); ``discharge: line``
uses the (a) line source; ``template_config`` (optional) names the config whose
``input_templates`` are rendered (default: the case's config), so the PD2 runs
(ruling R20) can have their own SIZE.h/build but share templates.
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
from .plume_source_sink import plume_source_sink
from ..theory.melt import MeltParams
from .sources import discharge_addmass, place_sources

PKG = Path(__file__).resolve().parents[1]
CFG_DIRS = {"a_": "a_line2d", "b_": "b_twoplume3d", "c_": "c_fjord3d"}


STAGE2_NAME = "data.stage2"
# (c) data.diagnostics averaging / snapshot periods [s]: mean3h, snap1h
C_DIAG_PERIODS = (10800.0, 3600.0)


def is_two_stage(d: Path) -> bool:
    """True if directory ``d`` (case input or run dir) holds a second-stage
    namelist: the single test shared by run_campaign.py and analysis/c_runs.py."""
    return (Path(d) / STAGE2_NAME).exists()


def load_cases(path: Path | None = None) -> dict:
    return yaml.safe_load(Path(path or PKG / "cases.yaml").read_text())


def _profile(strat: str, zc, p):
    """b_twoplume3d: Sermilik profile centred at H - H2 = 100 m (R11); else 200 m."""
    if strat == "twolayer":
        z_int = p["z_int"] if "z_int" in p else p["H"] - p["H2"]
        return profiles.two_layer(zc, z_int, profiles.T_PW, profiles.S_PW,
                                  profiles.T_AW, profiles.S_AW, p["dz_int"])
    if strat == "sermilik":
        z_c = p["H"] - p["H2"] if "H2" in p else 200.0
        return profiles.continuous_sermilik(zc, z_c=z_c)
    raise ValueError(f"unknown stratification {strat!r}")


def _c_forcing(p, bins, zc, T1d, S1d, hFacC, info) -> tuple[np.ndarray, float]:
    """Config (c) forcing files (Task 7 Step 4); returns (addMass [kg s^-1],
    externForcingCycle [s] = shelf period).

    Plume (Q0 > 0): ``plume_source_sink`` (Cowton et al. 2015 sinks/outflow,
    ``n_out_levels`` outflow cells, 1 in (c), ruling R18) in the front column
    i = front_i, shared equally by ``n_plume_cols`` (6) cross-fjord cells
    centred on j = ny//2 (lateral spreading of the outflow, R18).
    Known class-(b) differences (notes/pd1_c_fjord3d.md): the plume forcing
    is fixed from the initial profile (Cowton recompute it every step), and
    sum(q) = Q0 includes the plume melt volume, removed with the sinks.
    PD3 option (iii): addMass with local-water sinks/sources plus RBCS in
    the outflow cells, mask = q tau/V and targets T_out, S_out, so the RBCS
    tendency equals q (T_out - T)/V (notes/pd3_addmass_sink.md).
    OB T/S (all cases): ``OBEt/OBEs.bin`` records at t_n = (n + 1/2) dt_rec
    (centred records, eesupp/src/get_periodic_interval.F) over one shelf
    period; pycnocline centre at 200 m - A sin(2 pi t/P) (A = shelf_amp).
    Wind (wind_tau > 0): ``windx.bin``, tau for the first wind_hours, then 0.
    """
    nz, ny, nx = hFacC.shape
    rho0, fi, dz = p["rho0"], p["front_i"], p["dz"]
    V = p["dx"] * p["dy"] * dz
    add = np.zeros((nz, ny, nx))
    mask = np.zeros((nz, ny, nx))
    Tr = np.array(bins["T.init.bin"], dtype=float)
    Sr = np.array(bins["S.init.bin"], dtype=float)
    if p["Q0"] > 0:
        q, T_out, S_out, sol = plume_source_sink(
            zc, T1d, S1d, p["Q0"], dz, T0=p["T0_src"], S0=p["S0_src"],
            alpha=p["alpha"], n_out_levels=p["n_out_levels"], return_solution=True)
        n = p["n_plume_cols"]
        js = list(range(ny // 2 - n // 2, ny // 2 - n // 2 + n))
        for j in js:
            if hFacC[-1, j, fi] == 0:
                raise ValueError(f"plume column (j={j}, i={fi}) is land")
        out = np.flatnonzero(q > 0)
        for j in js:
            add[:, j, fi] = rho0 * q / n
            mask[out, j, fi] = (q[out] / n) * p["rbcs_tau"] / V
            Tr[out, j, fi], Sr[out, j, fi] = T_out[out], S_out[out]
        if mask.max() > 1.0:
            raise ValueError(f"RBCS mask {mask.max():.3g} > 1: lower rbcs_tau")
        info.update(plume_representation="bpt_sink_source_rbcs", plume_j=js,
                    outflow_k=[int(k) for k in out], T_out=float(T_out[out[0]]),
                    S_out=float(S_out[out[0]]), z_neutral=float(sol.z_neutral),
                    Q_outflow=float(q[out].sum()), Q_entrained=float(-q[q < 0].sum()),
                    net_addmass_kg_s=float(add.sum()),
                    net_addmass_over_rho0_Q0=float(add.sum() / (rho0 * p["Q0"])),
                    outflow_turnover_s=float(V / (q[out[0]] / n)),
                    melt_volume_in_sinks_m3_s=_melt_volume(sol, float(zc[out[-1]] - dz / 2)),
                    sum_rule_note=("sum(addMass)/rho0 = Q0 exactly; the sinks remove "
                                   "entrainment + plume melt volume (melt_volume_in_sinks_m3_s) "
                                   "from the ambient (class (b), pd1_c_fjord3d.md); plume forcing "
                                   "fixed from the initial profile"))
    else:
        info.update(plume_representation="none", net_addmass_kg_s=0.0)
    info["Q0_total"] = p["Q0"]

    per = p["externForcingPeriod"]
    cyc = p["shelf_period_days"] * 86400.0
    nrec = int(round(cyc / per))
    if abs(nrec * per - cyc) > 1e-6:
        raise ValueError("shelf period must be a multiple of externForcingPeriod")
    t = (np.arange(nrec) + 0.5) * per
    eta = p["shelf_amp"] * np.sin(2.0 * np.pi * t / cyc)   # pycnocline uplift [m]
    obT = np.empty((nrec, nz, ny))
    obS = np.empty((nrec, nz, ny))
    for r in range(nrec):
        Tz, Sz = profiles.continuous_sermilik(zc, z_c=200.0 - eta[r])
        obT[r], obS[r] = Tz[:, None], Sz[:, None]
    bins["OBEt.bin"], bins["OBEs.bin"] = obT, obS
    _shelf_relaxation_zone(p, bins, hFacC, mask, Tr, Sr, obT, obS, info)
    # Stevens east OB: depth-mean u_ob = Q0/A_ob, same nrec records as OBEt/OBEs
    _stevens_east(p, bins, hFacC, p["Q0"], info, nrec=nrec)
    if p["wind_tau"] > 0:
        n_on = int(round(p["wind_hours"] * 3600.0 / per))
        w = np.zeros((nrec, ny, nx))
        w[:n_on] = p["wind_tau"]          # down-fjord (+x, away from the ice)
        bins["windx.bin"] = w
    return add, cyc


def sponge_zone_mask(p, hFacC) -> np.ndarray:
    """RBCS mask (nz, ny, nx) that reproduces the former OBCS T/S sponge exactly.

    pkg/obcs/obcs_sponge.F:540-565 (east, useLinearSponge=.FALSE.), for
    isl = 1..L (L = spongeThickness), cell i = OB_Ie - isl:

        trelax = ((L-isl) T_ob + isl theta)/L ,  lam = 1/(((L-isl) tau_b + isl tau_i)/L)
        d theta/dt = -lam ((L-isl)/L) (theta - T_ob)                 [K s-1]

    RBCS (pkg/rbcs/rbcs_add_tendency.F:88-100) gives -(mask/tauRelaxT)(theta - relaxT),
    so with relaxT = T_ob and tauRelaxT = rbcs_tau:

        mask(isl) = rbcs_tau (L-isl) / (L tau(isl)),  tau(isl) = ((L-isl) tau_b + isl tau_i)/L

    (0 at isl = L).  OB_Ieast=-2 -> Fortran OB_Ie = Nx-1 (obcs_readparms.F:671-672),
    0-based tracer column nx-2, so the zone is 0-based columns nx-2-isl.  Applied
    where the cell is wet (the old tendency has maskC).  [s]/[s] -> dimensionless.
    """
    nz, ny, nx = hFacC.shape
    L, tb, ti = int(p["sponge_thickness"]), p["sponge_tau_bound"], p["sponge_tau_inner"]
    m = np.zeros((nz, ny, nx))
    for isl in range(1, L + 1):
        tau = ((L - isl) * tb + isl * ti) / L
        m[:, :, nx - 2 - isl] = np.where(hFacC[:, :, nx - 2 - isl] > 0,
                                         p["rbcs_tau"] * (L - isl) / (L * tau), 0.0)
    return m


def _shelf_relaxation_zone(p, bins, hFacC, plume_mask, Tr, Sr, obT, obS, info) -> None:
    """RBCS shelf relaxation zone + plume-outflow cells; writes rbcs_mask/T/S.bin.

    With Stevens, OBEt/OBEs hold tracer TENDENCIES (obcs_calc_stevens.F:744-752,
    added in obcs_apply_ts.F:150-160) but the OBCS T/S sponge reads them as
    TARGET values (obcs_sponge.F:545-550), so the T/S sponge is off and the same
    relaxation is done with RBCS: mask from ``sponge_zone_mask``, targets = the
    shelf T/S record of OBEt/OBEs at the same time (broadcast in x), records at
    the same times as OBEt (rbcsForcingOffset = 0 and the same
    GET_PERIODIC_INTERVAL, pkg/rbcs/rbcs_fields_load.F:72-76).  Cases without
    shelf forcing (shelf_amp = 0) still get the zone, relaxing to the constant
    initial profile (as the old sponge did).  Plume cells (mask q tau/V,
    targets T_out/S_out) are repeated in every record.
    """
    nz, ny, nx = hFacC.shape
    nrec = obT.shape[0]
    zone = sponge_zone_mask(p, hFacC)
    if np.any((zone > 0) & (plume_mask > 0)):
        raise ValueError("RBCS shelf relaxation zone overlaps the plume-outflow cells")
    mask = plume_mask + zone
    if mask.max() > 1.0:
        raise ValueError(f"RBCS mask {mask.max():.3g} > 1: lower rbcs_tau")
    cols = np.flatnonzero(zone.max(axis=(0, 1)) > 0)
    rT = np.broadcast_to(Tr, (nrec, nz, ny, nx)).copy()
    rS = np.broadcast_to(Sr, (nrec, nz, ny, nx)).copy()
    rT[:, :, :, cols] = obT[:, :, :, None]
    rS[:, :, :, cols] = obS[:, :, :, None]
    bins["rbcs_mask.bin"], bins["rbcs_T.bin"], bins["rbcs_S.bin"] = mask, rT, rS
    info.update(rbcs_zone_columns_0based=[int(c) for c in cols],
                rbcs_zone_mask_max=float(zone.max()))


def _check_stevens_sponge(rendered: dict) -> None:
    """Raise if a rendered data.obcs combines useStevens* with a T/S sponge
    (OBCSsponge_Theta/Salt): Stevens stores tendencies in OBEt/OBEs
    (obcs_calc_stevens.F:744-752), the sponge would read them as targets
    (obcs_sponge.F:545-550)."""
    txt = rendered.get("data.obcs")
    if txt is None:
        return
    live = "\n".join(re.split(r"[#!]", l)[0] for l in txt.splitlines())
    flag = lambda k: re.search(rf"\b{k}\s*=\s*\.?(t|true)\b", live, re.I) is not None
    if any(flag(f"useStevens{s}") for s in ("North", "South", "East", "West")) and \
            (flag("OBCSsponge_Theta") or flag("OBCSsponge_Salt")):
        raise ValueError("data.obcs: useStevens* with OBCSsponge_Theta/Salt=.TRUE. is invalid "
                         "(OBEt/OBEs hold tendencies under Stevens); use RBCS")


def _stevens_east(p, bins, hFacC, Q0_total, info, T1d=None, S1d=None, nrec=1) -> None:
    """Eastern open-boundary files for Stevens (1990), shared by (b) and (c).

    Stevens (pkg/obcs/obcs_calc_stevens.F:195-228) takes the depth-mean normal
    velocity from ``OBEu`` (area-weighted with drF*hFacW over the water column)
    and adds the baroclinic part of the previous step's boundary velocity; the
    tracers relax toward ``OBEt/OBEs`` (TrelaxStevens/SrelaxStevens).  The
    Orlanski OB left the depth-mean u free (it drew ~800 m3/s IN in the (b)
    gate run b_3d_twolayer_x0_40 and drained (c): eta -11 m in 24 h).  Here
    the net discharge is prescribed to leave:

        u_ob = Q0_total / A_ob                      [m s-1, + = eastward = out]
        A_ob = dy * sum_{j,k} hFacW(i_ob,j,k) drF_k [m2]

    Index: OB_Ieast=-2 -> Fortran i_ob = Nx-1 (pkg/obcs/obcs_readparms.F:671-672),
    the u-point on the WEST face of cell Nx-1, i.e. 0-based column nx-2, and
    OB_Ie = i_ob is the i used for hFacW and OBEu (obcs_calc_stevens.F:182-199).
    Its hFacW = min(hFacC(i_ob), hFacC(i_ob-1)) = columns nx-2 and nx-3
    (0-based).  All cells are full (flat bottom, hFacC in {0,1}); dry rows
    (side walls) get 0.  The grid is the one actually generated (land rows).

    Files, (nz, ny) per record: OBEu = u_ob on wet faces, OBEv = 0 (tangential
    u is zero).  ``nrec`` is the number of records of each file:

    * nrec = 1 (b, no periodicExternalForcing): files are read once, record 1,
      at nIter0 (pkg/obcs/obcs_fields_load.F:347-363); OBEt/OBEs = the ambient
      profile ``T1d/S1d`` broadcast over y (written here).
    * nrec > 1 (c, periodicExternalForcing): every OB file with a name is read
      at records iRec0/iRec1 of GET_PERIODIC_INTERVAL
      (pkg/obcs/obcs_fields_load.F:102-105, 159-173), so OBEu/OBEv must have
      the same nrec = cycle/period records as the time-varying OBEt/OBEs
      (written by the caller); the records of OBEu/OBEv are identical.
    """
    nz, ny, nx = hFacC.shape
    i_ob = nx - 2                                           # 0-based u-point column
    hFacW = np.minimum(hFacC[:, :, i_ob], hFacC[:, :, i_ob - 1])   # (nz, ny)
    area_col = hFacW * p["dz"]                              # hFacW drF [m]
    A_ob = float(area_col.sum() * p["dy"])                  # [m2]
    if A_ob <= 0.0:
        raise ValueError("no wet open-boundary face for the Stevens eastern OB")
    u_ob = Q0_total / A_ob                                  # [m s-1]
    u2d = np.where(hFacW > 0, u_ob, 0.0)
    if nrec == 1:
        bins["OBEu.bin"], bins["OBEv.bin"] = u2d, np.zeros((nz, ny))
    else:
        bins["OBEu.bin"] = np.broadcast_to(u2d, (nrec, nz, ny))
        bins["OBEv.bin"] = np.zeros((nrec, nz, ny))
    if T1d is not None:
        bins["OBEt.bin"] = np.broadcast_to(T1d[:, None], (nz, ny))
        bins["OBEs.bin"] = np.broadcast_to(S1d[:, None], (nz, ny))
    info.update(ob_u_m_per_s=u_ob, ob_area_m2=A_ob, ob_column_0based=i_ob)


def _melt_volume(sol, z_top: float) -> float:
    """Meltwater volume flux of the half cone below z_top [m^3 s^-1]:
    int 2 b m rho_i/rho0 dz (wall contact width 2b, theory/plume.py)."""
    ok = np.isfinite(sol.Q) & (sol.z <= z_top)
    mp = MeltParams()
    mw = 2.0 * sol.b[ok] * sol.melt[ok] * mp.rho_i / mp.rho0
    return float(np.trapezoid(mw, sol.z[ok]))


def _c_stages(p) -> list[dict]:
    """Template variables per run stage of (c). Wind cases follow the thesis
    two-stage protocol (verification/BayIW_Rectan/README.md, "Two-stage run
    strategy"): stage 1 with the wind file, pickup at its end; stage 2
    restarts from it without the wind file (otherwise the cyclic forcing
    would repeat the pulse every cycle)."""
    dt, n_tot = p["deltaT"], int(round(p["endTime"] / p["deltaT"]))
    off = {"windLine": "# zonalWindFile: none (no wind in this case/stage)"}
    if p["wind_tau"] <= 0:
        return [dict(off, nIter0=0, nTimeSteps=n_tot, pChkptFreq=432000.0,
                     stage_label="single stage")]
    for per, what in ((dt, "deltaT"), *((x, "diagnostics period") for x in C_DIAG_PERIODS)):
        r = p["wind_stage1_s"] / per
        if abs(r - round(r)) > 1e-9:
            raise ValueError(f"wind_stage1_s={p['wind_stage1_s']} is not a multiple of the "
                             f"{what} {per} s (stage boundary must align)")
    n1 = int(round(p["wind_stage1_s"] / dt))
    if n_tot <= n1:   # short run inside stage 1: one stage, wind file on
        return [dict(nIter0=0, nTimeSteps=n_tot, pChkptFreq=432000.0,
                     windLine=" zonalWindFile='windx.bin',",
                     stage_label="single stage (wind on; run not longer than stage 1)")]
    return [dict(nIter0=0, nTimeSteps=n1, pChkptFreq=p["wind_stage1_s"],
                 windLine=" zonalWindFile='windx.bin',",
                 stage_label="stage 1 of 2 (wind on)"),
            dict(off, nIter0=n1, nTimeSteps=n_tot - n1, pChkptFreq=432000.0,
                 stage_label="stage 2 of 2 (wind off, restart from pickup)")]


def _git_state() -> dict:
    """Commit of the repository that generated the case (provenance in case.json)."""
    import subprocess
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PKG, capture_output=True,
                              text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=PKG,
                                    capture_output=True, text=True, check=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        return {"git_commit": "unknown", "git_dirty": None}
    return {"git_commit": head, "git_dirty": dirty}


def make_case(name: str, out_dir: Path, template_dir: Path | None = None,
              cases_file: Path | None = None, overrides: dict | None = None) -> Path:
    """Write all inputs for case ``name`` into ``out_dir``; return ``out_dir``.

    ``overrides`` (optional) replaces merged case parameters after cases.yaml
    is read, e.g. ``{"endTime": 7200.0}`` for a timed short run (Task 4); the
    values used are recorded in case.json.
    """
    spec = load_cases(cases_file)
    if name not in spec["cases"]:
        raise KeyError(f"unknown case {name!r}")
    case = spec["cases"][name]
    cfg = case["config"]
    p = {**spec["configs"][cfg], **{k: v for k, v in case.items() if k != "config"}}
    p.update(overrides or {})
    p["case"], p["config"] = name, cfg
    out_dir = Path(out_dir)

    # ---- grid
    dz = p["dz"]
    nz = int(round(p["H"] / dz))
    if "dx0" in p:  # stretched line-plume grid: (a) and its PD2 variants (R20)
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

    depth = np.full((ny, nx), p["H"])
    for (jl, il) in p.get("land_columns", []):  # optional land columns [j, i]
        depth[jl, il] = 0.0
    for jl in p.get("land_j_rows", []):  # whole rows y = j are land (side walls)
        depth[jl, :] = 0.0
    for il in p.get("land_i_cols", []):  # whole columns x = i are land
        depth[:, il] = 0.0
    bins = {"delX.bin": dxv, "delY.bin": dyv, "delZ.bin": dzv, "topog.bin": -depth}
    T1d, S1d = _profile(p["strat"], zc, p)
    bins["T.init.bin"] = np.broadcast_to(T1d[:, None, None], (nz, ny, nx))
    bins["S.init.bin"] = np.broadcast_to(S1d[:, None, None], (nz, ny, nx))
    R, length = icefront_fields(nx, ny, nz, dxv, dy, front_i, p["H"])
    bins["R_icefront.bin"], bins["icefrontlength.bin"] = R, length

    # ---- discharge; hFacC (0/1, flat full cells) from the generated bathymetry
    hFacC = ((np.arange(nz) * dz)[:, None, None] < depth[None]).astype(float)
    rho0, kb = p["rho0"], nz - 1
    info: dict = {}
    if p.get("discharge") == "line":  # (a) and its PD2 variants
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
        _stevens_east(p, bins, hFacC, Q0, info, T1d=T1d, S1d=S1d)
    else:
        add, p["externForcingCycle"] = _c_forcing(p, bins, zc, T1d, S1d, hFacC, info)
    bins["addmass.bin"] = add

    # ---- templates
    if template_dir is None:
        # ``template_config`` lets several builds (SIZE.h) share one template set (R20)
        template_dir = PKG / "config" / p.get("template_config", cfg) / "input_templates"
    template_dir = Path(template_dir)
    if not template_dir.is_dir():
        raise FileNotFoundError(f"template dir {template_dir} does not exist")
    tfiles = [f for f in sorted(template_dir.iterdir())
              if f.is_file() and not f.name.startswith(".")]
    stages = _c_stages(p) if cfg == "c_fjord3d" else [{}]
    p.update(stages[0])
    rendered = {f.name: Template(f.read_text()).substitute(p) for f in tfiles}
    if len(stages) > 1:  # wind cases: stage 2 restarts from the stage-1 pickup
        p2 = {**p, **stages[1]}
        rendered.update({STAGE2_NAME if f.name == "data" else f"{f.name}.stage2": Template(f.read_text()).substitute(p2)
                         for f in tfiles if f.name == "data"})
        info["stages"] = stages
    _check_stevens_sponge(rendered)
    missing = []
    for fname, text in rendered.items():
        for line in text.splitlines():
            if line.lstrip().startswith(("#", "!")):
                continue
            for ref in re.findall(r"'([^']+\.bin)'", line):
                if ref not in bins:
                    missing.append(f"{fname}: {ref}")
    if missing:  # nothing has been written yet
        raise RuntimeError("namelists name files that were not written: "
                           + "; ".join(missing))

    out_dir.mkdir(parents=True, exist_ok=True)
    for fname, arr in bins.items():
        write_bin(out_dir / fname, arr)
    for fname, text in rendered.items():
        (out_dir / fname).write_text(text)

    meta = {"params": p, **_git_state(),
            "dx_min": float(dxv.min()), "dx_max": float(dxv.max()),
            "sum_addmass_kg_s": float(add.sum()),
            **info}
    (out_dir / "case.json").write_text(json.dumps(meta, indent=2, default=float))
    return out_dir
