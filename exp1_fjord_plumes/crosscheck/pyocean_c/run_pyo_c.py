#!/usr/bin/env python
"""pyOCEAN cross-check of the (c) fjord eastern open boundary.

Implements the configuration pre-registered in ``exp1_fjord_plumes/notes/pd1_c_fjord3d.md``, last section
("Verificación cruzada de la OB de (c) con pyOCEAN", 2026-11-06, approved by Esteban):

* pyOCEAN ``SigmaBaroclinic``, hydrostatic, ALE z* (``~/pyOCEAN/core/dynamics/sigma_baroclinic.py``);
* east OB: Flather (implicit theta form) on the depth mean, Orlanski + nudging on u', 10-cell sponge with the
  MITgcm OBCS profile, tau = (1200 s inner, 60 s boundary), tracers relaxed toward a time-varying exterior
  (``~/pyOCEAN/docs/manual/theory/open-boundaries.md``);
* grid and stratification of (c): 200 x 30 cells of 200 m (rows 0 and 29 land, the (c) side walls), 50 levels of
  10 m, f0 = 1.37e-4 s^-1, Sermilik profile, linear EOS with the (c) tAlpha/sBeta;
* Ah = Kh = 30 m^2 s^-1, Ar = Kr = 1e-5 m^2 s^-1;
* exterior = the (c) shelf series (pycnocline at 200 - 30 sin(2 pi t / 6 d) m), taken from
  ``exp1_fjord_plumes.gendata.make_case._c_forcing`` (the very OBEt/OBEs records MITgcm reads);
* source: Q = 100 m^3 s^-1 as surface runoff at the fjord head (``pyo_shelf_runoff``) or none (``pyo_shelf``).

Index order of pyOCEAN arrays: (k, j, i) with k = 0 at the BED; u-faces (N, ny, nx + 1).

Outputs (never committed) in ``--out`` (default /Volumes/LaCie/portfolio_runs/exp1/crosscheck_pyo/<case>/):
  day_XX.npz    daily snapshot: eta, u, T, S, h (layer thicknesses), time
  series.npz    hourly records (columns in ``SERIES_COLS``)
  status.json   finite / completed flags, cost, configuration actually used
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "VECLIB_NUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "5")

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

PYOCEAN = Path.home() / "pyOCEAN"
PORTFOLIO = Path(__file__).resolve().parents[3]
for _p in (str(PORTFOLIO), str(PYOCEAN)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from exp1_fjord_plumes.gendata import make_case  # noqa: E402
from exp1_fjord_plumes.crosscheck.pyocean_c.eval_pyo_c import T_MID, pycnocline_depth  # noqa: E402

RUNS = Path("/Volumes/LaCie/portfolio_runs/exp1/crosscheck_pyo")
CASES = {"pyo_shelf": 0.0, "pyo_shelf_runoff": 100.0}        # runoff Q [m^3 s^-1] (pre-registration)
MITGCM_ANALOGUE = "c_fjord_shelf"                            # shelf only, shelf_amp = 30 m (cases.yaml)
DT = 60.0                                                    # [s] = (c) deltaT (cases.yaml c_fjord3d)
SPONGE_CELLS = 10                                            # pre-registered
SPONGE_TAU = (1200.0, 60.0)                                  # (tau_inner, tau_bound) [s], pre-registered
AH = KH = 30.0                                               # [m^2 s^-1] pre-registered (= (c) viscAh, diffKhT/S)
AR = KR = 1.0e-5                                             # [m^2 s^-1] pre-registered (= (c) viscAr, diffKrT/S)
CD = 2.5e-3                                                  # (c) bottomDragQuadratic (data template), see README
S_REF = 34.0                                                 # [psu] (c) sRef (data template); tRef = 0
I_MID = 100                                                  # mid-fjord column x = nx/2 (criterion 4)
SERIES_COLS = ["t", "eta_mean", "eta_mean_absmax_hour", "Q_ob_out_hour", "cum_ob_out", "cum_runoff",
               "pyc_mid", "umax_ob3", "umax_int", "umax", "n_trunc"]


def case_params(Q: float = 0.0) -> dict:
    """Merged (c) parameters of ``c_fjord_shelf`` (cases.yaml config + case), as make_case.make_case merges them."""
    spec = make_case.load_cases()
    case = spec["cases"][MITGCM_ANALOGUE]
    p = {**spec["configs"][case["config"]], **{k: v for k, v in case.items() if k != "config"}}
    p["Q0"] = float(Q)
    return p


def c_inputs(Q: float = 0.0) -> dict:
    """Grid, initial profile and shelf OB records of (c), computed by the portfolio's own make_case code.

    Same lines as make_case.make_case (grid, depth with land rows/columns, hFacC) and then ``_c_forcing``, which
    writes OBEt/OBEs: records at t_n = (n + 1/2) dt_rec over one 6-day cycle (eesupp/src/get_periodic_interval.F),
    pycnocline at 200 m - A sin(2 pi t_n / P) (make_case.py, ``_c_forcing``).  With Q > 0 the plume columns of (c)
    (``info['plume_j']``, 6 cross-fjord cells centred on ny/2, ruling R18) are returned for the runoff.
    """
    p = case_params(Q)
    dz, nx, ny = p["dz"], p["nx"], p["ny"]
    nz = int(round(p["H"] / dz))
    zc = -(np.arange(nz) + 0.5) * dz                          # [m] level centres, surface first
    depth = np.full((ny, nx), p["H"])
    for jl in p.get("land_j_rows", []):
        depth[jl, :] = 0.0
    for il in p.get("land_i_cols", []):
        depth[:, il] = 0.0
    hFacC = ((np.arange(nz) * dz)[:, None, None] < depth[None]).astype(float)
    T1d, S1d = make_case._profile(p["strat"], zc, p)
    bins = {"T.init.bin": np.broadcast_to(T1d[:, None, None], (nz, ny, nx)),
            "S.init.bin": np.broadcast_to(S1d[:, None, None], (nz, ny, nx))}
    info: dict = {}
    make_case._c_forcing(p, bins, zc, T1d, S1d, hFacC, info)
    return dict(p=p, nz=nz, zc=zc, T1d=T1d, S1d=S1d, obT=np.asarray(bins["OBEt.bin"]),
                obS=np.asarray(bins["OBEs.bin"]), rec_dt=float(p["externForcingPeriod"]),
                cycle=float(p["shelf_period_days"] * 86400.0), plume_j=info.get("plume_j"))


def exterior(inp: dict, t_end: float):
    """t [s] -> {'east': OBData}: the (c) shelf records, linear in time and periodic as MITgcm's
    GET_PERIODIC_INTERVAL (records centred at (n + 1/2) dt_rec, pkg/obcs/obcs_fields_load.F:102-173), built as
    snapshots from t = -dt_rec/2 to t_end + dt_rec/2 for ``forcing.boundary_time.interpolated_boundary``
    (pyOCEAN, linear between snapshots).  u_normal = 0, eta/ubar = None (Flather toward a resting exterior,
    pyOCEAN open-boundaries.md "The exterior state").  zw: the (c) level interfaces, surface first [m]."""
    from boundaries.obcs import OBData
    from forcing.boundary_time import interpolated_boundary
    obT, obS, rdt = inp["obT"], inp["obS"], inp["rec_dt"]
    nrec, nz, ny = obT.shape
    zw = np.repeat((-np.arange(nz + 1) * inp["p"]["dz"])[:, None], ny, axis=1)
    u0 = np.zeros((nz, ny))
    n_hi = int(np.ceil(t_end / rdt)) + 1
    times, states = [], []
    for n in range(-1, n_hi + 1):
        r = n % nrec
        times.append((n + 0.5) * rdt)
        states.append({"east": OBData(side="east", T=obT[r].copy(), S=obS[r].copy(), u_normal=u0, zw=zw)})
    return interpolated_boundary(times, states)


def build(case: str, t_end: float, dt: float = DT):
    """pyOCEAN model, initial state and runoff Q of ``case`` (module docstring)."""
    from core.dynamics.sigma_baroclinic import SigmaBaroclinic
    from core.grid.horizontal import create_cartesian_grid
    from forcing.sources.glofas import RiverInput
    Q = CASES[case]
    inp = c_inputs(100.0)                     # plume_j of (c) (only the cells are used; the OB series has no Q)
    p, nz = inp["p"], inp["nz"]
    nx, ny = p["nx"], p["ny"]
    rho0 = p["rho0"]
    H = np.full((ny, nx), p["H"])
    for jl in p["land_j_rows"]:               # (c) side walls; pyOCEAN land = NaN depth
        H[jl, :] = np.nan
    # (c) land column nx-1 is NOT kept: pyOCEAN's open edge is the domain side (face i = nx); see README
    hg = create_cartesian_grid(nx=nx, ny=ny, dx=p["dx"], dy=p["dy"], f0=p["f0"])
    river = None
    if Q > 0:
        js = list(inp["plume_j"])
        origin = np.datetime64("2000-01-01T00:00", "ns")
        times = origin + np.array([0, int(np.ceil(t_end / 86400.0)) + 2], dtype="timedelta64[D]")
        river = RiverInput(cells=np.array([[j, p["front_i"]] for j in js]), times=times,
                           Q=np.full((2, len(js)), Q / len(js)), origin=origin)
    # Linear EOS: MITgcm rho = rhoNil (1 - tAlpha (T - tRef) + sBeta (S - sRef)) (model/src/find_rho_2d.F,
    # eosType LINEAR); pyOCEAN rho = rho0 - alpha (T - T0) + beta (S - S0) (Vallis 2006 eq. 1.156) with
    # alpha = rho0 tAlpha [kg m^-3 K^-1], beta = rho0 sBeta [kg m^-3 psu^-1].
    m = SigmaBaroclinic(
        hg, H, dt=dt, n_layers=nz, rho0=rho0, t_coef=rho0 * p["tAlpha"], t_ref=0.0,
        salinity=True, s_coef=rho0 * p["sBeta"], s_ref=S_REF,
        closure="constant", nu_const=AR, kv=KR, ah=AH, kh=KH, z0=None,
        bottom_drag="quadratic", cd=CD, hbbl=p["dz"],
        vertical_coordinate="ale", regrid="zstar", regrid_params=dict(dz_top=np.full(nz, p["dz"])),
        pgf="fv", rho_ref="auto",
        open_boundaries=("east",), sponge_cells=SPONGE_CELLS, sponge_tau=SPONGE_TAU,
        open_boundary_state=exterior(inp, t_end + 2 * dt), open_boundary_mode={"east": "radiating"},
        mix_diagnostic=False, runoff=river)
    T0 = np.broadcast_to(inp["T1d"][::-1, None, None], (nz, ny, nx)).copy()   # bed first
    S0 = np.broadcast_to(inp["S1d"][::-1, None, None], (nz, ny, nx)).copy()
    st = m.initial_state(np.zeros((ny, nx)), T0, S=S0)
    return m, st, Q, inp


def _wet(m) -> np.ndarray:
    return np.asarray(m.fs.active, bool)


def eta_mean(m, st) -> float:
    """Area-weighted domain-mean sea level over the wet cells [m]."""
    w = _wet(m)
    return float(np.sum((m.area * st.eta)[w]) / np.sum(m.area[w]))


def ob_out(m) -> float:
    """Net OUTWARD volume flux through the east open faces in the last step [m^3 s^-1]: the theta-weighted Flather
    flux the continuity equation used (``WetDryFreeSurface.last_open_Qu``, core/dynamics/wetting_drying.py:458-472;
    east = +x = model sign outward)."""
    return float(np.sum(m.fs.last_open_Qu[:, -1]))


def umax_cols(u: np.ndarray) -> tuple[float, float, float]:
    """max|u| [m s^-1] on the 3 u-face columns next to the OB (nx-2 .. nx, the last being the OB face), on the
    interior u-faces 10 .. 189, and overall (criterion 3; the MITgcm analogue is UVEL i = 196-198 vs 10-189,
    notes/pd1_b_twoplume3d.md)."""
    a = np.abs(np.nan_to_num(u))
    return float(a[:, :, -3:].max()), float(a[:, :, 10:190].max()), float(a.max())


def pyc_mid(m, st) -> float:
    """y-mean over the wet rows of the pycnocline depth (T = T_MID crossing, linear) at column I_MID [m]."""
    w = _wet(m)[:, I_MID]
    z = m.layer_heights(st.eta, st.h)[:, :, I_MID]           # (N, ny) layer centres [m], bed first
    d = [pycnocline_depth(st.T[:, j, I_MID], z[:, j]) for j in np.flatnonzero(w)]
    return float(np.mean(d))


def write_status(path: Path, **kw) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(kw, indent=1, default=float))
    tmp.replace(path)


def run(case: str, days: float, out: Path, threads: int = 5, dt: float = DT) -> dict:
    from core.model.threads import set_threads
    set_threads(threads)
    out.mkdir(parents=True, exist_ok=True)
    t_end = days * 86400.0
    n_total = int(round(t_end / dt))
    m, st, Q, inp = build(case, t_end, dt)
    steps_hour, steps_day = int(round(3600.0 / dt)), int(round(86400.0 / dt))
    cfg = dict(case=case, Q_runoff=Q, dt=dt, days=days, nx=m.H.shape[1], ny=m.H.shape[0], nz=m.n,
               sponge_cells=SPONGE_CELLS, sponge_tau_inner_bound=SPONGE_TAU, ah=AH, kh=KH, ar=AR, kr=KR, cd=CD,
               regrid="zstar", pgf="fv", rho_ref="auto", eos="linear", t_coef=inp["p"]["rho0"] * inp["p"]["tAlpha"],
               s_coef=inp["p"]["rho0"] * inp["p"]["sBeta"], s_ref=S_REF,
               runoff_cells=None if Q == 0 else [[int(j), 0] for j in inp["plume_j"]], threads=threads, T_mid=T_MID)
    status = dict(config=cfg, state="RUNNING", finite=True, completed=False, t_s=0.0, steps=0)
    spath = out / "status.json"
    write_status(spath, **status)
    np.savez(out / "day_00.npz", time=st.time, eta=st.eta, u=st.u, T=st.T, S=st.S, h=st.h)
    series = []
    cum_out = cum_rof = 0.0
    q_hour = 0.0
    eabs_hour = 0.0
    t0 = time.perf_counter()
    t_warm = None
    err = None
    n = 0
    try:
        for n in range(n_total):
            st = m.step(st)
            q = ob_out(m)
            cum_out += q * dt
            cum_rof += Q * dt
            q_hour += q * dt
            em = eta_mean(m, st)
            eabs_hour = max(eabs_hour, abs(em))
            if n == 2:
                t_warm = time.perf_counter()
            if (n + 1) % steps_hour == 0:
                fin = bool(np.isfinite(st.u).all() and np.isfinite(st.T).all() and np.isfinite(st.S).all()
                           and np.isfinite(st.eta).all())
                ob3, uin, uall = umax_cols(st.u)
                series.append([st.time, em, eabs_hour, q_hour / 3600.0, cum_out, cum_rof,
                               pyc_mid(m, st) if fin else np.nan, ob3, uin, uall, float(st.n_trunc)])
                q_hour, eabs_hour = 0.0, 0.0
                np.savez(out / "series.npz", data=np.array(series), columns=np.array(SERIES_COLS))
                if not fin:
                    raise FloatingPointError(f"non-finite fields at t = {st.time / 86400:.3f} d")
            if (n + 1) % steps_day == 0:
                d = (n + 1) // steps_day
                np.savez(out / f"day_{d:02d}.npz", time=st.time, eta=st.eta, u=st.u, T=st.T, S=st.S, h=st.h)
                el = time.perf_counter() - t0
                r = series[-1]
                print(f"{case} day {d:2d}  wall {el / 60:7.1f} min  eta_mean {r[1]: .3e} m  Q_ob {r[3]: 9.3f} "
                      f"m3/s  pyc_mid {r[6]:6.1f} m  umax ob3/int {r[7]:.3f}/{r[8]:.3f} m/s", flush=True)
                status.update(t_s=st.time, steps=n + 1, wall_s=el)
                write_status(spath, **status)
    except Exception as e:   # blow-up or a pyOCEAN abort (e.g. CFL truncations): recorded, not hidden
        err = f"{type(e).__name__}: {e}"
        traceback.print_exc()
    wall = time.perf_counter() - t0
    fin = bool(np.isfinite(st.u).all() and np.isfinite(st.T).all() and np.isfinite(st.S).all()
               and np.isfinite(st.eta).all()) and err is None
    done = err is None and n + 1 == n_total
    per = (time.perf_counter() - t_warm) / max(n - 2, 1) if t_warm else wall / max(n + 1, 1)
    status.update(state="OK" if done and fin else "FAILED", finite=fin, completed=done, t_s=st.time,
                  steps=(n + 1) if err is None else n, wall_s=wall, s_per_step=per, error=err,
                  n_trunc=int(st.n_trunc), n_limited=int(st.n_limited))
    write_status(spath, **status)
    if series:
        np.savez(out / "series.npz", data=np.array(series), columns=np.array(SERIES_COLS))
    print(f"{case}: {status['state']}  {st.time / 86400:.3f} d  wall {wall / 60:.1f} min  {per:.3f} s/step",
          flush=True)
    return status


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", choices=sorted(CASES), required=True)
    ap.add_argument("--days", type=float, default=30.0)
    ap.add_argument("--out", type=Path, default=None, help="output dir (default RUNS/<case>)")
    ap.add_argument("--threads", type=int, default=5)
    a = ap.parse_args()
    out = a.out or RUNS / a.case
    st = run(a.case, a.days, out, a.threads)
    raise SystemExit(0 if st["state"] == "OK" else 2)


if __name__ == "__main__":
    main()
