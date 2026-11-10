#!/opt/anaconda3/envs/pyocean/bin/python
r"""POST HOC comparison of our Ocean0 and Ocean1 at 2 km with MITgcm-JPL (Nakayama, Menemenlis)
ISOMIP+ output (Esteban's requests, 2026-11-09 and 2026-11-10). Not a pre-registered
criterion; it changes no verdict.

Usage::

    jpl_compare.py [--cases ocean0 ocean1] [--out exp2_isomip_cavity/figures/jpl_comparison.json]

The JSON holds one entry per case ("ocean0", "ocean1"). Windows (``CASES``):

* ocean0: ours = pd0 'last6' (monthly records 19-24 of ocean0_2km_r21); JPL's Ocean0 file
  holds one year of calendar-month means (12 records), records 7-12 are taken as their last
  6 months;
* ocean1: ours = pd0 'year20' (records 229-240 of ocean1_2km_r21); JPL's Ocean1 file holds 20
  years of calendar-month means (240 records, time 0-19.92 years), its last 12 records
  (time >= 19 years) are taken as year 20.

Regions, each model with its own draft and shelf mask: whole shelf; deep = shelf and ice
base < -300 m (AD16 Eq. 37). JPL's shelf = cells where meltRate is not masked (in all 12
records) and iceDraft < 0. Means are area-weighted (both grids are uniform 2 km; ours uses
rA).

Quantities (time means over the window, then area means):

* melt [m a^-1 water equivalent]: ours -SHIfwFlx / rho_fw (metrics.melt_rate_m_per_a);
  JPL meltRate [m s^-1 w.e.] x 365 x 86400;
* u* [m s^-1]: ours SHIuStar; JPL frictionVelocity;
* thermal driving [degC]: AD16 output definition Tw - T_zd (boundary-layer potential
  temperature minus interface freezing temperature). We have no Tw/Tb diagnostic, so ours is
  inverted from the heat balance AD16 Eq. 24 with kappa_i = 0 (exact for our code):
  Tw - Tb = rho_fw m L / (rho_sw cw gammaT) = -SHIfwFlx L / (rhoConst cw SHIgammT), applied to
  the monthly means: <fw> / <gammT> is the gammaT-weighted mean of Tw - Tb, not the plain
  mean that JPL's thermalDriving is (if they averaged it per time step);
* circulation (both models, time-mean fields of the window): max |psi_bt| in the cavity
  (y > 50 km, x < 640 km) and near the front (630 < x < 650 km), max |psi(x, z)| for
  x < 640 km, and max - min of each for x < 640 km. Ours at cell corners
  (metrics.barotropic_streamfunction / overturning_streamfunction); JPL's
  barotropicStreamfunction / overturningStreamfunction as submitted (AD16 asks for them
  interpolated to the standard grid: x, y centres, 5 m levels), so the extrema are not
  taken at the same points;
* spatial correlation of the time-mean melt maps on the common 2 km grid: cell centres
  match exactly (JPL x, y = our XC, YC; checked); our extra row y = 81 km (closed wall) is
  dropped; cells where either model has shelf, melt = 0 off-shelf.
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

from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402
from exp2_isomip_cavity.theory import melt as ml  # noqa: E402

JPL_DIR = REPO / "exp2_isomip_cavity" / "data_external" / "osf_jpl"
JPL_FILE = JPL_DIR / "Ocean0_yosi_COM_MITgcm_JPL_new.nc"
JPL_LAST6 = slice(6, 12)                # records 7-12 (1-based) of their 1-year file
RUNS = Path("/Volumes/LaCie/portfolio_runs/exp2")
CASES = {
    "ocean0": {"run": RUNS / "ocean0_2km_r21", "window": "last6",
               "jpl": JPL_FILE, "jpl_n": 12, "jpl_slice": JPL_LAST6},
    "ocean1": {"run": RUNS / "ocean1_2km_r21", "window": "year20",
               "jpl": JPL_DIR / "Ocean1_yosi_COM_MITgcm_JPL_new.nc", "jpl_n": 240,
               "jpl_slice": slice(228, 240)},
}
X_CAV, Y_CAV, X_FRONT = mt.X_CAVITY_MAX_M, mt.Y_CAVITY_MIN_M, mt.X_FRONT_BAND_M
JPL_GAMMA_T = 0.0325                    # MITgcm_JPL_Summary.pdf, "Parameter values"
SEC_PER_YEAR = ml.SECONDS_PER_YEAR


def _wmean(f, w, sel):
    return float((f * w)[sel].sum() / w[sel].sum())


def psi_stats(bt, moc, xb, yb, xm) -> dict:
    """Circulation extrema [Sv] of time-mean streamfunctions: bt (ny', nx') at (xb, yb),
    moc (nz', nx') at xm."""
    cav = (yb[:, None] > Y_CAV) & (xb[None, :] < X_CAV)
    front = (xb[None, :] > X_FRONT[0]) & (xb[None, :] < X_FRONT[1]) & np.ones_like(cav)
    west = np.broadcast_to(xb[None, :] < X_CAV, bt.shape)
    mw = np.broadcast_to(xm[None, :] < X_CAV, moc.shape)
    def ext(f, sel, how):            # None if the region holds no finite value
        v = f[sel]
        v = v[np.isfinite(v)]
        if not v.size:
            return None
        return float(np.abs(v).max()) if how == "abs" else float(v.max() - v.min())

    return {"psi_bt_max_abs_cavity_Sv": ext(bt, cav, "abs"),
            "psi_bt_max_abs_front_Sv": ext(bt, front, "abs"),
            "psi_moc_max_abs_Sv": ext(moc, mw, "abs"),
            "bt_strength_max_minus_min_Sv": ext(bt, west, "range"),
            "moc_strength_max_minus_min_Sv": ext(moc, mw, "range")}


def jpl_fields(path: Path, rec: slice = JPL_LAST6, n_expected: int | None = None) -> dict:
    import netCDF4
    with netCDF4.Dataset(path) as d:
        t = np.asarray(d["time"][:], float)
        if n_expected is not None and t.size != n_expected:
            raise ValueError(f"{path.name}: {t.size} records, expected {n_expected}")
        melt = d["meltRate"][:]
        mask_any = np.ma.getmaskarray(melt).any(axis=0)
        draft = np.asarray(d["iceDraft"][:], dtype=float)
        shelf = ~mask_any & (draft < 0)

        def tmean(name):
            return np.ma.filled(d[name][rec].mean(axis=0), np.nan).astype(float)

        x, y = np.asarray(d["x"][:], float), np.asarray(d["y"][:], float)
        psi = psi_stats(tmean("barotropicStreamfunction") / mt.SV,
                        tmean("overturningStreamfunction") / mt.SV, x, y, x)
        idx = list(range(t.size))[rec]
        return {"x": x, "y": y, "time_years_window": [t[idx[0]] / SEC_PER_YEAR,
                                                     t[idx[-1]] / SEC_PER_YEAR],
                "records": [k + 1 for k in idx], "draft": draft, "shelf": shelf,
                "melt": tmean("meltRate") * SEC_PER_YEAR,
                "ustar": tmean("frictionVelocity"), "td": tmean("thermalDriving"), "psi": psi}


def ours_fields(run_dir: Path, window: str = "last6") -> dict:
    rd = mt.load_run(run_dir)                     # require_ok + R21 melt_law guard
    win = mt.calendar_window(rd, window)
    fw = rd.mean_field(mt._SHI, "SHIfwFlx", win)
    gm = rd.mean_field(mt._SHI, "SHIgammT", win)
    rho = float(rd.params["rho_ref"])
    with np.errstate(divide="ignore", invalid="ignore"):
        td = np.where(gm > 0, -fw * ml.LATENT_HEAT / (rho * ml.CW * gm), np.nan)
    return {"rd": rd, "x": rd.xc[0, :], "y": rd.yc[:, 0], "draft": rd.draft,
            "shelf": rd.shelf, "w": rd.rA, "melt": mt.melt_rate_m_per_a(fw),
            "ustar": rd.mean_field(mt._SHI, "SHIuStar", win), "td": td,
            "gammaT": float(rd.params["gammaT"]), "records": [k + 1 for k in win],
            "psi": psi_stats(mt.barotropic_streamfunction(rd, win),
                             mt.overturning_streamfunction(rd, win), rd.xg[0],
                             mt._corner_y(rd), rd.xg[0])}


def region_stats(f: dict, w) -> dict:
    out = {}
    for reg, sel in (("shelf", f["shelf"]), ("deep300", f["shelf"] & (f["draft"] < -300.0))):
        out[reg] = {"n_cells": int(sel.sum()),
                    "melt_m_per_a": _wmean(f["melt"], w, sel),
                    "ustar_m_s": _wmean(f["ustar"], w, sel),
                    "thermal_driving_degC": _wmean(f["td"], w, sel)}
    return out


def compare(run_dir: Path, jpl_path: Path, window: str = "last6", jpl_slice: slice = JPL_LAST6,
            jpl_n: int | None = None) -> dict:
    o, j = ours_fields(run_dir, window), jpl_fields(jpl_path, jpl_slice, jpl_n)
    ny = len(j["y"])
    if not (np.array_equal(o["x"], j["x"]) and np.array_equal(o["y"][:ny], j["y"])):
        raise ValueError("grids do not align: our XC/YC differ from JPL x/y")
    if o["shelf"][ny:].any():
        raise ValueError("our extra rows beyond JPL's y range hold shelf cells")
    ours = region_stats(o, o["w"])
    jpl = region_stats(j, np.ones_like(j["melt"]))
    so, sj = o["shelf"][:ny], j["shelf"]
    either = so | sj
    mo = np.where(so, o["melt"][:ny], 0.0)
    mj = np.where(sj, j["melt"], 0.0)
    r = float(np.corrcoef(mo[either], mj[either])[0, 1])
    both = so & sj
    r_both = float(np.corrcoef(o["melt"][:ny][both], j["melt"][both])[0, 1])
    rel = {reg: {k: (ours[reg][k] - jpl[reg][k]) / jpl[reg][k] for k in ours[reg]
                 if k != "n_cells"} for reg in ours}
    return {
        "label": "POST HOC (Esteban's requests 2026-11-09/10); not a pd0 criterion",
        "ours": {"run": str(run_dir), "window": window, "records_1based": o["records"],
                 "gammaT": o["gammaT"], **ours, "circulation": o["psi"]},
        "jpl": {"file": jpl_path.name, "records_1based": j["records"],
                "time_years_window_start": j["time_years_window"], "gammaT": JPL_GAMMA_T,
                **jpl, "circulation": j["psi"]},
        "rel_diff_circulation_ours_vs_jpl": {k: (o["psi"][k] - j["psi"][k]) / j["psi"][k]
                                             for k in o["psi"]
                                             if o["psi"][k] is not None and j["psi"][k]},
        "rel_diff_ours_vs_jpl": rel,
        "gammaT_ratio_ours_over_jpl": o["gammaT"] / JPL_GAMMA_T,
        "melt_map_correlation": {
            "r_either_shelf_zero_off": r, "n_cells_either": int(either.sum()),
            "r_both_shelf": r_both, "n_cells_both": int(both.sum()),
            "n_shelf_only_ours": int((so & ~sj).sum()), "n_shelf_only_jpl": int((sj & ~so).sum()),
            "grid": "identical 2 km cell centres (JPL x, y == our XC, YC); our row y = 81 km "
                    "(wall) dropped"},
        "definitions": {
            "thermal_driving_ours": "-<SHIfwFlx> L / (rhoConst cw <SHIgammT>) (AD16 Eq. 24, "
                                    "kappa_i = 0): gammaT-weighted Tw - Tb",
            "thermal_driving_jpl": "thermalDriving: Tw - T_zd (AD16 output spec), Tw over a "
                                   "20 m layer below the draft",
            "ustar_jpl": "u from a 20 m layer and the 4 horizontal neighbours, no-normal-flow "
                         "velocities at vertical ice faces excluded (summary PDF)",
            "ustar_ours": "SHELFICEboundaryLayer: top cell plus the part of the next one that "
                          "completes one cell thickness (metrics.boundary_layer_velocity_sq)"},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="post hoc comparison with MITgcm-JPL")
    ap.add_argument("--cases", nargs="+", default=list(CASES), choices=list(CASES))
    ap.add_argument("--out", type=Path,
                    default=REPO / "exp2_isomip_cavity" / "figures" / "jpl_comparison.json")
    a = ap.parse_args(argv)
    res = json.loads(a.out.read_text()) if a.out.exists() else {}
    if "ours" in res:                              # pre-2026-11-10 layout: Ocean0 only
        res = {"ocean0": res}
    for c in a.cases:
        cs = CASES[c]
        res[c] = compare(cs["run"], cs["jpl"], cs["window"], cs["jpl_slice"], cs["jpl_n"])
        print(c, json.dumps({k: res[c][k] for k in ("ours", "jpl", "melt_map_correlation")},
                            indent=1, default=float))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=2, default=float) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
