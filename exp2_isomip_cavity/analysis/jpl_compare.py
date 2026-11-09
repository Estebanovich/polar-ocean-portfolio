#!/opt/anaconda3/envs/pyocean/bin/python
r"""POST HOC comparison of our Ocean0 at 2 km with MITgcm-JPL (Nakayama, Menemenlis) ISOMIP+
Ocean0 output (Esteban's request, 2026-11-09). Not a pre-registered criterion; it changes no
verdict.

Usage::

    jpl_compare.py [--run /Volumes/LaCie/portfolio_runs/exp2/ocean0_2km_r21]
                   [--jpl exp2_isomip_cavity/data_external/osf_jpl/Ocean0_yosi_COM_MITgcm_JPL_new.nc]
                   [--out exp2_isomip_cavity/figures/jpl_comparison.json]

Windows: ours = pd0 'last6' (monthly records 19-24 of the 2-year run); JPL's file holds one
year of monthly means (12 records, ``time`` = start of each month), records 7-12 are taken
as their last 6 months.

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

JPL_FILE = REPO / "exp2_isomip_cavity" / "data_external" / "osf_jpl" / \
    "Ocean0_yosi_COM_MITgcm_JPL_new.nc"
JPL_LAST6 = slice(6, 12)                # records 7-12 (1-based) of their 1-year file
JPL_GAMMA_T = 0.0325                    # MITgcm_JPL_Summary.pdf, "Parameter values"
SEC_PER_YEAR = ml.SECONDS_PER_YEAR


def _wmean(f, w, sel):
    return float((f * w)[sel].sum() / w[sel].sum())


def jpl_fields(path: Path) -> dict:
    import netCDF4
    with netCDF4.Dataset(path) as d:
        melt = d["meltRate"][:]
        mask_any = np.ma.getmaskarray(melt).any(axis=0)
        draft = np.asarray(d["iceDraft"][:], dtype=float)
        shelf = ~mask_any & (draft < 0)

        def tmean(name):
            return np.ma.filled(d[name][JPL_LAST6].mean(axis=0), np.nan).astype(float)

        return {"x": np.asarray(d["x"][:], float), "y": np.asarray(d["y"][:], float),
                "time_s": np.asarray(d["time"][:], float).tolist(),
                "draft": draft, "shelf": shelf,
                "melt": tmean("meltRate") * SEC_PER_YEAR,
                "ustar": tmean("frictionVelocity"), "td": tmean("thermalDriving")}


def ours_fields(run_dir: Path) -> dict:
    rd = mt.load_run(run_dir)                     # require_ok + R21 melt_law guard
    win = mt.calendar_window(rd, "last6")
    fw = rd.mean_field(mt._SHI, "SHIfwFlx", win)
    gm = rd.mean_field(mt._SHI, "SHIgammT", win)
    rho = float(rd.params["rho_ref"])
    with np.errstate(divide="ignore", invalid="ignore"):
        td = np.where(gm > 0, -fw * ml.LATENT_HEAT / (rho * ml.CW * gm), np.nan)
    return {"rd": rd, "x": rd.xc[0, :], "y": rd.yc[:, 0], "draft": rd.draft,
            "shelf": rd.shelf, "w": rd.rA, "melt": mt.melt_rate_m_per_a(fw),
            "ustar": rd.mean_field(mt._SHI, "SHIuStar", win), "td": td,
            "gammaT": float(rd.params["gammaT"]), "records": [k + 1 for k in win]}


def region_stats(f: dict, w) -> dict:
    out = {}
    for reg, sel in (("shelf", f["shelf"]), ("deep300", f["shelf"] & (f["draft"] < -300.0))):
        out[reg] = {"n_cells": int(sel.sum()),
                    "melt_m_per_a": _wmean(f["melt"], w, sel),
                    "ustar_m_s": _wmean(f["ustar"], w, sel),
                    "thermal_driving_degC": _wmean(f["td"], w, sel)}
    return out


def compare(run_dir: Path, jpl_path: Path) -> dict:
    o, j = ours_fields(run_dir), jpl_fields(jpl_path)
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
        "label": "POST HOC (Esteban's request 2026-11-09); not a pd0 criterion",
        "ours": {"run": str(run_dir), "records_1based": o["records"], "gammaT": o["gammaT"],
                 **ours},
        "jpl": {"file": jpl_path.name, "records_1based": list(range(7, 13)),
                "gammaT": JPL_GAMMA_T, **jpl},
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
    ap = argparse.ArgumentParser(description="post hoc comparison with MITgcm-JPL Ocean0")
    ap.add_argument("--run", type=Path,
                    default=Path("/Volumes/LaCie/portfolio_runs/exp2/ocean0_2km_r21"))
    ap.add_argument("--jpl", type=Path, default=JPL_FILE)
    ap.add_argument("--out", type=Path,
                    default=REPO / "exp2_isomip_cavity" / "figures" / "jpl_comparison.json")
    a = ap.parse_args(argv)
    res = compare(a.run, a.jpl)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=2, default=float) + "\n")
    print(json.dumps({k: res[k] for k in ("ours", "jpl", "melt_map_correlation")},
                     indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
