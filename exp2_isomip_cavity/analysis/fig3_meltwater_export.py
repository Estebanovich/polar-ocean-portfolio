#!/opt/anaconda3/envs/pyocean/bin/python
r"""Figure 3 of Exp. 2: meltwater export route of the 4 km ensemble, and the pd0 H2 verdict.

Usage::

    fig3_meltwater_export.py [--table figures/ensemble_table.json | --runs-root R]
                             [--run-suffix _r21] [--out-dir exp2_isomip_cavity/figures]

Panels (last 6 months, records 19-24; TRAC01 = meltwater volume fraction C, pd1 sect. 4):

  (a) y-z section of C at x_front (pd0 H2: per row, the open-ocean column just past the
      shelf, ``metrics.x_front_index``) for a cold (Tbot = -0.5) and a warm (Tbot = 1.5)
      core run with z_c = 400 m; dry cells blank; the run's D_out as a dashed line;
  (b) D_out vs Tbot per z_c (core runs, s = 1), with Spearman rho of each fully evaluated
      z_c; Tbot = 2.0 (Exp. 3 test data, not in H2) open markers; Ocean0 black star;
  (c) map of the vertically integrated tracer, sum_k C hFacC drF [m of meltwater] (static
      cell thickness; the r* factor, of order eta / H, is left out), of the warm z_c = 400 run,
      with the x_front columns marked: the export route.

The barotropic streamfunction is not used (its cavity maximum at 4 km is 2.2x that at 2 km,
pd1_calibracion.md R21 resolution test).

Outputs: fig3_meltwater_export.{png,pdf}, fig3_h2.json (verdict, gaps, N1-N3, caption).
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

from exp2_isomip_cavity.analysis import ensemble as en  # noqa: E402
from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402
from exp2_isomip_cavity.analysis.fig2_melt_sensitivity import ZC_COLOURS  # noqa: E402
from pfcommon import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

FIG_DIR = REPO / "exp2_isomip_cavity" / "figures"
SECTION_CASES = {"cold": "ens_T-0.5_zc400_s1.0", "warm": "ens_T1.5_zc400_s1.0"}
MAP_CASE = "ens_T1.5_zc400_s1.0"
CAPTION = ("Ensemble at 4 km horizontal resolution; last 6 months (records 19-24) of each "
           "2-year run. C = meltwater volume fraction (passive tracer, source -SHIfwFlx/rho_fw "
           "at the ice base, restored to 0 at x > 790 km). D_out = depth of the outgoing "
           "(u C > 0) tracer transport at x_front, transport weighted (pd0 H2). Open "
           "markers: Tbot = 2.0, Exp. 3 test data, not used in H2.")


def front_section(rd: mt.RunData, records) -> dict:
    """Time-mean C at x_front (pd0 H2 column per row), shape (nz, ny); dry cells masked."""
    C = rd.mean_field("dynDiag", "TRAC01", records)
    cols = mt.x_front_index(rd)
    sec = np.full((rd.nz, rd.ny), np.nan)
    dry = np.ones((rd.nz, rd.ny), bool)
    for j, i in enumerate(cols):
        if i < 0:
            continue
        sec[:, j] = C[:, j, i]
        dry[:, j] = rd.hfac[:, j, i] == 0
    return {"C": np.ma.masked_array(sec, dry | ~np.isfinite(sec)), "i_front": cols,
            "y_km": rd.yc[:, 0] / 1e3, "z_m": rd.rc, "yg_km": rd.yg[:, 0] / 1e3,
            "zf_m": rd.rf}


def column_tracer(rd: mt.RunData, records) -> np.ma.MaskedArray:
    """sum_k C hFacC drF [m of meltwater] of the time mean over ``records``; dry masked."""
    C = rd.mean_field("dynDiag", "TRAC01", records)
    col = (C * rd.hfac * rd.drF[:, None, None]).sum(axis=0)
    return np.ma.masked_array(col, ~rd.wet)


def _ok(r):
    return r["status"] == "OK" and r["D_out"] is not None


def make_figure(table: dict, h2: dict, sections: dict, map_run, out_stem: Path) -> dict:
    rs = en.rows(table)
    plotted = {"a": [], "b": {"core": [], "test": [], "reference": []}, "c": None}
    style.apply()
    fig = plt.figure(figsize=(7.0, 6.8))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.05], hspace=0.45, wspace=0.42,
                          left=0.09, right=0.93, top=0.95, bottom=0.14)
    secs = {k: front_section(rd, mt.LAST6) for k, (_, rd) in sections.items()}
    vmax = max(float(s["C"].max()) for s in secs.values()) if secs else 1.0
    cmap = plt.get_cmap("cmo.dense").copy()
    cmap.set_bad("0.8")
    by_case = {r["case"]: r for r in rs}
    for n, (k, (case, rd)) in enumerate(sections.items()):
        ax = fig.add_subplot(gs[0, n])
        s = secs[k]
        pm = ax.pcolormesh(s["yg_km"].tolist() + [s["yg_km"][-1] + np.diff(s["yg_km"])[-1]],
                           s["zf_m"], s["C"], cmap=cmap, vmin=0, vmax=vmax, shading="flat")
        r = by_case.get(case)
        if r is not None and _ok(r):
            ax.axhline(-r["D_out"], color="w", ls="--", lw=1.0)
            ax.text(0.02, 0.04, f"D_out = {r['D_out']:.0f} m", transform=ax.transAxes,
                    fontsize=7, color="k")
        T = rd.actual.get("Tbot")
        ax.set_ylim(float(s["zf_m"][-1]), 0.0)
        ax.set(xlabel="y (km)", ylabel="z (m)" if n == 0 else None,
               title=f"(a{n + 1}) C at x_front, Tbot = {T:g} °C")
        plotted["a"].append(case)
    if secs:
        fig.colorbar(pm, cax=fig.add_subplot(gs[0, 2]),
                     label="C (meltwater fraction), z_c = 400 m")

    # (b)
    axB = fig.add_subplot(gs[1, 0])
    for z in en.CORE_ZC:
        col = ZC_COLOURS[z]
        core = [r for r in en.core_rows(table, z) if _ok(r)]
        v = h2["per_z_c"][f"{z:.0f}"]
        lab = f"z_c = {z:.0f} m: " + (f"rho = {v['rho']:.2f}" if v.get("rho") is not None
                                       else "not evaluable" if v["verdict"] ==
                                       "NOT_EVALUABLE" else "rho undefined")
        axB.plot([r["Tbot"] for r in core], [r["D_out"] for r in core], "o-", color=col,
                 ms=4, lw=0.9, label=lab)
        plotted["b"]["core"] += [r["case"] for r in core]
        test = [r for r in rs if r["split"] == "test" and r["z_c"] == z and _ok(r)]
        axB.plot([r["Tbot"] for r in test], [r["D_out"] for r in test], "o", mfc="none",
                 mec=col, ms=5.5, mew=1.2)
        plotted["b"]["test"] += [r["case"] for r in test]
    ref = [r for r in rs if r["group"] == "reference" and _ok(r)]
    if ref:
        axB.plot([r["Tbot"] for r in ref], [r["D_out"] for r in ref], "*", color="k", ms=9,
                 label="Ocean0")
        plotted["b"]["reference"] = [r["case"] for r in ref]
    axB.axvspan(1.75, 2.25, color="0.9", lw=0, zorder=0)
    axB.invert_yaxis()
    axB.set(xlabel="Tbot (°C)", ylabel="D_out (m)", title="(b) outflow depth (pd0 H2)")
    axB.legend(fontsize=6.5, frameon=False, loc="best")

    # (c)
    axC = fig.add_subplot(gs[1, 1])
    if map_run is not None:
        case, rd = map_run
        col = column_tracer(rd, mt.LAST6)
        xg = np.append(rd.xg[0], rd.xg[0, -1] + rd.dxG[0, -1]) / 1e3
        yg = np.append(rd.yg[:, 0], rd.yg[-1, 0] + rd.dyG[-1, 0]) / 1e3
        cm = plt.get_cmap("cmo.dense").copy()
        cm.set_bad("0.8")
        pc = axC.pcolormesh(xg, yg, col, cmap=cm, shading="flat")
        fig.colorbar(pc, cax=fig.add_subplot(gs[1, 2]), label=r"$\int C\,dz$ (m)")
        cols = mt.x_front_index(rd)
        jj = np.nonzero(cols >= 0)[0]
        axC.plot(rd.xc[0, cols[jj]] / 1e3, rd.yc[jj, 0] / 1e3, "w.", ms=2.5)
        axC.contour(rd.xc[0] / 1e3, rd.yc[:, 0] / 1e3, rd.shelf.astype(float), [0.5],
                    colors="w", linewidths=0.6)
        T = rd.actual.get("Tbot")
        axC.set(xlabel="x (km)", ylabel="y (km)",
                title=f"(c) column meltwater, Tbot = {T:g} °C, z_c = 400 m")
        plotted["c"] = case
    fig.text(0.09, 0.01, CAPTION, ha="left", va="bottom", fontsize=6.5, wrap=True,
             color="0.25")
    paths = style.save(fig, out_stem)
    plt.close(fig)
    return {"paths": [str(p) for p in paths], "plotted": plotted}


def run(table: dict, out_dir: Path, *, sections: dict, map_run) -> dict:
    out_dir = Path(out_dir)
    h2 = en.h2_verdicts(table)
    fig = make_figure(table, h2, sections, map_run, out_dir / "fig3_meltwater_export")
    out = {"H2": h2, "N1_N3": en.n_summary(table), "caption": CAPTION,
           "window": "last 6 months, records 19-24 (pd0)",
           "panels": {"a": [c for c, _ in sections.values()],
                      "c": None if map_run is None else map_run[0]},
           "figure": fig["paths"]}
    (out_dir / "fig3_h2.json").write_text(json.dumps(out, indent=2) + "\n")
    return {**out, "plotted": fig["plotted"]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Exp. 2 Figure 3 and pd0 H2")
    ap.add_argument("--table", type=Path, default=None, help="ensemble table JSON")
    ap.add_argument("--runs-root", type=Path, default=en.DEFAULT_RUNS)
    ap.add_argument("--run-suffix", default="_r21")
    ap.add_argument("--out-dir", type=Path, default=FIG_DIR)
    a = ap.parse_args(argv)
    if a.table is not None and a.table.exists():
        tab = en.read_table(a.table)
    else:
        tab = en.collect(a.runs_root, a.run_suffix)
        if a.table is not None:
            en.write_table(tab, a.table)

    def load(case):
        try:
            return case, mt.load_run(a.runs_root / f"{case}{a.run_suffix}")
        except Exception as e:                    # noqa: BLE001 - a gap: panel left out
            print(f"WARNING: {case} not plotted: {e}", file=sys.stderr)
            return None

    sections = {k: v for k, v in ((k, load(c)) for k, c in SECTION_CASES.items()) if v}
    res = run(tab, a.out_dir, sections=sections, map_run=load(MAP_CASE))
    print(json.dumps({k: v["verdict"] for k, v in res["H2"]["per_z_c"].items()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
