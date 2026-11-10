#!/opt/anaconda3/envs/pyocean/bin/python
r"""Figure 2 of Exp. 2: melt sensitivity of the 4 km ensemble, and the pd0 H1 verdict.

Usage::

    fig2_melt_sensitivity.py [--melt-key {mean_melt_shelf,mean_melt_deep300}]
                             [--table figures/ensemble_table.json | --runs-root R --run-suffix _r21]
                             [--out-dir exp2_isomip_cavity/figures]

The H1 verdict (fig2_h1.json, key H1) always uses the whole-shelf mean melt (ruling R22,
2026-11-09, before any ensemble output was evaluated; pd0 addendum). The deep (< -300 m)
fit is written as H1_secondary_deep300, labelled secondary. ``--melt-key`` (default
mean_melt_shelf) only selects the <m> plotted in the panels.

Panels (last 6 months, records 19-24, of every run):

  (a) <m> vs Tbot, one colour per z_c (core runs, s = 1), with the pd0 H1 quadratic fit of
      each fully evaluated z_c; the Tbot = 2.0 runs (Exp. 3 test data, not fitted) as open
      markers; the Ocean0 reference (linear profile, Tbot = 1.0) as a black star;
  (b) heat map of <m>(Tbot, z_c) over the core grid plus the Tbot = 2.0 test column
      (hatched, outlined); gaps are left blank and labelled;
  (c) effect of the ice-draft scale s: <m> vs Tbot coloured by s (0.8, 1.0, 1.2), the
      geometry runs labelled with their z_c. DESCRIPTIVE: the Latin hypercube does not allow
      a clean fit of s.

Outputs: fig2_melt_sensitivity.{png,pdf}, fig2_h1.json (verdict, gaps, N1-N3, caption).
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
from pfcommon import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

FIG_DIR = REPO / "exp2_isomip_cavity" / "figures"
RES_JSON = FIG_DIR / "pd1_calibration_eval_r21.json"
ZC_COLOURS = {250.0: style.OKABE_ITO[0], 400.0: style.OKABE_ITO[1], 550.0: style.OKABE_ITO[2]}
S_COLOURS = {0.8: style.OKABE_ITO[4], 1.0: style.OKABE_ITO[7], 1.2: style.OKABE_ITO[5]}
MELT_LABEL = {"mean_melt_shelf": "whole-shelf mean melt",
              "mean_melt_deep300": "mean melt where ice base < -300 m"}


def _ok(r, key) -> bool:
    return r["status"] == "OK" and r[key] is not None


def caption(resolution: dict | None) -> str:
    """Caption. The resolution sentence states that the Ocean0 melt agreement between 4 and
    2 km is largely the calibration target (GammaT calibrated separately at each resolution
    to the same Ocean0 deep melt, pd0 K1), and cites the |psi_bt| and melt-map numbers."""
    res = resolution or {}
    rd = res.get("shelf_mean_melt", {}).get("rel_diff_4_vs_2")
    psi = res.get("psi_bt_max_cavity", {}).get("rel_diff_4_vs_2")
    r = res.get("melt_map_correlation", {}).get("r")
    parts = []
    if rd is not None:
        parts.append(f"Ocean0 whole-shelf melt differs by {100 * rd:+.2f} % (4 vs 2 km), but "
                     "this agreement is largely the calibration target (GammaT calibrated "
                     "separately at each resolution to the same Ocean0 deep melt, pd0 K1), "
                     "not independent evidence of convergence")
    if psi is not None:
        parts.append(f"the cavity |psi_bt| maximum differs by {100 * psi:+.0f} %")
    if r is not None:
        parts.append(f"the melt maps correlate at r = {r:.2f}")
    res_txt = (": " + "; ".join(parts)) if parts else ""
    return ("Ensemble at 4 km horizontal resolution (GammaT_cal,4km). 2 km vs 4 km "
            f"resolution test (notes/pd1_calibracion.md){res_txt}. pd0 H1 requires BOTH a "
            "quadratic fit with R^2 >= 0.98 AND a power law <m> = k (Tbot - T_f,GL)^n with "
            "n in [1.3, 2.0]; the curves in (a) are the quadratic fits, and n and the H1 "
            "verdict of each z_c are in its legend. Last 6 months (records 19-24) of each "
            "2-year run. Open markers / hatched column: Tbot = 2.0 runs, Exp. 3 test data, "
            "not used in H1. Panel (c) is descriptive.")


def h1_label(z: float, v: dict) -> str:
    """Legend entry of one z_c: R^2, n and the pd0 H1 verdict with its reason."""
    head = f"z_c = {z:.0f} m"
    if v["verdict"] == "NOT_EVALUABLE":
        return f"{head}: H1 not evaluable (gap)"
    q, p = v["quadratic"], v["power_law"]
    n = p["n"]
    txt = f"{head}: R² = {q['R2']:.3f}, " + ("n undefined" if n is None else f"n = {n:.2f}")
    if v["verdict"] == "PASS":
        return f"{txt}, H1 PASS"
    why = []
    if not q["pass"]:
        why.append(f"R² < {en.H1_R2_MIN}")
    if n is None:
        why.append("n undefined")
    elif n > en.H1_N_RANGE[1]:
        why.append(f"n > {en.H1_N_RANGE[1]}")
    elif n < en.H1_N_RANGE[0]:
        why.append(f"n < {en.H1_N_RANGE[0]}")
    return f"{txt}, H1 FAIL ({'; '.join(why)})"


def make_figure(table: dict, h1: dict, out_stem: Path, cap: str) -> dict:
    key = h1["melt_key"]
    rs = en.rows(table)
    plotted = {"a": {"core": [], "test": [], "reference": [], "fit_curves": [],
                     "legend": []},
               "b": {}, "c": {}}
    style.apply()
    fig = plt.figure(figsize=(7.0, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.42, wspace=0.32,
                          left=0.09, right=0.97, top=0.95, bottom=0.22)
    axA, axB, axC = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

    # (a)
    for z in en.CORE_ZC:
        col = ZC_COLOURS[z]
        core = [r for r in en.core_rows(table, z) if _ok(r, key)]
        v = h1["per_z_c"][f"{z:.0f}"]
        lab = h1_label(z, v)
        axA.plot([r["Tbot"] for r in core], [r[key] for r in core], "o", color=col, ms=4.5,
                 label=lab)
        plotted["a"]["legend"].append(lab)
        plotted["a"]["core"] += [r["case"] for r in core]
        if v["verdict"] != "NOT_EVALUABLE":
            q = v["quadratic"]
            t = np.linspace(min(en.CORE_TBOT), max(en.CORE_TBOT), 100)
            axA.plot(t, q["a"] * t**2 + q["b"] * t + q["c"], "-", color=col, lw=1.0)
            plotted["a"]["fit_curves"].append(f"{z:.0f}")
        test = [r for r in rs if r["split"] == "test" and r["z_c"] == z and _ok(r, key)]
        axA.plot([r["Tbot"] for r in test], [r[key] for r in test], "o", mfc="none", mec=col,
                 ms=5.5, mew=1.2)
        plotted["a"]["test"] += [r["case"] for r in test]
    ref = [r for r in rs if r["group"] == "reference" and _ok(r, key)]
    if ref:
        axA.plot([r["Tbot"] for r in ref], [r[key] for r in ref], "*", color="k", ms=9,
                 label="Ocean0 (linear profile)")
        plotted["a"]["reference"] = [r["case"] for r in ref]
    axA.plot([], [], "o", mfc="none", mec="0.3", label="Tbot = 2.0 (Exp. 3 test, not fitted)")
    axA.plot([], [], "-", color="0.3", lw=1.0, label="quadratic fit (H1 part 1 only)")
    axA.axvspan(1.75, 2.25, color="0.9", lw=0, zorder=0)
    axA.set(xlabel="Tbot (°C)", ylabel=r"$\langle m \rangle$ (m a$^{-1}$)",
            title=f"(a) {MELT_LABEL[key]}, s = 1"
                  + ("" if h1.get("role") == "criterion" else " [secondary]"))
    lo, hi = axA.get_ylim()
    axA.set_ylim(lo, lo + 1.55 * (hi - lo))              # room for the legend above the data
    axA.legend(fontsize=5.8, loc="upper left", frameon=False)

    # (b)
    Ts = list(en.CORE_TBOT) + [en.TEST_TBOT]
    grid = np.full((len(en.CORE_ZC), len(Ts)), np.nan)
    status = [["missing"] * len(Ts) for _ in en.CORE_ZC]
    for r in rs:
        if r["group"] not in ("factorial", "extrapolation") or r["z_c"] not in en.CORE_ZC:
            continue
        i = en.CORE_ZC.index(r["z_c"])
        j = [abs(r["Tbot"] - t) < 1e-9 for t in Ts].index(True)
        if _ok(r, key):
            grid[i, j] = r[key]
            status[i][j] = "ok"
        else:
            status[i][j] = "gap"
    xe = np.arange(len(Ts) + 1) - 0.5
    ye = np.arange(len(en.CORE_ZC) + 1) - 0.5
    pm = axB.pcolormesh(xe, ye, np.ma.masked_invalid(grid), cmap="cmo.amp", shading="flat")
    fig.colorbar(pm, ax=axB, label=r"$\langle m \rangle$ (m a$^{-1}$)", fraction=0.05)
    for i in range(len(en.CORE_ZC)):
        for j in range(len(Ts)):
            txt = "gap" if status[i][j] != "ok" else f"{grid[i, j]:.1f}"
            axB.text(j, i, txt, ha="center", va="center", fontsize=6.5,
                     color="k" if status[i][j] != "ok" else "w"
                     if grid[i, j] > np.nanmean(grid) else "k")
    axB.add_patch(plt.Rectangle((len(Ts) - 1.5, -0.5), 1, len(en.CORE_ZC), fill=False,
                                hatch="///", ec="0.3", lw=1.0))
    axB.set_xticks(range(len(Ts)), [f"{t:g}" for t in Ts])
    axB.set_yticks(range(len(en.CORE_ZC)), [f"{z:.0f}" for z in en.CORE_ZC])
    axB.set(xlabel="Tbot (°C)  [2.0: Exp. 3 test]", ylabel="z_c (m)",
            title="(b) ensemble grid, s = 1")
    plotted["b"] = {"grid": grid.tolist(), "Tbot": Ts, "z_c": list(en.CORE_ZC),
                    "test_column_Tbot": en.TEST_TBOT}

    # (c)
    for s, col in S_COLOURS.items():
        sel = [r for r in rs if r["split"] == "train" and r["group"] != "reference"
               and r["s"] is not None and abs(r["s"] - s) < 1e-9 and _ok(r, key)]
        if not sel:
            continue
        axC.plot([r["Tbot"] for r in sel], [r[key] for r in sel],
                 "o" if s == 1.0 else "D", color=col, ms=4 if s == 1.0 else 5.5,
                 alpha=0.6 if s == 1.0 else 1.0,
                 label=f"s = {s:.1f}" + (" (core)" if s == 1.0 else " (Latin hypercube)"))
        if s != 1.0:
            for r in sel:
                axC.annotate(f"{r['z_c']:.0f}", (r["Tbot"], r[key]),
                             xytext=(4, 3) if s < 1 else (4, -9),
                             textcoords="offset points", fontsize=6, color=col)
        plotted["c"][f"{s:.1f}"] = [r["case"] for r in sel]
    axC.set(xlabel="Tbot (°C)", ylabel=r"$\langle m \rangle$ (m a$^{-1}$)",
            title="(c) effect of the ice-draft scale s, DESCRIPTIVE (no fit; labels: z_c in m)")
    axC.legend(fontsize=7, frameon=False)
    fig.text(0.09, 0.01, cap, ha="left", va="bottom", fontsize=6.5, wrap=True, color="0.25")
    paths = style.save(fig, out_stem)
    plt.close(fig)
    return {"paths": [str(p) for p in paths], "plotted": plotted}


def run(table: dict, out_dir: Path, *, melt_key: str = en.H1_MELT_KEY,
        resolution: dict | None = None) -> dict:
    out_dir = Path(out_dir)
    h1 = en.h1_verdicts(table)                                     # criterion (R22)
    sec = en.h1_verdicts(table, melt_key="mean_melt_deep300")      # secondary
    cap = caption(resolution)
    shown = h1 if melt_key == en.H1_MELT_KEY else sec
    fig = make_figure(table, shown, out_dir / "fig2_melt_sensitivity", cap)
    out = {"H1": h1, "H1_secondary_deep300": sec, "N1_N3": en.n_summary(table),
           "caption": cap, "panel_c_label": "descriptive", "plotted_melt_key": melt_key,
           "melt_key_note": "R22 (2026-11-09 16:39 CST, before any ensemble output was "
                            "evaluated; pd0 addendum): H1 <m> = whole-shelf mean melt; the "
                            "deep (< -300 m) fit is secondary, not a criterion",
           "window": "last 6 months, records 19-24 (pd0)", "figure": fig["paths"]}
    (out_dir / "fig2_h1.json").write_text(json.dumps(out, indent=2) + "\n")
    return {**out, "plotted": fig["plotted"]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Exp. 2 Figure 2 and pd0 H1")
    ap.add_argument("--melt-key", default=en.H1_MELT_KEY, choices=en.MELT_KEYS,
                    help="<m> plotted in the panels (the H1 verdict is always whole-shelf)")
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
    res_json = json.loads(RES_JSON.read_text()).get("resolution_test") if RES_JSON.exists() \
        else None
    res = run(tab, a.out_dir, melt_key=a.melt_key, resolution=res_json)
    print(json.dumps({k: v["verdict"] for k, v in res["H1"]["per_z_c"].items()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
