#!/opt/anaconda3/envs/pyocean/bin/python
r"""Figures 1, 2 and 4 of Exp. 3 and the hypothesis verdicts (figures/hypotheses.json).

Usage::

    figures.py [--fig-dir exp3_melt_emulator/figures]

* Figure 1 ``fig1_method_errors``: error by method (rows) in the four groups of test runs,
  in range (loro, 25 main folds), z_c out of range (zc_out, 25), warm out of range
  (warm_out, 5) and the final test (Tbot = 2.0, 3). Top: rmse_map (Burgard et al. 2022
  Eq. 32); bottom: |integrated_rel_error| (B22 Eq. 33-34, one shelf). One point per run
  (open = run failing N1), the bar is the median. The grey band is the noise floor
  (figures/noise_floor.json; pd0: 0.157 m a-1 and 1.3 %), below which no difference is
  interpreted. Main variant only.
* Figure 2 ``fig2_maps_T2``: melt maps of the final test run at z_c = 400 m (Tbot = 2.0):
  MITgcm, the best of M1-M3 on that run (lowest rmse_map in final_test.json), M6 and M7 on
  one colour scale, and each prediction minus MITgcm. The maps are read from
  final_maps.npz (sha256 checked against final_test.json) by ``main`` and passed in.
* Figure 4 ``fig4_spread_error``: |integrated_rel_error| against the relative ensemble
  spread for M6 and M7 over the 33 extrapolation cases of H3c (zc_out, warm_out, final),
  with the H3c box (spread < 5 % and error > 20 %) and the 10 % spread threshold suggested
  by Ockenden et al. (2026, GMD 19:7979, conclusions; literatura/Ockenden_2026*.md L731).
  The H3c logic is analysis/hypotheses.py.

Every figure writes PNG + PDF through ``pfcommon.style.save`` and its numbers to a JSON
of the same stem. Figure functions take the JSON dicts; only Figure 2 takes maps.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp3_melt_emulator.analysis import hypotheses as hy  # noqa: E402
from pfcommon import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

FIG_DIR = REPO / "exp3_melt_emulator" / "figures"
METHODS = hy.METHODS
METHOD_LABELS = {
    "M0": "M0 mean map", "M1": "M1 linear, local", "M2": "M2 quadratic, local",
    "M3": "M3 quadratic, nonlocal", "M4": "M4 per-cell regression", "M5": "M5 GP on EOFs",
    "M6": "M6 deep ensemble", "M7": "M7 hybrid",
}
METHOD_COLOURS = {
    "M0": "#888888", "M1": style.OKABE_ITO[4], "M2": style.OKABE_ITO[0],
    "M3": style.OKABE_ITO[5], "M4": style.OKABE_ITO[2], "M5": style.OKABE_ITO[3],
    "M6": style.OKABE_ITO[1], "M7": style.OKABE_ITO[7],
}
FIG1_KINDS = {"loro": "in range (LORO)", "zc_out": "$z_c$ out of range",
              "warm_out": "warm out of range", "final": "final test ($T_{bot}$ = 2.0)"}
FIG1_XLIM = {"rmse_map": (1e-2, 2e2), "abs_ire": (1e-4, 3e1)}   # one x range per row
KIND_MARKERS = {"zc_out": "o", "warm_out": "s", "final": "*"}


def _write_json(rec: dict, stem: Path) -> None:
    Path(stem).with_suffix(".json").write_text(json.dumps(rec, indent=1) + "\n")


def _floors(noise: dict) -> dict:
    return {"rmse_map": float(noise["pooled_median_rms"]),
            "abs_ire": float(noise["median_abs_integrated_rel_diff"])}


# ---------------------------------------------------------------------------- Figure 1
def fig1_method_errors(cv: dict, final: dict, noise: dict, out_dir: Path) -> dict:
    style.apply()
    entries = hy.main_entries(cv, final)
    floors = _floors(noise)
    metrics = (("rmse_map", "RMSE of the melt map (m a$^{-1}$)"),
               ("abs_ire", "|integrated melt error| (fraction)"))
    fig, axes = plt.subplots(2, len(FIG1_KINDS), figsize=(10.0, 6.2), sharey=True)
    panels = {}
    for j, kind in enumerate(FIG1_KINDS):
        panels[kind] = {}
        for i, m in enumerate(METHODS):
            sel = sorted([e for e in entries if hy.kind_of(e) == kind and e["method"] == m],
                         key=lambda e: e["case"])
            if not sel:
                continue
            vals = {"rmse_map": np.array([e["metrics"]["rmse_map"] for e in sel]),
                    "abs_ire": np.array([abs(e["metrics"]["integrated_rel_error"])
                                         for e in sel])}
            n1 = np.array([bool(e.get("n1_pass", True)) for e in sel])
            panels[kind][m] = {
                "n_runs": len(sel),
                "median_rmse_map": float(np.median(vals["rmse_map"])),
                "median_abs_ire": float(np.median(vals["abs_ire"])),
                "per_run": [{"case": e["case"], "rmse_map": float(r), "abs_ire": float(a),
                             "n1_pass": bool(p)} for e, r, a, p in
                            zip(sel, vals["rmse_map"], vals["abs_ire"], n1)]}
            jit = np.linspace(-0.28, 0.28, len(sel)) if len(sel) > 1 else np.zeros(1)
            for k, (key, _) in enumerate(metrics):
                ax = axes[k, j]
                v = np.maximum(vals[key], 1e-4)          # log axis: |error| >= 1e-4 shown at 1e-4
                col = METHOD_COLOURS[m]
                ax.scatter(v[n1], i + jit[n1], s=11, color=col, alpha=0.75, lw=0, zorder=3)
                ax.scatter(v[~n1], i + jit[~n1], s=13, facecolor="none", edgecolor=col,
                           lw=0.9, zorder=3)
                med = np.median(v)
                ax.plot([med, med], [i - 0.38, i + 0.38], color="k", lw=1.8, zorder=4)
    for k, (key, lab) in enumerate(metrics):
        for j, kind in enumerate(FIG1_KINDS):
            ax = axes[k, j]
            ax.set_xscale("log")
            lo, hi = FIG1_XLIM[key]
            ax.axvspan(lo, floors[key], color="0.85", zorder=0)
            ax.set_xlim(lo, hi)
            ax.set_xlabel(lab, fontsize=8)
            ax.grid(axis="x", which="major", color="0.9", lw=0.5)
            if k == 0:
                ax.set_title(FIG1_KINDS[kind])
            ax.text(0.02, 0.985, f"({'abcdefgh'[k * len(FIG1_KINDS) + j]})",
                    transform=ax.transAxes, va="top", fontweight="bold")
    axes[0, 0].set_yticks(range(len(METHODS)), [METHOD_LABELS[m] for m in METHODS])
    axes[0, 0].set_ylim(len(METHODS) - 0.5, -0.5)
    handles = [Line2D([], [], marker="o", ls="", color="0.4", ms=4, label="one test run"),
               Line2D([], [], marker="o", ls="", mfc="none", mec="0.4", ms=4,
                      label="run failing N1"),
               Line2D([], [], color="k", lw=1.8, label="median"),
               plt.Rectangle((0, 0), 1, 1, color="0.85",
                             label="noise floor (0.157 m a$^{-1}$; 1.3 %)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    stem = Path(out_dir) / "fig1_method_errors"
    style.save(fig, stem)
    plt.close(fig)
    rec = {"figure": "fig1_method_errors", "variant": "main",
           "noise_floor": floors,
           "note": "|errors| below 1e-4 are drawn at 1e-4 (log axis); values in the JSON "
                   "are exact",
           "panels": panels}
    _write_json(rec, stem)
    return rec


# ---------------------------------------------------------------------------- Figure 2
def fig2_maps(final: dict, maps: dict, run, out_dir: Path) -> dict:
    """Maps of ``run`` (a final-test dataset Run) from ``maps`` (final_maps.npz layout)."""
    style.apply()
    import cmocean
    ent = {e["method"]: e for e in final["entries"] if e["case"] == run.case}
    phys = {m: ent[m]["metrics"]["rmse_map"] for m in ("M1", "M2", "M3")}
    best = min(phys, key=phys.get)
    show = [best, "M6", "M7"]
    sh = run.shelf_mask
    truth = np.where(sh, run.melt, np.nan)
    pred = {m: np.where(sh, np.asarray(maps[f"{m}__{run.case}__mean"], float), np.nan)
            for m in show}
    allv = np.concatenate([truth[sh]] + [p[sh] for p in pred.values()])
    vmin, vmax = float(min(0.0, allv.min())), float(allv.max())
    dmax = float(max(np.nanmax(np.abs(p - truth)) for p in pred.values()))
    rows, cols = np.nonzero(sh)
    r0, r1 = max(rows.min() - 1, 0), min(rows.max() + 2, sh.shape[0])
    c0, c1 = max(cols.min() - 1, 0), min(cols.max() + 2, sh.shape[1])
    dx, dy = float(run.x[1] - run.x[0]), float(run.y[1] - run.y[0])
    xe = (np.r_[run.x[c0:c1] - dx / 2, run.x[c1 - 1] + dx / 2]) / 1e3
    ye = (np.r_[run.y[r0:r1] - dy / 2, run.y[r1 - 1] + dy / 2]) / 1e3
    crop = lambda a: np.ma.masked_invalid(a[r0:r1, c0:c1])

    fig, axes = plt.subplots(4, 2, figsize=(7.6, 7.2), sharex=True, sharey=True)
    titles = ["MITgcm"] + [METHOD_LABELS[m] for m in show]
    fields = [truth] + [pred[m] for m in show]
    stats = {}
    for i, (t, f) in enumerate(zip(titles, fields)):
        ax = axes[i, 0]
        pm = ax.pcolormesh(xe, ye, crop(f), cmap=cmocean.cm.thermal, vmin=vmin, vmax=vmax)
        if i == 0:
            ax.set_title(f"(a) {t}\nshelf mean {np.nanmean(truth[sh]):.2f} m a$^{{-1}}$",
                         fontsize=8, loc="left")
            axes[0, 1].axis("off")
            axes[0, 1].text(0.0, 0.5,
                            f"{run.case}\nTbot = {run.Tbot:.1f} °C, z$_c$ = {run.z_c:.0f} m, "
                            f"s = {run.s:.1f}\nfinal test (evaluated once)\n"
                            f"best of M1–M3 by RMSE on this run: {best}",
                            transform=axes[0, 1].transAxes, va="center", fontsize=8.5)
            continue
        m = show[i - 1]
        mt = ent[m]["metrics"]
        stats[m] = {"rmse_map": mt["rmse_map"], "integrated_rel_error":
                    mt["integrated_rel_error"], "rmse_gl": mt["rmse_gl"]}
        ax.set_title(f"({'abcdefg'[2 * i - 1]}) {t}\nRMSE {mt['rmse_map']:.1f} m a$^{{-1}}$, "
                     f"integrated melt {100 * mt['integrated_rel_error']:+.1f} %",
                     fontsize=8, loc="left")
        axd = axes[i, 1]
        pd_ = axd.pcolormesh(xe, ye, crop(f - truth), cmap=cmocean.cm.balance,
                             vmin=-dmax, vmax=dmax)
        axd.set_title(f"({'abcdefg'[2 * i]}) {m} − MITgcm", fontsize=8, loc="left")
    for ax in axes.flat:
        ax.set_aspect("equal")
    for ax in axes[:, 0]:
        ax.set_ylabel("y (km)")
    for ax in axes[-1, :]:
        ax.set_xlabel("x (km)")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    cax1 = fig.add_axes([0.08, 0.035, 0.38, 0.018])
    fig.colorbar(pm, cax=cax1, orientation="horizontal",
                 label="melt rate (m a$^{-1}$ freshwater; one scale)")
    cax2 = fig.add_axes([0.57, 0.035, 0.38, 0.018])
    fig.colorbar(pd_, cax=cax2, orientation="horizontal",
                 label="prediction − MITgcm (m a$^{-1}$)")
    stem = Path(out_dir) / "fig2_maps_T2"
    style.save(fig, stem)
    plt.close(fig)
    rec = {"figure": "fig2_maps_T2", "case": run.case, "Tbot": run.Tbot, "z_c": run.z_c,
           "best_physics": best, "best_physics_rule": "lowest rmse_map on this run among "
           "M1-M3 (final_test.json)", "physics_rmse_map": phys, "methods": stats,
           "colour_scale": {"vmin": vmin, "vmax": vmax, "diff_absmax": dmax},
           "mitgcm_shelf_mean": float(np.nanmean(truth[sh]))}
    _write_json(rec, stem)
    return rec


# ---------------------------------------------------------------------------- Figure 4
def fig4_spread_error(cv: dict, final: dict, out_dir: Path,
                      expected: int | None = hy.EXPECTED["h3c"]) -> dict:
    style.apply()
    res = hy.h3c(hy.main_entries(cv, final), expected)
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0), sharex=True, sharey=True)
    allc = res["M6"]["cases"] + res["M7"]["cases"]
    xs = np.array([c["spread_rel"] for c in allc])
    ys = np.maximum(np.array([c["abs_ire"] for c in allc]), 1e-4)
    xlo, xhi = min(xs.min() / 1.5, 0.01), max(xs.max() * 1.5, 0.5)
    ylo, yhi = min(ys.min() / 2.0, 1e-3), max(ys.max() * 1.5, 1.0)
    for ax, m, lab in zip(axes, ("M6", "M7"), "ab"):
        r = res[m]
        ax.fill_between([xlo, hy.H3C_SPREAD], hy.H3C_ERROR, yhi, color="#f4c7c3", zorder=0,
                        label="H3c violation: spread < 5 %, error > 20 %")
        ax.axvline(hy.O26_SPREAD, color="0.3", ls="--", lw=1.0,
                   label="spread = 10 % (Ockenden et al. 2026)")
        ax.axhline(hy.H3C_ERROR, color="0.3", ls=":", lw=1.0, label="error = 20 %")
        for kind, mk in KIND_MARKERS.items():
            cs = [c for c in r["cases"] if c["kind"] == kind]
            ax.scatter([c["spread_rel"] for c in cs],
                       [max(c["abs_ire"], 1e-4) for c in cs], marker=mk,
                       s=60 if mk == "*" else 22, color=METHOD_COLOURS[m],
                       edgecolor="k", lw=0.4, zorder=3, label=f"{kind} ({len(cs)})")
        for q, v in enumerate(r["violations"]):          # alternate above / below the point
            ax.annotate(v["case"].removeprefix("ens_"), (v["spread_rel"], v["abs_ire"]),
                        xytext=(-4, 9) if q % 2 == 0 else (6, -12),
                        textcoords="offset points", fontsize=6.5,
                        ha="right" if q % 2 == 0 else "left",
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="0.3"))
        c = r["ockenden_10pct_confusion"]
        n_hi = c["spread_gt_10pct_error_gt_20pct"] + c["spread_le_10pct_error_gt_20pct"]
        ax.set_title(f"({lab}) {METHOD_LABELS[m]}: H3c {r['verdict']} "
                     f"({r['n_violations']} of {r['n_cases']} in the box)\n"
                     f"spread > 10 % in {c['spread_gt_10pct_error_gt_20pct']} of the {n_hi} "
                     "cases with error > 20 %", loc="left", fontsize=8.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
        ax.set_xlabel("relative ensemble spread")
    axes[0].set_ylabel("|integrated melt error| (fraction)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, -0.02, "spread = std (ddof = 0) of the 5 member integrated melts / "
             "|integrated melt of the ensemble mean|; main folds; H3c rule from the PD0 "
             "preregistration", ha="center", fontsize=7, color="0.3")
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    stem = Path(out_dir) / "fig4_spread_error"
    style.save(fig, stem)
    plt.close(fig)
    rec = {"figure": "fig4_spread_error", "variant": "main", **res}
    _write_json(rec, stem)
    return rec


# ---------------------------------------------------------------------------- main
def load_final_maps(final: dict) -> dict:
    path = Path(final["maps"])
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    if sha != final["maps_sha256"]:
        raise ValueError(f"{path}: sha256 {sha} != final_test.json {final['maps_sha256']}")
    with np.load(path, allow_pickle=False) as d:
        return {k: d[k] for k in d.files}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fig-dir", type=Path, default=FIG_DIR)
    a = ap.parse_args(argv)
    from exp3_melt_emulator.data import load_runs
    cv = json.loads((a.fig_dir / "cv_results.json").read_text())
    final = json.loads((a.fig_dir / "final_test.json").read_text())
    noise = json.loads((a.fig_dir / "noise_floor.json").read_text())
    rec = hy.write(a.fig_dir / "cv_results.json", a.fig_dir / "final_test.json",
                   a.fig_dir / "hypotheses.json")
    print("hypotheses:", json.dumps(rec["verdicts"]))
    fig1_method_errors(cv, final, noise, a.fig_dir)
    runs = {r.case: r for r in load_runs()}
    maps = load_final_maps(final)
    case = next(c for c in final["test"] if runs[c].z_c == 400.0)
    fig2_maps(final, maps, runs[case], a.fig_dir)
    fig4_spread_error(cv, final, a.fig_dir)
    print("wrote fig1, fig2, fig4 in", a.fig_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
