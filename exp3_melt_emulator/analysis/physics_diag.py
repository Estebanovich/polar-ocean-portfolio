#!/opt/anaconda3/envs/pyocean/bin/python
r"""Figure 3 of Exp. 3: physics diagnostics (descriptive, spec Sec. 7.3).

Usage::

    physics_diag.py [--fig-dir exp3_melt_emulator/figures]

Runs: the factorial runs with s = 1 (z_c = 250, 400, 550 m; Tbot = -0.5 ... 1.5 for
training, Tbot = 2.0 for the final test). Melt in m a-1 of freshwater, months 19-24.

(a) Integrated melt against Tbot per z_c: MITgcm against M1-M3. Integrated melt
    M = rho_fw sum_j m_j a_j (Burgard et al. 2022 Eq. 34), in Gt a-1 with rho_fw = 1000
    kg m-3 (AD16 Table 4, exp2 theory/melt.py RHO_FW) and a_j = dx dy. The method values
    are held-out predictions: the main LORO fold of each training run (cv_results.json)
    and the final test (final_test.json) at Tbot = 2.0. integrated_true of every entry
    is checked against the dataset.
(b) POST HOC relative to Exp. 2 (its pre-registered H1 used Tbot - T_f,GL and Tbot <= 1.5):
    power laws <m> = k F^n of the MITgcm whole-shelf mean melt, by ordinary least squares of
    ln <m> on ln F (as exp2 analysis/ensemble.py power_law_fit), with
      (b1) F = Tbot - T_f,GL, T_f,GL = lambda1 S + lambda2 + lambda3 p (AD16 Eq. 25,
           Table 4, as data.freezing_temperature) at S = 34.7 (Exp. 2 pd0 H1) and the
           deepest ice base of the run's shelf;
      (b2) F = the shelf-mean thermal forcing <TF> of data.cell_features (TF_mean), the
           nonlocal forcing of Jourdain et al. 2020 Eq. 1 / Favier et al. 2019 Eq. 5.
    Fitted per z_c over Tbot <= 1.5 (the Exp. 2 range) and over all Tbot including 2.0;
    (b2) also pooled over the 28 non-reference runs (s = 0.8, 1.0, 1.2).
(c) Melt sensitivities of Lambert & Burgard (2025, Eq. 1-3) on the final test, per method
    and z_c, with REF = Tbot 0.0, +1 = 1.0 and +2 = 2.0 (the test run). One model per method
    is used, the final-fold fit (trained on all 25 runs): its in-sample maps at Tbot 0.0
    and 1.0 and its test map at 2.0 (final_maps.npz). M is the shelf-mean melt;
      Eq. 1 dM/dT|(REF,+1) and dM/dT|(+1,+2), Eq. 2 d2M/dT2, Eq. 3 deep amplification
      DA = (dM10/dT)/(dM/dT) over (REF,+1); the same ratio over (+1,+2) is also given
      (an extension of Eq. 3, labelled DA_p1_p2). M10 = mean over the deepest 10 % of the
      draft (metrics.deep_amplification). M0 (a constant mean map) has no sensitivity and
      is left out. T here is Tbot, the temperature below the thermocline, which LB25 apply
      as their warming (+1, +2 degC at depth).
    REF and +1 are IN-SAMPLE (the final-fold model trained on those runs); only +2 is held
    out; the figure titles and fig3_physics.json ("caveat") say so. DA values with
    |dM/dT| < 0.5 m a-1 degC-1 are flagged (small denominator, post hoc display rule) and
    the (c2) axis is clipped at 12 with the clipped values annotated.

Writes fig3_physics.{png,pdf} (pfcommon.style.save) and fig3_physics.json.
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

from exp2_isomip_cavity.theory.melt import RHO_FW  # noqa: E402
from exp3_melt_emulator.analysis.figures import (  # noqa: E402
    FIG_DIR, METHOD_COLOURS, METHOD_LABELS, load_final_maps)
from exp3_melt_emulator.analysis.hypotheses import kind_of  # noqa: E402
from exp3_melt_emulator.data import cell_features, freezing_temperature  # noqa: E402
from exp3_melt_emulator.metrics import deep_amplification, lambert_sensitivities  # noqa: E402
from pfcommon import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

S_GL = 34.7                     # Exp. 2 pd0 H1: Sbot = 34.7 for T_f,GL (exp2 metrics.S_BOT_TF)
PHYS = ("M1", "M2", "M3")
LAMBERT_METHODS = tuple(f"M{i}" for i in range(1, 8))
LAMBERT_T = (0.0, 1.0, 2.0)     # REF, +1, +2 (Lambert & Burgard 2025, Sect. 2)
ZC_COLOURS = {250.0: style.OKABE_ITO[0], 400.0: style.OKABE_ITO[1], 550.0: style.OKABE_ITO[2]}
POST_HOC = "POST HOC relative to Exp. 2"
# Display caveats, chosen after seeing the values (labelled post hoc in the JSON):
DA_DENOM_MIN = 0.5      # m a-1 degC-1: |dM/dT| below this makes DA = (dM10/dT)/(dM/dT) a ratio
                        # of small numbers (~3x the 0.157 m a-1 noise floor of the melt maps)
DA_AXIS_MAX = 12.0      # (c2) y-axis top; larger DA values are drawn at the top and annotated
FIT_MIN_RANGE = 1.5     # a power-law fit whose forcing spans < 1.5x (max/min) is flagged
                        # poorly constrained
GT = 1e12                       # kg per Gt


def powerlaw_fit(forcing, melt, labels=None) -> dict:
    """m = k F^n by OLS of ln m on ln F (as exp2 analysis/ensemble.py power_law_fit);
    runs with F <= 0 or m <= 0 are excluded and listed (exp2 refuses instead). R^2 in log
    space. An empirical fit, not a theory: the form of Exp. 2 pd0 H1, as Holland et al. 2008
    Sect. 6 fit their modelled melt (exponent 1.61; literatura/Holland_2008*.md L923-925).
    ``forcing_ratio`` = max/min of the used forcing; ``poorly_constrained`` if < 1.5."""
    x, m = np.asarray(forcing, float), np.asarray(melt, float)
    labels = list(labels) if labels is not None else list(range(len(x)))
    ok = (x > 0) & (m > 0)
    out = {"n_runs": int(ok.sum()), "excluded": [l for l, o in zip(labels, ok) if not o]}
    if ok.sum() < 3 or np.ptp(x[ok]) == 0:
        return {**out, "n": None, "k": None, "R2_log": None,
                "reason": "fewer than 3 usable runs or constant forcing"}
    lx, ly = np.log(x[ok]), np.log(m[ok])
    n, lnk = np.polyfit(lx, ly, 1)
    res = ly - (n * lx + lnk)
    r2 = 1.0 - float((res**2).sum()) / float(((ly - ly.mean()) ** 2).sum())
    ratio = float(x[ok].max() / x[ok].min())
    return {**out, "n": float(n), "k": float(np.exp(lnk)), "R2_log": r2,
            "forcing_ratio": ratio, "poorly_constrained": ratio < FIT_MIN_RANGE}


def tf_gl(run) -> float:
    """T_f,GL (degC): AD16 Eq. 25 liquidus at S = 34.7 and the deepest shelf ice base."""
    return float(freezing_temperature(S_GL, np.nanmin(run.draft[run.shelf_mask])))


def shelf_mean(run, field) -> float:
    return float(np.mean(np.asarray(field, float)[run.shelf_mask]))


def to_gt(run, total) -> float:
    """sum_j m_j (m a-1, uniform cells) -> rho_fw a sum_j m_j in Gt a-1 (B22 Eq. 34)."""
    a = abs(float(run.x[1] - run.x[0]) * float(run.y[1] - run.y[0]))
    return float(total) * a * RHO_FW / GT


def _factorial(runs):
    return [r for r in runs if r.group != "reference" and r.s == 1.0]


def _zcs(runs):
    return sorted({r.z_c for r in _factorial(runs)})


def _preds(cv, final):
    """Held-out integrated predictions: main loro (training runs) and final (test runs)."""
    out = {}
    for e in cv["entries"]:
        if kind_of(e) == "loro" and e.get("variant", "main") == "main":
            out[(e["method"], e["case"])] = e
    for e in final["entries"]:
        out[(e["method"], e["case"])] = e
    return out


def part_a(runs, cv, final) -> dict:
    pr = _preds(cv, final)
    res = {}
    for zc in _zcs(runs):
        rs = sorted([r for r in _factorial(runs) if r.z_c == zc], key=lambda r: r.Tbot)
        row = {"Tbot": [r.Tbot for r in rs], "cases": [r.case for r in rs],
               "MITgcm": [to_gt(r, r.melt[r.shelf_mask].sum()) for r in rs]}
        for r in rs:                       # every method's entry against the dataset
            tot = float(r.melt[r.shelf_mask].sum())
            for (m, case), e in pr.items():
                if case == r.case and not np.isclose(e["metrics"]["integrated_true"], tot,
                                                     rtol=1e-9):
                    raise ValueError(f"{m} {r.case}: integrated_true "
                                     f"{e['metrics']['integrated_true']} != dataset {tot}")
        for m in PHYS:
            row[m] = [to_gt(r, pr[(m, r.case)]["metrics"]["integrated_pred"]) for r in rs]
        row["source"] = ["final" if r.split == "test" else "loro" for r in rs]
        res[f"{zc:.0f}"] = row
    res["units"] = "Gt a-1 (freshwater, rho_fw = 1000 kg m-3)"
    return res


def part_b(runs) -> dict:
    allr = [r for r in runs if r.group != "reference"]
    info = {r.case: {"Tbot": r.Tbot, "z_c": r.z_c, "s": r.s, "split": r.split,
                     "mean_melt": shelf_mean(r, r.melt), "T_f_GL": tf_gl(r),
                     "TF_shelf_mean": float(cell_features(r)["TF_mean"][r.shelf_mask][0])}
            for r in allr}
    fits = {"Tbot_minus_Tf_GL": {}, "TF_shelf_mean": {}}
    for zc in _zcs(runs):
        rs = sorted([r for r in _factorial(runs) if r.z_c == zc], key=lambda r: r.Tbot)
        for sub, sel in (("Tbot_le_1.5", [r for r in rs if r.Tbot <= 1.5]), ("all_Tbot", rs)):
            cs = [r.case for r in sel]
            m = [info[c]["mean_melt"] for c in cs]
            fits["Tbot_minus_Tf_GL"].setdefault(f"{zc:.0f}", {})[sub] = powerlaw_fit(
                [info[c]["Tbot"] - info[c]["T_f_GL"] for c in cs], m, cs)
            fits["TF_shelf_mean"].setdefault(f"{zc:.0f}", {})[sub] = powerlaw_fit(
                [info[c]["TF_shelf_mean"] for c in cs], m, cs)
    cs = sorted(info)
    fits["TF_shelf_mean"]["pooled_all_runs"] = powerlaw_fit(
        [info[c]["TF_shelf_mean"] for c in cs], [info[c]["mean_melt"] for c in cs], cs)
    return {"label": POST_HOC,
            "note": "descriptive; Exp. 2 H1 (pre-registered) fitted Tbot - T_f,GL over "
                    "Tbot <= 1.5 only; the shelf-mean TF forcing and the Tbot = 2.0 runs are "
                    "added here after the fact",
            "rule": "<m> = k F^n, OLS of ln <m> on ln F; <m> = MITgcm whole-shelf mean melt "
                    "(m a-1), months 19-24",
            "poorly_constrained_rule": "post hoc display flag: forcing max/min < 1.5 over the "
                                       "fitted runs (e.g. <TF> at z_c = 550 spans 0.18-0.26 "
                                       "degC), so n is weakly determined",
            "runs": info, "fits": fits}


def part_c(runs, maps) -> dict:
    res = {}
    for zc in _zcs(runs):
        by_t = {r.Tbot: r for r in _factorial(runs) if r.z_c == zc}
        if not all(t in by_t for t in LAMBERT_T):
            continue
        rr = [by_t[t] for t in LAMBERT_T]

        def field(m, r):
            if m == "MITgcm":
                return r.melt
            return maps[f"{m}__{r.case}__" + ("mean" if r.split == "test" else
                                              "mean_in_sample")]

        row = {}
        for m in ("MITgcm",) + LAMBERT_METHODS:
            f = [field(m, r) for r in rr]
            M = {t: shelf_mean(r, x) for t, r, x in zip(LAMBERT_T, rr, f)}
            s = lambert_sensitivities(M)                          # LB25 Eq. 1-2
            row[m] = {"M": [M[t] for t in LAMBERT_T], **s,
                      "DA_ref_p1": deep_amplification(f[0], f[1], rr[0], rr[1]),   # Eq. 3
                      "DA_p1_p2": deep_amplification(f[1], f[2], rr[1], rr[2])}
            row[m]["DA_ref_p1_small_denominator"] = abs(s["dMdT_ref_p1"]) < DA_DENOM_MIN
            row[m]["DA_p1_p2_small_denominator"] = abs(s["dMdT_p1_p2"]) < DA_DENOM_MIN
        row["cases"] = [r.case for r in rr]
        res[f"{zc:.0f}"] = row
    res["definition"] = ("Lambert & Burgard 2025 Eq. 1-3 with REF/+1/+2 = Tbot 0.0/1.0/2.0; "
                         "M = shelf-mean melt (m a-1); final-fold model of each method "
                         "(in-sample at 0.0 and 1.0, test at 2.0); DA_p1_p2 is Eq. 3 applied "
                         "to the (+1, +2) interval (extension)")
    res["caveat"] = (
        "REF (Tbot 0.0) and +1 (1.0) are IN-SAMPLE maps of the final-fold model (trained on "
        "those runs); only +2 (Tbot 2.0) is held out. For the flexible methods (M4-M7) the "
        "REF/+1 values largely reproduce the training data, so dM/dT|(REF,+1) and DA_ref_p1 "
        "are not an independent test of learned sensitivity, and dM/dT|(+1,+2) restates the "
        "final-test error. DA_p1_p2 is an extension of LB25 Eq. 3 (defined over REF,+1 only). "
        "DA flagged *_small_denominator where |dM/dT| < 0.5 m a-1 degC-1 (post hoc display "
        "rule): a ratio of small numbers, e.g. M4 at z_c = 550, DA_ref_p1 = 25.4 with "
        "dM/dT = 0.36, a small-denominator artefact.")
    res["da_small_denominator_threshold"] = DA_DENOM_MIN
    return res


def _figure(a, b, c, out_dir) -> None:
    style.apply()
    fig = plt.figure(figsize=(10.5, 6.8))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.40)
    zcs = [k for k in a if k != "units"]
    # (a)
    ax = fig.add_subplot(gs[0, 0])
    styles = {"M1": ":", "M2": "--", "M3": "-."}
    for k in zcs:
        col = ZC_COLOURS.get(float(k), "k")
        ax.plot(a[k]["Tbot"], a[k]["MITgcm"], "o-", color=col, lw=1.6, ms=4,
                label=f"MITgcm, z$_c$ = {k} m")
        for m in PHYS:
            ax.plot(a[k]["Tbot"], a[k][m], styles[m], color=col, lw=0.9)
    for m in PHYS:
        ax.plot([], [], styles[m], color="0.3", lw=0.9, label=METHOD_LABELS[m])
    ax.axvspan(1.75, 2.25, color="0.92", zorder=0)
    ax.text(2.0, 0.62, "final test", transform=ax.get_xaxis_transform(), ha="center",
            va="center", rotation=90, fontsize=7)
    ax.set_xlabel("$T_{bot}$ (°C)")
    ax.set_ylabel("integrated melt (Gt a$^{-1}$)")
    ax.set_title("(a) integrated melt, held-out predictions", loc="left", fontsize=8.5)
    top = max(max(a[k][m]) for k in zcs for m in ("MITgcm",) + PHYS)
    ax.set_ylim(top=1.75 * top)
    ax.legend(fontsize=6.0, frameon=False, loc="upper left", ncol=2, columnspacing=0.6,
              handlelength=1.8)
    # (b1), (b2)
    for j, (key, xl) in enumerate((("Tbot_minus_Tf_GL", "$T_{bot} - T_{f,GL}$ (°C)"),
                                   ("TF_shelf_mean", "shelf-mean TF $\\langle TF\\rangle$ (°C)"))):
        ax = fig.add_subplot(gs[0, 1 + j])
        for k in zcs:
            col = ZC_COLOURS.get(float(k), "k")
            cs = a[k]["cases"]
            x = np.array([b["runs"][cc]["Tbot"] - b["runs"][cc]["T_f_GL"] if j == 0
                          else b["runs"][cc]["TF_shelf_mean"] for cc in cs])
            y = np.array([b["runs"][cc]["mean_melt"] for cc in cs])
            test = np.array([b["runs"][cc]["split"] == "test" for cc in cs])
            ax.scatter(x[~test], y[~test], color=col, s=16, zorder=3)
            ax.scatter(x[test], y[test], facecolor="none", edgecolor=col, s=26, zorder=3)
            f = b["fits"][key][k]["all_Tbot"]
            f5 = b["fits"][key][k]["Tbot_le_1.5"]
            if f["n"] is not None and (x > 0).all():
                xx = np.geomspace(x.min(), x.max(), 30)
                ax.plot(xx, f["k"] * xx ** f["n"], color=col, lw=0.9)
                n5 = f"{f5['n']:.2f}" if f5["n"] is not None else "–"
                flag = " *" if f["poorly_constrained"] else ""
                ax.plot([], [], color=col, lw=0.9,
                        label=f"z$_c$ {k}: n = {f['n']:.2f} (≤ 1.5: {n5}){flag}")
        if j == 1:
            p = b["fits"][key]["pooled_all_runs"]
            others = [v for v in b["runs"].values() if v["s"] != 1.0]
            ax.scatter([v["TF_shelf_mean"] for v in others], [v["mean_melt"] for v in others],
                       marker="^", color="0.5", s=14, zorder=2, label="geometry runs (s ≠ 1)")
            if p["n"] is not None:
                allx = [v["TF_shelf_mean"] for v in b["runs"].values()]
                xx = np.geomspace(min(allx), max(allx), 30)
                ax.plot(xx, p["k"] * xx ** p["n"], color="0.5", lw=0.9, ls="--",
                        label=f"pooled, 28 runs: n = {p['n']:.2f}")
            if any(b["fits"][key][k]["all_Tbot"]["poorly_constrained"] for k in zcs):
                ax.plot([], [], " ", label="* poorly constrained (F spans < 1.5×)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xl)
        ax.set_ylabel("MITgcm shelf-mean melt (m a$^{-1}$)")
        ax.set_title(f"(b{j + 1}) $\\langle m\\rangle = k F^n$ (POST HOC vs Exp. 2)",
                     loc="left", fontsize=8.5)
        if j == 1:                                   # room above the data for the legend
            lo_, hi_ = ax.get_ylim()
            ax.set_ylim(lo_, hi_ * 25)
            from matplotlib.ticker import NullFormatter
            ticks = [0.2, 0.3, 0.5, 1.0, 2.0]
            ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
            ax.xaxis.set_minor_formatter(NullFormatter())
        ax.legend(fontsize=6.3, frameon=False, loc="upper left")
    # (c1) dM/dT, (c2) DA
    names = ("MITgcm",) + LAMBERT_METHODS
    zc_c = [k for k in c if isinstance(c[k], dict)]
    off = np.linspace(-0.22, 0.22, max(len(zc_c), 1))
    ax = fig.add_subplot(gs[1, 0:2])
    ax2 = fig.add_subplot(gs[1, 2])
    for o, k in zip(off, zc_c):
        col = ZC_COLOURS.get(float(k), "k")
        xs = np.arange(len(names)) + o
        ax.scatter(xs, [c[k][n]["dMdT_ref_p1"] for n in names], marker="o", facecolor="none",
                   edgecolor=col, s=20, zorder=3)
        ax.scatter(xs, [c[k][n]["dMdT_p1_p2"] for n in names], marker="o", color=col, s=20,
                   zorder=3, label=f"z$_c$ = {k} m")
        for key_, filled in (("DA_ref_p1", False), ("DA_p1_p2", True)):
            v = np.array([c[k][n][key_] for n in names])
            small = np.array([c[k][n][key_ + "_small_denominator"] for n in names])
            vis = np.minimum(v, DA_AXIS_MAX)
            kw = dict(color=col) if filled else dict(facecolor="none", edgecolor=col)
            inb = v <= DA_AXIS_MAX
            ax2.scatter(xs[inb], vis[inb], marker="o", s=20, zorder=3, **kw)
            ax2.scatter(xs[~inb], vis[~inb], marker="^", s=34, zorder=3, **kw)
            ax2.scatter(xs[small], vis[small], marker="x", color="0.25", s=14, lw=0.8,
                        zorder=4)
            for xi, vi, n in zip(xs[~inb], v[~inb], np.array(names)[~inb]):
                d = c[k][n]["dMdT_ref_p1" if key_ == "DA_ref_p1" else "dMdT_p1_p2"]
                ax2.annotate(f"{n}, z$_c$ {k}: DA = {vi:.1f} (clipped);\n"
                             f"dM/dT = {d:.2f}: small-denominator artefact",
                             (xi, DA_AXIS_MAX), xytext=(0, 7), textcoords="offset points",
                             ha="center", va="bottom", fontsize=6.2)
    ax.scatter([], [], marker="o", facecolor="none", edgecolor="0.3", s=20,
               label="(REF, +1): 0 → 1 °C, both in-sample")
    ax.scatter([], [], marker="o", color="0.3", s=20,
               label="(+1, +2): 1 → 2 °C, +2 held out")
    for a_ in (ax, ax2):
        a_.set_xticks(range(len(names)), names, rotation=45 if a_ is ax2 else 0)
        a_.axvspan(-0.5, 0.5, color="0.93", zorder=0)
        a_.set_xlim(-0.5, len(names) - 0.5)
    ax.set_ylabel("dM/dT (m a$^{-1}$ °C$^{-1}$)")
    ax.set_title("(c1) dM/dT, Lambert & Burgard 2025 Eq. 1; final-fold models: "
                 "REF, +1 in-sample, +2 held out", loc="left", fontsize=8.5)
    lo_, hi_ = ax.get_ylim()
    ax.set_ylim(lo_, hi_ + 0.3 * (hi_ - lo_))
    ax.legend(fontsize=7, frameon=False, ncol=3, loc="upper left", columnspacing=1.0,
              handletextpad=0.2)
    ax2.axhline(1.0, color="0.5", lw=0.7, ls=":")
    ax2.set_ylabel("deep amplification DA")
    ax2.set_ylim(0, DA_AXIS_MAX * 1.3)
    ax2.axhline(DA_AXIS_MAX, color="0.6", lw=0.6, ls="--")
    ax2.set_title("(c2) DA, LB25 Eq. 3: open = (REF, +1), in-sample;\n"
                  "filled = (+1, +2), extension of Eq. 3", loc="left", fontsize=8)
    ax2.set_xlabel("× |dM/dT| < 0.5 m a$^{-1}$ °C$^{-1}$: ill-conditioned;"
                   " ▲ clipped at 12", fontsize=6.5)
    style.save(fig, Path(out_dir) / "fig3_physics")
    plt.close(fig)


def make(runs, cv: dict, final: dict, maps: dict, out_dir: Path) -> dict:
    """Figure 3 and fig3_physics.json; ``maps`` in the final_maps.npz layout."""
    a, b, c = part_a(runs, cv, final), part_b(runs), part_c(runs, maps)
    _figure(a, b, c, out_dir)
    rec = {"figure": "fig3_physics", "a_integrated_vs_Tbot": a, "b_power_law": b,
           "c_lambert": c}
    (Path(out_dir) / "fig3_physics.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fig-dir", type=Path, default=FIG_DIR)
    a = ap.parse_args(argv)
    from exp3_melt_emulator.data import load_runs
    cv = json.loads((a.fig_dir / "cv_results.json").read_text())
    final = json.loads((a.fig_dir / "final_test.json").read_text())
    rec = make(load_runs(), cv, final, load_final_maps(final), a.fig_dir)
    fb = rec["b_power_law"]["fits"]
    print(json.dumps({k: {z: {s: v[s]["n"] for s in v} for z, v in d.items()
                          if z != "pooled_all_runs"} for k, d in fb.items()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
