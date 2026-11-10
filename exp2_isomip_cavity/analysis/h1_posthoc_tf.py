#!/opt/anaconda3/envs/pyocean/bin/python
r"""POST HOC (ruling R23, approved by Esteban 2026-11-10): H1 refitted against the forcing at
the ice base instead of the far-field Tbot. It does NOT change the pre-registered H1 verdict
(FAIL at z_c = 250, 400, 550; figures/fig2_h1.json).

Hypothesis tested: n > 2 against (Tbot - T_f,GL) arises because Tbot is not the temperature
the ice sees (the ramp's thermocline puts shallow ice in colder water).

Usage::

    h1_posthoc_tf.py [--table figures/ensemble_table.json] [--runs-root R]
                     [--run-suffix _r21] [--out-dir exp2_isomip_cavity/figures]

Forcings per core run (group factorial, s = 1, Tbot <= 1.5; status OK):

* TF_eff (a priori, inputs only): area-weighted mean over wet ice-shelf columns of
  T_ramp(z_d) - T_f(S_ramp(z_d), p(z_d)), z_d = the column's ice-base depth (icetopo.bin, after
  hFac snapping), T_ramp, S_ramp = gendata.profiles.ramp with the run's case.json Tbot, z_c,
  d; T_f = lam1 S + lam2 + lam3 p (AD16 Eq. 25, Table 4), p = rho_ref g |z_d|. Variants:
  whole shelf and deep region (base < -300 m, AD16 Eq. 37);
* TD_diag (model output): the thermal-driving estimate of analysis/jpl_compare.py
  (``ours_fields`` / ``region_stats``, imported, not reimplemented):
  Tw - Tb = -<SHIfwFlx> L / (rho_ref cw <SHIgammT>) (AD16 Eq. 24, kappa_i = 0) of the
  last-6-month means. Same caveat as there: <fw>/<gammT> is the gammaT-weighted mean of
  Tw - Tb, not its plain mean;
* reference: Tbot - T_f,GL, the pre-registered H1 variable.

Fits: <m>_shelf = k TF^n by OLS of ln m on ln TF, per z_c and pooled over the 15 runs; R^2
in log space; leave-one-out range of n. Runs with TF <= 0 are excluded and listed.
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
from exp2_isomip_cavity.analysis import jpl_compare as jc  # noqa: E402
from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402
from exp2_isomip_cavity.analysis.fig2_melt_sensitivity import ZC_COLOURS  # noqa: E402
from exp2_isomip_cavity.gendata import profiles as pr  # noqa: E402
from exp2_isomip_cavity.theory import melt as ml  # noqa: E402
from pfcommon import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

LABEL = "POST HOC (R23)"
H08_N = 1.61        # Holland et al. 2008 sect. 6 (literatura/Holland_2008*.md L924-925):
                    # power law fitted to their model melt, exponent 1.61 (r^2 = 0.994)
FIG_DIR = REPO / "exp2_isomip_cavity" / "figures"
MELT_KEY = "mean_melt_shelf"
FORCINGS = ("TF_eff_shelf", "TF_eff_deep300", "TD_diag_shelf", "TD_diag_deep300",
            "Tbot_minus_Tf_GL")


# ---- forcings -------------------------------------------------------------------------------

def tf_local(rd: mt.RunData, region: str = "shelf") -> tuple[np.ndarray, np.ndarray]:
    """A priori far-field thermal forcing at each column's ice base [degC] and the column
    areas [m^2] over ``region`` (see module doc). Inputs only."""
    prof = rd.actual.get("init_profile", rd.params.get("init_profile"))
    if prof != "ramp":
        raise ValueError(f"{rd.run_dir}: init_profile {prof!r}, TF_eff needs a ramp run")
    Tbot, z_c = float(rd.actual["Tbot"]), float(rd.actual["z_c"])
    d = float(rd.actual.get("d", rd.params.get("d")))
    sel = mt._mask(rd, region)                           # shelf | deep300 (AD16 Eq. 37)
    zd = rd.draft[sel]                                   # ice-base depth, negative [m]
    T, S = pr.ramp(zd, Tbot, z_c, d)
    p = float(rd.params["rho_ref"]) * float(rd.params["gravity"]) * np.abs(zd)
    tf = T - (ml.LAMBDA1 * S + ml.LAMBDA2 + ml.LAMBDA3 * p)   # AD16 Eq. 25, Table 4
    return tf, rd.rA[sel]


def tf_eff_a_priori(rd: mt.RunData, region: str = "shelf") -> float:
    """Area-weighted mean of ``tf_local`` [degC] (TF_eff, R23)."""
    tf, w = tf_local(rd, region)
    return float((tf * w).sum() / w.sum())


def td_diagnosed(run_dir: Path) -> dict:
    """jpl_compare thermal-driving estimate [degC], area mean over shelf and deep300."""
    f = jc.ours_fields(Path(run_dir), "last6")
    st = jc.region_stats(f, f["w"])
    return {reg: st[reg]["thermal_driving_degC"] for reg in ("shelf", "deep300")}


# ---- fits -----------------------------------------------------------------------------------

def _ols(x, m):
    n, lnk = np.polyfit(np.log(x), np.log(m), 1)
    y = np.log(m)
    res = y - (n * np.log(x) + lnk)
    r2 = 1.0 - float((res**2).sum()) / float(((y - y.mean()) ** 2).sum())
    return float(n), float(np.exp(lnk)), r2


def power_fit(x, m, labels) -> dict:
    """m = k x^n by OLS in log-log; x <= 0 (or m <= 0) excluded and listed; LOO range of n."""
    x, m = np.asarray(x, float), np.asarray(m, float)
    ok = (x > 0) & (m > 0)
    excl = [l for l, o in zip(labels, ok) if not o]
    xs, ms = x[ok], m[ok]
    if len(xs) < 3 or np.ptp(xs) == 0:
        return {"n": None, "k": None, "R2_log": None, "n_runs": int(len(xs)),
                "loo_n_range": None, "excluded_nonpositive_forcing": excl,
                "reason": "fewer than 3 usable runs or constant forcing"}
    n, k, r2 = _ols(xs, ms)
    loo = [_ols(np.delete(xs, i), np.delete(ms, i))[0] for i in range(len(xs))]
    return {"n": n, "k": k, "R2_log": r2, "n_runs": int(len(xs)),
            "loo_n_range": [min(loo), max(loo)], "excluded_nonpositive_forcing": excl}


# ---- analysis -------------------------------------------------------------------------------

def analyse(table: dict, runs_root, suffix: str = "_r21") -> dict:
    runs = []
    for z in en.CORE_ZC:
        for r in en.core_rows(table, z):
            if r["status"] != "OK" or r[MELT_KEY] is None:
                continue
            run_dir = Path(runs_root) / f"{r['case']}{suffix}"
            rd = mt.load_run(run_dir)
            td = td_diagnosed(run_dir)
            runs.append({"case": r["case"], "Tbot": r["Tbot"], "z_c": z, "melt": r[MELT_KEY],
                         "TF_eff_shelf": tf_eff_a_priori(rd, "shelf"),
                         "TF_eff_deep300": tf_eff_a_priori(rd, "deep300"),
                         "TD_diag_shelf": td["shelf"], "TD_diag_deep300": td["deep300"],
                         "Tbot_minus_Tf_GL": r["Tbot"] - r["T_f_GL"]})
    fits = {}
    for f in FORCINGS:
        fits[f] = {}
        for key, sel in [(f"{z:.0f}", [u for u in runs if u["z_c"] == z]) for z in en.CORE_ZC] \
                + [("pooled", runs)]:
            fits[f][key] = power_fit([u[f] for u in sel], [u["melt"] for u in sel],
                                     [u["case"] for u in sel])
    return {"label": LABEL, "melt_key": MELT_KEY, "runs": runs, "fits": fits,
            "H1_preregistered_verdict_unchanged": "FAIL",
            "note": "post hoc; does not change the pre-registered H1 verdict (FAIL at every "
                    "z_c, figures/fig2_h1.json)",
            "definitions": {
                "TF_eff": "area mean over wet ice-shelf columns of T_ramp(z_d) - T_f(S_ramp(z_d), "
                          "rho_ref g |z_d|), AD16 Eq. 25; inputs only (icetopo.bin, case.json)",
                "TD_diag": "jpl_compare: -<SHIfwFlx> L / (rho_ref cw <SHIgammT>), AD16 Eq. 24 "
                           "kappa_i = 0, last 6 months; gammaT-weighted mean of Tw - Tb",
                "Tbot_minus_Tf_GL": "pre-registered H1 variable (reference)",
                "TD_diag_caveat": "partly circular: TD_diag is computed from the melt itself "
                                  "(m = c gammaT TD_diag), so its exponent is 1 + d ln<gammaT>"
                                  " / d ln TD_diag, i.e. it measures the u* response, not an "
                                  "independent forcing-melt relation",
                "fit": "OLS of ln <m>_shelf on ln forcing; R2 in log space; LOO = leave one "
                       "run out"}}


# ---- R24 part 2: quadratic predictors of Favier et al. (2019) -----------------------------

R24_LABEL = "POST HOC (R24)"
QUAD = ("P_nonlocal", "P_local")


def quad_predictors(rd: mt.RunData) -> dict:
    """Shelf-area means of the quadratic melt predictors [degC^2], TF = ``tf_local``:

    * P_nonlocal = < TF_loc <TF> >: Favier et al. 2019 Eq. 5 (literatura/Favier_2019*.md
      L491-497, "quadratic, local and non-local"; Burgard et al. 2022 Eq. 18, L687-702,
      "quadratic-semilocal"). Since <TF> is one number per run, P_nonlocal = <TF>^2
      exactly, i.e. TF_eff^2 (R23): a power law in it has q = n(TF_eff) / 2;
    * P_local = < TF_loc |TF_loc| >: Favier 2019 Eq. 4 (L447-455, "quadratic, local");
      Burgard 2022 Eq. 14 (L648-654, "quadratic-local"). Not reducible to TF_eff.
    Melt ~ k P with q = 1 is what the quadratic family predicts."""
    tf, w = tf_local(rd, "shelf")
    mean = float((tf * w).sum() / w.sum())
    return {"P_nonlocal": float((tf * mean * w).sum() / w.sum()),
            "P_local": float((tf * np.abs(tf) * w).sum() / w.sum()), "TF_eff_shelf": mean}


def analyse_quad(table: dict, runs_root, suffix: str = "_r21") -> dict:
    runs = []
    for z in en.CORE_ZC:
        for r in en.core_rows(table, z):
            if r["status"] != "OK" or r[MELT_KEY] is None:
                continue
            rd = mt.load_run(Path(runs_root) / f"{r['case']}{suffix}")
            runs.append({"case": r["case"], "Tbot": r["Tbot"], "z_c": z, "melt": r[MELT_KEY],
                         **quad_predictors(rd)})
    fits = {}
    for f in QUAD:
        fits[f] = {}
        for key, sel in [(f"{z:.0f}", [u for u in runs if u["z_c"] == z]) for z in en.CORE_ZC] \
                + [("pooled", runs)]:
            fits[f][key] = power_fit([u[f] for u in sel], [u["melt"] for u in sel],
                                     [u["case"] for u in sel])
    return {"label": R24_LABEL, "melt_key": MELT_KEY, "runs": runs, "fits": fits,
            "exponent_name": "q (fit m = k P^q; the key 'n' of each fit holds q)",
            "quadratic_family_predicts_q": 1.0,
            "H1_preregistered_verdict_unchanged": "FAIL",
            "identity": "P_nonlocal = <TF_loc <TF>> = <TF>^2 = TF_eff_shelf^2 exactly, so "
                        "q(P_nonlocal) = n(TF_eff_shelf) / 2 (no information beyond R23)",
            "sources": {"P_nonlocal": "Favier et al. 2019 Eq. 5 (L491-497); Burgard et al. "
                                      "2022 Eq. 18 (L687-702)",
                        "P_local": "Favier et al. 2019 Eq. 4 (L447-455); Burgard et al. 2022 "
                                   "Eq. 14 (L648-654)"},
            "note": "post hoc; far-field (a priori) ramp T, S at each column's ice base, AD16 "
                    "liquidus; does not change the pre-registered H1 verdict"}


def write_quad_outputs(res: dict, out_dir: Path) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    title = "POST HOC (R24): whole-shelf melt against quadratic melt predictors"
    style.apply()
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.9))
    for ax, f, xl in ((axs[0], "P_nonlocal", r"$\langle TF_{loc}\langle TF\rangle\rangle$ "
                                             r"(°C$^2$)"),
                      (axs[1], "P_local", r"$\langle TF_{loc}|TF_{loc}|\rangle$ (°C$^2$)")):
        for z in en.CORE_ZC:
            sel = [u for u in res["runs"] if u["z_c"] == z and u[f] > 0]
            fz = res["fits"][f][f"{z:.0f}"]
            lab = f"z_c = {z:.0f} m" + (f": q = {fz['n']:.2f}" if fz["n"] is not None else "")
            ax.loglog([u[f] for u in sel], [u["melt"] for u in sel], "o", color=ZC_COLOURS[z],
                      ms=4.5, label=lab)
        fp = res["fits"][f]["pooled"]
        xs = np.array([u[f] for u in res["runs"] if u[f] > 0])
        if fp["n"] is not None and xs.size:
            t = np.geomspace(xs.min(), xs.max(), 50)
            ax.loglog(t, fp["k"] * t ** fp["n"], "k-", lw=1.0,
                      label=f"pooled fit: q = {fp['n']:.2f}, R²(log) = {fp['R2_log']:.3f}")
            ms = np.array([u["melt"] for u in res["runs"] if u[f] > 0])
            x0, m0 = np.exp(np.log(xs).mean()), np.exp(np.log(ms).mean())
            ax.loglog(t, m0 * (t / x0), "k:", lw=1.0, label="guide: q = 1 (quadratic family)")
        ax.set(xlabel=xl, ylabel=r"$\langle m \rangle$ whole shelf (m a$^{-1}$)")
        ax.legend(fontsize=6, frameon=False)
    axs[0].set_title("(a) local x non-local (Favier 2019 Eq. 5)", fontsize=8)
    axs[1].set_title("(b) local quadratic (Favier 2019 Eq. 4)", fontsize=8)
    fig.suptitle(title, fontsize=9)
    fig.text(0.01, -0.02, "POST HOC, not a pre-registered test; the pd0 H1 verdict stays FAIL. "
             "TF_loc = far-field ramp T - T_f(S, p) at each column's ice base (a priori). "
             "(a) equals <TF>^2 = TF_eff^2 exactly,\nso its q is n(TF_eff)/2 (R23). 4 km core "
             "runs (s = 1, Tbot -0.5..1.5), last 6 months.", fontsize=6.5, color="0.25",
             ha="left")
    fig.tight_layout()
    paths = style.save(fig, out_dir / "fig2_posthoc_quad")
    plt.close(fig)
    js = {**res, "title": title, "figure": [str(p) for p in paths]}
    (out_dir / "fig2_posthoc_quad.json").write_text(json.dumps(js, indent=2) + "\n")
    return js


def write_outputs(res: dict, out_dir: Path) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    title = "POST HOC (R23): whole-shelf melt against the forcing at the ice base"
    style.apply()
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.9))
    for ax, f, xl in ((axs[0], "TF_eff_shelf", r"$TF_{eff}$, a priori (°C)"),
                      (axs[1], "TD_diag_shelf", r"$T_w - T_b$, diagnosed (°C)")):
        for z in en.CORE_ZC:
            sel = [u for u in res["runs"] if u["z_c"] == z and u[f] is not None and u[f] > 0]
            fz = res["fits"][f][f"{z:.0f}"]
            lab = f"z_c = {z:.0f} m" + (f": n = {fz['n']:.2f}" if fz["n"] is not None else "")
            ax.loglog([u[f] for u in sel], [u["melt"] for u in sel], "o", color=ZC_COLOURS[z],
                      ms=4.5, label=lab)
        fp = res["fits"][f]["pooled"]
        xs = np.array([u[f] for u in res["runs"] if u[f] is not None and u[f] > 0])
        if fp["n"] is not None and xs.size:
            t = np.geomspace(xs.min(), xs.max(), 50)
            ax.loglog(t, fp["k"] * t ** fp["n"], "k-", lw=1.0,
                      label=f"pooled fit: n = {fp['n']:.2f}, R²(log) = {fp['R2_log']:.3f}")
            ms = np.array([u["melt"] for u in res["runs"] if u[f] is not None and u[f] > 0])
            x0, m0 = np.exp(np.log(xs).mean()), np.exp(np.log(ms).mean())
            ax.loglog(t, m0 * (t / x0) ** H08_N, "k:", lw=1.0,
                      label=f"guide: slope {H08_N} (Holland et al. 2008, not a fit)")
        ax.set(xlabel=xl, ylabel=r"$\langle m \rangle$ whole shelf (m a$^{-1}$)")
        ax.legend(fontsize=6, frameon=False)
    axs[0].set_title("(a) a priori ramp forcing at the ice base", fontsize=8)
    axs[1].set_title("(b) diagnosed thermal driving (jpl_compare)", fontsize=8)
    fig.suptitle(title, fontsize=9)
    fig.text(0.01, -0.02, "POST HOC, not a pre-registered test; the pd0 H1 verdict stays FAIL. "
             "4 km core runs (s = 1, Tbot -0.5..1.5), last 6 months. (b): gammaT-weighted "
             "Tw - Tb inverted from AD16 Eq. 24.\n(b) is partly circular: Tw - Tb = m L / (cw gammaT) "
             "is computed from m itself, so m = c gammaT (Tw - Tb) and its exponent is "
             "1 + d ln<gammaT> / d ln(Tw - Tb) (the u* response), not an independent test.",
             fontsize=6.5, color="0.25", ha="left")
    fig.tight_layout()
    paths = style.save(fig, out_dir / "fig2_posthoc_tf")
    plt.close(fig)
    js = {**res, "holland2008_guide_n": H08_N, "title": title,
          "figure": [str(p) for p in paths]}
    (out_dir / "fig2_posthoc_tf.json").write_text(json.dumps(js, indent=2) + "\n")
    return js


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="POST HOC H1 against the ice-base forcing (R23)")
    ap.add_argument("--table", type=Path, default=FIG_DIR / "ensemble_table.json")
    ap.add_argument("--runs-root", type=Path, default=en.DEFAULT_RUNS)
    ap.add_argument("--run-suffix", default="_r21")
    ap.add_argument("--out-dir", type=Path, default=FIG_DIR)
    ap.add_argument("--quad", action="store_true",
                    help="R24 part 2: quadratic local / non-local predictors instead")
    a = ap.parse_args(argv)
    if a.quad:
        res = analyse_quad(en.read_table(a.table), a.runs_root, a.run_suffix)
        write_quad_outputs(res, a.out_dir)
        for f in QUAD:
            print(f, {k: (None if v["n"] is None else round(v["n"], 3), v["R2_log"] and
                          round(v["R2_log"], 4), v["loo_n_range"] and
                          [round(x, 2) for x in v["loo_n_range"]],
                          v["excluded_nonpositive_forcing"]) for k, v in res["fits"][f].items()})
        return 0
    res = analyse(en.read_table(a.table), a.runs_root, a.run_suffix)
    write_outputs(res, a.out_dir)
    for f in FORCINGS:
        print(f, {k: (None if v["n"] is None else round(v["n"], 3), v["R2_log"] and
                      round(v["R2_log"], 4), v["loo_n_range"] and
                      [round(x, 2) for x in v["loo_n_range"]],
                      v["excluded_nonpositive_forcing"]) for k, v in res["fits"][f].items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
