"""Render exp2_isomip_cavity/README.md from README.template.md and the figures/*.json metrics.

Same mechanism as Experiment 1 (exp1_fjord_plumes/scripts/render_readme.py): every result number
in the README comes from a metrics file, never typed by hand. Placeholders are ``{{name}}``; an
unknown placeholder raises. Values whose source file does not exist render as "[pending]" only
with ``--allow-pending``. The void pre-R21 records in figures/void_pre_r21/ are never read.

Usage: render_readme.py [--allow-pending] [--check]   (--check: exit 1 if README.md would change)
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

EXP = Path(__file__).resolve().parents[1]
FIG = EXP / "figures"
TEMPLATE = EXP / "README.template.md"
OUT = EXP / "README.md"
PENDING = "[pending]"
ZC = ("250", "400", "550")


def _load(name: str):
    """A metrics file of figures/ (top level only: never figures/void_pre_r21/)."""
    if "/" in name or "void" in name:
        raise ValueError(f"{name}: only top-level figures/*.json are read, never the void records")
    p = FIG / name
    return json.loads(p.read_text()) if p.exists() else None


def f(x, nd: int = 2) -> str:
    return f"{x:.{nd}f}"


def pct(x, nd: int = 1, sign: bool = False) -> str:
    """Fraction -> percent string (no % sign)."""
    return f"{100 * x:+.{nd}f}" if sign else f"{100 * x:.{nd}f}"


def rng(vals, nd: int = 2) -> str:
    return f"{min(vals):.{nd}f}–{max(vals):.{nd}f}"


def _case_label(c: str) -> str:
    """ens_T-0.5_zc250_s1.0 -> Tbot −0.5 °C, z_c 250 m, s 1.0 (typographic minus later)."""
    m = re.match(r"ens_T(-?[\d.]+)_zc(\d+)_s([\d.]+)", c)
    return f"Tbot = {m[1]} °C, z_c = {m[2]} m, s = {m[3]}" if m else c


def values() -> dict[str, str]:
    v: dict[str, str] = {}
    # ---- calibration (K1, V1), resolution test, Ocean0 supercooling ------------------------
    ce = _load("pd1_calibration_eval_r21.json")
    if ce:
        for res in ("2km", "4km"):
            k = ce["K1"][res]
            v[f"K1_{res}_deep"] = f(k["deep_melt_last6_m_per_a"])
            v[f"gT_cal_{res}"] = f"{k['gammaT_cal']:.6f}"
        v["V1_lo"], v["V1_hi"] = (str(x) for x in ce["V1"]["range"])
        v["V1_verdict"] = "PASS" if ce["V1"]["pass"] else "FAIL"
        rt = ce["resolution_test"]
        v["res_shelf_2km"] = f(rt["shelf_mean_melt"]["2km"])
        v["res_shelf_4km"] = f(rt["shelf_mean_melt"]["4km"])
        v["res_shelf_pct"] = pct(rt["shelf_mean_melt"]["rel_diff_4_vs_2"], 2)
        v["res_moc_pct"] = pct(rt["overturning_max"]["rel_diff_4_vs_2"], 1)
        v["res_psibt_2km"] = f(rt["psi_bt_max_cavity"]["2km"])
        v["res_psibt_4km"] = f(rt["psi_bt_max_cavity"]["4km"])
        v["res_psibt_pct"] = pct(rt["psi_bt_max_cavity"]["rel_diff_4_vs_2"], 0, sign=True)
        v["res_psibt_ratio"] = f(rt["psi_bt_max_cavity"]["4km"] / rt["psi_bt_max_cavity"]["2km"], 1)
        v["res_map_r"] = f(rt["melt_map_correlation"]["r"])
        v["res_map_n"] = str(rt["melt_map_correlation"]["n_cells"])
        runs = ce["runs"]
        v["cal_n_runs"] = str(len(runs))
        v["cal_cfl_max"] = f(max(r["R17_cfl"]["max"] for r in runs.values()), 3)
        v["cal_N1_max"] = pct(max(r["N1"]["value"] for r in runs.values()), 1)
        sc = runs["ocean0_4km_r21"]["R15_supercooling_end_state"]
        v["sc4_count"], v["sc4_nwet"] = str(sc["count"]), f"{sc['n_wet']:,}"
        v["sc4_min"] = f(sc["min_dT"], 2)
        v["sc4_top"] = str(sc["top_cells_any_hfac"]["count"])
    for res in ("2km", "4km"):
        cal = _load(f"calibration_{res}_r21.json")
        if cal:
            for p in cal["points"]:
                v[f"cal_{res}_{str(p['gammaT']).replace('.', 'p')}"] = f(p["deep_melt_last6_m_per_a"], 1)
            v[f"cal_{res}_seg"] = "–".join(str(x) for x in cal["segment"])
    # ---- Ocean1 (V2, V3, N1-N3, CFL, supercooling) -------------------------------------------
    f1 = _load("fig1_ocean.json")
    if f1:
        vd = f1["verdicts"]
        v["V2a"] = f(vd["V2a"]["value"])
        for k in ("V2b", "V2c", "V2d"):
            v[k] = f(vd[k]["value"], 3)
        for k in ("V2a", "V2b", "V2c", "V2d"):
            lo, hi = vd[k]["range"]
            v[f"{k}_range"] = f"{lo:g}–{hi:g}" if k == "V2a" else f"{lo:.2f}–{hi:.2f}"
            v[f"{k}_verdict"] = vd[k]["verdict"]
        ph = vd["V2b_posthoc_full_domain"]
        v["V2b_ph"], v["V2b_ph_jpl"] = f(ph["value"], 3), f(ph["jpl_value"], 3)
        v["V2b_ph_x"], v["V2b_ph_z"] = f"{ph['x_km']:.0f}", f"{ph['z_m']:.0f}"
        v["V2b_below_pct"] = pct(1 - vd["V2b"]["value"] / vd["V2b"]["range"][0], 0)
        for kind in ("bt", "moc"):
            x = vd["V3"][kind]
            v[f"V3_{kind}_C"] = f(x["C_bar"], 3)
            v[f"V3_{kind}_band"] = rng(x["band"], 3)
            v[f"V3_{kind}_verdict"] = x["verdict"]
            v[f"V3_{kind}_fit"] = (f"y = {x['fit_240_months']['slope']:.4f} x + "
                                   f"{x['fit_240_months']['intercept']:.4f}, "
                                   f"r = {x['fit_240_months']['r']:.3f}")
            v[f"V3_{kind}_fit_pub"] = (f"y = {x['fit_published']['slope']:g} x + "
                                       f"{x['fit_published']['intercept']:g}, "
                                       f"r = {x['fit_published']['r']:g}")
        v["V3_mbar"] = f(vd["V3"]["moc"]["m_bar"])
        be, bh = vd["V3"]["bt"]["margin_extended"], vd["V3"]["bt"]["margin_held_at_anchors"]
        assert be["side"] == bh["side"] == "inside", "V3 barotropic wording assumes an inside margin"
        v["V3_bt_ext_pct"], v["V3_bt_held_pct"] = pct(be["rel"], 1), pct(bh["rel"], 1)
        v["V3_bt_ext_lo"] = f(vd["V3"]["bt"]["band"][0], 3)
        v["V3_bt_held_lo"] = f(vd["V3"]["bt"]["band_held_at_anchors"][0], 2)
        me, mh = vd["V3"]["moc"]["margin_extended"], vd["V3"]["moc"]["margin_held_at_anchors"]
        v["V3_moc_ext_pct"], v["V3_moc_ext_sv"] = pct(me["rel"], 1), f(me["abs_Sv"], 4)
        v["V3_moc_held_pct"], v["V3_moc_held_sv"] = pct(mh["rel"], 1), f(mh["abs_Sv"], 4)
        v["V3_moc_ext_lo"] = f(vd["V3"]["moc"]["band"][0], 3)
        v["V3_moc_held_lo"] = f(vd["V3"]["moc"]["band_held_at_anchors"][0], 2)
        v["o1_N1"] = pct(vd["N1"]["value"], 2)
        v["o1_N2"] = f"{vd['N2']['value']:.1e}"
        v["o1_cfl"] = f(vd["R17_cfl_whole_run"]["max"], 3)
        v["o1_deep"] = f(vd["mean_melt_deep300_year20"])
        s1 = vd["R15_supercooling_end_state"]
        v["sc1_count"], v["sc1_nwet"] = str(s1["count"]), f"{s1['n_wet']:,}"
        v["sc1_min"], v["sc1_below01"] = f(s1["min_dT"], 2), str(s1["n_below_minus_0.1"])
        v["sc1_top"] = str(s1["top_cells_any_hfac"]["count"])
        v["sc1_wall"] = str(s1["wall_rows"]["count"])
    # ---- ensemble: N1-N3, H1 ------------------------------------------------------------------
    tab = _load("ensemble_table.json")
    if tab:
        n = len(tab["case"])
        v["ens_n"] = str(n)
        v["ens_n3_pass"] = str(sum(bool(x) for x in tab["N3_pass"]))
        v["ens_n2_pass"] = str(sum(bool(x) for x in tab["N2_pass"]))
        v["ens_n2_max"] = f"{max(abs(x) for x in tab['N2']):.1e}"
        fails = [(c, x) for c, x, p in zip(tab["case"], tab["N1"], tab["N1_pass"]) if not p]
        v["ens_n1_nfail"] = str(len(fails))
        v["ens_n1_fail_list"] = "; ".join(f"{_case_label(c)} ({pct(x, 1)} %)" for c, x in fails)
        i0 = tab["case"].index("ocean0_4km")
        v["ocean0_4km_melt"] = f(tab["mean_melt_shelf"][i0])
        v["ocean0_4km_Dout"] = f"{tab['D_out'][i0]:.0f}"
        test = [(c, d) for c, d, s in zip(tab["case"], tab["D_out"], tab["split"]) if s == "test"]
        v["test_Dout"] = ", ".join(f"{d:.0f}" for _, d in test)
        melt = dict(zip(tab["case"], tab["mean_melt_shelf"]))
        v["test_melt"] = ", ".join(f(melt[c], 1) for c, _ in test)
    h1 = _load("fig2_h1.json")
    if h1:
        v["T_f_GL"] = f(h1["H1"]["T_f_GL_degC"], 3)
        lo, hi = 1.3, 2.0
        for z in ZC:
            x = h1["H1"]["per_z_c"][z]
            v[f"H1_R2_{z}"] = f"{x['quadratic']['R2']:.4f}"
            v[f"H1_n_{z}"] = f(x["power_law"]["n"])
            v[f"H1_k_{z}"] = f"{x['power_law']['k']:.3g}"
            v[f"H1_verdict_{z}"] = x["verdict"]
            v[f"H1_melt_{z}"] = ", ".join(f(m, 1) for m in x["melt"])
            v[f"H1_ratio_{z}"] = f(x["melt"][-1] / x["melt"][0], 1)
            sx = h1["H1_secondary_deep300"]["per_z_c"][z]
            v[f"H1deep_n_{z}"] = f(sx["power_law"]["n"])
        v["H1_band"] = f"[{lo}, {hi}]"
        v["H1_R2_min"] = f"{min(h1['H1']['per_z_c'][z]['quadratic']['R2'] for z in ZC):.4f}"
    # ---- H2 -------------------------------------------------------------------------------
    h2 = _load("fig3_h2.json")
    if h2:
        for z in ZC:
            x = h2["H2"]["per_z_c"][z]
            v[f"H2_Dout_{z}"] = ", ".join(f"{d:.0f}" for d in x["D_out"])
            v[f"H2_rho_{z}"] = f"{x['rho']:.1f}"
            v[f"H2_verdict_{z}"] = x["verdict"]
    # ---- post hoc R23 / R24 ---------------------------------------------------------------
    tf = _load("fig2_posthoc_tf.json")
    if tf:
        for var in ("TF_eff_shelf", "TF_eff_deep300", "TD_diag_shelf"):
            for z in ZC + ("pooled",):
                x = tf["fits"][var][z]
                v[f"{var}_{z}"] = f(x["n"])
                v[f"{var}_{z}_loo"] = rng(x["loo_n_range"])
                v[f"{var}_{z}_R2"] = f(x["R2_log"], 2 if x["R2_log"] < 0.99 else 4)
        v["holland_n"] = f"{tf['holland2008_guide_n']:.2f}"
        r550 = [r["TF_eff_shelf"] for r in tf["runs"] if r["z_c"] == 550.0]
        v["TF_eff_550_range"] = rng(r550)
    qd = _load("fig2_posthoc_quad.json")
    if qd:
        for z in ZC + ("pooled",):
            x = qd["fits"]["P_local"][z]
            v[f"Ploc_{z}"] = f(x["n"])
            v[f"Ploc_{z}_loo"] = rng(x["loo_n_range"])
            v[f"Ploc_{z}_R2"] = f(x["R2_log"])
        v["Pnonloc_pooled"] = f(qd["fits"]["P_nonlocal"]["pooled"]["n"])
    r24 = _load("fig2_posthoc_r24.json")
    if r24:
        v["r24_n"] = f(r24["power_law"]["n"])
        v["r24_loo"] = rng(r24["loo_n_range"])
        v["r24_R2q"] = f"{r24['quadratic']['R2']:.7f}"
        v["r24_melt"] = ", ".join(f(m, 1) for m in r24["melt"])
        v["r24_verdict"] = {"no apoya": "does not support", "apoya": "supports"}[r24["verdict"]]
        v["r24_TF_n"] = f(r24["context_fits_hol"]["TF_eff_shelf"]["n"])
        v["r24_Ploc_q"] = f(r24["context_fits_hol"]["P_local"]["n"])
        v["r24_N1_max"] = pct(max(r24["N"]["N1"].values()), 2)
        v["r24_n_runs"] = str(len(r24["cases"]))
        loo = [x for z in ZC for x in r24["core_fits_Tbot"][z]["loo_n_range"]]
        v["H1_loo_min"], v["H1_loo_max"] = f(min(loo)), f(max(loo))
    # ---- post hoc comparison with MITgcm-JPL ------------------------------------------------
    jp = _load("jpl_comparison.json")
    if jp:
        for case in ("ocean0", "ocean1"):
            c = jp[case]
            for side in ("ours", "jpl"):
                for reg in ("shelf", "deep300"):
                    r = c[side][reg]
                    v[f"j{case}_{side}_{reg}_melt"] = f(r["melt_m_per_a"])
                    v[f"j{case}_{side}_{reg}_ustar"] = f"{r['ustar_m_s']:.5f}"
                    v[f"j{case}_{side}_{reg}_td"] = f(r["thermal_driving_degC"], 3)
                circ = c[side]["circulation"]
                v[f"j{case}_{side}_moc"] = f(circ["psi_moc_max_abs_Sv"], 3)
                v[f"j{case}_{side}_moc_full"] = f(circ["psi_moc_max_abs_full_domain_Sv"], 3)
                v[f"j{case}_{side}_bt_front"] = f(circ["psi_bt_max_abs_front_Sv"], 3)
                v[f"j{case}_{side}_bt_cav"] = f(circ["psi_bt_max_abs_cavity_Sv"], 3)
            v[f"j{case}_jpl_gT"] = f"{c['jpl']['gammaT']:g}"
            v[f"j{case}_r"] = f(c["melt_map_correlation"]["r_either_shelf_zero_off"])
            v[f"j{case}_r_both"] = f(c["melt_map_correlation"]["r_both_shelf"])
            v[f"j{case}_n_jpl_only"] = str(c["melt_map_correlation"]["n_shelf_only_jpl"])
            v[f"j{case}_n_ours_only"] = str(c["melt_map_correlation"]["n_shelf_only_ours"])
            v[f"j{case}_n_jpl"] = str(c["jpl"]["shelf"]["n_cells"])
            v[f"j{case}_n_ours"] = str(c["ours"]["shelf"]["n_cells"])
        v["jocean0_jpl_deep_ustar_3"] = f"{jp['ocean0']['jpl']['deep300']['ustar_m_s']:.4f}"
    # ---- Exp. 3 export ------------------------------------------------------------------------
    e3 = _load("exp3_dataset.json")
    if e3:
        v["e3_n"], v["e3_train"], v["e3_test"] = str(e3["n_runs"]), str(e3["n_train"]), str(e3["n_test"])
        v["e3_excluded"] = str(len(e3["excluded"]))
        v["e3_n1"] = str(e3["n_N1_fail"])
        v["e3_sha12"] = e3["sha256"][:12]
        v["e3_mb"] = f(e3["bytes"] / 1e6, 1)
    return v


def render(allow_pending: bool) -> str:
    text = TEMPLATE.read_text()
    v = values()
    missing = sorted({k for k in re.findall(r"\{\{(\w+)\}\}", text) if k not in v})
    if missing and not allow_pending:
        raise KeyError(f"no value for placeholders: {missing}")
    return re.sub(r"\{\{(\w+)\}\}", lambda m: typographic_minus(v.get(m.group(1), PENDING)), text)


def typographic_minus(s: str) -> str:
    """Hyphen-minus before a digit -> U+2212 minus ("-1.5" -> "−1.5"); hyphens inside ranges
    such as "0.50-0.61" and inside case names (ens_T-0.5) are left alone."""
    return re.sub(r"(?<![\w.])-(?=\d)", "−", s)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--allow-pending", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    new = render(a.allow_pending)
    if a.check:
        ok = OUT.exists() and OUT.read_text() == new
        print("README.md is up to date" if ok else "README.md differs from the rendered template")
        sys.exit(0 if ok else 1)
    OUT.write_text(new)
    n = new.count(PENDING)
    print(f"wrote {OUT} ({len(new.split())} words, {n} pending values)")
