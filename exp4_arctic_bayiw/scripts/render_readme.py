#!/opt/anaconda3/envs/pyocean/bin/python
"""Render exp4_arctic_bayiw/README.md from README.template.md and the committed metrics files.

Same mechanism as Experiments 1-3 (exp3_melt_emulator/scripts/render_readme.py): every number in
the README comes from a committed file, never typed by hand. The sources are figures/*.json and,
for the configuration (grid, time step, wind pulse, ice), the case registry cases.yaml, the file
from which every run's inputs are generated. Placeholders are ``{{name}}``; an unknown placeholder
raises. Values whose source file does not exist render as "[pending]" only with
``--allow-pending``. A test (tests/test_render_readme.py) fails if the template holds a digit
that is not a whitelisted label (NUMBER_WHITELIST).

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
BANDS = ("kelvin_subinertial", "near_inertial", "superinertial", "seiche")
BAND_SHORT = {"low": "low", "kelvin_subinertial": "K", "near_inertial": "NI", "superinertial": "SI",
              "seiche": "S"}

# Digits allowed in the template outside placeholders: labels and citation years, never results.
CITATION_YEARS = ("2027", "1982", "1986", "1994", "1997", "2000", "2001", "2005", "2006", "2012", "2014",
                  "2016", "2017", "2018", "2020", "2022", "2023", "2024", "2026")
NUMBER_WHITELIST = (
    r"\b(?:R1k|R1b|R[1-3]|F1c|F1g|F[1-5]|V[1-4][ab]?)\b",           # run and criterion names
    r"\bH4[a-d]?\b", r"\bR_?1\b", r"\bc_?1\b", r"\bR₁", r"\bc₁", r"\bc_n\b",
    r"\bExp(?:eriment|\.)? [1-5]\b", r"\bPD[0-3]\b", r"^#+ [0-9]+(?:\.[0-9]+)*\.?",
    r"\b[Ee]qs?\. [0-9]+(?:\.[0-9]+)*(?:[–-][0-9]+)?", r"§[0-9]+(?:\.[0-9]+)*",
    r"\b(?:Figure|Fig\.|Figures) [0-9](?:[a-e])?(?:[–-][0-9])?", r"\bAnimations? [1-6]\b", r"\bAmendments? [1-3](?:[–-][1-3])?",
    r"\bSentinel-2\b", r"\bMODIS\b", r"\bIBCAO\b", r"\bS100\b", r"\bTEOS-10\b", r"\bGGL90\b",
    r"\bM2\b", r"\bMP4\b", r"\bk = 1\b", r"\b(?:class )?\([a-d]\)", r"\(i{1,3}\)", r"\([1-5]\)", r"\bsha256\b", r"\b(?:i|ii)\b(?=\))",
    r"\bv5\.2\b", r"\bcheckpoint69j\b", r"\bUTM 33N\b", r"\bΣ₅",
    r"\b(?:" + "|".join(CITATION_YEARS) + r")[a-b]?\b",
)


def strip_allowed(text: str) -> str:
    """The template minus placeholders, code spans, URLs, comments, the References section and the
    whitelisted tokens: what is left must contain no digit."""
    text = text.split("\n## References", 1)[0]
    text = re.sub(r"<!--.*?-->", " ", text, flags=re.S)
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    text = re.sub(r"\{\{\w+\}\}", " ", text)
    text = re.sub(r"`[^`\n]*`", " ", text)
    text = re.sub(r"\]\([^)]*\)", "]", text)
    text = re.sub(r"https?://\S+", " ", text)
    for pat in NUMBER_WHITELIST:
        text = re.sub(pat, " ", text, flags=re.M)
    return text


WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")


def word(n: int, cap: bool = False) -> str:
    """Count -> English word (counts in the text come from the JSON too)."""
    w = WORDS[n] if 0 <= n < len(WORDS) else str(n)
    return w.capitalize() if cap else w


def _load(name: str):
    p = FIG / name
    return json.loads(p.read_text()) if p.exists() else None


def _cases():
    import yaml
    return yaml.safe_load((EXP / "cases.yaml").read_text())


def f(x, nd: int = 2) -> str:
    return f"{x:.{nd}f}"


def gj(x, nd: int = 2, sign: bool = False) -> str:
    """J -> GJ string."""
    return f"{x / 1e9:+.{nd}f}" if sign else f"{x / 1e9:.{nd}f}"


def pct(x, nd: int = 1, sign: bool = False) -> str:
    """Fraction -> percent string (no % sign)."""
    return f"{100 * x:+.{nd}f}" if sign else f"{100 * x:.{nd}f}"


def rng(vals, nd: int = 2, sep: str = "–") -> str:
    return f"{min(vals):.{nd}f}{sep}{max(vals):.{nd}f}"


def typographic_minus(s: str) -> str:
    """Hyphen-minus before a digit -> U+2212 minus; hyphens in ranges ("1.7-6.6") are left alone."""
    return re.sub(r"(?<![\w.])-(?=\d)", "−", s)


def values() -> dict[str, str]:
    v: dict[str, str] = {}
    cs = _cases()
    fj = cs["configs"]["fjord_vmf"]
    x0, x1, y0, y1 = fj["bbox_utm"]
    v["box_km"] = f"{(x1 - x0) / 1e3:.0f} × {(y1 - y0) / 1e3:.0f} km"
    v["dx_m"], v["nz"] = f"{fj['dx_m']:.0f}", str(fj["nz"])
    v["dx400_m"] = f"{cs['configs']['fjord_vmf_400']['dx_m']:.0f}"
    v["n_sponge"] = str(fj["n_sponge"])
    v["dz_top"] = f(fj["vgrid"]["dz_top"], 1)
    v["H_grid"] = f"{fj['vgrid']['H_m']:.0f}"
    v["deltaT"] = f"{fj['deltaT']:.0f}"
    v["forced_h"], v["t_end_h"] = f"{fj['forced_h']:.0f}", f"{fj['t_end_h']:.0f}"
    w = fj["wind"]
    v["tau_max"], v["sigma_h"], v["center_h"] = f(w["tau_max"], 1), f(w["sigma_h"], 0), f(w["center_h"], 1)
    v["ice_h"], v["ice_A"] = f(fj["ice"]["h_m"], 0), f(fj["ice"]["A"], 0)
    v["pulse_end_h"] = f(w["center_h"] + 3 * w["sigma_h"], 1)
    # fjord runs whose wind blows at another angle than the config (along-axis) one: cross-axis runs
    cross = [n for n, c in cs["cases"].items()
             if str(c.get("config", "")).startswith("fjord_vmf") and c.get("wind_angle_deg") is not None
             and abs(float(c["wind_angle_deg"]) - float(fj["wind_angle_deg"])) > 1e-6]
    assert cross, "README: the cross-axis fjord run is named in §2.3"
    v["cross_axis_runs"] = " and ".join(cross)
    r1k = cs["cases"]["R1k"]
    v["R1k_t_end"] = f"{r1k['t_end_h']:.0f}"
    bay = cs["configs"]["bay_thesis"]                       # the bay runs keep the thesis set-up
    v["bay_dt"], v["bay_forced"] = f"{bay['deltaT']:.0f}", f"{bay['forced_h']:.0f}"
    bw = bay["wind"]
    v["bay_tau"], v["bay_sigma"], v["bay_center"] = f(bw["max_tau"], 1), f(bw["sigma_t_h"], 0), f(bw["center_t_h"], 1)
    v["bay_patch_r"] = str(bw["radius"])
    v["R1_t_end"] = f"{cs['cases']['R1'].get('t_end_h', bay['t_end_h']):.0f}"
    v["R2_t_end"] = f"{cs['cases']['R2']['t_end_h']:.0f}"
    assert cs["cases"]["R3"]["t_end_h"] == cs["cases"]["R2"]["t_end_h"], "README: 'R2 and R3' share a duration"
    mods = sorted(d.parent.name for d in (EXP / "config").glob("*/code")
                  if all((d / x).exists() for x in ("obcs_check.F", "seaice_dynsolver.F", "obcs_add_tides.F")))
    v["mod_configs"] = ", ".join(f"`{m}`" for m in mods)

    th = _load("theory.json")
    if th:
        c = th["cases"]
        sm, sb, wm, wb = (c["vmf_summer"], c["vmf_summer_Hmean"], c["vmf_winter"], c["vmf_winter_Hmean"])
        bay, b78 = c["thesis_bay_thesis_f"], c["thesis_bay_78N_vmf_summer"]
        v["lat"] = f"{th['lat_fjord_deg']:.0f}"
        v["f78"] = f"{sm['f_s-1']:.3e}".replace("e-04", " × 10⁻⁴")
        v["Ti78"] = f(sm["inertial_period_h"], 1)
        v["f_thesis"] = f"{bay['f_s-1']:.2e}".replace("e-05", " × 10⁻⁵")
        v["Ti_thesis"] = f(bay["inertial_period_h"], 1)
        v["c1_summer"] = rng([sb["c_n"][0], sm["c_n"][0]])
        v["c1_winter"] = rng([wb["c_n"][0], wm["c_n"][0]])
        v["R1_summer"] = rng([sb["rossby_radius_m"][0] / 1e3, sm["rossby_radius_m"][0] / 1e3], 1)
        v["R1_winter"] = rng([wb["rossby_radius_m"][0] / 1e3, wm["rossby_radius_m"][0] / 1e3], 2)
        res = th["resolution"]
        v["R1dx_summer"] = rng([res["vmf_summer_Hmean"]["R1_over_dx_200m"], res["vmf_summer"]["R1_over_dx_200m"]], 0)
        v["R1dx_winter"] = rng([res["vmf_winter_Hmean"]["R1_over_dx_200m"], res["vmf_winter"]["R1_over_dx_200m"]], 1)
        v["R1dx_summer_400"] = rng([res["vmf_summer_Hmean"]["R1_over_dx_400m"], res["vmf_summer"]["R1_over_dx_400m"]], 1)
        v["TK_summer"] = rng([sm["kelvin_circuit_period_h"][0], sb["kelvin_circuit_period_h"][0]], 0)
        v["poincare1"] = rng([sb["poincare_cross_1_over_f"][0], sm["poincare_cross_1_over_f"][0]])
        v["poincare_bay"] = rng([b78["poincare_cross_1_over_f"][0], bay["poincare_cross_1_over_f"][0]])
        g = th["fjord_geometry"]
        v["fj_area"] = f(g["area_m2"] / 1e6, 0)
        v["fj_L"], v["fj_W"] = f(g["along_axis_length_m"] / 1e3, 1), f(g["mean_width_m"] / 1e3, 1)
        v["fj_Hmax"], v["fj_Hmean"] = f(g["max_depth_m"], 0), f(g["mean_depth_m"], 1)
        v["fj_P"] = f(g["contour_outer_m"] / 1e3, 0)
        v["fj_axis"] = f(g["axis_angle_deg"], 1)
        v["wind_from"] = f"{th['wind']['down_fjord_from_direction_deg_true_utm']:03.0f}"
        bg = th["thesis_bay_geometry"]
        v["bay_dx"] = f(sum(bg["dx_bay_m"]) / len(bg["dx_bay_m"]) / 1e3, 1)
        v["bay_L"], v["bay_W"], v["bay_H"] = f(bg["length_m"] / 1e3, 0), f(bg["width_m"] / 1e3, 1), f(bg["max_depth_m"], 0)
        v["c1_bay"], v["c1_bay78"] = f(bay["c_n"][0]), f(b78["c_n"][0])
        v["R1_bay"], v["R1_bay78"] = f(bay["rossby_radius_m"][0] / 1e3, 1), f(b78["rossby_radius_m"][0] / 1e3, 1)
        e = th["eos"]["fitted_vmf"]
        v["eos_ta"] = f"{e['tAlpha']:.4e}".replace("e-05", " × 10⁻⁵")
        v["eos_sb"] = f"{e['sBeta']:.4e}".replace("e-04", " × 10⁻⁴")
        v["eos_err"] = f(th["eos"]["check_on_final_profiles"]["max_abs_err_kg_m3"], 3)
        ps = th["profiles_summary"]
        v["n_summer"], v["n_winter"] = str(ps["summer"]["n_qc_pass"]), str(ps["winter"]["n_qc_pass"])
        v["n_winter_years"] = str(len(ps["winter"]["years_present"]))
        v["v2b_range"] = rng(th["v2b_range_m_s-1"]["fitted_vmf"])
        b = {x["name"]: x for x in th["bands"]["fjord_78N_280h"]["bands"]}
        v["low_period"] = f"{b['low']['period_h'][0]:.0f}"
        v["ni_period"] = rng(b["near_inertial"]["period_h"], 1)
        v["k_period"] = f"{b['kelvin_subinertial']['period_h'][0]:.1f}–{b['kelvin_subinertial']['period_h'][1]:.0f}"
        v["si_period"] = rng(b["superinertial"]["period_h"], 1)
        v["seiche_period"] = f(b["seiche"]["period_h"][1], 1)
        v["ni_bins"] = f(b["near_inertial"]["n_bins"], 1)
        v["n_resolved"] = word(sum(x["name"] != "low" for x in th["bands"]["fjord_78N_280h"]["bands"]))
        v["ni_lo"], v["ni_hi"] = (f(x, 2) for x in th["bands"]["fjord_78N_280h"]["ni_factors"])
    bt = _load("bathy_report.json")
    if bt:
        s = bt["sounds"]
        v["aks_w"] = f(s["Akselsundet"]["width_npi_m"] / 1e3, 2)
        v["mar_w"] = f(max(s["Mariasundet"]["passages_npi_m"]) / 1e3, 2)
        v["mar_grid_w"] = f(s["Mariasundet"]["grid"]["wet_length_m"] / 1e3, 2)
        v["aks_sill"], v["mar_sill"] = f(s["Akselsundet"]["sill_depth_grid_m"], 0), f(s["Mariasundet"]["sill_depth_grid_m"], 0)
        tid = bt["tid_fraction"]["cells_inner"]
        v["tid_measured"], v["tid_interp"] = pct(tid["measured"]), pct(tid["indirect"])
        v["tid_land"] = pct(tid["land"])
    v1 = _load("v1_thesis.json")
    if v1:
        v["V1_n"] = str(v1["n_common_eta_records"])
        v["V1_deta"] = f"{v1['eta_max_abs_diff_m']:.1f}"
        v["V1_tol"] = _sci(v1["tolerance_eta_m"])
        v["V1_verdict"] = "PASS" if v1["V1_exacto_pass"] else "FAIL"
    fx = _load("flux_vs_thesis.json")
    bx = _load("budget_vs_thesis.json")
    if fx and bx:
        v["nash_rel"] = _sci(fx["max_rel_diff"])
        v["kf_rel"] = _sci(bx["max_rel_diff"])
        v["contrast_tol"] = pct(fx["tolerance_rel"], 0)

    val = _load("validation.json")
    if val:
        a = val["V2a"]
        v["V2a_peak"] = f(a["f_peak_hz"] / a["f_hz"], 3)
        v["V2a_off"] = f(a["offset_df"], 2)
        v["V2a_verdict"] = "PASS" if a["pass"] else "FAIL"
        ct = a["classification"]["column_test"]
        v["col_peak"], v["col_off"] = f(ct["peak_over_f"], 3), f(ct["offset_df"], 2)
        b2 = val["V2b"]
        v["V2b_c"] = f(b2["c_m_s"])
        v["V2b_n"] = str(b2["n_stations"])
        v["V2b_verdict"] = "PASS" if b2["pass"] else "FAIL"
        fw = b2["descriptive_free_wave"]
        v["fw_c"], v["fw_se"] = f(fw["c_m_s"]), f(fw["c_stderr_m_s"])
        arr = [x / 3600.0 for x in b2["arrival_s"]]
        end = w["center_h"] + 3 * w["sigma_h"]
        v["V2b_early"] = str(sum(x <= end for x in arr))
        v3 = val["V3"]
        v["V3_decay"], v["V3_c"] = f(v3["presence"]["decay_scale_km"], 1), f(v3["speed"]["c_m_s"])
        v["V3_verdict"] = "PASS" if v3["pass"] else "FAIL"
        t = val["thresholds"]
        v["V3_c_range"] = rng(t["V3"]["c_range_m_s"])
        v["V3_decay_range"] = rng(t["V3"]["decay_scale_km"], 0)
        v["V4a_tol"] = pct(t["V4a"]["max_rel_diff"], 0)
        v["V4b_tol"] = f(t["V4b"]["max_abs_diff"], 2)
        v["factor"] = f(t["factor_min"], 1)
        v["factor_inv"] = f"1/{t['factor_min']:.1f}"
        v["H4a_low_max"] = f(t["H4a"]["low_band_max_frac"], 2)
        v["H4a_phi"] = f(t["H4a"]["phi_ni_support"], 1)
        v["H4b_max"] = f(t["H4b"]["energy_ratio_max"], 1)
        v4 = val["V4"]
        v["V4_F1"], v["V4_V4"] = gj(v4["V4a"]["F1_out_J"], 2, True), gj(v4["V4a"]["V4_out_J"], 2, True)
        v["V4a_rel"] = pct(v4["V4a"]["rel_diff"], 0)
        v["V4b_max"] = f(max(v4["V4b"]["abs_diff"].values()), 2)
        v["V4b_min"] = f(min(v4["V4b"]["abs_diff"].values()), 2)
        v["V4_verdict"] = "PASS" if v4["pass"] else "FAIL"
        nf = len(val["summary"]["failed"])
        v["n_failed"], v["n_failed_Word"] = word(nf), word(nf, True)
    d2 = _load("pd2_v2_diagnostics.json")
    if d2:
        pc = d2["V2a_peak_map"]["peak_over_f_percentiles"]
        v["v2a_med"], v["v2a_p5"], v["v2a_p95"] = f(pc["50"]), f(pc["5"]), f(pc["95"])
        v["v2a_frac"] = pct(d2["V2a_peak_map"]["fraction_within_1df_of_f"], 0)
    d4 = _load("pd2_v4_diagnostics.json")
    if d4:
        g1 = sum(x["gross_J"] for x in d4["F1"].values())
        v["v4_aks_ratio"] = f(d4["V4"]["Akselsundet"]["gross_J"] / d4["F1"]["Akselsundet"]["gross_J"], 2)
        v["v4_mar_ratio"] = f(d4["V4"]["Mariasundet"]["gross_J"] / d4["F1"]["Mariasundet"]["gross_J"], 2)
        v["aks_gross_share"] = pct(d4["F1"]["Akselsundet"]["gross_J"] / g1, 0)
    ice = _load("ice_ob_check.json")
    if ice:
        m = [x["min_cells"] for x in ice.values()]
        v["ice_min_cells"] = f"{min(m):.0f}"

    h = _load("hypotheses.json")
    if h:
        nh = sum(bool(re.fullmatch(r"H4[a-d]", k)) for k in h)
        v["n_hyp"], v["n_hyp_Word"] = word(nh), word(nh, True)
        a = h["H4a"]
        v["H4a_low_share"] = f(a["low_band_share"], 3)
        v["H4a_low"] = gj(a["bands_out_J"]["low"], 3, True)
        v["H4a_res"] = gj(a["resolved_sum_out_J"], 3, True)
        for k in BANDS:
            v[f"H4a_{BAND_SHORT[k]}"] = gj(a["bands_out_J"][k], 3, True)
            v[f"H4a_phi_{BAND_SHORT[k]}"] = f(a["phi"][k], 2)
        v["H4a_outcome"] = _outcome(a["outcome"])
        v["H4a_FA_share"] = f(a["descriptive_F_plus_A"]["low_band_share"], 3)
        b4 = h["H4b"]
        v["H4b_EF2"], v["H4b_EF3"] = gj(b4["i"]["E_in_F2_J"], 3), gj(b4["i"]["E_in_F3_J"], 3)
        v["H4b_ratio"] = f(b4["i"]["ratio"], 3)
        v["H4b_i"] = "met" if b4["i"]["met"] else "not met"
        v["H4b_ke_ratio"] = f"{b4['i']['post_hoc_ke_bc_only']['ratio']:.4f}"
        v["H4b_Fbc"], v["H4b_Fbt"] = gj(b4["ii"]["F_bc_in_J"], 4), gj(b4["ii"]["F_bt_in_J"], 4)
        v["H4b_ii"] = "met" if b4["ii"]["met"] else "not met"
        c4 = h["H4c"]
        v["H4c_i"] = f(c4["i"]["ratio"], 2)
        v["H4c_WF1"], v["H4c_WF3"] = gj(c4["i"]["W_NI_F1_J"], 2), gj(c4["i"]["W_NI_F3_J"], 2)
        v["H4c_ii"] = f(c4["ii"]["ratio"], 1)
        v["H4c_F1"], v["H4c_F3"] = gj(c4["ii"]["F_res_out_F1_J"], 4, True), gj(c4["ii"]["F_res_out_F3_J"], 5, True)
        v["H4c_outcome"] = _outcome(c4["outcome"])
        d = h["H4d"]["R1k_vs_R1"]
        a1 = d["post_hoc_amendment1_numerator"]
        v["a1_R_ratio"] = f(a1["ratio"], 2)
        v["a1_R1k_eta"], v["a1_R1_eta"] = pct(a1["eta_kpp"], 2), pct(a1["eta_closure"], 2)
        v["a1_R_outcome"] = PD0_WORD[a1["literal_rule_would_read"]]
        v["R1_window"] = f"{d['window_h'][1]:.0f}"
        v["H4d_R1k_eta"], v["H4d_R1_eta"] = pct(d["eta_kpp"], 2), pct(d["eta_closure"], 3)
        v["H4d_R_ratio"] = f(d["ratio"], 2)
        v["H4d_R_outcome"] = _outcome(d["outcome"])
        v["H4d_WR1k"], v["H4d_WR1"] = gj(d["W_wind_J"]["kpp"], 1), gj(d["W_wind_J"]["closure"], 1)
        v["H4d_W_ratio"] = f(d["W_wind_J"]["closure"] / d["W_wind_J"]["kpp"], 1)
        v["H4d_num_ratio"] = f(d["numerator_J"]["kpp"] / d["numerator_J"]["closure"], 2)
        d = h["H4d"]["F1_vs_F1c"]
        v["H4d_F1_eta"], v["H4d_F1c_eta"] = pct(d["eta_kpp"], 2), pct(d["eta_closure"], 2)
        v["H4d_F_ratio"] = f(d["ratio"], 1)
        v["H4d_F_outcome"] = _outcome(d["outcome"])
        v["H4d_F1_num"] = gj(d["numerator_J"]["kpp"], 3, True)
        v["H4d_F1c_num"] = gj(d["numerator_J"]["closure"], 3, True)
        r = h["post_hoc_h4d_resolved_bands"]
        v["ph_eta_F1"], v["ph_eta_F1c"] = pct(r["F1"]["eta_resolved"], 3), pct(r["F1c"]["eta_resolved"], 3)
        v["ph_eta_ratio"] = f(r["ratio_F1_over_F1c"], 2)
        br = h["post_hoc"]["bridge"]
        for run in ("R1k", "R2", "F1"):
            s = br[run]["flux_share_signed"]
            v[f"br_{run}_low"], v[f"br_{run}_ni"] = f(s["low"], 2), f(s["near_inertial"], 2)
        if "Seiche" not in br["R2"]["station_fraction"]:
            raise KeyError("hypotheses.json post_hoc.bridge.R2.station_fraction has no 'Seiche' station")
        v["br_R2_station_ni"] = f(br["R2"]["station_fraction"]["Seiche"]["near_inertial"], 3)
        cm = br["F1"]["station_fraction"]["coast_south_mean"]
        v["br_F1_coast_k"], v["br_F1_coast_ni"] = f(cm["kelvin_subinertial"], 2), f(cm["near_inertial"], 2)
        al = h["post_hoc"]["arneborg_liljebladh_net_in_J"]
        v["al_F1"] = gj(al["F1"]["total"]["F_bc_in_J"], 3, True)
        v["al_F1c"] = gj(al["F1c"]["total"]["F_bc_in_J"], 3, True)
        v["al_R1"] = gj(al["R1"]["total"]["F_bc_in_J"], 3, True)
        v["al_R2"] = gj(al["R2"]["total"]["F_bc_in_J"], 3, True)
        sl = h["post_hoc"]["strait_limitation"]
        v["mar_wR1"] = rng([sl["Mariasundet"]["width_over_R1"]["vmf_summer"],
                            sl["Mariasundet"]["width_over_R1"]["vmf_summer_Hmean"]])
        v["aks_wR1"] = rng([sl["Akselsundet"]["width_over_R1"]["vmf_summer"],
                            sl["Akselsundet"]["width_over_R1"]["vmf_summer_Hmean"]])
        bpf = h["post_hoc"]["band_flux_per_strait"]["F1c"]["total"]
        v["F1c_K"] = gj(bpf["resolved_out_J"]["kelvin_subinertial"], 3, True)
    bud = _load("budget.json")
    if bud:
        r = bud["runs"]
        v["ape_F1_end"] = gj(r["F1"]["energy_J"]["ape"]["end"], 1)
        v["ape_R2_end"] = gj(r["R2"]["energy_J"]["ape"]["end"], 0)
        v["ape_F2_end"] = gj(r["F2"]["energy_J"]["ape"]["end"], 2)
        v["ke_F2_end"] = gj(r["F2"]["energy_J"]["ke_bc"]["end"], 3)
        v["P_F1"] = gj(r["F1"]["P_J"], 1)
        v["r1_F_out"] = gj(r["R1"]["F_bc_pressure_out_J"], 2)
        v["r1_A_in"] = gj(r["R1"]["totals_in_J"]["A_bc"], 2)
        _expect(r["R1"]["F_bc_pressure_out_J"] > 0 and r["R1"]["totals_in_J"]["A_bc"] > 0,
                "R1: pressure work out of the bay and advection into it")
    f0 = _load("fig0_schematic.json")
    if f0:
        v["fig0_date"], v["fig0_year"] = f0["sentinel2"]["date_text"], f0["sentinel2"]["year"]
    mr = _load("posthoc_moving_reference.json")
    if mr:
        R = mr["runs"]
        v["mr_F1_low"] = gj(R["F1"]["low_band_out_J"], 3, True)
        v["mr_F1_res"] = gj(sum(R["F1"]["total"]["band_out_J"][k] for k in BANDS), 3, True)
        v["mr_F1_eta"] = pct(R["F1"]["eta"], 2)
        v["mr_F1_eta_res"] = pct(R["F1"]["eta_resolved"], 3)
        tb = R["F1"]["total"]["band_out_J"]
        v["mr_F1_share"] = f(abs(tb["low"]) / sum(abs(x) for x in tb.values()), 2)
        v["mr_F1_K"] = gj(tb["kelvin_subinertial"], 3, True)
        v["mr_F1c_low"] = gj(R["F1c"]["low_band_out_J"], 3, True)
        v["mr_F1c_eta"] = pct(R["F1c"]["eta"], 2)
        v["mr_ratio_R"] = f(R["R1k_0_120h"]["eta"] / R["R1"]["eta"], 2)
        if h:   # share of the fixed-reference low band that the moving reference removes
            v["mr_drift_frac"] = f(1 - R["F1"]["low_band_out_J"] / h["H4a"]["bands_out_J"]["low"], 2)
        pd_ = mr["profile_drift"]
        v["drift_F1c"] = f(pd_["F1c"]["max_abs_change_kg_m-3"], 3)
        v["drift_F1"] = f(pd_["F1"]["max_abs_change_kg_m-3"], 3)
        v["drift_F2"] = f(pd_["F2"]["max_abs_change_kg_m-3"], 3)
        v["drift_bay"] = rng([pd_["R1"]["max_abs_change_kg_m-3"], pd_["R1k"]["max_abs_change_kg_m-3"]], 3)
    # the Blender animations of §3.6 (analysis/blender_render.py movie writes these summaries)
    ba, bb = _load("blender_shotA_F1.json"), _load("blender_shotB_F3F2.json")
    if ba:
        v["blA_sigma"], v["blA_z"] = f(ba["sigma_kg_m3"], 3), f(ba["sigma_depth_t0_m"], 0)
        v["blA_lim"] = f"{ba['colour_limit_m']:g}"
        v["bl_ve_sea"], v["bl_ve_land"] = f"{ba['ve_sea']:g}", f"{ba['ve_land']:g}"
        sa = ba["saturation"]["F1"]
        v["blA_sat"], v["blA_sat_t"] = pct(sa["worst"]), f(sa["worst_t_h"], 1)
    if bb:
        assert bb["ve_sea"] == ba["ve_sea"] and bb["ve_land"] == ba["ve_land"], "README: one VE for both shots"
        v["blB_sigma"], v["blB_z"] = f(bb["sigma_kg_m3"], 3), f(bb["sigma_depth_t0_m"], 0)
        v["blB_step"] = f(bb["step_kg_m3"], 4)
        v["blB_lim"] = f"{bb['colour_limit_m']:g}"
        v["blB_sponge"] = word(int(bb["sponge_cells"]))
        v["blB_share"] = f"{100 * bb['sigma_min_share']:.0f}"
        for r in ("F3", "F2"):
            sr = bb["saturation"][r]
            v[f"blB_sat_{r}"], v[f"blB_sat_{r}_t"] = pct(sr["worst"]), f(sr["worst_t_h"], 1)
    if val and h:
        check_prose(val, h, mr)
    return v


PD0_WORD = {"inconclusive_T_fin": "inconclusive with T_fin", "not_distinguishable": "not distinguishable",
            "supported": "supported", "refuted": "refuted", "alternative": "alternative (Kelvin)",
            "mixed": "mixed", "undefined": "undefined (zero denominator)"}


def _outcome(o: str) -> str:
    return PD0_WORD[o]


def _expect(cond: bool, what: str) -> None:
    """The README prose states this in fixed words; refuse to render if the JSON says otherwise."""
    if not cond:
        raise ValueError(f"README prose no longer matches the metrics: expected {what}")


def check_prose(val: dict, h: dict, mr: dict | None) -> None:
    """Signs and outcomes behind every fixed word of README.template.md (final review, Minor 11)."""
    _expect(not val["V2a"]["pass"] and not val["V2b"]["pass"] and val["V3"]["pass"] and not val["V4"]["pass"],
            "V2a, V2b and V4 fail and V3 passes")
    _expect(val["V3"]["sense"]["direction"] == "east", "V3 travelling east along the south coast (cyclonic)")
    v4 = val["V4"]["V4a"]
    _expect(v4["F1_out_J"] < 0 < v4["V4_out_J"], "V4: inward at 200 m, outward at 400 m")
    a = h["H4a"]
    _expect(a["outcome"] == "inconclusive_T_fin" and a["bands_out_J"]["low"] < 0 and a["resolved_sum_out_J"] < 0,
            "H4a inconclusive, low band and resolved sum inward")
    b = h["H4b"]
    _expect(b["i"]["met"] and not b["ii"]["met"], "H4b(i) met and H4b(ii) not met")
    c = h["H4c"]
    _expect(c["outcome"] == "not_distinguishable" and c["ii"]["F_res_out_F1_J"] < 0 and c["ii"]["F_res_out_F3_J"] < 0,
            "H4c not distinguishable with both resolved fluxes inward")
    thr = c["threshold"]
    _expect(1 / thr < c["i"]["ratio"] < thr, "H4c(i) in neither tail")
    dr, df = h["H4d"]["R1k_vs_R1"], h["H4d"]["F1_vs_F1c"]
    _expect(dr["outcome"] == "supported" and df["outcome"] == "refuted" and df["numerator_J"]["kpp"] < 0,
            "H4d supported in the bay and refuted in the fjord with an inward F1 numerator")
    _expect(dr["W_wind_J"]["kpp"] < dr["W_wind_J"]["closure"] and dr["numerator_J"]["kpp"] > dr["numerator_J"]["closure"],
            "H4d bay: less wind work and a larger radiated flux with KPP")
    _expect(dr["post_hoc_amendment1_numerator"]["literal_rule_would_read"] != "supported",
            "the Amendment 1 numerator would not support H4d in the bay")
    k = h["post_hoc"]["band_flux_per_strait"]["F1c"]["total"]["resolved_out_J"]["kelvin_subinertial"]
    _expect(k > 0, "F1c: the Kelvin band carries the net outward")
    al = h["post_hoc"]["arneborg_liljebladh_net_in_J"]
    _expect(al["F1"]["total"]["F_bc_in_J"] > 0 and al["F1c"]["total"]["F_bc_in_J"] < 0
            and al["R1"]["total"]["F_bc_in_J"] < 0 and al["R2"]["total"]["F_bc_in_J"] < 0,
            "A&L band: enters in F1, leaves in F1c, R1 and R2")
    br = h["post_hoc"]["bridge"]
    for run in ("R1k", "R2"):
        s_ = br[run]["flux_share_signed"]
        _expect(s_["near_inertial"] > s_["low"] > 0, f"{run}: outward and mostly near-inertial")
    _expect(br["F1"]["flux_share_signed"]["low"] < 0, "F1: low-band share inward")
    if mr:
        R = mr["runs"]
        tb = R["F1"]["total"]["band_out_J"]
        _expect(R["F1"]["low_band_out_J"] < 0 and R["F1"]["eta"] < 0, "moving reference: F1 low band and eta still negative")
        _expect(abs(tb["low"]) / sum(abs(x) for x in tb.values()) > h["general"]["thresholds"]["H4a"]["low_band_max_frac"]
                and sum(tb[x] for x in BANDS) > 0 > a["resolved_sum_out_J"],
                "moving reference: low band still dominant, F1 resolved sum changes sign")
        _expect(R["F1c"]["low_band_out_J"] > 0 > h["post_hoc_h4d_resolved_bands"]["F1c"]["low_out_J"],
                "moving reference: the F1c low band turns outward")


def _sci(x: float) -> str:
    """1.09e-09 -> '1.1 × 10⁻⁹' (superscript exponent)."""
    m, e = f"{x:.1e}".split("e")
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{m} × 10{str(int(e)).translate(sup)}"


def render(allow_pending: bool) -> str:
    text = TEMPLATE.read_text()
    v = values()
    missing = sorted({k for k in re.findall(r"\{\{(\w+)\}\}", text) if k not in v})
    if missing and not allow_pending:
        raise KeyError(f"no value for placeholders: {missing}")
    return re.sub(r"\{\{(\w+)\}\}", lambda m: typographic_minus(v.get(m.group(1), PENDING)), text)


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
    print(f"wrote {OUT} ({len(new.split())} words, {new.count(PENDING)} pending values)")
