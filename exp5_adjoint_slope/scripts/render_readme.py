#!/opt/anaconda3/envs/pyocean/bin/python
"""Render exp5_adjoint_slope/README.md from README.template.md and the committed metrics files.

Same mechanism as Experiments 1-4 (exp4_arctic_bayiw/scripts/render_readme.py): every number in
the README comes from a committed file (figures/*.json, or cases.yaml for the configuration),
never typed by hand. Placeholders are ``{{name}}``. An unknown placeholder raises KeyError; a
metrics file that does not exist raises FileNotFoundError, unless ``--allow-pending`` (its values
then render as "[pending]"). tests/test_render_readme.py fails if the template holds a digit
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
SOURCES = ("validation", "sensitivity", "hypotheses", "mean_state", "diagnostics", "crosscheck")

CITATION_YEARS = ("1963", "1991", "1995", "1999", "2020", "2000", "2004", "2006", "2007", "2010", "2011", "2012",
                  "2014", "2015", "2016", "2018", "2019", "2022", "2023", "2024", "2026", "2027")
NUMBER_WHITELIST = (
    r"\b(?:G0|GA|GB|GM30|D90|E0z|E[0-4]|EA2[1-9]|EA30|EA21–EA30|F0b?|FC[1-3]|FP[1-3]|GS|GSU|GV1|GV2d|V3E|H5A|H5BC)\b",
    r"\bV[1-5][de]?\b", r"\bH5[a-c]?\b", r"\bExp(?:eriment|\.)? [1-5]\b", r"\bPD[0-2]b?\b",
    r"^#+ [0-9]+(?:\.[0-9]+)*\.?", r"\b[Ee]qs?\. [0-9]+", r"§[0-9]+(?:\.[0-9]+)*",
    r"\b(?:Figure|Fig\.|Figures) [0-9](?:[a-d])?(?:[–-][0-9])?", r"\bTables? [1-3]\b",
    r"^[1-9]\. ", r"\(i{1,2}\)", r"\(\w\)", r"\b(?:class )?\([a-d]\)", r"\bθ₀\b", r"\bκᵥ\b",
    r"\b[23]-D\b", r"< (?:1|30) %", r"> 50 %", r"2σ_Δ", r"\b30 days\b", r"< 0\.1\b", r"±2 s(?:tandard errors|\.e\.)", r"\bJ > 0\b", r"\bp\. [0-9]+(?:[–-][0-9]+)?", r"\bv5\.2\b",
    r"\b(?:" + "|".join(CITATION_YEARS) + r")[a-b]?\b",
)


def strip_allowed(text: str) -> str:
    """The template minus placeholders, code, URLs, comments, the References section and the
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


def _load(name: str, allow_pending: bool):
    p = FIG / f"{name}.json"
    if p.exists():
        return json.loads(p.read_text())
    if allow_pending:
        return None
    raise FileNotFoundError(f"metrics file missing: {p}")


def _sci(x: float, nd: int = 1) -> str:
    m, e = f"{x:.{nd}e}".split("e")
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{m} × 10{str(int(e)).translate(sup)}"


def gw(x_tw: float, nd: int = 2, sign: bool = False) -> str:
    """TW -> GW string."""
    return f"{1e3 * x_tw:+.{nd}f}" if sign else f"{1e3 * x_tw:.{nd}f}"


def pct(x: float, nd: int = 1) -> str:
    """Fraction -> percent (no % sign); scientific below 0.1 %."""
    return f"{100 * x:.{nd}f}" if 100 * x >= 0.1 else _sci(100 * x)


WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")


def word(n: int, cap: bool = False) -> str:
    w = WORDS[n] if 0 <= n < len(WORDS) else str(n)
    return w.capitalize() if cap else w


def verdict(ok: bool) -> str:
    return "passes" if ok else "fails"


def typographic_minus(s: str) -> str:
    return re.sub(r"(?<![\w.])-(?=\d)", "−", s)


def values(allow_pending: bool = False) -> dict[str, str]:
    d = {n: _load(n, allow_pending or n == "crosscheck") for n in SOURCES}
    v: dict[str, str] = {}
    import yaml
    reg = yaml.safe_load((EXP / "cases.yaml").read_text())
    c10, c2 = reg["configs"]["chan10"], reg["configs"]["chan2p5"]
    v["dx_coarse_km"], v["dx_fine_km"] = f"{c10['dx_m'] / 1e3:g}", f"{c2['dx_m'] / 1e3:g}"
    v["Lx_km"], v["Ly_km"], v["H_m"] = f"{c10['Lx_m'] / 1e3:.0f}", f"{c10['Ly_m'] / 1e3:.0f}", f"{c10['depth_total_m']:.0f}"
    v["nr"] = str(c10["nr"])
    dt_c = float(re.search(r"deltaT=([0-9.]+)", (EXP / "config/chan10/input/data").read_text())[1])
    v["g0_years"] = f"{reg['cases']['G0']['n_timesteps'] * dt_c / (365 * 86400):.0f}"
    cs = reg["cases"]
    for n in ("GA", "GB", "GM30", "EA21", "H5BC_21_p", "FC1"):
        if n in cs and "window_days" in cs[n]:
            v[f"win_{n.split('_')[0]}"] = f"{cs[n]['window_days']:g}"
    dt_f = float(re.search(r"deltaT=([0-9.]+)", (EXP / "config/chan2p5/input/data").read_text())[1])
    v["win_FC"] = f"{cs['FC1']['n_timesteps'] * dt_f / 86400:g}"
    v["win_GV1"] = f"{cs['GV1']['window_days']:g}"
    from exp5_adjoint_slope.analysis import h5a as _h
    v["h5a_ymax_km"], v["h5a_y0_km"] = f"{_h.Y_MAX / 1e3:g}", f"{_h.Y0 / 1e3:g}"
    v["w1_km"], v["w2_km"], v["w3_km"] = (f"{w / 1e3:g}" for w in _h.WIDTHS)
    if d["sensitivity"]:
        v["fp_amp"] = f"{max(abs(x) for x in d['sensitivity']['fp_pattern']):g}"
    if d["validation"]:
        A, B = d["validation"]["A"], d["validation"]["B"]
        v1 = A["V1d"]
        v["v1_n"], v["v1_worst"] = str(v1["n_points"]), _sci(v1["value"])
        tl = [p["rel_err_tlm"] for p in v1["points"] if p["rel_err_tlm"] is not None]
        v["v1_tlm_worst"] = _sci(max(tl)) if tl else PENDING
        v["v1_verdict"] = verdict(v1["pass"])
        for u, name in (("fu", "tau"), ("theta", "theta"), ("diffkr", "kappa")):
            x = A[f"V2d_{u}"]["value"]
            v[f"v2d_{name}"] = f"{x:.2f}" if x >= 0.01 else _sci(x)
            v[f"v2d_{name}_verdict"] = verdict(A[f"V2d_{u}"]["pass"])
            v[f"v3d_{name}"] = pct(A[f"V3d_{u}"]["value"])
            v[f"v3d_{name}_verdict"] = verdict(A[f"V3d_{u}"]["pass"])
        v["v2d_tau_plus"], v["v2d_tau_minus"] = gw(A["V2d_fu"]["dj_plus"], 2, True), gw(A["V2d_fu"]["dj_minus"], 2, True)
        v["v3d_tau_pred"] = gw(A["V3d_fu"]["pred"])
        v["v2e"], v["v2e_verdict"] = f"{B['V2e']['value']:.3f}", verdict(B["V2e"]["pass"])
        e = B["V3e"]
        v["v3e"], v["v3e_verdict"] = pct(e["value"]), verdict(e["pass"])
        v["v3e_pred"], v["v3e_fwd"], v["v3e_se"] = gw(e["pred"]), gw(e["fwd"]), gw(e["se"])
        h = B["H5a"]
        for k, key in (("37.5", "w1"), ("18.75", "w2"), ("75", "w3")):
            w = h["widths"][k]
            v[f"h5a_{key}_i"], v[f"h5a_{key}_ii"] = gw(abs(w["dj_uniform"])), gw(abs(w["dj_dipole"]))
            v[f"h5a_{key}_verdict"] = "holds" if w["pass"] else "does not hold"
        w = h["widths"]["37.5"]
        v["h5a_paired"], v["h5a_paired_se"] = gw(w["paired_diff"]), gw(w["paired_se"])
        v["h5a_robust"] = "robust" if h["robust"] else "not robust"
        dr = h["direct_37.5km_years26-30"]
        v["h5a_fwd_i"], v["h5a_fwd_ii"] = gw(abs(dr["i"]["fwd"][0])), gw(abs(dr["ii"]["fwd"][0]))
        v["h5a_v2e_i"] = f"{dr['i']['v2e']['value']:.2f}"
        v["h5a_low_ii"] = f"{100 * (1 - abs(dr['ii']['pred']) / abs(dr['ii']['fwd'][0])):.0f}"
        v["h5a_low_i"] = f"{100 * (1 - abs(dr['i']['pred']) / abs(dr['i']['fwd'][0])):.0f}"
        c = B["H5c"]
        v["h5c_verdict"] = c["verdict"]
        lead = dict(zip(c["leads_days"], c["sens_mean"]))
        v["h5c_s0"], v["h5c_s30"], v["h5c_s60"] = (gw(lead[k], 0) for k in (0.0, 30.0, 60.0))
        p = B["H5b_coarse_prediction"]
        v["h5b_pred"], v["h5b_pred_se"] = gw(p["dj_mean"]), gw(p["dj_se"])
        v["h5b_pred_n_pos"] = str(sum(x > 0 for x in p["members"]))
        v["n_members"] = str(len(p["members"]))
    if d["mean_state"]:
        J = d["mean_state"]["J"]
        v["j_total"], v["j_resolved"] = gw(J["J_TW"], 2, True), gw(J["J_resolved_TW"], 2, True)
        v["j_bolus"], v["j_redi"] = gw(J["J_bolus_TW"], 2, True), gw(J["J_redi_TW"], 2, True)
        v["shelf_theta"] = f"{J['shelf_theta_degC']:.2f}"
        v["y_qb_km"] = f"{d['mean_state']['y_qb_m'] / 1e3:.0f}"
    if d["diagnostics"]:
        g = d["diagnostics"]["adjoint_growth"]
        for k in ("30", "90", "180", "360"):
            v[f"grow_{k}"] = _sci(g[k]["max_abs_dJdtheta0"])
        v["n_E"] = word(len(d["diagnostics"]["ensemble_J_TW"]), cap=True)
        v["n_E_word"] = word(len(d["diagnostics"]["ensemble_J_TW"]))
        v["theta_noise"] = _sci(d["diagnostics"]["theta_noise_degC"], 0)
        s = d["diagnostics"]["ensemble_std_TW"]
        v["ens_std_m1"], v["ens_std_m2"], v["ens_std_m3"] = _sci(s[0] * 1e3), gw(s[1]), gw(s[2])
        fd = d["diagnostics"]["d90_fd"]
        v["d90_ref"], v["d90_plus"], v["d90_minus"] = (gw(fd[k], 2, True) for k in ("fc_ref", "fc_plus", "fc_minus"))
    if d["crosscheck"]:
        x = d["crosscheck"]
        h = x["H5b"]
        n = word(h["n"])
        if h["indicative"]:
            v["h5b_verdict"] = f"indicative only (N = {h['n']})"
        elif not h["detectable"]:
            v["h5b_verdict"] = f"not detectable with N = {h['n']}"
        else:
            v["h5b_verdict"] = "holds" if h["pass"] else "fails"
        v["h5b_n"] = n
        v["h5b_dj"], v["h5b_2sd"] = gw(h["dj"]), gw(2 * h["sigma_delta"])
        v["h5b_ratio"] = f"{h['ratio']:.2f}"
        jc, jp = list(x["J_control_TW"].values()), list(x["J_pert_TW"].values())
        v["h5b_jc_min"], v["h5b_jc_max"] = gw(min(jc)), gw(max(jc))
        v["h5b_jp_min"], v["h5b_jp_max"] = gw(min(jp)), gw(max(jp))
        f = list(x["V5_eddy_fraction"].values())
        v["v5_min"], v["v5_max"] = (f"{round(100 * y, -1):.0f}" if abs(y) >= 1 else f"{round(100 * y):.0f}"
                                    for y in (min(f), max(f)))
        v["v5_verdict"] = verdict(min(f) > 0.5)
        if x.get("F1_equilibrium"):
            v["eke_peak"] = _sci(max(x["EKE_spinup_m2s2"]["F1"]))
            v["eke_final"] = _sci(x["F1b_equilibrium"]["last"])
            v["eke_drift_f1"] = f"{100 * x['F1_equilibrium']['rel_change']:.0f}"
            v["eke_drift_f1b"] = f"{100 * x['F1b_equilibrium']['rel_change']:.0f}"
    return v


def render(allow_pending: bool) -> str:
    text = TEMPLATE.read_text()
    v = values(allow_pending)
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
