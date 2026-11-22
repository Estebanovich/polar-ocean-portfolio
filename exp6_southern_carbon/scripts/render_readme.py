#!/opt/anaconda3/envs/pyocean/bin/python
"""Render exp6_southern_carbon/README.md from README.template.md and the committed metrics files.

Same mechanism as Experiments 1-5 (exp5_adjoint_slope/scripts/render_readme.py): every number in
the README comes from a committed file (figures/*.json, or cases.yaml and the analysis constants for
the configuration), never typed by hand. Placeholders are ``{{name}}``. An unknown placeholder
raises KeyError; a metrics file that does not exist raises FileNotFoundError, unless
``--allow-pending`` (its values then render as "[pending]"). tests/test_render_readme.py fails if
the template holds a digit that is not a whitelisted label (NUMBER_WHITELIST).

Usage: render_readme.py [--allow-pending] [--check]   (--check: exit 1 if README.md would change)
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

EXP = Path(__file__).resolve().parents[1]
REPO = EXP.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
FIG = EXP / "figures"
TEMPLATE = EXP / "README.template.md"
OUT = EXP / "README.md"
PENDING = "[pending]"
SOURCES = ("spinup", "validation", "hypotheses", "figdata", "crosscheck")
RHO0 = 1035.0
MONTHS = ("January", "February", "March", "April", "May", "June", "July", "August", "September",
          "October", "November", "December")

CITATION_YEARS = ("1982", "1984", "1990", "2001", "2005", "2006", "2007", "2008", "2011", "2013", "2014", "2015", "2016", "2017", "2018",
                  "2019", "2020", "2023", "2024", "2025", "2026")
NUMBER_WHITELIST = (
    r"\b(?:G0s?|GC|GT|F0|FT|T10|T2|R_CP)\b", r"\bH6[a-c]\b", r"\bV6[a-e]\b", r"\bPD0\b",
    r"\bExp(?:eriment|\.)? [1-6]\b", r"\bCO₂\b", r"\bPO₄\b", r"\bpCO₂\b", r"\bG₂₀₀\b", r"\bΔDIC₁₀₀₀\b",
    r"\bSSP3-7\.0\b", r"\bCMIP6\b", r"\b(?:UKESM1|ACCESS-OM2)\b", r"\b(?:SST|DIC)\b",
    r"^#+ [0-9]+(?:\.[0-9]+)*\.?", r"\b[Ee]qs?\. [0-9]+(?:(?:, | and |–)[0-9]+)*", r"§[0-9]+(?:\.[0-9]+)*",
    r"\b(?:Figure|Fig\.|Figures|Figs\.) [0-9](?:[–-][0-9])?", r"\bTables? [1-3]\b",
    r"^[1-9]\. ", r"\(i{1,3}\)", r"\(\w\)", r"\b(?:class )?\([a-d]\)", r"±25 %", r"\b2σ\b", r"\[0\.5, 2\]",
    r"\bupper-200 m\b", r"\b(?:200|1000) m\b", r"\b1 January\b",   # definitions (Jarnikova et al. 2025 eqs. 1, 4)
    r"< 1 %", r"> 0\b", r"\bS\(t_b\)", r"\b(?:m|s)⁻[¹²³]", r"\bkg⁻¹", r"\byr⁻¹", r"\bm⁻²", r"\bpkg/\w+",
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


def sci(x: float, nd: int = 1) -> str:
    m, e = f"{x:.{nd}e}".split("e")
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{m} × 10{str(int(e)).translate(sup)}"


def num(x: float, nd: int = 3, sign: bool = False) -> str:
    """Fixed point when readable, scientific below 10^-nd."""
    if x != 0 and abs(x) < 10 ** -nd:
        return ("+" if sign and x > 0 else "") + sci(x)
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def verdict(ok: bool) -> str:
    return "passes" if ok else "fails"


def word(n: int) -> str:
    w = ("zero", "one", "two", "three", "four", "five", "six")
    return w[n] if 0 <= n < len(w) else str(n)


def typographic_minus(s: str) -> str:
    return re.sub(r"(?<![\w.])-(?=\d)", "−", s)


def config_values() -> dict[str, str]:
    import yaml
    from exp6_southern_carbon.analysis import collect
    from exp6_southern_carbon.analysis.hypotheses import BRANCH_YEARS
    reg = yaml.safe_load((EXP / "cases.yaml").read_text())
    g, f = reg["configs"]["g50"], reg["configs"]["f10"]
    cs = reg["cases"]
    v = {"dx_coarse_km": f"{g['dx_m'] / 1e3:g}", "dx_fine_km": f"{f['dx_m'] / 1e3:g}",
         "Lx_km": f"{g['nx'] * g['dx_m'] / 1e3:g}", "Ly_km": f"{g['ny'] * g['dx_m'] / 1e3:g}", "nr": str(g["nr"]),
         "dt_s": f"{g['deltaT']:g}", "accel": f"{cs['G0']['deltaT_tracer'] / g['deltaT']:g}",
         "g0s_years": str(cs["G0s"]["n_years"]), "gt_years": str(cs["GT"]["n_years"]),
         "branch_len": str(cs["GT_w1900"]["n_years"]), "fine_len": str(cs["FT_1900_c1"]["n_years"]),
         "wind_factor": f"{cs['GT_w1900']['wind_factor']:g}",
         "win_lo": str(collect.WINDOW[0]), "win_hi": str(collect.WINDOW[1]),
         "R_rows": str(collect.J_EDGE - 1), "R_edge_km": f"{collect.J_EDGE * g['dx_m'] / 1e3:g}",
         "r_cp": f"{collect.R_CP:g}", "n_members": word(len([c for c in cs if c.startswith("FT_1900_c")]))}
    for k, y in enumerate(BRANCH_YEARS, start=1):
        v[f"y{k}"] = str(y)
    v["y_first"], v["y_last"] = "1850", "2100"
    return v


def values(allow_pending: bool = False) -> dict[str, str]:
    d = {n: _load(n, allow_pending or n == "crosscheck") for n in SOURCES}
    v = config_values()
    yrs = [v[f"y{k}"] for k in range(1, 5)]
    if d["spinup"]:
        s = d["spinup"]["G0"]
        a = s["v6a"]
        v["g0_years"] = f"{s['years']:.0f}"
        v["v6a_dflux"] = sci(abs(a["d_flux_century"]))
        v["v6a_ddic"] = sci(abs(a["d_deep_dic_century"]) * 1e6 / RHO0)
        v["v6a_verdict"] = verdict(a["pass"])
        v["moc_R"] = f"{s['moc_R_sv']:.2f}"
        v["moc_edge_frac"] = f"{100 * s['moc_R_sv'] / s['moc_max_sv']:.0f}"
        v["upw_frac"] = f"{100 * s['upwelling_fraction_R']:.0f}"
        v["mld_R"] = f"{s['mld_winter_R_m']:.0f}"
        v["mld_month"] = MONTHS[int(s["mld_month_max"]) - 1]
        v["flux_g0"] = f"{s['flux_last_mol_m2_yr']:.2f}"
        if "G0s" in d["spinup"]:
            b = d["spinup"]["G0s"]["v6a"]
            v["g0s_dflux"], v["g0s_ddic"] = sci(abs(b["d_flux_century"])), sci(abs(b["d_deep_dic_century"]) * 1e6 / RHO0)
    if d["validation"]:
        x = d["validation"]["V6b"]
        v["v6b_max"], v["v6b_verdict"] = sci(x["max_frac"]), verdict(x["pass"])
    if d["figdata"]:
        gR = d["figdata"]["sections"]["g200_R"]
        v["g200_first"], v["g200_last"] = num(gR[yrs[0]]), num(gR[yrs[-1]])
        v["g200_drop"] = f"{100 * (1 - gR[yrs[-1]] / gR[yrs[0]]):.0f}"
    if d["hypotheses"]:
        h = d["hypotheses"]
        a = h["h6a"]
        v["h6a_verdict"] = verdict(a["pass"])
        v["dS"], v["drift_gc"] = num(a["dS"], sign=True), num(a["drift_gc"])
        v["dS_pct"] = f"{a['dS_pct_of_S1900']:.0f}"
        v["zero_cross"] = (f"changes sign in the {a['zero_crossing_year']} branch" if a["zero_crossing_year"]
                           else "keeps its sign in every branch")
        b = h["h6b"]
        v["h6b_verdict"], v["h6b_n"] = verdict(b["pass"]), word(b["n_within"])
        for k, y in enumerate(yrs, start=1):
            t = h["terms"][y]
            v[f"s_gt_{k}"], v[f"s_nat_{k}"] = num(h["s_gt"][y], sign=True), num(a["s_nat"][y], sign=True)
            v[f"s_ant_{k}"] = num(a["s_ant"][y], sign=True)
            v[f"g200_{k}"] = num(t["parent"]["g200"])
            v[f"ratio_{k}"] = f"{b['ratio'][y]:.2f}"
            for key in ("mld", "moc", "bio", "temp"):
                v[f"{key}_{k}"] = f"{t[key] * 1e-9:+.1f}" if key in t else PENDING
            v[f"meas_{k}"] = f"{h['meas'][y] * 1e-9:+.1f}"
            if "meas_parts" in t:                                       # post hoc (22 nov): term by term
                m, off = t["meas_parts"], t["mld"] + t["moc"]
                v[f"phys_off_{k}"], v[f"phys_meas_{k}"] = f"{off * 1e-9:+.1f}", f"{(m['adv'] + m['mix']) * 1e-9:+.1f}"
                v[f"phys_ratio_{k}"] = f"{off / (m['adv'] + m['mix']):.2f}"
                v[f"mix_meas_{k}"], v[f"bio_meas_{k}"] = f"{m['mix'] * 1e-9:+.1f}", f"{m['bio'] * 1e-9:+.1f}"
    if d["crosscheck"]:
        c = d["crosscheck"]["H6c"]
        n = int(c["n"])
        if c["indicative"]:
            v["h6c_verdict"] = f"is indicative only (N = {n})"
        elif not c["detectable"]:
            v["h6c_verdict"] = f"fails: the change is not detectable with N = {n}"
        else:
            v["h6c_verdict"] = verdict(c["pass"])
        v["h6c_dS"], v["h6c_2sd"] = num(c["dS"], sign=True), num(2 * c["sigma"])
        v["h6c_ratio"], v["h6c_n"] = f"{c['ratio']:.2f}", word(n)
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
