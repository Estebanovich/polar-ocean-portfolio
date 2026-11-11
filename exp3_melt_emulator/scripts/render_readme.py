"""Render exp3_melt_emulator/README.md from README.template.md and the figures/*.json metrics.

Same mechanism as Experiments 1 and 2 (exp2_isomip_cavity/scripts/render_readme.py): every
result number in the README comes from a metrics file in figures/, never typed by hand.
Placeholders are ``{{name}}``; an unknown placeholder raises. Values whose source file does
not exist render as "[pending]" only with ``--allow-pending``. Standard library only (runs
under pyocean and anemoi).

Usage: render_readme.py [--allow-pending] [--check]   (--check: exit 1 if README.md would change)
"""
from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from pathlib import Path

EXP = Path(__file__).resolve().parents[1]
FIG = EXP / "figures"
TEMPLATE = EXP / "README.template.md"
OUT = EXP / "README.md"
PENDING = "[pending]"
METHODS = ("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
ZC = ("250", "400", "550")
KINDS = ("loro", "zc_out", "warm_out", "final")

# Digits allowed in the template outside placeholders (test_every_number_from_json). They are
# labels, citation years and two design constants of the Exp. 2 ensemble, never results.
CITATION_YEARS = ("2027", "2003", "2008", "2009", "2010", "2016", "2017", "2018", "2019", "2020",
                  "2021", "2022", "2023", "2025", "2026")
DESIGN_CONSTANTS = ("4 km", "d = 300 m")   # Exp. 2 grid spacing and ramp thickness
NUMBER_WHITELIST = (
    r"\bM[0-7]\b", r"\bH[1-3][a-c]?\b", r"\bN[1-3]\b", r"\bR[0-9]{1,2}\b",
    r"\bExp(?:eriment|\.)? [1-5]\b", r"\bOcean[01]\b", r"\bocean0_4km\b", r"\bISMIP6\b",
    r"\bsha256\b", r"^#+ [0-9]+(?:\.[0-9]+)*\.?",       # section headings
    r"\bEqs?\. [0-9]+(?:[–-][0-9]+)?(?:(?:, | and )[0-9]+(?:[–-][0-9]+)?)*",
    r"§[0-9]+(?:\.[0-9]+)*", r"\b(?:Figure|Fig\.|Figures) [0-9][a-d]?(?:[–-][0-9])?",
    r"\bTable [0-9]+", r"\bAmendment [0-9]\b", r"\([a-d][12]\)",
    r"\b[a-d][12]\b(?=\)| and)",                       # panel labels "(b1) and (b2)"
    r"\([1-4]\) (?=[A-Z])",                            # numbered panels of Figure 0
    r"max\(TF, 0\)",                                   # TF+ = max(TF, 0)
    r"\b1 − error/error",                              # skill = 1 - err/err_M0
    r"(?<![0-9.])\+[12](?![0-9.%])",                  # the +1 / +2 scenarios of LB25
    r"\b(?:" + "|".join(CITATION_YEARS) + r")[a-b]?\b",
)


def strip_allowed(text: str) -> str:
    """The template minus placeholders, code spans, URLs, comments, the References section,
    the design constants and the whitelisted tokens: what is left must contain no digit."""
    text = text.split("\n## References", 1)[0]
    text = re.sub(r"<!--.*?-->", " ", text, flags=re.S)
    text = re.sub(r"```.*?```", " ", text, flags=re.S)   # fenced code (commands)
    text = re.sub(r"\{\{\w+\}\}", " ", text)
    text = re.sub(r"`[^`\n]*`", " ", text)
    text = re.sub(r"\]\([^)]*\)", "]", text)            # link targets (figure paths, URLs)
    text = re.sub(r"https?://\S+", " ", text)
    for c in DESIGN_CONSTANTS:
        text = text.replace(c, " ")
    for pat in NUMBER_WHITELIST:
        text = re.sub(pat, " ", text, flags=re.M)
    return text


def _load(name: str):
    p = FIG / name
    return json.loads(p.read_text()) if p.exists() else None


def f(x, nd: int = 2) -> str:
    return f"{x:.{nd}f}"


def pct(x, nd: int = 1, sign: bool = False) -> str:
    """Fraction -> percent string (no % sign)."""
    return f"{100 * x:+.{nd}f}" if sign else f"{100 * x:.{nd}f}"


def g(x) -> str:
    """Compact form of a hyperparameter (1e-08 -> 1e-8, 0.001 -> 1e-3, 1000.0 -> 1e3)."""
    if x == int(x) and abs(x) < 1e3:
        return str(int(x))
    m, e = f"{x:.0e}".split("e")
    return f"{m}e{int(e)}"


def typographic_minus(s: str) -> str:
    """Hyphen-minus before a digit -> U+2212 minus; hyphens inside case names (ens_T-0.5)
    and ranges written with an en dash are left alone."""
    return re.sub(r"(?<![\w.])-(?=\d)", "−", s)


def _case(c: str) -> str:
    """ens_T2.0_zc400_s1.0 -> Tbot = 2.0 °C, z_c = 400 m, s = 1.0."""
    m = re.match(r"ens_T(-?[\d.]+)_zc(\d+)_s([\d.]+)", c)
    return f"Tbot = {m[1]} °C, z_c = {m[2]} m, s = {m[3]}" if m else c


def _md_table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def values() -> dict[str, str]:
    v: dict[str, str] = {}
    # ---- frozen configuration: data, folds, methods, thresholds ----------------------------
    cfg = _load("protocol_config.json")
    if cfg:
        folds = {fo["name"]: fo for fo in cfg["folds"]}
        by_kind: dict[str, list] = {}
        for fo in cfg["folds"]:
            by_kind.setdefault(fo["kind"], []).append(fo)
        final, prof = folds["final"], folds["profile_test"]
        train = final["train"]
        v["n_train"] = str(len(train))
        v["n_test"] = str(len(final["test"]))
        v["n_runs_total"] = str(len(train) + len(final["test"]) + len(prof["test"]))
        v["n_factorial"] = str(sum(c.endswith("_s1.0") for c in train))
        v["n_geometry"] = str(sum(not c.endswith("_s1.0") for c in train))
        v["n_n1"] = str(len(cfg["n1_fail_cases"]))
        v["n1_list"] = "; ".join(_case(c) for c in cfg["n1_fail_cases"])
        tb = sorted({float(re.match(r"ens_T(-?[\d.]+)_", c)[1]) for c in train if c.endswith("_s1.0")})
        v["fact_Tbot"] = ", ".join(f"{t:g}" for t in tb)
        zc = sorted({int(re.search(r"_zc(\d+)_", c)[1]) for c in train if c.endswith("_s1.0")})
        v["fact_zc"] = ", ".join(str(z) for z in zc)
        v["zc_a"], v["zc_b"], v["zc_c"] = (str(z) for z in zc)
        v["s_values"] = " and ".join(sorted({re.search(r"_s([\d.]+)$", c)[1] for c in train
                                             if not c.endswith("_s1.0")}))
        v["test_Tbot"] = re.match(r"ens_T(-?[\d.]+)_", final["test"][0])[1]
        lo = by_kind["loro"]
        v["loro_n"] = str(len(lo))
        v["loro_train"] = str(len(lo[0]["train"]))
        zo = sorted(by_kind["zc_out"], key=lambda fo: fo["name"])
        v["zc_out_n"] = str(len(zo))
        v["zc_out_train"] = "/".join(str(len(fo["train"])) for fo in zo)
        v["zc_out_test"] = "/".join(str(len(fo["test"])) for fo in zo)
        zt = [[int(re.search(r"_zc(\d+)_", c)[1]) for c in fo["test"]] for fo in zo]
        v["zc_out_ranges"] = ", ".join(f"{min(z)}–{max(z)} m" for z in zt)
        wo = folds["warm_out"]
        tw = lambda cs: [float(re.match(r"ens_T(-?[\d.]+)_", c)[1]) for c in cs]  # noqa: E731
        v["warm_train_max"] = f"{max(tw(wo['train'])):g}"
        v["warm_test_min"] = f"{min(tw(wo['test'])):g}"
        v["warm_n_train"], v["warm_n_test"] = str(len(wo["train"])), str(len(wo["test"]))
        v["n_folds_cv"] = str(len(cfg["folds"]) - 1)
        v["profile_methods"] = ", ".join(cfg["profile_test_methods"])
        th = cfg["thresholds"]
        v["gl_frac"] = pct(th["gl_frac"], 0)
        v["H3a_min_better"] = str(th["H3a"]["min_runs_better"])
        v["H3b_ratio"] = f(th["H3b"]["max_ratio_M7_over_M6"])
        v["H3c_spread"] = pct(th["H3c"]["spread_max"], 0)
        v["H3c_err"] = pct(th["H3c"]["error_min"], 0)
        v["ock_spread"] = pct(th["H3c"]["ockenden_flag"], 0)
        m = cfg["methods"]
        v["c_const"] = f"{m['M1']['constants']['c']:.4f}"
        v["M4_degrees"] = " or ".join(str(d) for d in m["M4"]["degrees"])
        k = len(m["M4"]["input_features"])
        v["M4_n_terms"] = str((k + 1) * (k + 2) // 2)
        m5 = m["M5"]
        v["M5_var"] = pct(m5["var_fraction"], 0)
        v["M5_const_b"] = f"{g(m5['const_bounds'][0])}–{g(m5['const_bounds'][1])}"
        v["M5_rbf_b"] = f"{g(m5['rbf_bounds'][0])}–{g(m5['rbf_bounds'][1])}"
        v["M5_white0"] = g(m5["white_init"])
        v["M5_white_b"] = f"[{g(m5['white_bounds'][0])}, {g(m5['white_bounds'][1])}]"
        v["M5_restarts"] = str(m5["n_restarts"])
        v["M5_kcap"] = m5["k_cap"].replace(" - ", " − ")
        m6, m7 = m["M6"], m["M7"]
        v["nn_members"] = str(m6["n_members"])
        v["nn_seeds"] = f"{m6['seeds'][0]}–{m6['seeds'][-1]}"
        v["nn_hidden"] = "-".join(str(h) for h in m6["hidden"])
        v["nn_lr"] = g(m6["lr"])
        v["nn_batch"] = str(m6["batch_size"])
        v["nn_max_epochs"] = str(m6["max_epochs"])
        v["nn_patience"] = str(m6["patience"])
        v["nn_val"] = str(m6["n_val"])
        v["nn_var_floor"] = g(m6["var_floor"])
        v["M6_n_feat"] = str(len(m6["input_features"]))
        v["M6_feats"] = ", ".join(m6["input_features"])
        v["M7_n_feat"] = str(len(m7["input_features"]))
        v["M7_feats"] = ", ".join(m7["input_features"])
        v["M7_n0"] = f"{m7['n0']:g}"
        assert m7["tf_kind"] == "auto" and m7["tf_kind_folds"] is None
    # ---- noise floor --------------------------------------------------------------------------
    nf = _load("noise_floor.json")
    if nf:
        v["nf_rms"] = f(nf["pooled_median_rms"], 3)
        v["nf_ire"] = pct(nf["median_abs_integrated_rel_diff"], 1)
        v["nf_min"], v["nf_max"] = f(nf["min_rms"], 3), f(nf["max_rms"], 3)
        v["nf_n"] = str(nf["n_runs"])
        v["nf_months"] = re.search(r"months (\d+)-(\d+)", nf["windows"]["last6"]).expand(r"\1–\2")
    # ---- hypotheses and ranking ---------------------------------------------------------------
    hy = _load("hypotheses.json")
    if hy:
        for k_, verdict in hy["verdicts"].items():
            v[f"verdict_{k_}"] = verdict
        a = hy["H3a"]
        v["H3a_M6"], v["H3a_M7"] = pct(a["mean_abs_ire"]["M6"]), pct(a["mean_abs_ire"]["M7"])
        v["H3a_nbetter"] = str(a["n_runs_M7_better"])
        v["H3a_nruns"] = str(len(a["cases"]))
        v["H3a_better"] = ", ".join(_case(c) for c in a["runs_M7_better"])
        v["H3a_floor"] = pct(a["floor"], 1)
        b = hy["H3b"]
        v["H3b_M6"], v["H3b_M7"] = f(b["mean_rmse_map"]["M6"]), f(b["mean_rmse_map"]["M7"])
        v["H3b_ratio_obs"] = f(b["ratio_M7_over_M6"], 3)
        v["H3b_nfolds"] = str(b["n_folds"])
        c = hy["H3c"]
        for mm in ("M6", "M7"):
            x = c[mm]
            v[f"H3c_{mm}_ncases"] = str(x["n_cases"])
            v[f"H3c_{mm}_nviol"] = str(x["n_violations"])
            cf = x["ockenden_10pct_confusion"]
            for key, short in (("spread_gt_10pct_error_gt_20pct", "hh"),
                               ("spread_gt_10pct_error_le_20pct", "hl"),
                               ("spread_le_10pct_error_gt_20pct", "lh"),
                               ("spread_le_10pct_error_le_20pct", "ll")):
                v[f"ock_{mm}_{short}"] = str(cf[key])
        v["H3c_M7_viol_list"] = "; ".join(
            f"{_case(x['case'])} ({x['kind']}{'' if x['n1_pass'] else ', fails N1'}): spread "
            f"{pct(x['spread_rel'])} %, |integrated error| {pct(x['abs_ire'])} %"
            for x in c["M7"]["violations"])
        v["H3c_M6_lowhigh"] = ", ".join(_case(x) for x in c["M6"]["ockenden_10pct_confusion"]["cases_low_spread_high_error"])
        rk = hy["ranking"]["main"]
        for kind in KINDS:
            r = rk[kind]
            v[f"n_{kind}_runs"] = str(r["methods"]["M0"]["n_runs"])
            v[f"order_rmse_{kind}"] = " < ".join(r["order_by_mean_rmse_map"])
            v[f"order_ire_{kind}"] = " < ".join(r["order_by_mean_abs_ire"])
            v[f"order_gl_{kind}"] = " < ".join(r["order_by_rms_rmse_gl"])
            v[f"floor_best_{kind}"] = ", ".join(r["within_rmse_floor_of_best"])
            for mm in METHODS:
                x = r["methods"][mm]
                v[f"rmse_{kind}_{mm}"] = f(x["mean_rmse_map"])
                v[f"ire_{kind}_{mm}"] = pct(x["mean_abs_ire"], 1)
                v[f"gl_{kind}_{mm}"] = f(x["rms_rmse_gl"])
        # the wording of the abstract and of section 3.1 assumes this order
        best_rmse = {k_: rk[k_]["order_by_mean_rmse_map"][0] for k_ in KINDS}
        best_ire = {k_: rk[k_]["order_by_mean_abs_ire"][0] for k_ in KINDS}
        assert best_rmse == {"loro": "M4", "zc_out": "M2", "warm_out": "M4", "final": "M4"}, best_rmse
        assert best_ire == {"loro": "M6", "zc_out": "M2", "warm_out": "M6", "final": "M6"}, best_ire
        assert hy["verdicts"] == {"H3a": "FAIL", "H3b": "PASS", "H3c_M6": "PASS",
                                  "H3c_M7": "FAIL"}, hy["verdicts"]
        n1 = hy["ranking"]["n1"]["kinds"]
        n1_best = {k_: n1[k_]["order_by_mean_rmse_map"][0] for k_ in ("loro", "zc_out", "warm_out")}
        assert n1_best == {"loro": "M6", "zc_out": "M2", "warm_out": "M5"}, n1_best   # abstract wording
        for k_, short in (("loro", "loro"), ("warm_out", "warm"), ("zc_out", "zc")):
            v[f"n1_{short}_best"] = n1_best[k_]
            v[f"n1_{short}_best_rmse"] = f(n1[k_]["methods"][n1_best[k_]]["mean_rmse_map"])
            v[f"n1_{short}_M4_rank"] = str(n1[k_]["order_by_mean_rmse_map"].index("M4") + 1)
        # prose claims of section 3 (cheap self-checks: a change in the JSON fails the render)
        for k_ in KINDS:
            for key in ("mean_rmse_map", "mean_abs_ire"):
                assert rk[k_]["methods"]["M3"][key] > rk[k_]["methods"]["M2"][key], (k_, key)
            assert all(rk[k_][o][0] != "M7" for o in
                       ("order_by_mean_rmse_map", "order_by_mean_abs_ire", "order_by_rms_rmse_gl"))
        assert all(rk["zc_out"][o][0] == "M2" for o in
                   ("order_by_mean_rmse_map", "order_by_mean_abs_ire", "order_by_rms_rmse_gl"))
        assert rk["zc_out"]["order_by_mean_rmse_map"][-1] == "M4"
        assert {x["case"] for x in hy["H3c"]["M7"]["violations"]} <= {
            c for c in (x["case"] for x in hy["H3c"]["M7"]["cases"]) if "_zc550_" in c}
        v["n1_loro_M4"] = f(n1["loro"]["methods"]["M4"]["mean_rmse_map"])
        v["n1_zc_M4"] = f(n1["zc_out"]["methods"]["M4"]["mean_rmse_map"])
        v["n1_warm_M4"] = f(n1["warm_out"]["methods"]["M4"]["mean_rmse_map"])
        am = hy["amendments"][0]
        v["am_frozen_h"] = f"{am['frozen_wall_limit_hours']:g}"
        v["am_new_h"] = f"{am['wall_limit_hours']:g}"
    # ---- bootstrap intervals over runs (R17) ------------------------------------------------
    bs = _load("bootstrap.json")
    if bs:
        v["boot_B"], v["boot_seed"] = str(bs["n_boot"]), str(bs["seed"])
        v["boot_level"] = pct(bs["level"], 0)
        for variant, pre in (("main", ""), ("n1", "n1_")):
            for kind, d in bs["variants"][variant].items():
                for mm, x in d.items():
                    ci = lambda st, sc, nd: "–".join(  # noqa: E731
                        f"{sc * y:.{nd}f}" for y in x[st]["ci95"])
                    v[f"{pre}rmse_ci_{kind}_{mm}"] = ci("rmse_map", 1, 2)
                    v[f"{pre}ire_ci_{kind}_{mm}"] = ci("abs_integrated_rel_error", 100, 1)
                    v[f"{pre}gl_ci_{kind}_{mm}"] = ci("rmse_gl_rms", 1, 2)
                    if hy and variant == "main" and kind in KINDS:   # same runs as the ranking
                        r = hy["ranking"]["main"][kind]["methods"][mm]
                        assert abs(r["mean_rmse_map"] - x["rmse_map"]["mean"]) < 1e-9
                        assert abs(r["mean_abs_ire"] - x["abs_integrated_rel_error"]["mean"]) < 1e-9
                        assert abs(r["rms_rmse_gl"] - x["rmse_gl_rms"]["mean"]) < 1e-9
        mains = bs["variants"]["main"]
        for kind in ("loro", "warm_out", "final", "zc_out"):           # wording of 3.1 and abstract
            iv = [mains[kind][mm]["rmse_map"]["ci95"] for mm in ("M4", "M5", "M6", "M7")]
            common = max(x[0] for x in iv) <= min(x[1] for x in iv)
            assert common == (kind != "zc_out"), kind
        z = mains["zc_out"]
        assert z["M2"]["rmse_map"]["ci95"][1] < z["M4"]["rmse_map"]["ci95"][0]
        assert all(z["M4"]["rmse_map"]["ci95"][0] > z[mm]["rmse_map"]["ci95"][1] for mm in ("M5", "M6", "M7"))
    # ---- POST HOC paired differences (R18): which differences are resolved ---------------------
    if bs and "paired_differences_post_hoc" in bs:
        pdv = bs["paired_differences_post_hoc"]["variants"]
        short = {"rmse_map": "rmse", "abs_integrated_rel_error": "ire", "rmse_gl_rms": "gl"}
        for variant, pre in (("main", ""), ("n1", "n1_")):
            for kind, d in pdv[variant].items():
                for st, sh in short.items():
                    v[f"{pre}lead_{sh}_{kind}"] = ", ".join(sorted(d[st]["leaders"]))
        L = lambda variant, kind, st: set(pdv[variant][kind][st]["leaders"])  # noqa: E731
        P = lambda kind, pair, st="rmse_map": pdv["main"][kind][st]["pairs"][pair]  # noqa: E731
        # every sentence of the abstract, 3.1 and the conclusions that rests on these
        assert all(L("main", "zc_out", st) == {"M2"} for st in short)
        assert all(L("n1", "zc_out", st) == {"M2"} for st in short)
        assert P("zc_out", "M4-M0")["resolved"] and P("zc_out", "M4-M0")["mean_diff"] > 0
        assert {"M2", "M4", "M5", "M6", "M7"} <= L("main", "loro", "rmse_map")
        assert L("main", "loro", "abs_integrated_rel_error") == {"M6"}
        assert "M2" in L("n1", "loro", "abs_integrated_rel_error")
        assert {"M1", "M2", "M4", "M5", "M6", "M7"} <= L("main", "warm_out", "rmse_map")
        assert L("main", "warm_out", "abs_integrated_rel_error") == {"M4", "M5", "M6", "M7"}
        assert L("main", "final", "rmse_map") == {"M4", "M5", "M6"}
        assert all(P("final", f"M7-{m}")["resolved"] and P("final", f"M7-{m}")["mean_diff"] > 0
                   for m in ("M4", "M5", "M6"))
        assert L("main", "final", "abs_integrated_rel_error") == {"M4", "M5", "M6", "M7"}
        assert P("loro", "M7-M6")["resolved"] and P("loro", "M7-M6")["mean_diff"] > 0
        for variant in ("main", "n1"):          # a method without a network always leads the map
            assert all(L(variant, k, "rmse_map") - {"M6", "M7"} for k in pdv[variant])
        assert pdv["main"]["warm_out"]["rmse_map"]["best"] == "M4"
        oth = sorted(L("main", "warm_out", "rmse_map") - {"M4"})
        v["lead_rmse_warm_out_others"] = ", ".join(oth[:-1]) + " or " + oth[-1]
        lo = [bs["variants"]["main"]["loro"][m]["rmse_map"]["mean"] for m in L("main", "loro", "rmse_map")]
        v["lead_rmse_loro_range"] = f"{min(lo):.2f}–{max(lo):.2f}"
        wo = [bs["variants"]["main"]["warm_out"][m]["rmse_map"]["mean"] for m in L("main", "warm_out", "rmse_map")]
        v["lead_rmse_warm_range"] = f"{min(wo):.2f}–{max(wo):.2f}"
        x = P("loro", "M7-M6")
        v["loro_M7_M6_diff"] = f(x["mean_diff"])
        v["loro_M7_M6_ci"] = f"{x['ci95'][0]:.2f}–{x['ci95'][1]:.2f}"
        x = P("zc_out", "M7-M2")
        v["zc_M7_M2_ci"] = f"{x['ci95'][0]:.2f}–{x['ci95'][1]:.2f}"
        x = P("loro", "M6-M4")
        v["loro_M6_M4_ci"] = typographic_minus(f"{x['ci95'][0]:.2f} to {x['ci95'][1]:.2f}")
        for a_, b_ in (("M1", "M4"), ("M1", "M5"), ("M1", "M6"), ("M2", "M4"), ("M2", "M5"),
                       ("M2", "M6"), ("M3", "M4"), ("M3", "M5"), ("M3", "M6")):
            y = P("final", f"{b_}-{a_}")
            assert y["resolved"] and y["mean_diff"] < 0, (a_, b_)   # "each resolved better than M1-M3"
        assert L("main", "warm_out", "rmse_map") | {"M0", "M3"} == set(METHODS)
        def better(kind, a_, b_, st="rmse_map", variant="main"):     # a_ resolved better than b_
            pr = pdv[variant][kind][st]["pairs"]
            y = pr.get(f"{a_}-{b_}") or {**pr[f"{b_}-{a_}"], "mean_diff": -pr[f"{b_}-{a_}"]["mean_diff"]}
            return y["resolved"] and not y["below_floor"] and y["mean_diff"] < 0
        for nn in ("M6", "M7"):
            assert not any(better(k, nn, "M5") for k in pdv["main"])         # never better than M5
            for o in ("M1", "M2", "M3"):
                assert better("loro", nn, o) and better("final", nn, o)
                assert better("warm_out", "M6", o)
            assert better("zc_out", "M2", nn)
        among = [(a_, b_) for a_ in ("M4", "M5", "M6", "M7") for b_ in ("M4", "M5", "M6", "M7")
                 if a_ < b_ and (better("loro", a_, b_) or better("loro", b_, a_))]
        assert among == [("M6", "M7")], among                                  # only M7 worse than M6
        assert better("final", "M5", "M6") and not better("final", "M4", "M6")
        assert all(better("loro", "M6", o, "abs_integrated_rel_error") for o in METHODS if o != "M6")
        rows = []
        n1v, n1p = hy["ranking"]["n1"]["kinds"], pdv["n1"]
        for kind in ("loro", "zc_out", "warm_out"):
            rows.append([f"`{kind}`", n1v[kind]["order_by_mean_rmse_map"][0],
                         ", ".join(sorted(n1p[kind]["rmse_map"]["leaders"])),
                         n1v[kind]["order_by_mean_abs_ire"][0],
                         ", ".join(sorted(n1p[kind]["abs_integrated_rel_error"]["leaders"]))])
        v["table_n1"] = _md_table(["Held-out set", "lowest mean rmse_map", "leading (rmse_map)",
                                   "lowest mean \\|integrated error\\|", "leading (\\|integrated error\\|)"], rows)
    pb = _load("probabilistic.json")
    if pb:
        d = pb["variants"]["main"]
        kinds = ("loro", "zc_out", "warm_out", "final", "profile_test")
        rows = []
        for mm in ("M5", "M6", "M7"):
            rows.append([mm] + [f"{f(d[k][mm]['coverage_68'])} / {f(d[k][mm]['coverage_95'])} / "
                                f"{f(d[k][mm]['crps'])}" if mm in d[k] else "—" for k in kinds])
        v["table_prob"] = _md_table(["Method"] + [f"`{k}`" for k in kinds], rows)
        rows = []
        for mm in METHODS[1:]:
            rows.append([mm] + [f"{f(d[k][mm]['skill_rmse_map'])} / {f(d[k][mm]['skill_abs_integrated_rel_error'])}"
                                if mm in d[k] else "—" for k in kinds])
        v["table_skill"] = _md_table(["Method"] + [f"`{k}`" for k in kinds], rows)
        cov = sorted(int(k.split("_")[1]) for k in d["loro"]["M6"] if k.startswith("coverage_"))
        v["cov_lo"], v["cov_hi"] = str(cov[0]), str(cov[1])
        v["prof_M7_c95"] = f(d["profile_test"]["M7"]["coverage_95"])
        v["prof_M7_crps"] = f(d["profile_test"]["M7"]["crps"], 1)
        v["warm_M7_c95"] = f(d["warm_out"]["M7"]["coverage_95"])
        v["final_M7_c95"] = f(d["final"]["M7"]["coverage_95"])
        v["loro_M6_c95"], v["loro_M6_c68"] = f(d["loro"]["M6"]["coverage_95"]), f(d["loro"]["M6"]["coverage_68"])
        assert d["warm_out"]["M7"]["coverage_95"] < min(d["warm_out"][m]["coverage_95"] for m in ("M5", "M6"))
        assert d["final"]["M7"]["coverage_95"] < min(d["final"][m]["coverage_95"] for m in ("M5", "M6"))
    # ---- cross-validation record --------------------------------------------------------------
    cv = _load("cv_results.json")
    if cv:
        v["cv_units"] = str(len(cv["units"]))
        v["cv_entries"] = str(len(cv["entries"]))
        v["cv_wall_h"] = f(cv["session"]["wall_s"] / 3600, 1)
        v["cv_cpu_h"] = f"{cv['cpu_s_total'] / 3600:.0f}"
        v["cv_workers"] = str(cv["session"]["workers"])
        v["cv_resumed"] = str(cv["session"]["units_resumed"])
        m7 = [u for u in cv["units"] if u["method"] == "M7" and u["variant"] == "main"]
        v["M7_main_units"] = str(len(m7))
        v["M7_main_local"] = str(sum(u["info"]["tf_kind"] == "local" for u in m7))
        ns = [mb["n"] for u in m7 for mb in u["info"]["members"]]
        v["M7_n_range"] = f"{min(ns):.2f}–{max(ns):.2f}"
        v["M7_n_median"] = f(statistics.median(ns))
        nn = [u for u in cv["units"] if u["method"] in ("M6", "M7")]
        v["nn_hit_max"] = str(sum(u["info"]["n_hit_max_epochs"] for u in nn))
        v["nn_units"] = str(len(nn))
        v["nn_members_trained"] = str(sum(len(u["info"]["members"]) for u in nn))
        assert all(not mb["hit_max_epochs"] for u in nn for mb in u["info"]["members"])
        m4 = [u for u in cv["units"] if u["method"] == "M4" and u["variant"] == "main"]
        v["M4_main_units"] = str(len(m4))
        v["M4_main_quad"] = str(sum(u["info"]["degree"] == 2 for u in m4))
        m5 = [u["info"]["k"] for u in cv["units"] if u["method"] == "M5" and u["variant"] == "main"]
        v["M5_k_range"] = f"{min(m5)}–{max(m5)}"
        prof = {e["method"]: e for e in cv["entries"]
                if e["fold"] == "profile_test" and e["variant"] == "main"}
        for mm, e in prof.items():
            v[f"prof_{mm}_ire"] = pct(e["metrics"]["integrated_rel_error"], 1 if mm == "M2" else 0, sign=True)
            v[f"prof_{mm}_rmse"] = f(e["metrics"]["rmse_map"], 1)
        v["prof_M6_spread"] = pct(prof["M6"]["members"]["spread_rel"], 1)
        v["prof_M7_spread"] = pct(prof["M7"]["members"]["spread_rel"], 0)
        p1 = {e["method"]: e for e in cv["entries"]
              if e["fold"] == "profile_test" and e["variant"] == "n1"}
        v["prof_n1_M7_ire"] = pct(p1["M7"]["metrics"]["integrated_rel_error"], 0, sign=True)
    # ---- final test ---------------------------------------------------------------------------
    ft = _load("final_test.json")
    if ft:
        v["final_attempts"] = str(len(ft["attempts"]))
        v["final_status"] = ft["attempts"][-1]["status"]
        v["final_wall_min"] = f"{ft['wall_s'] / 60:.0f}"
        v["final_created"] = ft["created"][:10]
        v["final_maps_sha12"] = ft["maps_sha256"][:12]
        v["dataset_sha12"] = ft["dataset_sha256"][:12]
        v["config_sha12"] = ft["config_sha256"][:12]
        ent = {(e["case"], e["method"]): e for e in ft["entries"]}
        rows = []
        for mm in METHODS:
            row = [mm]
            for c in ft["test"]:
                e = ent[(c, mm)]
                sp = (e.get("members") or {}).get("spread_rel")
                row.append(f"{f(e['metrics']['rmse_map'], 1)} / {pct(e['metrics']['integrated_rel_error'], 1, sign=True)}"
                           + (f" ({pct(sp, 1)})" if sp is not None else ""))
            rows.append(row)
        hdr = ["Method"] + ["z_c = " + re.search(r"_zc(\d+)_", c)[1] + " m" for c in ft["test"]]
        v["table_final"] = _md_table(hdr, rows)
        for c in ft["test"]:
            z = re.search(r"_zc(\d+)_", c)[1]
            v[f"final_true_{z}"] = f(ent[(c, "M0")]["metrics"]["integrated_true"], 0)
    # ---- figures 2 and 3 ----------------------------------------------------------------------
    f2 = _load("fig2_maps_T2.json")
    if f2:
        v["f2_zc"] = f"{f2['z_c']:.0f}"
        v["f2_best"] = f2["best_physics"]
        v["f2_mean"] = f(f2["mitgcm_shelf_mean"], 1)
        for mm, x in f2["methods"].items():
            v[f"f2_{mm}_rmse"] = f(x["rmse_map"], 1)
            v[f"f2_{mm}_ire"] = pct(x["integrated_rel_error"], 0, sign=True)
            v[f"f2_{mm}_gl"] = f(x["rmse_gl"], 2)
    f3 = _load("fig3_physics.json")
    if f3:
        fb = f3["b_power_law"]["fits"]
        for z in ZC:
            for var, key in (("Tbot_minus_Tf_GL", "dT"), ("TF_shelf_mean", "TF")):
                v[f"pl_{key}_{z}"] = f(fb[var][z]["all_Tbot"]["n"])
                v[f"pl_{key}_{z}_le"] = f(fb[var][z]["Tbot_le_1.5"]["n"])
                v[f"pl_{key}_{z}_weak"] = "yes" if fb[var][z]["all_Tbot"]["poorly_constrained"] else "no"
        pooled = fb["TF_shelf_mean"]["pooled_all_runs"]
        v["pl_pooled_n"], v["pl_pooled_runs"] = f(pooled["n"]), str(pooled["n_runs"])
        v["pl_pooled_R2"] = f(pooled["R2_log"])
        tf550 = [r["TF_shelf_mean"] for r in f3["b_power_law"]["runs"].values()
                 if r["z_c"] == 550.0 and r["s"] == 1.0]
        v["TF550_range"] = f"{min(tf550):.2f}–{max(tf550):.2f}"
        cl = f3["c_lambert"]
        v["da_thr"] = f(cl["da_small_denominator_threshold"], 1)
        tl = [float(re.match(r"ens_T(-?[\d.]+)_", c)[1]) for c in cl["250"]["cases"]]
        v["lb_T"] = f"{tl[0]:.1f}, {tl[1]:.1f} and {tl[2]:.1f}"
        le = [k_ for k_ in fb["Tbot_minus_Tf_GL"]["250"] if k_.startswith("Tbot_le_")]
        v["exp2_Tmax"] = le[0].removeprefix("Tbot_le_")
        rows = []
        for mm in ("MITgcm",) + METHODS[1:]:
            row = [mm]
            for z in ZC:
                x = cl[z][mm]
                row.append(f"{f(x['dMdT_ref_p1'])} / {f(x['dMdT_p1_p2'])} ({f(x['d2MdT2'])})")
            for z in ZC:
                x = cl[z][mm]
                flag = lambda b: " ×" if b else ""  # noqa: E731
                row.append(f"{f(x['DA_ref_p1'], 1)}{flag(x['DA_ref_p1_small_denominator'])} / "
                           f"{f(x['DA_p1_p2'], 1)}{flag(x['DA_p1_p2_small_denominator'])}")
            rows.append(row)
        v["table_lambert"] = _md_table(
            ["Model"] + [f"dM/dT, z_c = {z} m" for z in ZC] + [f"DA, z_c = {z} m" for z in ZC], rows)
        x = cl["550"]["M4"]
        v["M4_550_DA"], v["M4_550_dMdT"] = f(x["DA_ref_p1"], 1), f(x["dMdT_ref_p1"])
        for z in ZC:
            v[f"da_mit_{z}"] = f(cl[z]["MITgcm"]["DA_ref_p1"], 1)
    # ---- post hoc profile-test diagnosis ------------------------------------------------------
    pd = _load("profile_test_diag.json")
    if pd:
        c = pd["consistency"]
        v["prof_d"], v["prof_zc"] = f"{c['dataset_params']['d']:.0f}", f"{c['dataset_params']['z_c']:.0f}"
        v["prof_Tbot"] = f"{c['dataset_params']['Tbot']:g}"
        v["prof_cells"] = str(c["n_shelf_cells"])
        jt = pd["joint_TF_vs_draft_minus_zc"]
        v["prof_cells_above"] = str(jt["ref_cells_TF_above_train_max_in_bin"])
        v["prof_frac_above"] = pct(jt["ref_frac_cells_TF_above_train_max_in_bin"], 0)
        mr = pd["marginal_ranges"]
        v["prof_marg_out"] = pct(max(x["ref_frac_outside"] for x in mr.values()), 0)
        v["prof_TF_ref"] = f"{mr['TF']['ref_min'] + 0.0:.2f}".replace("-0.00", "0.00") + f" to {mr['TF']['ref_max']:.2f}"
        v["prof_TF_train"] = f"{mr['TF']['train_min']:.2f} to {mr['TF']['train_max']:.2f}"
        dec = pd["M7_refit_decomposition"]
        v["prof_excess_share"] = pct(dec["excess_share_in_cells_TF_above_train_max_in_bin"], 0)
        v["prof_g1_ire"] = pct(dec["integrated_rel_error_P_only_g1_mean"], 0, sign=True)
        K = dec["members_K"]
        v["prof_K_range"] = f"{min(K):.0f}–{max(K):.0f}"
        v["prof_K0_M2"] = f(pd["M2_check"]["K_equiv_at_n2"], 1)
        rows = []
        jb = {tuple(b["bin_m"]): b for b in jt["bins"]}
        for b in dec["bins"]:
            j = jb[tuple(b["bin_m"])]
            rows.append([f"[{b['bin_m'][0]}, {b['bin_m'][1]})", str(b["n_cells"]),
                         f(b["median_TF"]), f(j["train_TF"]["max"]),
                         f"{b['sum_true']:.0f}", f"{b['sum_pred']:.0f}",
                         f(b["median_g"]), f(b["median_g_required"])])
        v["table_profile"] = _md_table(
            ["draft − z_c band (m)", "cells", "TF, Ocean0 (median, °C)", "TF, training (max, °C)",
             "Σ melt, MITgcm", "Σ melt, M7", "g_θ learned (median)", "g needed (median)"], rows)
        band = max(dec["bins"], key=lambda b: b["excess"])     # the band holding most excess
        v["band_label"] = f"{-band['bin_m'][1]}–{-band['bin_m'][0]}"
        v["band_g"], v["band_g_need"] = f(band["median_g"], 1), f(band["median_g_required"], 2)
        bins = sorted(dec["bins"], key=lambda b: b["bin_m"][0])
        last = max(i for i, b in enumerate(bins) if b["excess"] >= 0.01 * dec["excess_total"])
        deep = bins[last + 1:]                              # below the last band with excess
        v["deep_lo"] = str(-deep[0]["bin_m"][0])
        errs = [b["sum_pred"] / b["sum_true"] - 1 for b in deep]
        v["prof_deep_err"] = f"{pct(min(errs), 0, sign=True)} to {pct(max(errs), 0, sign=True)}"
    return v


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
    n = new.count(PENDING)
    print(f"wrote {OUT} ({len(new.split())} words, {n} pending values)")
