"""Render exp1_fjord_plumes/README.md from README.template.md and the figures/*.json metrics.

Every number in the README comes from a metrics file (plan Task 8, Step 1: no hand-typed numbers).
Placeholders are ``{{name}}``. An unknown placeholder raises. Values whose source file does not exist
yet (e.g. the 5 m (b) subset) render as "[pending]" only with ``--allow-pending``; otherwise they raise.

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


def _load(name: str):
    p = FIG / name
    return json.loads(p.read_text()) if p.exists() else None


def _rng(vals, fmt="{:.2f}") -> str:
    return f"{fmt.format(min(vals))}–{fmt.format(max(vals))}"


def values() -> dict[str, str]:
    v: dict[str, str] = {}
    # ---- (a) Figure 1 and the convergence study --------------------------------------------
    f1 = _load("fig1_metrics.json")
    if f1:
        m = {n: c["modes"]["mitgcm"] for n, c in f1["cases"].items()}
        ps = {n: c["pass"] for n, c in f1["cases"].items()}
        v["a_C2_npass"] = str(sum(p["C2"] for p in ps.values()))
        v["a_C2pp_npass"] = str(sum(p["C2pp"] for p in ps.values()))
        v["a_C2p_npass"] = str(sum(p["C2p"] for p in ps.values()))
        v["a_ncases"] = str(len(ps))
        fails = [n for n, p in ps.items() if not p["C2"]]
        v["a_C2_fail_cases"] = ", ".join(n.replace("a_line2d_", "") for n in fails) or "none"
        v["a_C2_fail_diff"] = ", ".join(f"{m[n]['z_detach_diff_over_H']:.3f}" for n in fails) or "–"
        v["a_C2_diff_range"] = _rng([x["z_detach_diff_over_H"] for x in m.values()], "{:.3f}")
        v["a_C2p_ratio_range"] = _rng([x["u_ratio"] for x in m.values()])
        v["a_C2pp_ratio_range"] = _rng([x["melt_ratio"] for x in m.values()])
    cv = _load("pd2_convergence.json")
    if cv:
        r = cv["runs"]
        for k, tag in (("conv_dx10", "10"), ("conv_dx5", "5"), ("conv_dx2p5", "2p5")):
            v[f"a_RQ_{tag}"] = f"{r[k]['R_Q']:.3f}"
        v["a_conv_class"] = cv["convergence"]["class"]
        v["a_RQ_ext"] = f"{cv['convergence']['richardson']['R_ext']:.2f}"
    # ---- (b) Figure 2, D1, D2, clean-background excess ---------------------------------------
    f2 = _load("fig2_metrics.json")
    if f2:
        for s in ("twolayer", "sermilik"):
            c = f2["C4"][s]
            v[f"b_C4_{s}_model"] = f"{c['z_detach_model']:.1f}"
            v[f"b_C4_{s}_theory"] = f"{c['z_detach_theory']:.1f}"
            v[f"b_C4_{s}_diff"] = f"{c['diff_over_H']:.3f}"
    ex = _load("fig2_excess_clean.json")
    if ex:
        for s in ("twolayer", "sermilik"):
            e = ex[s]
            v[f"b_bg_frac_{s}"] = f"{e['Mbg_over_Msingle']:.2f}"
            v[f"b_eps_ratio_{s}"] = f"{e['E_single'] / e['eps']:.0f}"
            re_ = {int(k): x["R_E"] for k, x in e["R_E"].items()}
            for x0, val in re_.items():
                v[f"b_RE_{s}_{x0}"] = f"{val:.2f}"
            v[f"b_RE_{s}_sep_range"] = _rng([val for x0, val in re_.items() if x0 > 0])
    d2 = _load("d2_resolution.json")
    if d2:
        for k in ("dx10", "dx5"):
            v[f"b_D2_{k}_model"] = f"{d2[k]['z_detach_model']:.1f}"
            v[f"b_D2_{k}_theory"] = f"{d2[k]['z_detach_theory']:.1f}"
            v[f"b_D2_{k}_diff"] = f"{d2[k]['diff_over_H']:.3f}"
        v["b_D2_delta"] = f"{d2['dx5']['diff_over_H'] - d2['dx10']['diff_over_H']:.3f}"
    b5 = _load("fig2_excess_dx5.json")
    if b5:
        for x0, x in b5["R_E"].items():
            v[f"b5_RE_{x0}"] = f"{x['R_E']:.2f}"
            v[f"b5_dRE_{x0}"] = f"{x['diff_5m_minus_10m']:+.2f}"
            v[f"b5_zdet_{x0}"] = f"{x['c4']['z_detach_model']:.0f}"
        v["b5_Esingle_change_pct"] = f"{100 * (b5['E_single'] / ex['twolayer']['E_single'] - 1):.0f}"
    # ---- (c) Figure 3, control, amplitude, pyOCEAN -------------------------------------------
    f3 = _load("fig3_metrics.json")
    if f3:
        v["c_W_over_Rd"] = f"{f3['W_over_Rd']:.2f}"
        v["c_ctrl_Qex"] = f"{f3['control']['Q_ex_mean_m3_s']:.0f}"
        for c, x in f3["cases"].items():
            k = c.replace("c_fjord_", "")
            v[f"c_{k}_Qex"] = f"{x['Q_ex_mean_m3_s']:,.0f}"
            v[f"c_{k}_Qexmax"] = f"{x['Q_ex_max_m3_s']:,.0f}"
            v[f"c_{k}_xctrl"] = f"{x['Q_ex_over_control']:.0f}"
            v[f"c_{k}_heat"] = f"{x['heat_toward_front_mean_W'] / 1e9:.1f}"
        v["c_C5p_all"] = "all four cases" if f3["C5p_all_pass"] else "NOT all cases"
        h = {c.replace("c_fjord_", ""): x["heat_toward_front_mean_W"] for c, x in f3["cases"].items()}
        v["c_heat_sum_parts"] = f"{(h['plume'] + h['wind'] + h['shelf']) / 1e9:.1f}"
        q = {c.replace("c_fjord_", ""): x["Q_ex_mean_m3_s"] for c, x in f3["cases"].items()}
        v["c_shelf_over_plume_Qex"] = f"{q['shelf'] / q['plume']:.0f}"
    am = _load("c_shelf_amplitude.json")
    if am:
        for c in ("shelf", "all"):
            for col, tag in (("x20km_mid", "20"), ("x5km_near_face", "5")):
                f = am["cases"][f"c_fjord_{c}"][col]
                v[f"c_amp_{c}_{tag}"] = f"{f['amplitude']:.1f}"
                v[f"c_ampr_{c}_{tag}"] = f"{f['amplitude_over_nominal']:.2f}"
                v[f"c_lag_{c}_{tag}"] = f"{f['lag_days']:.2f}"
    py = _load("crosscheck_pyo_c.json")
    if py:
        rs = {c["case"]: c for c in py["cases"]}
        for case in ("pyo_shelf", "pyo_shelf_runoff"):
            r = rs[case]
            v[f"py_{case}_amp"] = f"{r['c4']['amplitude']:.1f}"
            v[f"py_{case}_lag"] = f"{r['c4']['lag_days']:.2f}"
            v[f"py_{case}_c3"] = f"{r['c3']['max_ratio']:.2f}"
        rr = rs["pyo_shelf_runoff"]["c2"]
        v["py_runoff_Qout"] = f"{rr['Q_out_mean_d20_30']:.2f}"
        v["py_runoff_rel"] = f"{100 * rr['rel_error']:.1f}"
        v["py_runoff_drift"] = f"{rr['eta_drift_d20_30_m']:.4f}"
    return v


def render(allow_pending: bool) -> str:
    text = TEMPLATE.read_text()
    v = values()
    missing = sorted({k for k in re.findall(r"\{\{(\w+)\}\}", text) if k not in v})
    if missing and not allow_pending:
        raise KeyError(f"no value for placeholders: {missing}")
    return re.sub(r"\{\{(\w+)\}\}", lambda m: v.get(m.group(1), PENDING), text)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--allow-pending", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    new = render(a.allow_pending)
    if a.check:
        sys.exit(0 if OUT.exists() and OUT.read_text() == new else 1)
    OUT.write_text(new)
    n = new.count(PENDING)
    print(f"wrote {OUT} ({len(new.split())} words, {n} pending values)")
