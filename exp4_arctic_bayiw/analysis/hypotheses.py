"""Task 7, Step 2: H4a-H4d evaluated mechanically (figures/hypotheses.json).

Thresholds: ``validate.THRESHOLDS`` (the frozen PD0 section 5 numbers; tests/test_validate.py
parses them from the PD0 text). Operational details: task-7 rulings (R25), fixed before any H4
number was computed. Nothing here redefines a threshold. Inputs: figures/partition.json and
figures/budget.json (analysis/partition.py); no model output is read here.

Flux (PD0 Enmienda 3, 2026-11-13, superseding ruling R25): the baroclinic flux of H4a, H4b(ii)
and H4c(ii) is the Kang & Fringer baroclinic pressure work alone, F_bc = sum u_n' p' dz (as H4d,
Enmienda 2, and V4); H4b(ii)'s barotropic counterpart is F_bt = H u_n,bar p_bt. F + A is
descriptive everywhere.

* **H4a** (F1): phi_b = F_b / sum of the 4 resolved bands of the net outward KF baroclinic
  pressure work F_bc through both straits; low-band share |F_low| / sum_5 |F_b|. If the share > 0.25:
  "no concluyente con T_fin"; else phi_NI > 0.5 "apoyada", phi_K > 0.5 "alternativa", else
  "mixta". Flag ``resolved_net_inward`` if sum_resolved <= 0 (the rule does not change). The
  F_bc + A_bc split is descriptive.
* **H4b** (F2 vs F3): (i) E_in(F2)/E_in(F3) <= 0.5, E_in the time-mean baroclinic energy
  (KE_bc + exact APE) of the inner fjord, fields high-passed at 2 df; (ii) F_bc,in >= |F_bt,in|
  with the baroclinic and barotropic pressure work, net inward through both straits, 0-T_fin,
  all bands. Reported separately, with the budget residuals.
* **H4c** (F1 vs F3): (i) W_NI(F1)/W_NI(F3), (ii) resolved outward F_bc, F1/F3: both
  >= 1.5 "apoyada", both <= 1/1.5 "refutada", else "no distinguible"; with the PD0 winter warning.
* **H4d** (R1k/R1 over 0-120 h, F1/F1c over 0-T_fin): eta = F_bc_pressure_out / (P_bt + P_bc)
  (Enmienda 2); ratio = eta(KPP) / eta(thesis closure): >= 1.5 "apoyada", <= 1/1.5 "refutada",
  else "no distinguible"; eta' = F_bc_pressure_out / P_bc and F_bc + A_bc out are descriptive.

A non-positive term in a ratio is reported in ``flags``; the rule is applied to the number as
computed (rulings: "report it and flag it; do not reinterpret").

Usage: hypotheses.py [--partition figures/partition.json] [--budget figures/budget.json]
                     [--out figures/hypotheses.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp4_arctic_bayiw.analysis import validate as va  # noqa: E402

T = va.THRESHOLDS
RESOLVED = va.RESOLVED
PKG = Path(__file__).resolve().parents[1]
FIG = PKG / "figures"
PD0_LABEL = {"supported": "apoyada", "alternative": "alternativa", "mixed": "mixta",
             "inconclusive_T_fin": "no concluyente con T_fin", "refuted": "refutada",
             "not_distinguishable": "no distinguible"}
WINTER_WARNING = ("PD0 5 (H4c): the winter mode 1 is poorly resolved at 200 m (R1/dx = 1.7-6.6, "
                  "PD0 2.3 point 1; winter at 200 m is a declared limitation, PD1). A supported "
                  "H4c may be partly a resolution effect.")
LOW_BAND_CAVEAT = (
    "The KF pressure work uses p' from rho' relative to the FIXED initial rho_b (Enmienda 1). The "
    "background profile drifts during the runs (the exact APE grows even without wind, F2), so the "
    "low band of the strait pressure work shares the fixed-reference drift: a term int p'_drift(z) "
    "Q'(z) dz wherever there is a mean baroclinic exchange Q'. The low band (which makes H4a "
    "inconclusive and sets the sign of eta(F1)) cannot be attributed to physical exchange without "
    "the POST HOC moving-reference check (figures/posthoc_moving_reference.json).")
WINDOWS_NOTE = ("Windows: fjord runs 0-279.67 h (840 records every 1200 s; the state after the last "
                "step is not a diagnostic record); R1 and R1k_0_120h 0-120 h (481 records every 900 s); "
                "R2/R3 0-280 h (1121 records every 900 s).")
FACTOR_NOTE = ("PD0 5: the factor-1.5 rule (H4c, H4d) is kept as frozen. It rests on the V4a "
               "resolution bound (25 %); V4 FAILED (figures/validation.json: V4a and V4b), so that "
               "bound is not demonstrated. No threshold is changed.")


def _phi(b: dict) -> tuple[dict, float]:
    s = sum(b[k] for k in RESOLVED)
    return {k: (b[k] / s if s != 0 else None) for k in RESOLVED}, s


def h4a(bands_out: dict, f_plus_a: dict | None = None) -> dict:
    """H4a from the band split [J, + out] of the baroclinic pressure work F_bc through both straits
    (F1; Enmienda 3); ``f_plus_a``: the F_bc + A_bc split, descriptive."""
    t = T["H4a"]
    phi, s = _phi(bands_out)
    low = abs(bands_out["low"]) / sum(abs(bands_out[k]) for k in ("low",) + RESOLVED)
    if low > t["low_band_max_frac"]:
        out = "inconclusive_T_fin"
    elif phi["near_inertial"] is not None and phi["near_inertial"] > t["phi_ni_support"]:
        out = "supported"
    elif phi["kelvin_subinertial"] is not None and phi["kelvin_subinertial"] > t["phi_k_alternative"]:
        out = "alternative"
    else:
        out = "mixed"
    r = {"bands_out_J": dict(bands_out), "resolved_sum_out_J": s, "phi": phi,
         "low_band_share": low, "resolved_net_inward": bool(s <= 0), "outcome": out,
         "outcome_pd0": PD0_LABEL[out],
         "rule": "low share > 0.25 -> inconclusive; else phi_NI > 0.5 supported, phi_K > 0.5 "
                 "alternative, else mixed (PD0 5; rulings: share = |F_low| / sum of the 5 |F_b|)"}
    if f_plus_a is not None:
        p, ps = _phi(f_plus_a)
        r["descriptive_F_plus_A"] = {
            "bands_out_J": dict(f_plus_a), "phi": p, "resolved_sum_out_J": ps,
            "low_band_share": abs(f_plus_a["low"]) / sum(abs(f_plus_a[k]) for k in ("low",) + RESOLVED)}
    return r


def h4b(e_f2: float, e_f3: float, f_bc_in: float, f_bt_in: float, fa_bc_in: float | None = None,
        fa_bt_in: float | None = None, residuals: dict | None = None, ke_only: dict | None = None) -> dict:
    """H4b (i) energy ratio and (ii) the source through the straits (F2): baroclinic against
    barotropic pressure work (Enmienda 3); F + A descriptive. J."""
    ratio = e_f2 / e_f3
    r = {"i": {"E_in_F2_J": e_f2, "E_in_F3_J": e_f3, "ratio": ratio,
               "threshold": T["H4b"]["energy_ratio_max"],
               "met": bool(ratio <= T["H4b"]["energy_ratio_max"])},
         "ii": {"F_bc_in_J": f_bc_in, "F_bt_in_J": f_bt_in,
                "met": bool(f_bc_in >= abs(f_bt_in))}}
    if fa_bc_in is not None:
        r["ii"]["descriptive_F_plus_A"] = {
            "FA_bc_in_J": fa_bc_in, "FA_bt_in_J": fa_bt_in,
            "FA_bc_in_ge_abs_FA_bt_in": bool(fa_bc_in >= abs(fa_bt_in))}
    if residuals is not None:
        r["budget_residual_D_bc_rel"] = residuals
    if ke_only is not None:
        r["i"]["post_hoc_ke_bc_only"] = {
            "label": "POST HOC sensitivity (written after seeing the exact-APE drift of the winter "
                     "runs): ratio of the high-passed KE_bc alone; not the criterion",
            "KE_bc_hp_F2_J": ke_only["F2"], "KE_bc_hp_F3_J": ke_only["F3"],
            "ratio": ke_only["F2"] / ke_only["F3"]}
    return r


def _tie(x: float, thr: float) -> bool:
    """Exact-arithmetic ties (x = thr) count as on the threshold even when the floating-point
    quotient is one ulp off (0.02/0.03 > 1/1.5); relative tolerance 1e-12, rounding only."""
    import math
    return math.isclose(x, thr, rel_tol=1e-12, abs_tol=0.0)


def _ratio_rule(rs: list, thr: float) -> str:
    if all(x >= thr or _tie(x, thr) for x in rs):
        return "supported"
    if all(x <= 1.0 / thr or _tie(x, 1.0 / thr) for x in rs):
        return "refuted"
    return "not_distinguishable"


def h4c(w_ni_f1: float, w_ni_f3: float, f_res_f1: float, f_res_f3: float,
        fa_res: dict | None = None) -> dict:
    """H4c from the near-inertial wind work [J] and the resolved outward baroclinic pressure work
    F_bc [J] (Enmienda 3); ``fa_res`` {F1, F3}: the resolved F_bc + A_bc, descriptive."""
    thr = T["H4c"]["ratio_min"]
    flags = [f"{k} <= 0" for k, v in (("w_ni_f1", w_ni_f1), ("w_ni_f3", w_ni_f3),
                                       ("f_res_f1", f_res_f1), ("f_res_f3", f_res_f3)) if v <= 0]
    ri, rii = w_ni_f1 / w_ni_f3, f_res_f1 / f_res_f3
    out = _ratio_rule([ri, rii], thr)
    r = {"i": {"W_NI_F1_J": w_ni_f1, "W_NI_F3_J": w_ni_f3, "ratio": ri},
         "ii": {"F_res_out_F1_J": f_res_f1, "F_res_out_F3_J": f_res_f3, "ratio": rii},
         "threshold": thr, "outcome": out, "outcome_pd0": PD0_LABEL[out], "flags": flags,
         "warning": WINTER_WARNING}
    if fa_res is not None:
        r["ii"]["descriptive_F_plus_A"] = {**fa_res, "ratio": fa_res["F1"] / fa_res["F3"]}
    return r


def h4d_pair(kpp: dict, closure: dict) -> dict:
    """One H4d pair from budget.json entries of the KPP run and the thesis-closure run."""
    thr = T["H4d"]["ratio_min"]
    flags, eta = [], {}
    for tag, b in (("kpp", kpp), ("closure", closure)):
        num, den = b["F_bc_pressure_out_J"], b["P_bt_J"] + b["P_bc_J"]
        if num <= 0:
            flags.append(f"{tag}: numerator F_bc_pressure_out <= 0 ({num:.4g} J)")
        if den <= 0:
            flags.append(f"{tag}: denominator P_bt + P_bc <= 0 ({den:.4g} J)")
        eta[tag] = num / den
    ratio = eta["kpp"] / eta["closure"]
    out = _ratio_rule([ratio], thr)
    return {"eta_kpp": eta["kpp"], "eta_closure": eta["closure"], "ratio": ratio, "threshold": thr,
            "outcome": out, "outcome_pd0": PD0_LABEL[out], "flags": flags,
            "numerator_J": {"kpp": kpp["F_bc_pressure_out_J"], "closure": closure["F_bc_pressure_out_J"]},
            "W_wind_J": {"kpp": kpp["P_bt_J"] + kpp["P_bc_J"],
                         "closure": closure["P_bt_J"] + closure["P_bc_J"]},
            "eta_prime_kpp_descriptive": kpp["F_bc_pressure_out_J"] / kpp["P_bc_J"],
            "eta_prime_closure_descriptive": closure["F_bc_pressure_out_J"] / closure["P_bc_J"],
            "descriptive_FA_bc_out_J": {"kpp": -kpp["totals_in_J"]["FA_bc"],
                                        "closure": -closure["totals_in_J"]["FA_bc"]}}


# ---- POST HOC descriptive block (task-7 rulings; never a verdict) ----------------------------------

POST_HOC = "POST HOC descriptive (task-7 rulings, R24): not a verdict, no threshold"
JG06 = ("Johnson & Garrett (2006, JPO 36:945, doi:10.1175/JPO2896.1), as cited by Jackson et al. "
        "(2018, JPO; literatura/Jackson_2018 ... .md, l. 704-709): 'unless a strait's width is "
        "considerably narrower than the deformation radius, the majority of incident Kelvin wave "
        "energy will be transmitted through the strait'. J&G 2006 itself is not in the local "
        "libraries [verificar]")


def strait_limitation(bathy: dict, th: dict) -> dict:
    """Strait width against R1 (theory.json): Mariasundet = the wider of its two NPI passages
    (Task 1b), Akselsundet = its NPI width; the model grid width is reported beside."""
    cases = ("vmf_summer", "vmf_summer_Hmean", "vmf_winter", "vmf_winter_Hmean")
    out = {"label": POST_HOC, "reference": JG06}
    for name, sd in bathy["sounds"].items():
        w = float(max(sd["passages_npi_m"])) if "passages_npi_m" in sd else float(sd["width_npi_m"])
        out[name] = {"width_m": w, "grid_wet_length_m": float(sd["grid"]["wet_length_m"]),
                     "R1_m": {c: float(th["cases"][c]["rossby_radius_m"][0]) for c in cases},
                     "width_over_R1": {c: w / float(th["cases"][c]["rossby_radius_m"][0]) for c in cases}}
    return out


def post_hoc(P: dict, bathy: dict | None = None, th: dict | None = None) -> dict:
    """R24 table (band flux per strait, resolved bands and the low band apart, + out), the
    Arneborg & Liljebladh (2001) 0.01-0.07 cph net inward flux, the bridge summary, the trapped
    fraction and (with bathy, th) the strait limitation. Read from partition.json."""
    out = {"label": POST_HOC, "band_flux_per_strait": {}, "arneborg_liljebladh_net_in_J": {},
           "bridge": {}, "trapped_fraction": {}}
    for run, r in P.items():
        st = r.get("straits")
        if st:
            out["band_flux_per_strait"][run] = {
                g: {"resolved_out_J": {k: v["F_bc_out_J"][k] for k in RESOLVED},
                    "low_out_J": v["F_bc_out_J"]["low"],
                    "resolved_F_plus_A_out_J": {k: v["FA_bc_out_J"][k] for k in RESOLVED},
                    "low_F_plus_A_out_J": v["FA_bc_out_J"]["low"]} for g, v in st.items()}
            tot = st["total"]["F_bc_out_J"]
            sab = sum(abs(tot[k]) for k in ("low",) + RESOLVED)
            out["bridge"][run] = {"flux_share_signed": {k: tot[k] / sab for k in ("low",) + RESOLVED}}
        if r.get("arneborg_liljebladh_band_post_hoc"):
            out["arneborg_liljebladh_net_in_J"][run] = {g: v for g, v in r["arneborg_liljebladh_band_post_hoc"].items()
                                                        if isinstance(v, dict)}
        if r.get("stations"):
            out["bridge"].setdefault(run, {})["station_fraction"] = {
                s: v["fraction"] for s, v in r["stations"].items()
                if s in ("shelf", "coast_south_mean", "Seiche", "Seiche2", "Omega")}
        if r.get("energy", {}).get("trapped_fraction_ke_bc"):
            out["trapped_fraction"][run] = r["energy"]["trapped_fraction_ke_bc"]
    if bathy is not None and th is not None:
        out["strait_limitation"] = strait_limitation(bathy, th)
    return out


def post_hoc_h4d_resolved(P: dict, B: dict) -> dict:
    """PI decision (Esteban, 2026-11-13): eta with the resolved bands of the outward KF pressure
    work only, eta_res = sum_4 resolved F_bc_out / (P_bt + P_bc), for F1 and F1c and their ratio.
    POST HOC, never a verdict."""
    out = {"label": "POST HOC descriptive (PI decision 2026-11-13): eta from the 4 resolved bands "
                    "only; never a verdict, the frozen H4d outcome stands",
           "explanation": ("The literal H4d F1/F1c 'refutada' comes from the sign of the low-frequency "
                           "net exchange term (the unresolved > 140 h band of the strait pressure work, "
                           "which also shares the fixed-reference drift), not from KPP radiating less."),
           "runs": ["F1", "F1c"]}
    for r in ("F1", "F1c"):
        b = P[r]["straits"]["total"]["F_bc_out_J"]
        res = sum(b[k] for k in RESOLVED)
        out[r] = {"resolved_out_J": res, "low_out_J": b["low"], "P_J": B[r]["P_J"],
                  "eta_resolved": res / B[r]["P_J"]}
    out["ratio_F1_over_F1c"] = out["F1"]["eta_resolved"] / out["F1c"]["eta_resolved"]
    return out


def evaluate(partition: Path | None = None, budget: Path | None = None, out: Path | None = None) -> dict:
    """H4a-H4d from partition.json and budget.json; writes hypotheses.json (``out``)."""
    P = json.loads(Path(partition or FIG / "partition.json").read_text())["runs"]
    B = json.loads(Path(budget or FIG / "budget.json").read_text())["runs"]
    tot = lambda run, key: P[run]["straits"]["total"][key]                   # noqa: E731
    res = {}
    res["H4a"] = h4a(tot("F1", "F_bc_out_J"), tot("F1", "FA_bc_out_J"))
    res["H4a"]["run"] = "F1"
    rec = P["F1"].get("record", {})
    res["H4a"]["window"] = {"t_h": [0.0, rec.get("t_last_h")], "n_records": rec.get("n")}
    b2 = B["F2"]["totals_in_J"]
    res["H4b"] = h4b(P["F2"]["energy"]["totals"]["E_in_hp_J"], P["F3"]["energy"]["totals"]["E_in_hp_J"],
                     b2["F_bc"], b2["F_bt"], b2["FA_bc"], b2["FA_bt"],
                     residuals={r: B[r]["D_bc_over_gross_sources"] for r in ("F2", "F3")},
                     ke_only={r: P[r]["energy"]["totals"]["ke_bc_hp_J"] for r in ("F2", "F3")})
    fres = {r: sum(tot(r, "F_bc_out_J")[k] for k in RESOLVED) for r in ("F1", "F3")}
    fares = {r: sum(tot(r, "FA_bc_out_J")[k] for k in RESOLVED) for r in ("F1", "F3")}
    res["H4c"] = h4c(P["F1"]["wind_work_J"]["bands"]["near_inertial"],
                     P["F3"]["wind_work_J"]["bands"]["near_inertial"], fres["F1"], fres["F3"], fares)
    res["H4d"] = {"R1k_vs_R1": {**h4d_pair(B["R1k_0_120h"], B["R1"]),
                                "window_h": B["R1k_0_120h"].get("window_h")},
                  "F1_vs_F1c": {**h4d_pair(B["F1"], B["F1c"]), "window_h": B["F1"].get("window_h")}}
    for h in ("H4b", "H4c"):
        res[h]["window"] = res["H4a"]["window"]
    if all("P_J" in B[r] for r in ("F1", "F1c")) and "F1c" in P:
        res["post_hoc_h4d_resolved_bands"] = post_hoc_h4d_resolved(P, B)
    res["general"] = {"factor_rule_note": FACTOR_NOTE, "low_band_caveat": LOW_BAND_CAVEAT,
                      "windows": WINDOWS_NOTE,
                      "flux_definition": ("PD0 Enmienda 3 (2026-11-13): H4a, H4b(ii), H4c(ii) and H4d "
                                          "are decided by the KF pressure work (F_bc; F_bt in H4b(ii)); "
                                          "F + A is descriptive"),
                      "thresholds": {k: T[k] for k in ("factor_min", "H4a", "H4b", "H4c", "H4d")},
                      "sources": "PD0 5, Enmiendas 1-2; task-7 rulings (R25); partition.json, budget.json"}
    bathy = PKG / "figures" / "bathy_report.json"
    res["post_hoc"] = post_hoc(P, json.loads(bathy.read_text()) if bathy.exists() else None, va.load_theory())
    if out is not None:
        Path(out).write_text(json.dumps(res, indent=1, default=va.jsonable) + "\n")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--partition", type=Path, default=FIG / "partition.json")
    ap.add_argument("--budget", type=Path, default=FIG / "budget.json")
    ap.add_argument("--out", type=Path, default=FIG / "hypotheses.json")
    a = ap.parse_args(argv)
    r = evaluate(a.partition, a.budget, a.out)
    for h in ("H4a", "H4c"):
        print(h, r[h]["outcome"])
    print("H4b (i)", r["H4b"]["i"]["met"], "(ii)", r["H4b"]["ii"]["met"])
    for p, v in r["H4d"].items():
        print("H4d", p, v["outcome"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
