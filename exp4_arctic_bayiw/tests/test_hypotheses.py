"""Tests for exp4_arctic_bayiw.analysis.hypotheses (Task 7, Step 2): H4a-H4d evaluated
mechanically with the frozen PD0 thresholds (validate.THRESHOLDS, which test_validate checks
against the PD0 text) and the operational rules of the task-7 rulings, on synthetic JSONs: pass,
fail, tie at the threshold and inconclusive for every hypothesis."""
from __future__ import annotations

import json

import pytest

from exp4_arctic_bayiw.analysis import hypotheses as hy, validate as va

RES = ("kelvin_subinertial", "near_inertial", "superinertial", "seiche")


def _bands(low, k, ni, si, se):
    return {"low": low, "kelvin_subinertial": k, "near_inertial": ni, "superinertial": si, "seiche": se}


def test_thresholds_are_the_validated_pd0_ones():
    T = va.THRESHOLDS
    assert hy.T is T
    assert (T["H4a"]["phi_ni_support"], T["H4a"]["low_band_max_frac"]) == (0.5, 0.25)
    assert T["H4b"]["energy_ratio_max"] == 0.5 and T["H4c"]["ratio_min"] == T["H4d"]["ratio_min"] == 1.5


# ---- H4a ------------------------------------------------------------------------------------------

@pytest.mark.parametrize("b, outcome", [
    (_bands(0.0, 1.0, 6.0, 2.0, 1.0), "supported"),               # phi_NI = 0.6
    (_bands(0.0, 6.0, 1.0, 2.0, 1.0), "alternative"),             # phi_K = 0.6
    (_bands(0.0, 4.0, 4.0, 1.0, 1.0), "mixed"),                   # 0.4 / 0.4
    (_bands(0.0, 2.0, 5.0, 2.0, 1.0), "mixed"),                   # phi_NI = 0.5 exactly: not > 0.5
    (_bands(4.0, 1.0, 9.0, 1.0, 1.0), "inconclusive_T_fin"),      # low share 4/16 = 0.25: not > 0.25
])
def test_h4a_outcomes(b, outcome):
    r = hy.h4a(b)
    if b["low"] == 4.0:
        assert r["low_band_share"] == pytest.approx(0.25)
        assert r["outcome"] == "supported"                         # tie at 0.25 is not inconclusive
    else:
        assert r["outcome"] == outcome
    assert set(r["phi"]) == set(RES) and not r["resolved_net_inward"]


def test_h4a_low_band_over_a_quarter_is_inconclusive_whatever_the_phis():
    r = hy.h4a(_bands(-5.0, 1.0, 9.0, 1.0, 1.0))                  # |low| 5 / 17 > 0.25
    assert r["outcome"] == "inconclusive_T_fin" and r["outcome_pd0"] == "no concluyente con T_fin"
    assert r["phi"]["near_inertial"] == pytest.approx(0.75)          # every phi still reported


def test_h4a_flags_a_net_inward_resolved_flux_without_changing_the_rule():
    r = hy.h4a(_bands(0.0, -3.0, 1.0, 0.5, 0.5))                  # resolved sum -1
    assert r["resolved_net_inward"] is True
    assert r["phi"]["kelvin_subinertial"] == pytest.approx(3.0)      # -3 / -1 > 0.5
    assert r["outcome"] == "alternative"


def test_h4a_reports_the_f_plus_a_split_as_descriptive():
    """PD0 Enmienda 3: the pressure work decides; F_bc + A_bc is descriptive."""
    r = hy.h4a(_bands(0.0, 1.0, 6.0, 2.0, 1.0), f_plus_a=_bands(0.0, 5.0, 1.0, 1.0, 1.0))
    assert r["outcome"] == "supported"
    assert r["descriptive_F_plus_A"]["phi"]["kelvin_subinertial"] == pytest.approx(0.625)
    assert "outcome" not in r["descriptive_F_plus_A"]


# ---- H4b ------------------------------------------------------------------------------------------

@pytest.mark.parametrize("e2, e3, met", [(4.0, 10.0, True), (6.0, 10.0, False), (5.0, 10.0, True)])
def test_h4b_i_energy_ratio(e2, e3, met):
    r = hy.h4b(e_f2=e2, e_f3=e3, f_bc_in=1.0, f_bt_in=0.5)
    assert r["i"]["ratio"] == pytest.approx(e2 / e3) and r["i"]["met"] is met


@pytest.mark.parametrize("bc, bt, met", [(2.0, -1.0, True), (1.0, 1.0, True), (0.5, -1.0, False),
                                         (-2.0, 0.1, False)])
def test_h4b_ii_enters_through_the_sounds(bc, bt, met):
    r = hy.h4b(e_f2=1.0, e_f3=10.0, f_bc_in=bc, f_bt_in=bt, fa_bc_in=0.3, fa_bt_in=0.2)
    assert r["ii"]["met"] is met and r["ii"]["F_bc_in_J"] == bc
    assert r["ii"]["descriptive_F_plus_A"] == {"FA_bc_in_J": 0.3, "FA_bt_in_J": 0.2,
                                               "FA_bc_in_ge_abs_FA_bt_in": True}


def test_h4b_reports_the_budget_residuals():
    r = hy.h4b(e_f2=1.0, e_f3=10.0, f_bc_in=1.0, f_bt_in=0.5, residuals={"F2": 0.1, "F3": -0.2})
    assert r["budget_residual_D_bc_rel"] == {"F2": 0.1, "F3": -0.2}


# ---- H4c ------------------------------------------------------------------------------------------

@pytest.mark.parametrize("w1, w3, f1, f3, outcome", [
    (3.0, 1.0, 2.0, 1.0, "supported"),
    (1.5, 1.0, 1.5, 1.0, "supported"),                              # tie at 1.5: >= 1.5
    (1.0, 3.0, 1.0, 2.0, "refuted"),
    (1.0, 1.5, 1.0, 1.5, "refuted"),                                # tie at 1/1.5
    (3.0, 1.0, 1.0, 1.0, "not_distinguishable"),
    (3.0, 1.0, 1.0, 3.0, "not_distinguishable"),
])
def test_h4c_outcomes(w1, w3, f1, f3, outcome):
    r = hy.h4c(w_ni_f1=w1, w_ni_f3=w3, f_res_f1=f1, f_res_f3=f3)
    assert r["outcome"] == outcome
    assert "winter" in r["warning"] and not r["flags"]


def test_h4c_flags_non_positive_terms_and_applies_the_rule_literally():
    r = hy.h4c(w_ni_f1=3.0, w_ni_f3=1.0, f_res_f1=-1.0, f_res_f3=1.0)
    assert r["ii"]["ratio"] == pytest.approx(-1.0)
    assert "f_res_f1 <= 0" in r["flags"] and r["outcome"] == "not_distinguishable"


# ---- H4d ------------------------------------------------------------------------------------------

def _pair(num_k, pbt_k, pbc_k, num_c, pbt_c, pbc_c):
    return ({"F_bc_pressure_out_J": num_k, "P_bt_J": pbt_k, "P_bc_J": pbc_k,
             "totals_in_J": {"FA_bc": -2 * num_k}},
            {"F_bc_pressure_out_J": num_c, "P_bt_J": pbt_c, "P_bc_J": pbc_c,
             "totals_in_J": {"FA_bc": -2 * num_c}})


@pytest.mark.parametrize("vals, outcome, ratio", [
    ((3.0, 50.0, 50.0, 1.0, 50.0, 50.0), "supported", 3.0),
    ((1.5, 50.0, 50.0, 1.0, 50.0, 50.0), "supported", 1.5),        # tie
    ((1.0, 50.0, 50.0, 3.0, 50.0, 50.0), "refuted", 1 / 3),
    ((1.2, 50.0, 50.0, 1.0, 50.0, 50.0), "not_distinguishable", 1.2),
])
def test_h4d_pair(vals, outcome, ratio):
    kpp, clo = _pair(*vals)
    r = hy.h4d_pair(kpp, clo)
    assert r["outcome"] == outcome and r["ratio"] == pytest.approx(ratio)
    assert r["eta_kpp"] == pytest.approx(vals[0] / 100.0)
    assert r["eta_prime_kpp_descriptive"] == pytest.approx(vals[0] / 50.0)
    assert r["descriptive_FA_bc_out_J"]["kpp"] == pytest.approx(2 * vals[0])
    assert not r["flags"]


def test_h4d_flags_a_non_positive_eta():
    kpp, clo = _pair(-1.0, 50.0, 50.0, 1.0, 50.0, 50.0)
    r = hy.h4d_pair(kpp, clo)
    assert any("numerator" in f for f in r["flags"]) and r["ratio"] == pytest.approx(-1.0)
    assert r["outcome"] == "refuted"                                # literal rule, flagged


# ---- end to end on synthetic JSONs ------------------------------------------------------------------

def _synthetic_jsons():
    part = {"runs": {
        "F1": {"straits": {"total": {"F_bc_out_J": _bands(0.1, 6.0, 1.0, 2.0, 1.0),
                                     "FA_bc_out_J": _bands(0.1, 1.0, 1.0, 1.0, 1.0)}},
               "wind_work_J": {"bands": {"near_inertial": 4.0}}},
        "F3": {"straits": {"total": {"F_bc_out_J": _bands(0.1, 2.0, 0.5, 0.5, 0.5),
                                     "FA_bc_out_J": _bands(0.1, 1.0, 1.0, 1.0, 1.0)}},
               "wind_work_J": {"bands": {"near_inertial": 2.0}},
               "energy": {"totals": {"E_in_hp_J": 10.0, "ke_bc_hp_J": 5.0}}},
        "F2": {"energy": {"totals": {"E_in_hp_J": 2.0, "ke_bc_hp_J": 0.5}}}}}

    def b(num, pbt, pbc, fa_bc=1.0, fa_bt=0.2, D=0.0):
        return {"F_bc_pressure_out_J": num, "P_bt_J": pbt, "P_bc_J": pbc,
                "totals_in_J": {"F_bc": fa_bc, "F_bt": fa_bt, "FA_bc": 0.5, "FA_bt": 0.1},
                "D_bc_over_gross_sources": D}
    bud = {"runs": {"F1": b(2.0, 50, 50), "F1c": b(1.0, 50, 50), "F2": b(0, 0, 0, D=0.05),
                    "F3": b(0, 1, 1, D=-0.1), "R1": b(1.0, 10, 10), "R1k_0_120h": b(1.0, 10, 10)}}
    return part, bud


def test_evaluate_on_synthetic_jsons(tmp_path):
    part, bud = _synthetic_jsons()
    (tmp_path / "partition.json").write_text(json.dumps(part))
    (tmp_path / "budget.json").write_text(json.dumps(bud))
    out = hy.evaluate(tmp_path / "partition.json", tmp_path / "budget.json", tmp_path / "hypotheses.json")
    assert out["H4a"]["outcome"] == "alternative"                   # phi_K = 0.6
    assert out["H4b"]["i"]["met"] and out["H4b"]["i"]["ratio"] == pytest.approx(0.2)
    assert out["H4b"]["ii"]["met"]                                   # F_bc 1.0 >= |F_bt| 0.2
    assert "Enmienda 3" in out["general"]["flux_definition"]
    assert out["H4c"]["outcome"] == "supported"                      # 2.0 and 3.0
    assert out["H4d"]["F1_vs_F1c"]["outcome"] == "supported"
    assert out["H4d"]["R1k_vs_R1"]["outcome"] == "not_distinguishable"
    assert "V4" in out["general"]["factor_rule_note"]
    saved = json.loads((tmp_path / "hypotheses.json").read_text())
    assert saved["H4a"]["outcome"] == "alternative"


def test_strait_limitation_compares_widths_with_r1():
    bathy = {"sounds": {"Mariasundet": {"passages_npi_m": [500.0, 480.0], "grid": {"wet_length_m": 1300.0}},
                        "Akselsundet": {"width_npi_m": 1100.0, "grid": {"wet_length_m": 1200.0}}}}
    th = {"cases": {c: {"rossby_radius_m": [r]} for c, r in
                    (("vmf_summer", 3600.0), ("vmf_summer_Hmean", 2500.0),
                     ("vmf_winter", 1300.0), ("vmf_winter_Hmean", 340.0))}}
    out = hy.strait_limitation(bathy, th)
    m = out["Mariasundet"]
    assert m["width_m"] == 500.0 and m["width_over_R1"]["vmf_summer"] == pytest.approx(500 / 3600)
    assert out["Akselsundet"]["width_m"] == 1100.0
    assert "Johnson & Garrett (2006" in out["reference"] and out["label"].startswith("POST HOC")


def test_post_hoc_block_is_labelled_and_never_a_verdict():
    part, bud = _synthetic_jsons()
    part["runs"]["F1"]["straits"]["Akselsundet"] = part["runs"]["F1"]["straits"]["total"]
    part["runs"]["F1"]["arneborg_liljebladh_band_post_hoc"] = {"total": {"FA_bc_in_J": 3.0, "F_bc_in_J": 1.0}}
    ph = hy.post_hoc(part["runs"])
    assert ph["label"].startswith("POST HOC")
    assert ph["band_flux_per_strait"]["F1"]["Akselsundet"]["low_out_J"] == pytest.approx(0.1)
    assert ph["band_flux_per_strait"]["F1"]["Akselsundet"]["resolved_out_J"]["kelvin_subinertial"] == 6.0
    assert ph["band_flux_per_strait"]["F1"]["Akselsundet"]["resolved_F_plus_A_out_J"]["kelvin_subinertial"] == 1.0
    assert ph["arneborg_liljebladh_net_in_J"]["F1"]["total"]["FA_bc_in_J"] == 3.0
    assert "outcome" not in json.dumps(ph)


def test_h4b_post_hoc_ke_only_ratio_is_descriptive_and_does_not_change_the_verdict():
    r = hy.h4b(e_f2=6.0, e_f3=10.0, f_bc_in=1.0, f_bt_in=0.5, ke_only={"F2": 0.1, "F3": 5.0})
    assert r["i"]["met"] is False                                   # decided by E_in only
    ph = r["i"]["post_hoc_ke_bc_only"]
    assert ph["ratio"] == pytest.approx(0.02) and ph["label"].startswith("POST HOC")
    assert "met" not in ph


# ---- fix round 1 -------------------------------------------------------------------------------------

def test_h4d_flags_a_non_positive_denominator():
    kpp, clo = _pair(1.0, -60.0, 50.0, 1.0, 50.0, 50.0)        # P_bt + P_bc = -10 (kpp)
    r = hy.h4d_pair(kpp, clo)
    assert any(f.startswith("kpp: denominator") for f in r["flags"])


def test_h4d_flags_the_closure_side():
    kpp, clo = _pair(1.0, 50.0, 50.0, -2.0, 50.0, 50.0)
    r = hy.h4d_pair(kpp, clo)
    assert any(f.startswith("closure: numerator") for f in r["flags"])
    assert not any(f.startswith("kpp") for f in r["flags"])


@pytest.mark.parametrize("kn, cn, outcome", [(2.0, 3.0, "refuted"), (3.0, 2.0, "supported"),
                                             (0.2, 0.3, "refuted"), (0.3, 0.2, "supported")])
def test_h4d_tie_at_exactly_the_factor_counts(kn, cn, outcome):
    """eta ratios that equal 1.5 or 1/1.5 in exact arithmetic are ties (>= / <=) even when the
    floating-point quotient is one ulp off (0.02/0.03 = 0.6666666666666667 > 1/1.5)."""
    kpp, clo = _pair(kn, 50.0, 50.0, cn, 50.0, 50.0)
    assert hy.h4d_pair(kpp, clo)["outcome"] == outcome


def test_h4c_tie_is_robust_to_rounding():
    assert hy.h4c(w_ni_f1=0.2, w_ni_f3=0.3, f_res_f1=0.2, f_res_f3=0.3)["outcome"] == "refuted"


def test_post_hoc_resolved_band_eta_for_h4d():
    part = {"F1": {"straits": {"total": {"F_bc_out_J": _bands(-10.0, 1.0, 2.0, 0.5, 0.5)}}},
            "F1c": {"straits": {"total": {"F_bc_out_J": _bands(0.1, 1.0, 0.5, 0.25, 0.25)}}}}
    bud = {"F1": {"P_J": 100.0}, "F1c": {"P_J": 200.0}}
    r = hy.post_hoc_h4d_resolved(part, bud)
    assert r["label"].startswith("POST HOC") and "outcome" not in r
    assert r["F1"]["eta_resolved"] == pytest.approx(4.0 / 100.0)
    assert r["F1c"]["eta_resolved"] == pytest.approx(2.0 / 200.0)
    assert r["ratio_F1_over_F1c"] == pytest.approx(4.0)
    assert "low-frequency" in r["explanation"]


def test_evaluate_reports_windows_and_the_low_band_caveat(tmp_path):
    part, bud = _synthetic_jsons()
    for r in part["runs"].values():
        r["record"] = {"n": 840, "dt_s": 1200.0, "t_last_h": 279.6666666666667}
    part["runs"]["F1c"] = {"straits": {"total": {"F_bc_out_J": _bands(0.1, 1.0, 1.0, 1.0, 1.0),
                                                 "FA_bc_out_J": _bands(0.1, 1.0, 1.0, 1.0, 1.0)}},
                           "record": {"n": 840, "dt_s": 1200.0, "t_last_h": 279.6666666666667}}
    for b in bud["runs"].values():
        b["P_J"] = b["P_bt_J"] + b["P_bc_J"]
        b["window_h"] = [0.0, 279.6666666666667]
    bud["runs"]["R1k_0_120h"]["window_h"] = [0.0, 120.0]
    (tmp_path / "partition.json").write_text(json.dumps(part))
    (tmp_path / "budget.json").write_text(json.dumps(bud))
    out = hy.evaluate(tmp_path / "partition.json", tmp_path / "budget.json")
    assert out["H4a"]["window"] == {"t_h": [0.0, 279.6666666666667], "n_records": 840}
    assert out["H4d"]["R1k_vs_R1"]["window_h"] == [0.0, 120.0]
    assert "fixed-reference drift" in out["general"]["low_band_caveat"]
    assert "post_hoc_h4d_resolved_bands" in out


# ---- final-review fixes -------------------------------------------------------------------------

def test_h4a_literal_rule_can_support_with_an_inward_resolved_flux():
    """Pinned (final review, Minor 3): the PD0 rule is applied literally; when the resolved sum is
    inward the phis change sign and an inward flux can come out 'supported', with the flag."""
    r = hy.h4a(_bands(0.0, 1.0, -4.0, 0.5, 0.5))                  # resolved sum -2, NI carries it
    assert r["resolved_net_inward"] is True
    assert r["phi"]["near_inertial"] == pytest.approx(2.0)
    assert r["outcome"] == "supported"


def test_zero_denominators_are_flagged_not_raised():
    r = hy.h4b(1.0, 0.0, 1.0, 1.0)
    assert r["i"]["ratio"] is None and r["i"]["met"] is False and r["i"]["flags"]
    c = hy.h4c(1.0, 0.0, 1.0, 1.0)
    assert c["i"]["ratio"] is None and c["outcome"] == "undefined" and any("== 0" in f for f in c["flags"])
    k, cl = _pair(1.0, 0.0, 0.0, 1.0, 50.0, 50.0)
    d = hy.h4d_pair(k, cl)
    assert d["eta_kpp"] is None and d["ratio"] is None and d["outcome"] == "undefined"
    assert any("denominator" in f for f in d["flags"])
    k, cl = _pair(1.0, 50.0, 50.0, 0.0, 50.0, 50.0)                # eta_closure = 0
    d = hy.h4d_pair(k, cl)
    assert d["ratio"] is None and d["outcome"] == "undefined"


def test_h4d_reports_the_amendment_1_numerator_as_post_hoc():
    k, cl = _pair(3.0, 50.0, 50.0, 1.0, 50.0, 50.0)
    k["totals_in_J"]["FA_bc"] = -2.0                                 # F + A out: +2 (KPP)
    cl["totals_in_J"]["FA_bc"] = 4.0                                 # F + A out: -4 (closure)
    d = hy.h4d_pair(k, cl)
    p = d["post_hoc_amendment1_numerator"]
    assert "POST HOC" in p["label"] and "not a verdict" in p["label"]
    assert p["eta_kpp"] == pytest.approx(0.02) and p["eta_closure"] == pytest.approx(-0.04)
    assert p["ratio"] == pytest.approx(-0.5) and p["literal_rule_would_read"] == "refuted"
    assert d["outcome"] == "supported"                               # the verdict is untouched
