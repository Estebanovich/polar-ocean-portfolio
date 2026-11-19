import pytest

import numpy as np

from exp6_southern_carbon.analysis import carbonate, decomp, hypotheses

YEAR_S = 31536000.0
BASE = {"mld_inv": 150.0, "area_R": 1e12, "g200": 0.3, "ddic1000": 100.0, "dpo4dz200": 4e-3,
        "dpo4_1000": 1.0, "moc_sv": 5.0, "r_cp": 117.0, "e": 0.3}


def test_no_change_gives_zero():
    t = decomp.terms(BASE, dict(BASE))
    assert t["mld"] == t["moc"] == t["bio"] == t["sum_circ_bio"] == 0.0


def test_moc_units():
    br = dict(BASE, moc_sv=BASE["moc_sv"] + 1.0)
    t = decomp.terms(dict(BASE, ddic1000=100.0), br)
    assert t["moc"] == pytest.approx(1e6 * 100e-6 * 1035.0 * YEAR_S)


def test_mld_units_and_bio_sign():
    br = dict(BASE, mld_inv=BASE["mld_inv"] + 10.0)
    t = decomp.terms(BASE, br)
    assert t["mld"] == pytest.approx(10.0 * 0.3 * 1035.0 * 1e12 * 1e-6)
    # extra PO4 supply -> more export -> DIC removed: bio < 0 (Jarnikova et al. 2025, eq. 3 sign)
    assert t["bio"] == pytest.approx(-117.0 * 10.0 * 4e-3 * 1035.0 * 1e12 * 1e-6 * 0.3)
    assert t["sum_circ_bio"] == pytest.approx(t["mld"] + t["moc"] + t["bio"])


def test_h6a_signed():
    r = hypotheses.h6a({1900: -1.0, 1980: -0.5, 2020: 0.1, 2060: 0.5}, {1900: 0.0, 1980: 0.0, 2020: 0.05, 2060: 0.1})
    assert r["dS"] == pytest.approx(1.5) and r["drift_gc"] == pytest.approx(0.1) and r["pass"]
    assert r["zero_crossing_year"] == 2020
    assert r["s_ant"][2060] == pytest.approx(0.4) and r["s_nat"][1900] == 0.0
    weak = hypotheses.h6a({1900: 0.0, 2060: 0.15}, {1900: 0.0, 2060: 0.1})
    assert not weak["pass"]
    wrong = hypotheses.h6a({1900: 0.5, 2060: -0.5}, {1900: 0.0, 2060: 0.0})
    assert not wrong["pass"] and wrong["zero_crossing_year"] == 2060


def test_h6b_all_branches():
    meas = {1900: 1.0, 1980: 2.0, 2020: -1.0, 2060: 4.0}
    ok = hypotheses.h6b({1900: 1.2, 1980: 1.6, 2020: -0.8, 2060: 4.9}, meas)
    assert ok["pass"] and ok["n_within"] == 4
    bad = hypotheses.h6b({1900: 1.3, 1980: 2.0, 2020: -1.0, 2060: 4.0}, meas)
    assert not bad["pass"] and bad["n_within"] == 3 and not bad["within"][1900]


# ---- temperature term (spec Sec. 7: reported, not in the H6b sum) ----

def test_conductance_is_ratio_of_area_weighted_sums():
    f = np.array([[2.0, 4.0]])                     # mol m-2 yr-1, + into the ocean
    dp = np.array([[1e-5, 3e-5]])                  # pCO2atm - pCO2ocean (atm)
    area = np.array([[1.0, 1.0]])
    assert decomp.conductance(f, dp, area) == pytest.approx(6.0 / 4e-5)


def test_temperature_term_warming_outgasses_by_about_4_percent_per_degree():
    dic, alk, po4 = np.full((1, 3), 2.10), np.full((1, 3), 2.30), np.full((1, 3), 1.5e-3)
    area = np.full((1, 3), 1e10)
    c = 1e5                                        # mol m-2 yr-1 atm-1
    t_p = np.full((1, 3), 5.0)
    term = decomp.temperature_term(dic, alk, po4, t_p, t_p + 1.0, area, c)
    p0 = carbonate.pco2(2.10, 2.30, 5.0, po4=1.5e-3)[0]
    assert term < 0.0                              # warmer branch -> higher pCO2 -> less uptake
    assert -term / (c * area.sum() * p0) == pytest.approx(0.0423, rel=0.15)   # Takahashi et al. 1993
    assert decomp.temperature_term(dic, alk, po4, t_p, t_p, area, c) == 0.0
