"""Tests for analysis/h1_posthoc_tf.py (POST HOC, ruling R23): H1 refitted against the
forcing at the ice base instead of Tbot.

Synthetic runs (test_metrics.make_run, 8 x 5 x 4 grid): the four shelf columns get ice-base
depths 550, 450, 350 and 250 m (the shelf/deep masks only use icetopo and wetness), so
TF_eff of a ramp profile has an independent analytic value computed here.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import ensemble as en
from exp2_isomip_cavity.analysis import h1_posthoc_tf as hp
from exp2_isomip_cavity.analysis import metrics as mt
from test_metrics import NX, NY, default_geometry, make_run

DRAFTS = np.array([-550.0, -450.0, -350.0, -250.0])      # shelf columns 0-3
L1, L2, L3 = -0.0573, 0.0832, -7.53e-8                  # AD16 Table 4
RHO, G = 1027.51, 9.81


def draft_field():
    d = default_geometry()[0].copy()
    d[:, :4] = DRAFTS[None, :]
    return d


def tf_expected(Tbot, z_c, d=300.0, region="shelf"):
    """Independent copy: ramp T, S at each column's base depth minus the AD16 liquidus."""
    zd = -DRAFTS
    r = np.clip((zd - (z_c - d / 2)) / d, 0, 1)
    T = -1.9 + (Tbot + 1.9) * r
    S = 33.8 + (34.7 - 33.8) * r
    tf = T - (L1 * S + L2 + L3 * RHO * G * zd)
    sel = np.ones(4, bool) if region == "shelf" else DRAFTS < -300
    return float(tf[sel].mean())                          # equal column areas


def write_ramp_run(path, Tbot, z_c, m=10.0, gamma=1e-4, name="ens_x"):
    fw = -m * 1000.0 / 31536000.0
    return make_run(path, nrec=24, draft=draft_field(),
                    shi=lambda k: {"SHIfwFlx": fw, "SHIgammT": gamma},
                    case_extra={"params": {"case": name, "init_profile": "ramp", "Tbot": Tbot,
                                           "z_c": z_c, "s": 1.0, "d": 300.0},
                                "actual": {"Tbot": Tbot, "z_c": z_c, "s": 1.0, "d": 300.0,
                                           "init_profile": "ramp"}})


def test_tf_eff_a_priori_analytic(tmp_path):
    rd = mt.load_run(write_ramp_run(tmp_path / "r", 0.5, 400.0))
    assert hp.tf_eff_a_priori(rd, "shelf") == pytest.approx(tf_expected(0.5, 400.0), rel=1e-12)
    assert hp.tf_eff_a_priori(rd, "deep300") == pytest.approx(
        tf_expected(0.5, 400.0, region="deep300"), rel=1e-12)


def test_tf_eff_refuses_non_ramp_run(tmp_path):
    rd = mt.load_run(make_run(tmp_path / "o", nrec=24))   # Ocean0 profile
    with pytest.raises(ValueError, match="ramp"):
        hp.tf_eff_a_priori(rd, "shelf")


def test_diagnosed_td_is_jpl_compare_definition(tmp_path):
    """Tw - Tb = -<SHIfwFlx> L / (rho_ref cw <SHIgammT>) (AD16 Eq. 24, kappa_i = 0)."""
    run = write_ramp_run(tmp_path / "r", 0.5, 400.0, m=10.0, gamma=2e-4)
    fw = -10.0 * 1000.0 / 31536000.0
    expect = -fw * 3.34e5 / (RHO * 3974.0 * 2e-4)
    res = hp.td_diagnosed(run)
    assert res["shelf"] == pytest.approx(expect, rel=1e-12)
    assert res["deep300"] == pytest.approx(expect, rel=1e-12)


def test_power_fit_recovers_n_and_loo_and_excludes_nonpositive():
    x = np.array([0.3, 0.6, 1.1, 1.7, 2.4, -0.1])
    m = 4.2 * np.abs(x) ** 1.73
    r = hp.power_fit(x, m, labels=list("abcdef"))
    assert r["n"] == pytest.approx(1.73, abs=1e-12) and r["k"] == pytest.approx(4.2, rel=1e-12)
    assert r["R2_log"] == pytest.approx(1.0, abs=1e-12)
    assert r["excluded_nonpositive_forcing"] == ["f"] and r["n_runs"] == 5
    assert r["loo_n_range"] == pytest.approx([1.73, 1.73], abs=1e-12)
    # LOO is computed: one off point moves the range
    m2 = m.copy(); m2[2] *= 1.5
    r2 = hp.power_fit(x, m2, labels=list("abcdef"))
    assert r2["loo_n_range"][0] < r2["loo_n_range"][1]


def test_analyse_pooled_pure_power_law_in_tf_eff(tmp_path):
    """15 core runs with m = 3 TF_eff^1.7: every per-z_c and the pooled fit give n = 1.7."""
    specs = []
    for z in en.CORE_ZC:
        for T in en.CORE_TBOT:
            name = f"ens_T{T}_zc{z:.0f}_s1.0"
            m = 3.0 * tf_expected(T, z) ** 1.7
            write_ramp_run(tmp_path / f"{name}_r21", T, z, m=m, name=name)
            specs.append({"case": name, "group": "factorial", "Tbot": T, "z_c": z, "s": 1.0})
    tab = en.collect(tmp_path, cases=specs)
    res = hp.analyse(tab, tmp_path)
    fx = res["fits"]["TF_eff_shelf"]
    for key in ("250", "400", "550", "pooled"):
        assert fx[key]["n"] == pytest.approx(1.7, abs=1e-10), key
    assert fx["pooled"]["n_runs"] == 15 and fx["pooled"]["R2_log"] == pytest.approx(1.0)
    assert res["label"] == "POST HOC (R23)"
    assert set(res["fits"]) >= {"TF_eff_shelf", "TF_eff_deep300", "TD_diag_shelf",
                                "TD_diag_deep300", "Tbot_minus_Tf_GL"}
    assert res["H1_preregistered_verdict_unchanged"] == "FAIL"
    # figure + json
    out = hp.write_outputs(res, tmp_path / "figs")
    for ext in ("png", "pdf", "json"):
        assert (tmp_path / "figs" / f"fig2_posthoc_tf.{ext}").stat().st_size > 0
    js = json.loads((tmp_path / "figs" / "fig2_posthoc_tf.json").read_text())
    assert js["label"] == "POST HOC (R23)" and js["holland2008_guide_n"] == 1.61
    assert "POST HOC" in out["title"]


# ---- R24 part 2: quadratic local / local x non-local predictors (Favier 2019) -------------

def tf_columns(Tbot, z_c, d=300.0):
    zd = -DRAFTS
    r = np.clip((zd - (z_c - d / 2)) / d, 0, 1)
    return -1.9 + (Tbot + 1.9) * r - (L1 * (33.8 + 0.9 * r) + L2 + L3 * RHO * G * zd)


def test_quadratic_predictors_analytic_and_identity(tmp_path):
    rd = mt.load_run(write_ramp_run(tmp_path / "r", 0.5, 400.0))
    q = hp.quad_predictors(rd)
    tf = tf_columns(0.5, 400.0)
    assert q["P_nonlocal"] == pytest.approx((tf * tf.mean()).mean(), rel=1e-12)
    assert q["P_nonlocal"] == pytest.approx(tf.mean() ** 2, rel=1e-12)    # identity
    assert q["P_local"] == pytest.approx((tf * np.abs(tf)).mean(), rel=1e-12)


def test_r24_pure_linear_in_p_recovers_q_one(tmp_path):
    """m = 7 P_nonlocal exactly: q = 1 per z_c and pooled (R2 = 1)."""
    specs = []
    for z in en.CORE_ZC:
        for T in en.CORE_TBOT:
            name = f"ens_T{T}_zc{z:.0f}_s1.0"
            m = 7.0 * tf_columns(T, z).mean() ** 2
            write_ramp_run(tmp_path / f"{name}_r21", T, z, m=m, name=name)
            specs.append({"case": name, "group": "factorial", "Tbot": T, "z_c": z, "s": 1.0})
    res = hp.analyse_quad(en.collect(tmp_path, cases=specs), tmp_path)
    for key in ("250", "400", "550", "pooled"):
        f = res["fits"]["P_nonlocal"][key]
        assert f["n"] == pytest.approx(1.0, abs=1e-10) and f["k"] == pytest.approx(7.0, rel=1e-9)
    assert res["fits"]["P_nonlocal"]["pooled"]["R2_log"] == pytest.approx(1.0, abs=1e-12)
    assert res["label"] == "POST HOC (R24)"
    out = hp.write_quad_outputs(res, tmp_path / "figs")
    for ext in ("png", "pdf", "json"):
        assert (tmp_path / "figs" / f"fig2_posthoc_quad.{ext}").stat().st_size > 0
    assert "POST HOC (R24)" in out["title"]
