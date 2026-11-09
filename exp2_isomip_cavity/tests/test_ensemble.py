"""Tests for analysis/ensemble.py: ensemble table, H1 (pd0) and H2 (pd0) verdicts.

Synthetic runs come from test_metrics.make_run (8 x 5 x 4 grid, shelf x < 640 km, draft
-40 m, x_front = column 4). Each run carries a uniform melt rate m (so mean_melt = m
exactly) and a tracer transport at x_front split between levels 1 (z = -30 m) and 2
(z = -50 m), so z_out = -30 w - 50 (1 - w) exactly. Expected values are analytic.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from exp2_isomip_cavity.analysis import ensemble as en
from exp2_isomip_cavity.theory import melt as ml
from test_metrics import NX, NY, NZ, make_run

D_GL = 572.1674970323206          # make_run's deepest_wet_ice_base_depth_m
TF_GL = ml.LAMBDA1 * 34.7 + ml.LAMBDA2 + ml.LAMBDA3 * 1027.51 * 9.81 * D_GL   # AD16 Eq. 25
CORE_T = [-0.5, 0.0, 0.5, 1.0, 1.5]
ZC = [250.0, 400.0, 550.0]


def _fw_for_melt(m):
    """Uniform SHIfwFlx [kg m^-2 s^-1] giving melt m [m/a] (m = -fw / 1000 * 31 536 000)."""
    return -m * 1000.0 / 31536000.0


def _ut(w):
    """UTRAC01 (u C hFacW) at every face: w at level 1, 1 - w at level 2 (both > 0)."""
    u = np.zeros((NZ, NY, NX))
    u[1], u[2] = w, 1.0 - w
    return u


def write_case(root, name, Tbot, z_c, s, m, w=0.5, group="factorial", status="OK"):
    fw = _fw_for_melt(m)
    make_run(root / f"{name}_r21", nrec=24, status=status,
             shi=lambda k: {"SHIfwFlx": fw}, dyn=lambda k: {"UTRAC01": _ut(w)},
             case_extra={"params": {"case": name, "Tbot": Tbot, "z_c": z_c, "s": s,
                                    "ensemble_group": group},
                         "actual": {"Tbot": Tbot, "z_c": z_c, "s": s}})
    return {"case": name, "group": group, "Tbot": Tbot, "z_c": z_c, "s": s}


def zout(w):
    return -30.0 * w - 50.0 * (1.0 - w)


def build_core(root, melt_fn, w_fn=lambda T, z: 0.5, skip=(), failed=()):
    specs = []
    for z in ZC:
        for T in CORE_T:
            name = f"ens_T{T}_zc{z:.0f}_s1.0"
            if (T, z) in skip:
                specs.append({"case": name, "group": "factorial", "Tbot": T, "z_c": z,
                              "s": 1.0})
                continue
            specs.append(write_case(root, name, T, z, 1.0, melt_fn(T, z), w_fn(T, z),
                                    status="FAILED" if (T, z) in failed else "OK"))
    return specs


# ---- collect ------------------------------------------------------------------------------

def test_collect_ok_run_values_and_failed_run_is_a_gap(tmp_path):
    specs = [write_case(tmp_path, "ens_a", 0.5, 400.0, 1.0, m=12.0, w=0.25),
             write_case(tmp_path, "ens_b", 1.0, 400.0, 1.0, m=20.0, status="FAILED")]
    tab = en.collect(tmp_path, cases=specs)
    assert tab["case"] == ["ens_a", "ens_b"]
    assert set(tab) >= {"case", "Tbot", "z_c", "s", "mean_melt_shelf", "mean_melt_deep300",
                        "total_flux", "N1", "N2", "z_out", "D_out", "split", "status"}
    assert len({len(v) for v in tab.values()}) == 1          # dict of equal-length lists
    i = 0
    assert tab["status"][i] == "OK" and tab["N3_pass"][i] is True
    assert tab["mean_melt_shelf"][i] == pytest.approx(12.0, rel=1e-12)
    assert tab["mean_melt_deep300"][i] is None               # no base < -300 m in the toy
    assert tab["z_out"][i] == pytest.approx(zout(0.25), rel=1e-12)
    assert tab["D_out"][i] == pytest.approx(-zout(0.25), rel=1e-12)
    assert tab["N1"][i] == pytest.approx(0.0, abs=1e-12) and tab["N1_pass"][i] is True
    assert tab["T_f_GL"][i] == pytest.approx(TF_GL, rel=1e-12)
    # total flux [kg/s] = m rho_fw / yr * shelf area (4 columns x 5 rows of 10 x 20 km)
    assert tab["total_flux"][i] == pytest.approx(12.0 * 1000 / 31536000 * 20 * 10e3 * 20e3,
                                                 rel=1e-12)
    j = 1                                                    # the FAILED run: a gap
    assert tab["status"][j].startswith("GAP_N3") and tab["N3_pass"][j] is False
    for k in ("mean_melt_shelf", "D_out", "z_out", "N1", "N2", "total_flux"):
        assert tab[k][j] is None                             # never imputed
    assert tab["Tbot"][j] == 1.0 and tab["z_c"][j] == 400.0  # but still listed


def test_collect_missing_run_is_a_gap_and_split(tmp_path):
    specs = [{"case": "ens_missing", "group": "factorial", "Tbot": 0.0, "z_c": 250.0,
              "s": 1.0},
             write_case(tmp_path, "ens_T2.0", 2.0, 250.0, 1.0, m=40.0, group="extrapolation")]
    tab = en.collect(tmp_path, cases=specs)
    assert tab["status"][0] == "GAP_MISSING" and tab["mean_melt_shelf"][0] is None
    assert tab["split"] == ["train", "test"]                 # Tbot = 2.0 is Exp. 3 test data


def test_collect_refuses_case_json_that_disagrees_with_registry(tmp_path):
    spec = write_case(tmp_path, "ens_a", 0.5, 400.0, 1.0, m=12.0)
    spec = {**spec, "Tbot": 0.75}
    tab = en.collect(tmp_path, cases=[spec])
    assert tab["status"][0].startswith("GAP_MISMATCH")
    assert tab["mean_melt_shelf"][0] is None


def test_registry_cases_have_the_28_members_and_the_reference():
    specs = en.registry_cases()
    assert specs[0]["case"] == "ocean0_4km" and specs[0]["group"] == "reference"
    groups = [s["group"] for s in specs[1:]]
    assert groups.count("factorial") == 15 and groups.count("geometry") == 10
    assert groups.count("extrapolation") == 3 and len(specs) == 29
    assert all(s["Tbot"] == 2.0 for s in specs if s["group"] == "extrapolation")


# ---- fits ---------------------------------------------------------------------------------

def test_power_law_recovers_n_and_k_exactly():
    T = np.array(CORE_T)
    m = 3.7 * (T - TF_GL) ** 1.65
    r = en.power_law_fit(T, m, TF_GL)
    assert r["n"] == pytest.approx(1.65, abs=1e-12)
    assert r["k"] == pytest.approx(3.7, rel=1e-12)


def test_quadratic_fit_exact_and_r2():
    T = np.array(CORE_T)
    m = 2.0 * T**2 - 1.0 * T + 5.0
    r = en.quadratic_fit(T, m)
    assert (r["a"], r["b"], r["c"]) == pytest.approx((2.0, -1.0, 5.0), abs=1e-10)
    assert r["R2"] == pytest.approx(1.0, abs=1e-12)
    # analytic R^2 for a known residual: y = T^2 + e, e orthogonal to {1, T, T^2}
    e = np.array([1.0, -4.0, 6.0, -4.0, 1.0])               # 4th difference: orthogonal
    y = T**2 + 0.1 * e
    r = en.quadratic_fit(T, y)
    ss_res = (0.1 * e) ** 2
    assert r["R2"] == pytest.approx(1 - ss_res.sum() / ((y - y.mean()) ** 2).sum(), rel=1e-10)


def test_power_law_undefined_when_melt_not_positive():
    T = np.array(CORE_T)
    m = np.array([-1.0, 1.0, 2.0, 3.0, 4.0])
    r = en.power_law_fit(T, m, TF_GL)
    assert r["n"] is None and "not positive" in r["reason"]


# ---- H1 -----------------------------------------------------------------------------------

def test_h1_pure_power_law_in_range_passes(tmp_path):
    specs = build_core(tmp_path, lambda T, z: 2.0 * (T - TF_GL) ** 1.6)
    v = en.h1_verdicts(en.collect(tmp_path, cases=specs), melt_key="mean_melt_shelf")
    for z in ZC:
        r = v["per_z_c"][f"{z:.0f}"]
        assert r["power_law"]["n"] == pytest.approx(1.6, abs=1e-10)
        assert r["quadratic"]["R2"] >= 0.98
        assert r["verdict"] == "PASS"
    assert v["T_f_GL_degC"] == pytest.approx(TF_GL, rel=1e-12)


def test_h1_linear_melt_fails_on_n_although_quadratic_fits(tmp_path):
    """m linear in Tbot: R^2 = 1 (quadratic passes) but n = 1 < 1.3, so H1 fails."""
    specs = build_core(tmp_path, lambda T, z: 10.0 * (T - TF_GL))
    v = en.h1_verdicts(en.collect(tmp_path, cases=specs), melt_key="mean_melt_shelf")
    r = v["per_z_c"]["400"]
    assert r["quadratic"]["R2"] == pytest.approx(1.0, abs=1e-12)
    assert r["quadratic"]["pass"] is True
    assert r["power_law"]["n"] == pytest.approx(1.0, abs=1e-10)
    assert r["power_law"]["pass"] is False and r["verdict"] == "FAIL"


def test_h1_gap_makes_that_z_c_not_evaluable_only(tmp_path):
    specs = build_core(tmp_path, lambda T, z: 2.0 * (T - TF_GL) ** 1.6,
                       failed={(1.0, 250.0)}, skip={(0.0, 550.0)})
    v = en.h1_verdicts(en.collect(tmp_path, cases=specs), melt_key="mean_melt_shelf")
    assert v["per_z_c"]["250"]["verdict"] == "NOT_EVALUABLE"
    assert v["per_z_c"]["250"]["gaps"] == ["ens_T1.0_zc250_s1.0"]
    assert v["per_z_c"]["550"]["verdict"] == "NOT_EVALUABLE"
    assert v["per_z_c"]["550"]["gaps"] == ["ens_T0.0_zc550_s1.0"]
    assert v["per_z_c"]["400"]["verdict"] == "PASS"


def test_h1_ignores_test_split_and_geometry_runs(tmp_path):
    specs = build_core(tmp_path, lambda T, z: 2.0 * (T - TF_GL) ** 1.6)
    specs.append(write_case(tmp_path, "ens_T2.0_zc400_s1.0", 2.0, 400.0, 1.0, m=1e4,
                            group="extrapolation"))
    specs.append(write_case(tmp_path, "ens_T0.3_zc400_s0.8", 0.3, 400.0, 0.8, m=1e4,
                            group="geometry"))
    v = en.h1_verdicts(en.collect(tmp_path, cases=specs), melt_key="mean_melt_shelf")
    r = v["per_z_c"]["400"]
    assert r["n_runs"] == 5 and r["power_law"]["n"] == pytest.approx(1.6, abs=1e-10)


def test_h1_requires_an_explicit_melt_key(tmp_path):
    with pytest.raises(TypeError):
        en.h1_verdicts({})                                   # no default mask (pd0 silent)
    with pytest.raises(ValueError):
        en.h1_verdicts({"case": []}, melt_key="mean_melt")


# ---- Spearman and H2 ----------------------------------------------------------------------

def test_spearman_known_values():
    x = [1, 2, 3, 4, 5]
    assert en.spearman(x, [10, 20, 30, 40, 50]) == pytest.approx(1.0)
    assert en.spearman(x, [5, 4, 3, 2, 1]) == pytest.approx(-1.0)
    # 1 - 6 sum d^2 / (n (n^2 - 1)): ranks [5, 3, 4, 1, 2] -> sum d^2 = 36 -> -0.8
    assert en.spearman(x, [50, 30, 40, 10, 20]) == pytest.approx(-0.8, abs=1e-12)
    # ties: average ranks, then Pearson on ranks
    rx = np.array([1, 2, 3, 4, 5.0]); ry = np.array([1, 2.5, 2.5, 4, 5])
    ref = np.corrcoef(rx, ry)[0, 1]
    assert en.spearman(x, [1, 7, 7, 9, 10]) == pytest.approx(ref, rel=1e-12)
    assert en.spearman(x, [3, 3, 3, 3, 3]) is None          # constant: undefined


def test_h2_monotone_passes_non_monotone_fails_boundary_passes(tmp_path):
    # D_out = 30 w + 50 (1 - w): w increasing in Tbot -> D_out decreasing
    order = {250.0: [0.1, 0.3, 0.5, 0.7, 0.9],              # rho = -1        -> PASS
             400.0: [0.1, 0.5, 0.3, 0.9, 0.7],              # D ranks 5,3,4,1,2: rho = -0.8 -> FAIL
             550.0: [0.1, 0.3, 0.7, 0.5, 0.9]}              # D ranks 5,4,2,3,1: rho = -0.9 -> PASS
    specs = build_core(tmp_path, lambda T, z: 10.0,
                       w_fn=lambda T, z: order[z][CORE_T.index(T)])
    v = en.h2_verdicts(en.collect(tmp_path, cases=specs))
    assert v["per_z_c"]["250"]["rho"] == pytest.approx(-1.0)
    assert v["per_z_c"]["250"]["verdict"] == "PASS"
    assert v["per_z_c"]["400"]["rho"] == pytest.approx(-0.8)
    assert v["per_z_c"]["400"]["verdict"] == "FAIL"
    assert v["per_z_c"]["550"]["rho"] == pytest.approx(-0.9)
    assert v["per_z_c"]["550"]["verdict"] == "PASS"         # pd0: rho <= -0.9
    assert v["per_z_c"]["250"]["D_out"] == pytest.approx([-zout(w) for w in order[250.0]])


def test_h2_gap_not_evaluable(tmp_path):
    specs = build_core(tmp_path, lambda T, z: 10.0, failed={(0.5, 400.0)})
    v = en.h2_verdicts(en.collect(tmp_path, cases=specs))
    assert v["per_z_c"]["400"]["verdict"] == "NOT_EVALUABLE"
    assert v["per_z_c"]["400"]["gaps"] == ["ens_T0.5_zc400_s1.0"]
    # constant D_out (w = 0.5 everywhere): rho undefined -> FAIL with a reason, not NaN
    assert v["per_z_c"]["250"]["rho"] is None and v["per_z_c"]["250"]["verdict"] == "FAIL"


def test_table_json_round_trip(tmp_path):
    specs = [write_case(tmp_path, "ens_a", 0.5, 400.0, 1.0, m=12.0)]
    tab = en.collect(tmp_path, cases=specs)
    p = en.write_table(tab, tmp_path / "t.json")
    assert en.read_table(p) == tab
    assert not math.isnan(tab["mean_melt_shelf"][0])
