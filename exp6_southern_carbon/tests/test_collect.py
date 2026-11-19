import numpy as np
import pytest

from exp6_southern_carbon.analysis import collect
from exp6_southern_carbon.gendata import co2

SPY = 31536
MON = SPY // 12


def test_annual_means_keyed_by_absolute_end_iteration():
    start = 5_000_000
    its = [start + MON * m for m in range(1, 25)]            # two model years of monthly records
    vals = [1.0] * 12 + [3.0] * 12
    a = collect.annual_means(its, vals, start, SPY)
    assert a == {start + SPY: 1.0, start + 2 * SPY: 3.0}


def test_annual_means_drops_incomplete_year():
    start = 0
    its = [MON * m for m in range(1, 19)]                     # 1.5 years
    a = collect.annual_means(its, np.ones(18), start, SPY)
    assert list(a) == [SPY]


def test_branch_sensitivity_years_6_to_15():
    n0 = 1_000_000
    fp = {n0 + k * SPY: 0.5 for k in range(1, 26)}
    fb = {n0 + k * SPY: 0.5 + (0.3 if 6 <= k <= 15 else 9.0) for k in range(1, 26)}
    assert collect.branch_sensitivity(fb, fp, n0, SPY) == pytest.approx(0.3)


def test_branch_sensitivity_missing_parent_year_raises():
    n0 = 0
    fb = {k * SPY: 1.0 for k in range(1, 26)}
    fp = {k * SPY: 1.0 for k in range(1, 10)}
    with pytest.raises(KeyError):
        collect.branch_sensitivity(fb, fp, n0, SPY)


def test_winter_mld_is_max_of_monthly_climatology():
    months = np.tile(np.array([50.0] * 5 + [400.0, 600.0, 650.0, 500.0] + [50.0] * 3), 3)
    months[19] = 950.0                                         # one August out of three years
    assert collect.winter_mld(months) == pytest.approx((650.0 * 2 + 950.0) / 3)


def test_export_fraction():
    assert collect.export_fraction({"adv": 3.0, "mix": 1.0, "gchem": -3.0}) == pytest.approx(0.75)


from pathlib import Path

T2 = Path("/Volumes/LaCie/portfolio_runs/exp6/T2")
T2_OK = (T2 / "run_status.json").exists() and '"OK"' in (T2 / "run_status.json").read_text()


@pytest.mark.skipif(not T2_OK, reason="T2 not available")
def test_readers_on_t2():
    c = collect.io.case(T2)
    n0, spy = int(c["nIter0"]), int(c["steps_per_year"])
    f = collect.flux_series(T2)
    assert list(f) == [n0 + spy, n0 + 2 * spy] and all(np.isfinite(list(f.values())))
    w = (2, 2)                                    # year 1 of a cold start has no initial pickup
    nb = n0                                       # a "branch" starting where T2 starts
    circ = collect.circulation(T2, nb, window=w)
    assert 0.0 < circ["mld_inv"] < 4000.0 and np.isfinite(circ["moc_sv"])
    assert collect.measured(T2, T2, nb, window=w) == 0.0
    assert collect.branch_sensitivity(f, f, nb, spy, window=w) == 0.0


def test_atm_pco2_fixed(tmp_path):
    (tmp_path / "data.dic").write_text(" &DIC_FORCING\n  dic_int1=1,\n  dic_pCO2=2.8431698608e-04,\n &\n")
    assert collect.atm_pco2(tmp_path, [10, 20]) == pytest.approx([2.8431698608e-04] * 2)


def test_atm_pco2_transient_matches_dic_atmos(tmp_path):
    (tmp_path / "data.dic").write_text(" &DIC_FORCING\n  dic_int1=2,\n  dic_int2=4,\n  dic_int3=100,\n"
                                       "  dic_int4=10,\n &\n")
    (tmp_path / "co2atmos.dat").write_text("2.8e-04\n2.8e-04\n3.0e-04\n3.0e-04\n")
    got = collect.atm_pco2(tmp_path, [105, 110])
    want = [co2.atm_pco2_at(i, [2.8e-4, 2.8e-4, 3.0e-4, 3.0e-4], 100, 10) for i in (105, 110)]
    assert got == pytest.approx(want)


@pytest.mark.skipif(not T2_OK, reason="T2 not available")
def test_temperature_term_on_t2_is_zero_against_itself():
    c = collect.io.case(T2)
    t = collect.temperature(T2, T2, int(c["nIter0"]), window=(2, 2))
    assert t["temp"] == 0.0 and t["conductance"] > 0.0 and np.isfinite(t["conductance"])
