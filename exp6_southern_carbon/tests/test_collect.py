import numpy as np
import pytest

from exp6_southern_carbon.analysis import collect

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
