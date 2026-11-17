"""analysis.h5a: the two tau_x(y) patterns of H5a (PD0 + PD1b scales) and the verdict."""
import numpy as np
import pytest

from exp5_adjoint_slope.analysis import h5a

YC = (np.arange(45) + 0.5) * 10e3


def test_patterns_have_equal_rms_over_shelf_and_slope():
    u, d = h5a.uniform(YC), h5a.dipole(YC, 37.5e3)
    m = YC < 275e3
    assert np.sqrt((u[m] ** 2).mean()) == pytest.approx(h5a.RMS)
    assert np.sqrt((d[m] ** 2).mean()) == pytest.approx(h5a.RMS)
    assert np.all(u[m] < 0) and np.all(u[~m] == 0)          # stronger easterlies: tau_x more negative
    # the dipole changes sign at the break (y0 = 125 km) and only the curl there
    assert d[YC < 115e3].max() <= 0 < d[(YC > 135e3) & (YC < 200e3)].min() or \
        d[YC < 115e3].min() >= 0 > d[(YC > 135e3) & (YC < 200e3)].max()


def test_dj_is_sens_times_pattern():
    s = np.zeros(45)
    s[10] = 2.0
    p = np.zeros(45)
    p[10] = 0.5
    assert h5a.dj(s, p) == pytest.approx(1.0)


def test_verdict_frozen_and_paired_significance():
    # members x rows: the dipole wins in every member -> pass, significant.
    rng = np.random.default_rng(0)
    d = h5a.dipole(YC, 37.5e3)
    sens = np.array([5 * d + 0.01 * rng.standard_normal(45) for _ in range(10)])
    r = h5a.evaluate(sens, YC)
    assert r["widths"]["37.5"]["pass"] and r["robust"] and r["widths"]["37.5"]["significant"]
    r2 = h5a.evaluate(-np.tile(h5a.uniform(YC), (10, 1)), YC)
    assert not r2["widths"]["37.5"]["pass"] and not r2["robust"]
