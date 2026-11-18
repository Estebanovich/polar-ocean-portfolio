import numpy as np
import pytest

from exp6_southern_carbon.analysis import drift, flux, gradients

YEAR_S = 31536000.0


def test_flux_into_ocean_units_and_sign():
    # DICCFLX = FluxCO2 > 0 when ff*pCO2atm > k0*fugf*pCO2ocean (pkg/dic/dic_surfforcing.F:261-265)
    f = flux.flux_into_ocean(np.array([1e-8, -2e-8]))
    np.testing.assert_allclose(f, [1e-8 * YEAR_S, -2e-8 * YEAR_S])


def test_channel_mean_area_weighted():
    f2d = np.array([[1.0, 3.0], [5.0, 100.0]])
    area = np.array([[1.0, 1.0], [2.0, 1.0]])
    wet = np.array([[True, True], [True, False]])
    assert flux.channel_mean(f2d, area, wet) == pytest.approx((1 + 3 + 10) / 4)


def test_sensitivity_window():
    parent = np.zeros(25)
    branch = parent.copy()
    branch[:5] += 100.0
    branch[5:15] += 2.0
    branch[15:] += -50.0
    assert flux.sensitivity(branch, parent, window=(6, 15)) == pytest.approx(2.0)


def test_g200_linear_profile():
    z = np.array([2.5, 50.0, 150.0, 250.0, 1200.0])
    dic = 2.0 + 1e-3 * z
    mask = np.ones((1, 1), bool)
    g = gradients.g200(dic[:, None, None], z, mask)
    # surface = level 1 (z = 2.5 m): [DIC(200) - DIC(2.5)] / 200 m, mol m-3 -> umol kg-1 (rho0 1035)
    assert g == pytest.approx(1e-3 * (200.0 - 2.5) / 200.0 * 1e6 / 1035.0)
    d = gradients.ddic1000(dic[:, None, None], z, mask)
    assert d == pytest.approx(1e-3 * (1000.0 - 2.5) * 1e6 / 1035.0)


def test_v6a_rule():
    yrs = np.arange(0, 1100, 10)
    flat = drift.v6a(yrs, np.zeros(yrs.size), np.zeros(yrs.size))
    assert flat["pass"]
    trend = drift.v6a(yrs, 0.02 * yrs / 100.0, np.zeros(yrs.size))
    assert not trend["pass"]
    small = drift.v6a(yrs, 0.005 * yrs / 100.0, 0.5e-3 * yrs / 100.0)   # deep DIC in mol m-3: 0.5e-3 < 1 umol/kg
    assert small["pass"]
    short = drift.v6a(np.arange(0, 900, 10), np.zeros(90), np.zeros(90))
    assert not short["pass"] and "1000" in short["reason"]


def test_v6a_deep_dic_threshold():
    yrs = np.arange(0, 1100, 10)
    big = drift.v6a(yrs, np.zeros(yrs.size), 1.2e-3 * yrs / 100.0)    # 1.2e-3 mol m-3 > 1.035e-3
    assert not big["pass"]
