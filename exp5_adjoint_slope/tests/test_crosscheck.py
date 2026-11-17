"""analysis.crosscheck (plan Task 8): H5b detectability and verdict, V5 eddy fraction."""
import numpy as np
import pytest

from exp5_adjoint_slope.analysis import crosscheck as cc


def test_sigma_delta():
    assert cc.sigma_delta([1, 2, 3], [4, 5, 6]) == pytest.approx(np.sqrt(2 / 3))


def test_h5b_needs_detectability():
    r = cc.h5b([0.0, 1.0, -1.0], [0.5, 1.5, -0.5], dj_pred=0.5)
    assert r["sign_ok"] and r["within_factor"] and not r["detectable"] and not r["pass"]
    r = cc.h5b([0.0, 0.01, -0.01], [1.0, 1.01, 0.99], dj_pred=0.8)
    assert r["detectable"] and r["sign_ok"] and r["within_factor"] and r["pass"]
    assert not cc.h5b([0.0, 0.01, -0.01], [1.0, 1.01, 0.99], dj_pred=-0.8)["pass"]
    assert not cc.h5b([0.0, 0.01, -0.01], [1.0, 1.01, 0.99], dj_pred=3.0)["within_factor"]


def test_h5b_indicative_with_two_members():
    assert cc.h5b([0.0, 0.01], [1.0, 1.01], dj_pred=1.0)["indicative"] is True
    assert cc.h5b([0.0, 0.01, 0.0], [1.0, 1.01, 1.0], dj_pred=1.0)["indicative"] is False


def _grid(nr=2, ny=4, nx=3):
    return {"hfac_s": np.ones((nr, ny, nx)), "drf": np.ones(nr), "dxg": np.ones(nx)}


def test_eddy_fraction_zero_without_variation():
    v = np.ones((5, 2, 4, 3))
    th = np.full((5, 2, 4, 3), 2.0)
    assert cc.eddy_fraction(v, th, _grid(), 2) == pytest.approx(0.0)


def test_eddy_fraction_one_for_pure_covariance():
    rng = np.random.default_rng(1)
    a = rng.standard_normal((50, 2, 4, 3))
    a -= a.mean(axis=(0, 3), keepdims=True)
    th = a + 0.0
    v = a.copy()
    # zero-mean v and theta (theta_f = 0): all the transport is v'theta'
    assert cc.eddy_fraction(v, th, _grid(), 2, theta_f=0.0) == pytest.approx(1.0)


def test_eke_removes_the_zonal_mean():
    # u = U(z,y) + a cos(kx): EKE = a^2/4 per component; a zonally uniform flow has none.
    nr, ny, nx = 2, 3, 8
    x = np.arange(nx)
    u = 0.3 + 0.1 * np.cos(2 * np.pi * x / nx) * np.ones((nr, ny, 1))
    v = np.full((nr, ny, nx), -0.2)
    wet = np.ones((nr, ny, nx), bool)
    e = cc.eke(u, v, u ** 2, v ** 2, wet)
    assert e == pytest.approx(0.5 * 0.1 ** 2 / 2)
    assert cc.eke(v, v, v ** 2, v ** 2, wet) == pytest.approx(0.0, abs=1e-15)


def test_eke_equilibrated_rule():
    # PD2 correction 2026-11-17: mean of the last 3 months within 10 % of the previous 3.
    assert cc.eke_equilibrated([1, 1, 1, 5, 5, 5, 9, 9, 9, 9.5, 9.5, 9.5])["pass"]
    r = cc.eke_equilibrated([1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4])
    assert not r["pass"] and r["rel_change"] == pytest.approx(1 / 3)
