"""analysis.ensemble: post hoc ensemble criteria of the PD2 amendment (notes/pd0_preregistro.md)."""
import numpy as np
import pytest

from exp5_adjoint_slope.analysis import ensemble as en


def test_v2e_ratio_of_member_means():
    dp, dm = np.array([1.0, 1.2, 0.8]), np.array([-0.9, -1.1, -0.9])
    r = en.v2e(dp, dm)
    assert r["value"] == pytest.approx(abs((dp + dm).mean()) / abs((dp - dm).mean()))
    assert r["pass"] is True
    assert en.v2e(dp, -0.5 * dp)["pass"] is False


def test_v3e_needs_agreement_and_detectability():
    dp, dm = np.array([1.0, 1.1, 0.9, 1.0]), -np.array([1.0, 1.1, 0.9, 1.0])
    r = en.v3e(1.1, dp, dm)
    assert r["fwd"] == pytest.approx(1.0) and r["detectable"] and r["pass"]
    assert r["value"] == pytest.approx(0.1)
    assert en.v3e(2.0, dp, dm)["pass"] is False
    noisy = en.v3e(0.1, np.array([1.0, -1.0, 1.2, -1.1]), np.array([-1.1, 1.0, -1.0, 1.1]))
    assert noisy["detectable"] is False and noisy["pass"] is False


def test_t_half_interpolates_after_the_maximum():
    leads = np.array([0.0, 30.0, 60.0])
    assert en.t_half(leads, np.array([10.0, 6.0, 2.0])) == pytest.approx(30 + 30 * (1 / 4))
    assert en.t_half(leads, np.array([10.0, 8.0, 7.0])) is None
    assert en.t_half(leads, np.array([4.0, 10.0, 5.0])) == pytest.approx(60.0)


def test_h5c_verdicts():
    leads = np.array([0.0, 30.0, 60.0])
    assert en.h5c(leads, np.array([10.0, 6.0, 2.0]))["verdict"] == "pass"
    assert en.h5c(leads, np.array([10.0, 2.0, 1.0]))["verdict"] == "fail"
    assert en.h5c(leads, np.array([10.0, 8.0, 7.0]))["verdict"] == "not testable beyond 60 days"
