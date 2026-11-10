"""Tests for exp2_isomip_cavity/gendata/profiles.py (Task 1)."""
from __future__ import annotations

import numpy as np
import pytest

from exp2_isomip_cavity.gendata import profiles as pr

Z = -(np.arange(36) + 0.5) * 20.0           # cell centres of the 36 x 20 m grid [m]


def test_ocean0_endpoints():
    # AD16 Eq. 21-22 with Table 6 (WARM): -1.9 degC / 33.8 at 0 m, 1.0 degC / 34.7 at 720 m
    T, S = pr.ocean0(np.array([0.0, -360.0, -720.0]))
    np.testing.assert_allclose(T, [-1.9, -0.45, 1.0], atol=1e-12)
    np.testing.assert_allclose(S, [33.8, 34.25, 34.7], atol=1e-12)


def test_cold_profile():
    # AD16 Table 5: T = -1.9 degC everywhere, S from 33.8 to 34.55
    T, S = pr.cold(np.array([0.0, -720.0]))
    np.testing.assert_allclose(T, [-1.9, -1.9], atol=1e-12)
    np.testing.assert_allclose(S, [33.8, 34.55], atol=1e-12)


def test_ramp_identity_with_ocean0():
    # N4: z_c = 360 m, d = 720 m makes r = z / zb,deep, i.e. AD16 Eq. 21-22
    z = np.linspace(0.0, -720.0, 1441)
    T, S = pr.ramp(z, 1.0, 360.0, d=720.0)
    T0, S0 = pr.ocean0(z)
    assert np.max(np.abs(T - T0)) < 1e-12
    assert np.max(np.abs(S - S0)) < 1e-12


@pytest.mark.parametrize("Tbot,z_c", [(-0.5, 250.0), (0.0, 400.0), (1.5, 550.0), (2.0, 250.0)])
def test_ramp_monotone_and_bounded(Tbot, z_c):
    z = np.linspace(0.0, -720.0, 2001)
    T, S = pr.ramp(z, Tbot, z_c)
    assert np.all(np.diff(T) >= 0) and np.all(np.diff(S) >= 0)   # z decreasing: T, S increase
    assert T.min() >= -1.9 - 1e-12 and T.max() <= Tbot + 1e-12
    assert S.min() >= 33.8 - 1e-12 and S.max() <= 34.7 + 1e-12
    # surface layer above z_c - d/2 and deep layer below z_c + d/2 are uniform
    assert T[0] == pytest.approx(-1.9) and T[-1] == pytest.approx(Tbot)
    T_mid, _ = pr.ramp(np.array([-z_c]), Tbot, z_c)
    assert T_mid[0] == pytest.approx(-1.9 + 0.5 * (Tbot + 1.9))


def test_linear_eos_reference_state():
    # AD16 Eq. 23: rho(Tref, Sref) = rho_ref
    assert pr.rho_linear(-1.0, 34.2) == pytest.approx(1027.51, abs=1e-12)
    assert pr.rho_linear(0.0, 34.2) == pytest.approx(1027.51 * (1 - 3.733e-5))


def test_load_anomaly_matches_trapezoidal_hydrostatic_integral():
    # Ocean0 profiles are linear in z, so the cell-centre density is exact for each cell
    # integral; compare with a fine trapezoidal integration of the continuous profile.
    dz = np.full(36, 20.0)
    T, S = pr.ocean0(Z)
    p = pr.load_anomaly(T, S, dz)
    assert p.shape == (37,) and p[0] == 0.0
    zf = np.concatenate([[0.0], -np.cumsum(dz)])
    zz = np.linspace(0.0, -720.0, 720 * 50 + 1)
    Tz, Sz = pr.ocean0(zz)
    integrand = pr.GRAVITY * (pr.rho_linear(Tz, Sz) - pr.RHO_REF)
    for k in range(1, 37):
        m = zz >= zf[k] - 1e-9
        ref = np.trapezoid(integrand[m], -zz[m])
        assert abs(p[k] - ref) <= 1e-8 * abs(ref)


def test_load_anomaly_at_draft_interpolates_within_cell():
    dz = np.full(36, 20.0)
    T, S = pr.ocean0(Z)
    p = pr.load_anomaly(T, S, dz)
    draft = np.array([[0.0, -20.0, -30.0, -720.0]])
    got = pr.load_anomaly_at_draft(p, dz, draft)
    np.testing.assert_allclose(got, [[0.0, p[1], 0.5 * (p[1] + p[2]), p[36]]], rtol=1e-14)
