"""Analytical tests for the buoyant plume models (criterion C1).

All tests except the volume-budget one switch melting off with
``MeltParams(melt_off=True)`` so that the classical Morton et al. (1956)
similarity solutions apply.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import cumulative_simpson

from exp1_fjord_plumes.theory.melt import MeltParams
from exp1_fjord_plumes.theory.plume import LinearEOS, line_plume, point_plume

EOS = LinearEOS()
G = EOS.g
T_AMB, S_AMB = 2.0, 34.0
T0, S0 = 0.0, 0.0
NOMELT = MeltParams(melt_off=True)


def _gprime0(Ta=T_AMB, Sa=S_AMB, T=T0, S=S0):
    # Jenkins (2011) Eq. 18 times g: g' = g [bS (Sa - S) - bT (Ta - T)]  [m s^-2]
    return G * (EOS.sBeta * (Sa - S) - EOS.tAlpha * (Ta - T))


def _strat_salinity(N, z_base, S_base=S_AMB):
    """Linear Sa(z) with constant N^2 = -g sBeta dSa/dz (T uniform)."""
    dSdz = -N**2 / (G * EOS.sBeta)
    return lambda z: S_base + dSdz * (np.asarray(z) - z_base)


def _above(sol, nb0):
    b0 = sol.b[0]
    return np.isfinite(sol.u) & (sol.z - sol.z[0] > nb0 * b0)


# ---------------------------------------------------------------------------
# Unstratified point source: Morton et al. (1956) Eq. 4 (top-hat form)
# b = (6 alpha / 5) z,  u ∝ z^(-1/3)
# ---------------------------------------------------------------------------

def _full_cone_uniform(alpha=0.1):
    z = np.linspace(-1000.0, 0.0, 4001)
    return point_plume(z, T_AMB, S_AMB, 1.0, T0, S0, NOMELT, "jenkins",
                       alpha=alpha, geometry="full_cone"), alpha


def test_full_cone_unstratified_width_slope():
    sol, alpha = _full_cone_uniform()
    sel = _above(sol, 20)
    dbdz = np.gradient(sol.b, sol.z)[sel]
    assert np.allclose(dbdz, 6 * alpha / 5, rtol=0.02)


def test_full_cone_unstratified_velocity_power():
    # The similarity law holds relative to the virtual origin z_v of the
    # finite source (Morton et al. 1956, "Determination of the virtual
    # source"), estimated here from a straight-line fit b = s (z - z_v).
    sol, _ = _full_cone_uniform()
    sel = _above(sol, 20)
    s, c = np.polyfit(sol.z[sel], sol.b[sel], 1)
    z_v = -c / s
    slope = np.polyfit(np.log(sol.z[sel] - z_v), np.log(sol.u[sel]), 1)[0]
    assert slope == pytest.approx(-1.0 / 3.0, abs=0.01)


def test_line_plume_unstratified_constant_velocity():
    alpha = 0.1
    z = np.linspace(-500.0, 0.0, 2001)
    sol = line_plume(z, T_AMB, S_AMB, 0.01, T0, S0, NOMELT, "jenkins",
                     alpha=alpha)
    sel = _above(sol, 20)
    u = sol.u[sel]
    assert np.all(np.abs(u / u.mean() - 1.0) < 0.02)
    dbdz = np.gradient(sol.b, sol.z)[sel]
    assert np.allclose(dbdz, alpha, rtol=0.02)


# ---------------------------------------------------------------------------
# Linearly stratified point source: maximum height
# ---------------------------------------------------------------------------
# Morton et al. (1956) Eq. 14: x = 0.410 alpha^(-1/2) F0^(1/4) G^(-3/8) x1,
# with G = N^2, F0 the source buoyancy flux, Gaussian profiles and
# alpha = 0.093; the vertical velocity vanishes at x1 ≈ 2.8 (p. 9). Hence
#     z_max N^(3/4) / F0^(1/4) = 0.410 * 2.8 * 0.093^(-1/2) = 3.76.
# The brief's 5.0 is NOT what the paper gives for F0 = g'0 Q0 (it equals
# 3.76 * pi^(1/4), i.e. F0 defined per pi [verificar]); see bitacora.md.
# Our equations are top-hat; written in fluxes (Q, M, F) the top-hat and
# Gaussian systems are identical when alpha_tophat = sqrt(2) alpha_gauss,
# so the paper's alpha = 0.093 is used as alpha = sqrt(2) * 0.093 here.
MTT_CONST = 0.410 * 2.8 * 0.093 ** -0.5
ALPHA_MTT_TOPHAT = math.sqrt(2.0) * 0.093


def _zmax_from_virtual_origin(Q0, N, alpha=ALPHA_MTT_TOPHAT):
    z_base = -500.0
    z = np.linspace(z_base, 0.0, 5001)
    sol = point_plume(z, T_AMB, _strat_salinity(N, z_base), Q0, T0, S0,
                      NOMELT, "jenkins", alpha=alpha, geometry="full_cone")
    assert not sol.reached_top
    # Pure-plume source (G0 = 1): virtual origin 5 b0 / (6 alpha) below the
    # source (Morton et al. 1956 Eq. 4 / Slater et al. 2016 Eq. 6).
    z_virtual = z_base - 5.0 * sol.b[0] / (6.0 * alpha)
    F0 = _gprime0() * Q0  # total source buoyancy flux [m^4 s^-3]
    return sol, sol.z_detach - z_virtual, F0


def test_full_cone_linear_strat_max_height():
    N = 1e-3
    sol, zmax, F0 = _zmax_from_virtual_origin(1e-3, N)
    assert zmax * N**0.75 / F0**0.25 == pytest.approx(MTT_CONST, abs=0.3)
    # neutral level below the top (MTT: x1 = 2.125 vs 2.8)
    assert sol.z_neutral < sol.z_detach


def test_point_plume_scaling_F0():
    N = 1e-3
    _, z1, F1 = _zmax_from_virtual_origin(1e-3, N)
    _, z2, F2 = _zmax_from_virtual_origin(2e-3, N)
    assert F2 / F1 == pytest.approx(2.0, rel=1e-12)
    assert z2 / z1 == pytest.approx(2.0 ** 0.25, rel=0.01)


# ---------------------------------------------------------------------------
# Geometry consistency and source conditions
# ---------------------------------------------------------------------------

def test_slater_source_conditions():
    # Slater et al. (2016) Eq. 5 with G0 = 1.
    alpha, Q0 = 0.1, 100.0
    z = np.linspace(-500.0, 0.0, 11)
    sol = point_plume(z, T_AMB, S_AMB, Q0, T0, S0, NOMELT, "jenkins",
                      alpha=alpha, geometry="half_cone")
    g0 = _gprime0()
    b0 = (32 * alpha * Q0**2 / (5 * math.pi**2 * g0)) ** 0.2
    u0 = (2 / math.pi) * (5 * math.pi**2 * g0 / (32 * alpha)) ** 0.4 * Q0**0.2
    assert sol.b[0] == pytest.approx(b0, rel=1e-12)
    assert sol.u[0] == pytest.approx(u0, rel=1e-12)


def test_half_cone_without_wall_terms_equals_full_cone():
    # Half a full cone: same b and u for half the volume flux when drag and
    # melt are off (Slater et al. 2016 Eqs. 1a-1d vs Morton et al. 1956).
    p = MeltParams(melt_off=True, Cd=0.0)
    z = np.linspace(-500.0, 0.0, 501)
    full = point_plume(z, T_AMB, S_AMB, 2.0, T0, S0, p, "jenkins",
                       geometry="full_cone")
    half = point_plume(z, T_AMB, S_AMB, 1.0, T0, S0, p, "jenkins",
                       geometry="half_cone")
    assert np.allclose(half.b, full.b, rtol=1e-6)
    assert np.allclose(half.u, full.u, rtol=1e-6)


def test_line_plume_stratified_detaches_below_surface():
    z_base = -500.0
    z = np.linspace(z_base, 0.0, 1001)
    sol = line_plume(z, T_AMB, _strat_salinity(5e-3, z_base), 1e-3, T0, S0,
                     MeltParams(), "mitgcm")
    assert not sol.reached_top
    assert z_base < sol.z_neutral < sol.z_detach < 0.0
    assert np.all(np.isnan(sol.u[sol.z > sol.z_detach]))


# ---------------------------------------------------------------------------
# Volume budget with melting (Jenkins 2011 Eq. 1)
# ---------------------------------------------------------------------------

def test_volume_budget_with_melt():
    alpha = 0.1
    p = MeltParams()
    z = np.linspace(-500.0, 0.0, 20001)
    sol = line_plume(z, 4.0, S_AMB, 0.01, T0, S0, p, "jenkins", alpha=alpha)
    ok = np.isfinite(sol.u)
    zz, Q = sol.z[ok], (sol.b * sol.u)[ok]
    entr = alpha * sol.u[ok]  # e = E0 U, Jenkins (2011) Eq. 6 [m s^-1]
    meltw = sol.melt[ok] * p.rho_i / p.rho0  # meltwater volume flux [m s^-1]
    budget = cumulative_simpson(entr + meltw, x=zz, initial=0.0)
    melt_part = cumulative_simpson(meltw, x=zz, initial=0.0)
    lhs = Q - Q[0]
    # melting must contribute far more than the tolerance for the test to bite
    assert melt_part[-1] / lhs[-1] > 1e-4
    rel = np.abs(lhs[1:] - budget[1:]) / np.abs(lhs[1:])
    assert rel.max() < 1e-6


# ---------------------------------------------------------------------------
# Follow-up R15: half cone with melt must terminate cleanly near u -> 0
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", ["mitgcm", "constant"])
def test_half_cone_with_melt_terminates_in_sermilik(mode):
    import warnings

    from exp1_fjord_plumes.gendata.profiles import continuous_sermilik

    H = 500.0
    z = np.linspace(-H, 0.0, 2001)
    Ta, Sa = continuous_sermilik(z)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)  # no overflow allowed
        sol = point_plume(z, Ta, Sa, 100.0, 0.0, 0.0, MeltParams(ice_depth=H),
                          mode, geometry="half_cone")
    assert np.isfinite(sol.z_detach)
    assert np.isfinite(sol.z_neutral)
    assert sol.z_neutral < sol.z_detach <= 0.0
    ok = np.isfinite(sol.Q)
    assert ok.sum() > 10 and np.all(sol.z[ok] <= sol.z_detach)
    assert np.all(np.diff(sol.Q[ok]) > 0.0)
    assert np.all(np.isfinite(sol.melt[ok])) and np.all(sol.melt[ok] > 0.0)
