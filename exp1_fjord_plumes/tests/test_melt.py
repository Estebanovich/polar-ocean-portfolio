"""Tests for the three-equation ice-ocean melt parameterisation.

The reference formulation is Holland & Jenkins (1999) as written in
Jenkins (2011, Eqs. 7-9) and Holland et al. (2008, Eqs. 2-4); the "mitgcm"
mode must reproduce ``pkg/icefront/icefront_thermodynamics.F`` exactly.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from exp1_fjord_plumes.theory.melt import (
    MeltParams,
    freezing_point,
    three_eq_melt,
    transfer_velocities,
)

S_AMB = 34.0  # psu
Z_FACE = -300.0  # m (negative downward)


@pytest.mark.parametrize("mode", ["mitgcm", "jenkins", "constant"])
def test_no_melt_at_freezing(mode):
    p = MeltParams()
    Tf = freezing_point(S_AMB, Z_FACE, p)
    m, Tb, Sb = three_eq_melt(Tf, S_AMB, Z_FACE, 0.5, p, mode)
    assert abs(m) < 1e-12
    assert Sb == pytest.approx(S_AMB, rel=1e-10)
    assert Tb == pytest.approx(Tf, abs=1e-10)


@pytest.mark.parametrize("mode", ["mitgcm", "jenkins", "constant"])
def test_melt_monotonic_in_thermal_driving(mode):
    p = MeltParams()
    Tf = freezing_point(S_AMB, Z_FACE, p)
    melts = [
        three_eq_melt(Tf + dT, S_AMB, Z_FACE, 0.5, p, mode)[0]
        for dT in (0.5, 1.0, 2.0, 4.0)
    ]
    assert melts[0] > 0.0
    assert all(b > a for a, b in zip(melts, melts[1:]))


def test_jenkins_mode_scales_with_velocity():
    # With constant Stanton numbers (Jenkins 2011, Table 1) gamma_S/gamma_T is
    # independent of u, so Tb and Sb are too, and m is exactly linear in u.
    # The brief's bound is [1.8, 2.0]; the upper bound gets a round-off margin
    # only (1e-12 relative), not a physical tolerance.
    p = MeltParams()
    Tf = freezing_point(S_AMB, Z_FACE, p)
    m1 = three_eq_melt(Tf + 4.0, S_AMB, Z_FACE, 0.3, p, "jenkins")[0]
    m2 = three_eq_melt(Tf + 4.0, S_AMB, Z_FACE, 0.6, p, "jenkins")[0]
    ratio = m2 / m1
    assert 1.8 <= ratio <= 2.0 * (1.0 + 1e-12)


def test_mitgcm_transfer_floor_and_linear_branch():
    # icefront_thermodynamics.F:91-97: gammaT = max(1e-2*sqrt(1.5e-3)*|w|, 1e-4)
    p = MeltParams()
    gT_slow, gS_slow = transfer_velocities(0.05, p, "mitgcm")
    assert gT_slow == pytest.approx(1e-4, rel=1e-14)
    assert gS_slow == pytest.approx(5.05e-3 * 1e-4, rel=1e-14)
    gT_fast, _ = transfer_velocities(1.0, p, "mitgcm")
    assert gT_fast == pytest.approx(1e-2 * math.sqrt(1.5e-3), rel=1e-14)
    # Below the crossover |w| = 1e-4 / (1e-2*sqrt(1.5e-3)) ~ 0.258 m/s the
    # melt rate does not depend on velocity.
    Tf = freezing_point(S_AMB, Z_FACE, p)
    ma = three_eq_melt(Tf + 2.0, S_AMB, Z_FACE, 0.05, p, "mitgcm")[0]
    mb = three_eq_melt(Tf + 2.0, S_AMB, Z_FACE, 0.20, p, "mitgcm")[0]
    assert ma == pytest.approx(mb, rel=1e-14)


def _fortran_icefront_fwflux(tLoc, sLoc, pLoc, w, thetaIce, rhoConst, cp_w,
                             latentHeat, cp_i):
    """Line-by-line transcription of icefront_thermodynamics.F:77-146.

    Returns FreshWaterFlux in kg m^-2 s^-1 (negative = melting).
    """
    a0, c0, b = -0.0575, 0.0901, -7.61e-4
    heatT = max(1.0e-2 * abs(w) * math.sqrt(1.5e-3), 1.0e-4)
    saltT = 5.05e-3 * heatT
    eps1 = rhoConst * cp_w * heatT
    eps2 = rhoConst * latentHeat * saltT
    eps3 = rhoConst * cp_i * saltT
    aqe = a0 * (-eps1 + eps3)
    recip_aqe = 0.5 / aqe
    eps4 = b * pLoc + c0
    eps6 = eps4 - tLoc
    eps7 = eps4 - thetaIce
    bqe = -eps1 * eps6 - sLoc * a0 * eps3 + eps3 * eps7 + eps2
    cqe = -(eps2 + eps3 * eps7) * sLoc
    discrim = bqe * bqe - 4.0 * aqe * cqe
    saltFreeze = (-bqe - math.sqrt(discrim)) * recip_aqe
    if saltFreeze < 0.0:
        saltFreeze = (-bqe + math.sqrt(discrim)) * recip_aqe
    thetaFreeze = a0 * saltFreeze + eps4
    fw = eps1 * (thetaFreeze - tLoc) / (latentHeat + cp_i * (thetaFreeze - thetaIce))
    return fw, thetaFreeze, saltFreeze


@pytest.mark.parametrize("w", [0.05, 0.8])
@pytest.mark.parametrize("dT", [0.3, 3.0])
def test_mitgcm_mode_matches_fortran_transcription(w, dT):
    p = MeltParams(ice_depth=500.0)
    z = -350.0
    depth = -z
    thetaIce = p.T_ice_surface * (p.ice_depth - depth) / p.ice_depth
    Tf = freezing_point(S_AMB, z, p)
    T = Tf + dT
    fw, TbF, SbF = _fortran_icefront_fwflux(T, S_AMB, depth, w, thetaIce,
                                            p.rho0, p.c_w, p.L, p.c_i)
    m, Tb, Sb = three_eq_melt(T, S_AMB, z, w, p, "mitgcm")
    assert Sb == pytest.approx(SbF, rel=1e-10)
    assert Tb == pytest.approx(TbF, rel=1e-10, abs=1e-12)
    # m is ice-thickness rate: m = -FreshWaterFlux / rho_ice
    assert m == pytest.approx(-fw / p.rho_i, rel=1e-10)
    assert fw < 0.0  # icefront sign convention: negative flux = melting


def test_melt_off_returns_zero():
    p = MeltParams(melt_off=True)
    m, Tb, Sb = three_eq_melt(5.0, S_AMB, Z_FACE, 0.5, p, "jenkins")
    assert m == 0.0
    assert np.isnan(Tb) and np.isnan(Sb)


def test_unknown_mode_raises():
    with pytest.raises(ValueError):
        three_eq_melt(1.0, S_AMB, Z_FACE, 0.5, MeltParams(), "bogus")
