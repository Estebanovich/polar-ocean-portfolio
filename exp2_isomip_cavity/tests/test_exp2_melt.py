"""Tests for exp2_isomip_cavity/theory/melt.py: the ISOMIP+ three-equation melt law
(AD16 sect. 3.1.8, Eq. 24-27, Table 4) as implemented in the Exp. 2 copy of
pkg/shelfice/shelfice_thermodynamics.F (config/ocean0/code/).

(File name test_exp2_melt.py: exp1_fjord_plumes/tests/test_melt.py already exists and the
test directories are not packages, so the same basename would collide.)
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy.optimize import brentq

from exp2_isomip_cavity.theory import melt as ml

PKG = Path(__file__).resolve().parents[1]
FORTRAN = PKG / "config" / "ocean0" / "code" / "shelfice_thermodynamics.F"

# AD16 Table 4
LAM1, LAM2, LAM3 = -0.0573, 0.0832, -7.53e-8      # degC/PSU, degC, degC/Pa
CW, LF, RHO_SW = 3974.0, 3.34e5, 1027.51           # J/kg/degC, J/kg, kg/m^3 (rhoConst)

# the hand case of the Task 2 brief
T, S, P, USTAR, GT = 1.0, 34.5, 5.0e6, 0.01, 0.05


def _hand_solution(T, S, p, ustar, gammaT):
    """Independent solution of AD16 Eq. 24-26 (bracketing root finder, not the quadratic):
    rho_fw m L = rho_sw cw u* GT (T - Tb); rho_fw m Sb = rho_sw u* GS (S - Sb);
    Tb = lam1 Sb + lam2 + lam3 p. Eliminating m gives f(Sb) = 0 on (0, inf)."""
    gT, gS = gammaT * ustar, gammaT / 35.0 * ustar

    def f(Sb):
        Tb = LAM1 * Sb + LAM2 + LAM3 * p
        return LF * gS * (S - Sb) - CW * gT * (T - Tb) * Sb

    Sb = brentq(f, 1e-6, 60.0, xtol=1e-15, rtol=1e-15, maxiter=500)
    Tb = LAM1 * Sb + LAM2 + LAM3 * p
    m = RHO_SW * gS * (S - Sb) / Sb                  # kg m^-2 s^-1, rho_fw m
    return m, Tb, Sb


def test_constants_are_ad16_table4():
    assert (ml.LAMBDA1, ml.LAMBDA2, ml.LAMBDA3) == (LAM1, LAM2, LAM3)
    assert (ml.CW, ml.LATENT_HEAT, ml.RHO_SW) == (CW, LF, RHO_SW)
    assert ml.RHO_FW == 1000.0
    assert ml.CD_TOP == 2.5e-3 and ml.U_TIDAL == 0.01
    assert ml.GAMMA_S_OVER_T == 1.0 / 35.0


def test_hand_case_on_liquidus():
    m, Tb, Sb = ml.three_equation(T, S, P, USTAR, GT)
    assert abs(Tb - (LAM1 * Sb + LAM2 + LAM3 * P)) < 1e-10        # AD16 Eq. 25


def test_hand_case_matches_independent_solution():
    m, Tb, Sb = ml.three_equation(T, S, P, USTAR, GT)
    m0, Tb0, Sb0 = _hand_solution(T, S, P, USTAR, GT)
    assert Sb == pytest.approx(Sb0, rel=1e-10)
    assert Tb == pytest.approx(Tb0, rel=1e-10, abs=1e-12)
    assert m == pytest.approx(m0, rel=1e-9)


def test_hand_case_heat_and_salt_balance_close():
    m, Tb, Sb = ml.three_equation(T, S, P, USTAR, GT)
    gT, gS = GT * USTAR, GT / 35.0 * USTAR
    heat_ocean = RHO_SW * CW * gT * (T - Tb)                      # W m^-2, AD16 Eq. 24 (RHS)
    heat_melt = m * LF                                            # W m^-2, AD16 Eq. 24 (LHS)
    assert abs(heat_melt - heat_ocean) / abs(heat_ocean) < 1e-10
    salt_ocean = RHO_SW * gS * (S - Sb)                           # AD16 Eq. 26
    assert abs(m * Sb - salt_ocean) / abs(salt_ocean) < 1e-10


@pytest.mark.parametrize("dT", [0.5, 2.0, 3.0])
def test_melting_above_freezing_point(dT):
    Tf = LAM1 * S + LAM2 + LAM3 * P
    m, Tb, Sb = ml.three_equation(Tf + dT, S, P, USTAR, GT)
    assert m > 0 and Sb < S and Tb < Tf + dT


@pytest.mark.parametrize("dT", [-0.05, -0.2])
def test_freezing_below_freezing_point(dT):
    Tf = LAM1 * S + LAM2 + LAM3 * P
    m, Tb, Sb = ml.three_equation(Tf + dT, S, P, USTAR, GT)
    assert m < 0 and Sb > S


def test_at_freezing_point_no_melt():
    Tf = LAM1 * S + LAM2 + LAM3 * P
    m, Tb, Sb = ml.three_equation(Tf, S, P, USTAR, GT)
    assert abs(m) < 1e-15 and Sb == pytest.approx(S, rel=1e-12)


def test_vectorised_matches_scalar():
    Ts = np.array([0.0, 1.0, -1.0])
    Ss = np.array([34.4, 34.5, 34.6])
    ps = np.array([3e6, 5e6, 4e6])
    us = np.array([0.01, 0.02, 0.012])
    m, Tb, Sb = ml.three_equation(Ts, Ss, ps, us, GT)
    for k in range(3):
        mk, Tbk, Sbk = ml.three_equation(Ts[k], Ss[k], ps[k], us[k], GT)
        assert (m[k], Tb[k], Sb[k]) == pytest.approx((mk, Tbk, Sbk), rel=1e-14)


def test_ustar_ad16_eq27():
    """u*^2 = Cd (u_w^2 + u_tidal^2), AD16 Eq. 27 as printed on PDF page 15 (ruling R21):
    Cd multiplies both terms, so at rest u* = sqrt(Cd) u_tidal, not u_tidal. Anchors computed
    by hand with Cd = 2.5e-3, u_tidal = 0.01 m/s."""
    # at rest: sqrt(2.5e-3 * 1e-4) = sqrt(2.5e-7) = 5.0e-4 m/s
    assert ml.friction_velocity(0.0) == pytest.approx(5.0e-4, rel=1e-12)
    # |u_w| = 0.01 m/s: sqrt(2.5e-3 * (1e-4 + 1e-4)) = sqrt(5.0e-7) = 7.0711e-4 m/s
    assert ml.friction_velocity(0.01 ** 2) == pytest.approx(7.0711e-4, rel=1e-4)
    # |u_w| = 0.2 m/s: sqrt(2.5e-3 * (0.04 + 1e-4)) = sqrt(1.0025e-4) = 1.00125e-2 m/s
    assert ml.friction_velocity(0.2 ** 2) == pytest.approx(1.00125e-2, rel=1e-5)


def test_mitgcm_freshwater_flux_sign_and_units():
    # SHIfwFlx = -rho_fw m (kg m^-2 s^-1, negative = melting); m in m/yr of freshwater
    m, Tb, Sb = ml.three_equation(T, S, P, USTAR, GT)
    fw = ml.shifwflx(T, S, P, USTAR, GT)
    assert fw == pytest.approx(-m, rel=1e-15)
    assert ml.melt_m_per_yr(fw) == pytest.approx(m / 1000.0 * 365 * 86400, rel=1e-15)
    assert ml.melt_m_per_yr(-1e-3) == pytest.approx(31.536, rel=1e-15)


def test_constants_agree_with_fortran_override_and_cases_yaml():
    """The Python law and the Fortran copy use the same ISOMIP+ constants; cases.yaml's
    Cd and u_tidal agree with them (data.shelfice takes Cd from cases.yaml; u_tidal is a
    PARAMETER in the Fortran)."""
    src = FORTRAN.read_text()

    def par(name):
        mt = re.search(rf"{name}\s*=\s*([-+0-9.]+)\s*_d\s*([-+0-9]+)", src)
        assert mt, f"{name} not found in {FORTRAN.name}"
        return float(mt.group(1)) * 10.0 ** int(mt.group(2))

    assert par("shiUTidal") == pytest.approx(ml.U_TIDAL, rel=1e-15)
    assert par("shiLam1") == pytest.approx(ml.LAMBDA1, rel=1e-15)
    assert par("shiLam2") == pytest.approx(ml.LAMBDA2, rel=1e-15)
    # MITgcm pressure is in dbar (pLoc = shelficeMass g SItodBar): lam3 [degC/dbar] = 1e4 lam3 [degC/Pa]
    assert par("shiLam3dBar") == pytest.approx(ml.LAMBDA3 * 1e4, rel=1e-13)
    cfg = yaml.safe_load((PKG / "cases.yaml").read_text())["configs"]["ocean0_2km"]
    assert cfg["u_tidal"] == ml.U_TIDAL and cfg["Cd"] == ml.CD_TOP


def test_fortran_ustar_has_cd_multiplying_both_terms():
    """R21: the Fortran override computes u*^2 = Cd (|u|^2 + u_tidal^2) (AD16 Eq. 27, PDF
    page 15), i.e. the drag coefficient multiplies the sum, not only |u|^2. Checked only
    inside the "#ifdef SHI_ISOMIP_PLUS_MELT" branches (the #else branch is MITgcm's own
    Holland and Jenkins form, which must stay as it is)."""
    lines = FORTRAN.read_text().splitlines()
    branches, cur = [], None
    for ln in lines:
        if ln.startswith("#ifdef SHI_ISOMIP_PLUS_MELT"):
            cur = []
        elif cur is not None and (ln.startswith("#else /* SHI_ISOMIP_PLUS_MELT")
                                  or ln.startswith("#endif /* SHI_ISOMIP_PLUS_MELT")):
            branches.append(cur)
            cur = None
        elif cur is not None:
            cur.append(ln)
    assert len(branches) == 3                                # declarations, init, u* law

    def flat(block):
        code = "\n".join(ln for ln in block if ln[:1] not in ("C", "c", "!", "#"))
        return re.sub(r"\s+", "", re.sub(r"\n\s{5}&", " ", code))   # join continuations

    law = [flat(b) for b in branches if "ustarSq=" in flat(b)]
    assert len(law) == 1
    assert law[0].count("ustarSq=") == 1
    assert "ustarSq=shiCdragfld(I,J,bi,bj)*(velSq(I,J)+shiUTidal*shiUTidal)" in law[0]
    assert "*velSq(I,J)+shiUTidal*shiUTidal" not in law[0]
