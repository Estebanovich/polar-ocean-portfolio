r"""ISOMIP+ three-equation melt law, as implemented in the Exp. 2 MITgcm configuration.

Formulation: Asay-Davis et al. (2016), hereafter AD16, sect. 3.1.8, Eq. 24-27, with the
Table 4 constants ("COM" melt parameterisation; Jenkins et al. 2010)::

    rho_fw m L     = -rho_sw cw u* GT (Tb - Tw)       heat     [W m^-2]        (AD16 Eq. 24)
    Tb             = lam1 Sb + lam2 + lam3 pb         liquidus [degC]          (AD16 Eq. 25)
    rho_fw m Sb    = -rho_sw u* GS (Sb - Sw)          salt     [PSU kg m^-2 s^-1] (Eq. 26)
    u*^2           = Cd (u_w^2 + u_tidal^2)           friction velocity [m^2 s^-2] (Eq. 27)

``rho_fw m`` is the meltwater mass flux [kg m^-2 s^-1], positive = melting; Tw, Sw are the
boundary-layer *potential* temperature and salinity (AD16's liquidus is the potential
freezing point, so no conversion to in situ temperature is made); the ice is perfectly
insulating (kappa_i = 0, AD16 Table 4). GS = GT / 35 (AD16 Table 4).

Eq. 27 is taken from the AD16 PDF (p. 15): Cd multiplies both terms, so the floor is
u* = sqrt(Cd) u_tidal = 5e-4 m s^-1 at rest (ruling R21; a text extraction had lost the
parentheses and the first implementation used Cd u_w^2 + u_tidal^2).

How the Exp. 2 MITgcm code computes it
--------------------------------------
``exp2_isomip_cavity/config/ocean0/code/shelfice_thermodynamics.F`` is a copy of
``pkg/shelfice/shelfice_thermodynamics.F`` with the changes marked ``ISOMIP+ (Exp. 2)``
(CPP option ``SHI_ISOMIP_PLUS_MELT``). It uses the BRIOS three-equation branch
(``shelfice_thermodynamics.F:560-683`` of the original) with SHELFICEconserve = T,
SHELFICEadvDiffHeatFlux = F (dFac = 0) and SHELFICEkappa = 0, which makes eps3 = eps8 = 0
(original lines 583-590). This module reproduces that branch line by line:

* eps1 = rhoConst cw gammaT, eps2 = rhoConst L gammaS                     (lines 570-573)
* eps4 = b pLoc + c0, eps6 = eps4 - tLoc                                   (lines 586-587)
* aqe Sb^2 + bqe Sb + cqe = 0, aqe = a0 eps1, bqe = eps1 eps6 - eps2,
  cqe = eps2 S                                                              (lines 591-603)
* Sb = (-bqe - sqrt(bqe^2 - 4 aqe cqe)) / (2 aqe), the other root if negative (636-638)
* Tb = a0 Sb + eps4                                                         (line 639)
* freshWaterFlux = rhoConst gammaS (1 - S / Sb)                             (lines 645-647)

with the ISOMIP+ changes: a0 = lam1, c0 = lam2, b = lam3 in degC/dbar (pLoc is in dbar,
``pLoc = shelficeMass g SItodBar``, original line 289); no SW_TEMP conversion (original
line 459); gammaT = GT u*, gammaS = GS u*, u*^2 = Cd (|u|^2 + u_tidal^2) (replacing the
Holland and Jenkins 1999 form of original lines 470-493).

Sign and units of the MITgcm diagnostic ``SHIfwFlx`` (``shelfice_thermodynamics.F:660``
and ``:701``; title "Ice shelf fresh water flux (positive upward)", units kg/m^2/s,
``pkg/shelfice/shelfice_init_fixed.F:250-252``): SHIfwFlx = freshWaterFlux
= rhoConst gammaS (Sb - S) / Sb = -rho_fw m, i.e. **negative = melting**. The protocol
melt rate (m/yr of freshwater, positive = melting) is m = -SHIfwFlx / rho_fw * 365 * 86400.

rho_sw: the code uses rUnit2mass = rhoConst = rho_ref = 1027.51 kg m^-3 for AD16's
rho_sw (1028); AD16 sect. 3.1.7 asks models that need rho_ref = rho_sw to use rho_ref.
"""
from __future__ import annotations

import numpy as np

# AD16 Table 4
LAMBDA1 = -0.0573          # liquidus slope, degC PSU^-1
LAMBDA2 = 0.0832           # liquidus intercept, degC
LAMBDA3 = -7.53e-8         # liquidus pressure coefficient, degC Pa^-1
CW = 3974.0                # seawater heat capacity, J degC^-1 kg^-1
LATENT_HEAT = 3.34e5       # latent heat of fusion, J kg^-1
RHO_FW = 1000.0            # freshwater density, kg m^-3
RHO_SW = 1027.51           # = rho_ref = MITgcm rhoConst (see module docstring), kg m^-3
CD_TOP = 2.5e-3            # top drag coefficient
U_TIDAL = 0.01             # RMS tidal velocity, m s^-1
GAMMA_S_OVER_T = 1.0 / 35.0
SECONDS_PER_YEAR = 365 * 86400.0      # 365-day calendar (context: diagnostics months)
PA_TO_DBAR = 1.0e-4                   # MITgcm SItodBar (model/src/ini_parms.F)


def friction_velocity(vel_sq, cd: float = CD_TOP, u_tidal: float = U_TIDAL):
    """u* = sqrt(Cd (|u|^2 + u_tidal^2)) [m s^-1] (AD16 Eq. 27, PDF p. 15); ``vel_sq`` =
    |u_w|^2 [m^2 s^-2]. At rest u* = sqrt(Cd) u_tidal = 5e-4 m s^-1."""
    return np.sqrt(cd * (np.asarray(vel_sq, dtype=float) + u_tidal ** 2))


def three_equation(T, S, p, ustar, gammaT, gammaS=None, *, rho_sw: float = RHO_SW,
                   cw: float = CW, latent_heat: float = LATENT_HEAT):
    """Solve AD16 Eq. 24-26 as the Exp. 2 shelfice_thermodynamics.F does.

    T, S: boundary-layer potential temperature [degC] and salinity [PSU]; p: pressure at
    the ice base [Pa]; ustar: friction velocity [m s^-1]; gammaT: nondimensional GT
    (gammaS defaults to GT / 35). Arrays broadcast.

    Returns (melt [kg m^-2 s^-1] = rho_fw m, positive = melting; Tb [degC]; Sb [PSU]).
    """
    T, S, p, ustar = (np.asarray(a, dtype=float) for a in (T, S, p, ustar))
    if gammaS is None:
        gammaS = gammaT * GAMMA_S_OVER_T
    a0, c0, b = LAMBDA1, LAMBDA2, LAMBDA3 / PA_TO_DBAR          # b in degC/dbar
    pLoc = p * PA_TO_DBAR                                         # dbar (line 289)
    gT = gammaT * ustar                                           # m s^-1, gamma_T = GT u*
    gS = gammaS * ustar
    eps1 = rho_sw * cw * gT                                       # lines 570-571
    eps2 = rho_sw * latent_heat * gS                              # lines 572-573
    eps4 = b * pLoc + c0                                          # line 586
    eps6 = eps4 - T                                               # line 587
    # kappa_i = 0 and dFac = 0: eps3 = eps3a = eps8 = 0 (lines 583-590)
    aqe = a0 * eps1                                               # line 591
    bqe = eps1 * eps6 - eps2                                      # lines 595-600
    cqe = eps2 * S                                                # lines 601-603
    disc = bqe * bqe - 4.0 * aqe * cqe                            # line 604
    recip_aqe = 0.5 / aqe                                         # line 593
    Sb = (-bqe - np.sqrt(disc)) * recip_aqe                       # line 636
    Sb = np.where(Sb < 0.0, (-bqe + np.sqrt(disc)) * recip_aqe, Sb)   # lines 637-638
    Tb = a0 * Sb + eps4                                           # line 639
    fw = rho_sw * gS * (1.0 - S / Sb)                             # lines 645-647 (SHIfwFlx)
    melt = -fw
    if melt.ndim == 0:
        return float(melt), float(Tb), float(Sb)
    return melt, Tb, Sb


def shifwflx(T, S, p, ustar, gammaT, **kw):
    """The MITgcm diagnostic SHIfwFlx [kg m^-2 s^-1] (negative = melting) for these inputs."""
    m = three_equation(T, S, p, ustar, gammaT, **kw)[0]
    return -m


def melt_m_per_yr(shifwflx_kg_m2_s):
    """Protocol melt rate [m/yr of freshwater, positive = melting] from SHIfwFlx."""
    return -np.asarray(shifwflx_kg_m2_s, dtype=float) / RHO_FW * SECONDS_PER_YEAR
