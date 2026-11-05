r"""Three-equation ice-ocean melt parameterisation.

Formulation (Holland & Jenkins 1999, as written in Jenkins 2011 Eqs. 7-9,
Holland et al. 2008 Eqs. 2-4 and Slater et al. 2016 Eqs. 2a-2c)::

    c_w gT (T - Tb) = m_w [L + c_i (Tb - Ti)]      heat      [W m^-2 / rho0]
    gS (S - Sb)     = m_w Sb                        salt      [psu m s^-1]
    Tb = lam1 Sb + lam2 + lam3 z                    liquidus  [degC]

with ``z`` the face elevation (m, negative downward), ``gT``/``gS`` the heat
and salt transfer velocities (m s^-1), ``m_w`` the meltwater volume flux per
unit area (m s^-1, Boussinesq, positive = melting) and ice salinity Si = 0.
Eliminating ``m_w`` and ``Tb`` gives a quadratic for ``Sb``; the physical
root is the positive one. The function returns the melt rate as ice
thickness per time, ``m = m_w rho0 / rho_i`` (m s^-1 of ice).

How MITgcm computes melt at a vertical ice face
-----------------------------------------------
(``/Volumes/ESTEBAN_MIT/MITgcm``, ``pkg/icefront``; read-only checkout)

* Transfer is velocity dependent with a floor, not constant:
  ``pkg/icefront/icefront_thermodynamics.F:91-97``::

      gammaT = max(1.0e-2 * |wVEL(k)| * sqrt(1.5e-3), 1.0e-4)   [m s^-1]
      gammaS = 5.05e-3 * gammaT

  i.e. Gamma_T = 1e-2, Cd = 1.5e-3, Gamma_S/Gamma_T = 5.05e-3 (so
  sqrt(Cd) Gamma_T = 3.87e-4, sqrt(Cd) Gamma_S = 1.96e-6), floor 1e-4 m/s.
  These are hard-wired, not namelist parameters. ``wVEL(I,J,K)`` is the
  vertical velocity at the top face of cell K in the same column. The floor
  wins for |w| < 0.258 m s^-1. NB: these differ from the Jenkins et al.
  (2010) / Jenkins (2011) Stanton numbers (1.1e-3, 3.1e-5) that Xu et al.
  (2012, Table 1) report for the same module.
* Liquidus: ``icefront_thermodynamics.F:77-79``: a0 = -0.0575 degC/psu,
  c0 = 0.0901 degC, b = -7.61e-4 degC/dbar, applied as
  ``Tb = a0 Sb + c0 + b pLoc`` with ``pLoc = |rC(k)|`` (dbar taken as m,
  line 109; half-cell mean at the ice base, line 119). Here lam3 = -b and
  ``z = -pLoc``.
* Ocean temperature is converted to in situ with ``SW_TEMP`` (line 114);
  this module takes T as given (potential = in situ is assumed).
* Ice temperature varies linearly from ``ICEFRONTthetaSurface`` at the
  surface to 0 degC at the ice-front base ``R_icefront``
  (``icefront_thermodynamics.F:121-123``); default thetaSurface = -20 degC
  (``icefront_readparms.F:88``).
* L = 334.0e3 J/kg and c_i = 2000 J/kg/K (``icefront_readparms.F:84-85``,
  defaults when SHELFICE is off), rho_ice = 917 kg/m^3 (line 86). c_w is the
  model ``HeatCapacity_Cp`` (default 3994, ``model/src/set_defaults.F:174``;
  ``verification/isomip/input/data:19`` sets 3974). The ocean density in
  the fluxes is ``rUnit2mass = rhoConst`` (``model/src/ini_parms.F:1574``).
* Freshwater-flux diagnostic ``ICFfwFlx`` (``icefront_thermodynamics.F:
  143-146, 188``)::

      FreshWaterFlux = rhoConst * c_w * gammaT * (Tb - T) / (L + c_i (Tb - Ti))

  in kg m^-2 s^-1, "positive increases ocean salinity" (line 50): **negative
  = melting**. Hence ``m_w = -ICFfwFlx / rhoConst`` and the ice-thickness
  melt rate returned here is ``m = -ICFfwFlx / rho_i``.

Modes
-----
``"mitgcm"``  transfer law above (exactly as the Fortran).
``"jenkins"`` gT = sqrt(Cd) GammaT |u|, gS = sqrt(Cd) GammaS |u| with
              Cd = 2.5e-3, sqrt(Cd) GammaT = 1.1e-3, sqrt(Cd) GammaS = 3.1e-5
              (Jenkins 2011, Table 1).
``"constant"`` gT = ``gammaT``, gS = ``gammaS`` (SHELFICE-style constant
              coefficients; ``pkg/shelfice/shelfice_readparms.F:118``).

The default ``MeltParams`` follows MITgcm icefront (for comparison runs);
``MeltParams.jenkins2011()`` returns the Jenkins (2011) Table 1 constants.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Literal

Mode = Literal["mitgcm", "jenkins", "constant"]
_MODES = ("mitgcm", "jenkins", "constant")


@dataclass(frozen=True)
class MeltParams:
    """Constants for the three-equation melt law. SI units throughout.

    Defaults reproduce ``pkg/icefront`` (see module docstring for
    ``file:line``) except ``c_w``, which must match the run's
    ``HeatCapacity_Cp`` (3974 here, as in ``verification/isomip``).
    """

    # --- "constant" mode transfer velocities [m s^-1]
    gammaT: float = 1.0e-4  # icefront_thermodynamics.F:95 floor; shelfice_readparms.F:118
    gammaS: float = 5.05e-7  # 5.05e-3 * gammaT (icefront_thermodynamics.F:96)
    # --- "mitgcm" mode (icefront_thermodynamics.F:91-97), hard-wired there
    mitgcm_GammaT: float = 1.0e-2  # [-]
    mitgcm_Cd: float = 1.5e-3  # [-]
    mitgcm_gammaS_over_gammaT: float = 5.05e-3  # [-]
    mitgcm_gammaT_min: float = 1.0e-4  # [m s^-1]
    # --- "jenkins" mode; Jenkins (2011) Table 1 lists the products
    #     sqrt(Cd) GammaT = 1.1e-3 and sqrt(Cd) GammaS = 3.1e-5.
    Cd: float = 2.5e-3  # drag coefficient [-]; also used for plume wall drag
    GammaT: float = 1.1e-3 / math.sqrt(2.5e-3)  # = 0.022 [-]
    GammaS: float = 3.1e-5 / math.sqrt(2.5e-3)  # = 6.2e-4 [-]
    # --- liquidus Tb = lam1 Sb + lam2 + lam3 z  (z < 0 below sea level)
    lam1: float = -0.0575  # degC psu^-1   icefront_thermodynamics.F:77
    lam2: float = 0.0901  # degC          icefront_thermodynamics.F:78
    lam3: float = 7.61e-4  # degC m^-1     icefront_thermodynamics.F:79 (b = -lam3)
    # --- thermodynamic constants
    L: float = 3.34e5  # J kg^-1        icefront_readparms.F:84
    c_w: float = 3974.0  # J kg^-1 K^-1  run's HeatCapacity_Cp (isomip/input/data:19)
    c_i: float = 2000.0  # J kg^-1 K^-1  icefront_readparms.F:85 (brief: 2009, Jenkins 2011)
    rho0: float = 1027.0  # kg m^-3  rhoConst of our runs [verificar]
    rho_i: float = 917.0  # kg m^-3  icefront_readparms.F:86
    # --- ice temperature: linear from T_ice_surface at z = 0 to 0 degC at
    #     z = -ice_depth (icefront_thermodynamics.F:121-123); uniform
    #     T_ice_surface if ice_depth is None.
    T_ice_surface: float = -20.0  # degC  icefront_readparms.F:88
    ice_depth: float | None = None  # m (positive)
    # --- Ruling R2: switch melting (and its heat/salt/volume fluxes) off
    melt_off: bool = False

    @classmethod
    def jenkins2011(cls, **overrides) -> "MeltParams":
        """Constants of Jenkins (2011), Table 1 (ice temperature not given
        there: T_ice_surface kept at -20 degC [verificar])."""
        base = cls(lam1=-5.73e-2, lam2=8.32e-2, lam3=7.61e-4, L=3.35e5,
                   c_w=3974.0, c_i=2009.0)
        return replace(base, **overrides)


def freezing_point(S: float, z: float, p: MeltParams) -> float:
    """Linearised liquidus Tf = lam1 S + lam2 + lam3 z [degC].

    Jenkins (2011) Eq. 9; z in m, negative downward.
    """
    return p.lam1 * S + p.lam2 + p.lam3 * z


def ice_temperature(z: float, p: MeltParams) -> float:
    """Ice temperature at the face [degC] (icefront_thermodynamics.F:121-123)."""
    if p.ice_depth is None:
        return p.T_ice_surface
    depth = -z
    return p.T_ice_surface * (p.ice_depth - depth) / p.ice_depth


def transfer_velocities(u: float, p: MeltParams, mode: Mode) -> tuple[float, float]:
    """Heat and salt transfer velocities (gT, gS) [m s^-1] for speed ``u``."""
    if mode == "mitgcm":
        # icefront_thermodynamics.F:91-97
        gT = max(p.mitgcm_GammaT * abs(u) * math.sqrt(p.mitgcm_Cd), p.mitgcm_gammaT_min)
        return gT, p.mitgcm_gammaS_over_gammaT * gT
    if mode == "jenkins":
        # Jenkins (2011) Eqs. 7-8: Stanton numbers sqrt(Cd) Gamma_{T,S} times U
        sCd = math.sqrt(p.Cd)
        return sCd * p.GammaT * abs(u), sCd * p.GammaS * abs(u)
    if mode == "constant":
        return p.gammaT, p.gammaS
    raise ValueError(f"unknown melt mode {mode!r}; expected one of {_MODES}")


def three_eq_melt(T: float, S: float, z: float, u: float, p: MeltParams,
                  mode: Mode) -> tuple[float, float, float]:
    """Solve the three-equation melt problem at one point of the ice face.

    Parameters
    ----------
    T, S : ocean (plume) temperature [degC] and salinity [psu]
    z : elevation of the face point [m], negative downward
    u : speed past the ice [m s^-1] (along-face velocity)
    p : constants; ``p.melt_off`` returns ``(0.0, nan, nan)``
    mode : ``"mitgcm"``, ``"jenkins"`` or ``"constant"`` (see module doc)

    Returns
    -------
    m : melt rate [m s^-1 of ice], positive = melting; m = m_w rho0 / rho_i
    Tb, Sb : interface temperature [degC] and salinity [psu]

    Quadratic for Sb (heat x salt balance, Holland & Jenkins 1999), written as
    A Sb^2 + B Sb + C = 0 with K = lam2 + lam3 z::

        A = lam1 (c_i gS - c_w gT)
        B = c_w gT (T - K) + gS [L + c_i (K - Ti)] - c_i lam1 gS S
        C = -gS S [L + c_i (K - Ti)]

    (identical to aqe, bqe, cqe of icefront_thermodynamics.F:105-132 divided
    by rho0). A > 0 and C < 0, so exactly one root is positive; it is taken
    in the cancellation-free form Sb = 2C / (-B - sqrt(B^2 - 4AC)).
    """
    if mode not in _MODES:
        raise ValueError(f"unknown melt mode {mode!r}; expected one of {_MODES}")
    if p.melt_off:
        return 0.0, math.nan, math.nan
    S = max(S, 0.0)  # icefront_thermodynamics.F:111
    gT, gS = transfer_velocities(u, p, mode)
    Ti = ice_temperature(z, p)
    K = p.lam2 + p.lam3 * z
    latent = p.L + p.c_i * (K - Ti)
    A = p.lam1 * (p.c_i * gS - p.c_w * gT)
    B = p.c_w * gT * (T - K) + gS * latent - p.c_i * p.lam1 * gS * S
    C = -gS * S * latent
    disc = B * B - 4.0 * A * C
    denom = -B - math.sqrt(disc)
    Sb = 2.0 * C / denom if denom != 0.0 else 0.0
    Tb = p.lam1 * Sb + K
    # heat balance, Jenkins (2011) Eq. 7 -> meltwater flux m_w [m s^-1]
    m_w = p.c_w * gT * (T - Tb) / (p.L + p.c_i * (Tb - Ti))
    return m_w * p.rho0 / p.rho_i, Tb, Sb

