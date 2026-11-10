"""Temperature / salinity profiles, linear equation of state and ice-shelf pressure load.

References: Asay-Davis et al. (2016), hereafter AD16. z is height [m], negative downwards.
Temperature is potential temperature [degC]; salinity is practical salinity [PSU].
"""
from __future__ import annotations

import numpy as np

ZB_DEEP = -720.0                     # AD16 Table 1, m
# AD16 Table 6 (WARM, Ocean0 initial and restoring state)
T0, S0 = -1.9, 33.8                  # surface, degC / PSU
TBOT_WARM, SBOT_WARM = 1.0, 34.7     # at 720 m
# AD16 Table 5 (COLD, Ocean1 initial state)
TBOT_COLD, SBOT_COLD = -1.9, 34.55
# AD16 Table 4 (linear EOS, Eq. 23)
RHO_REF = 1027.51                    # kg m^-3
T_REF, S_REF = -1.0, 34.2            # degC, PSU
ALPHA = 3.733e-5                     # degC^-1
BETA = 7.843e-4                      # PSU^-1
GRAVITY = 9.81                       # m s^-2, AD16 Table 1; MITgcm default (model/src/set_defaults.F:101)


def _linear(z, Tbot, Sbot, T0_=T0, S0_=S0):
    """AD16 Eq. 21-22. Eq. 22 as printed has (Sbot - T0); the intended (Sbot - S0) is used."""
    r = np.asarray(z, dtype=float) / ZB_DEEP
    return T0_ + (Tbot - T0_) * r, S0_ + (Sbot - S0_) * r


def ocean0(z_m):
    """WARM profile (AD16 Eq. 21-22, Table 6): (T [degC], S [PSU]) at heights z_m [m]."""
    return _linear(z_m, TBOT_WARM, SBOT_WARM)


def cold(z_m):
    """COLD profile (AD16 Eq. 21-22, Table 5): T = -1.9 degC, S 33.8 -> 34.55 PSU."""
    return _linear(z_m, TBOT_COLD, SBOT_COLD)


def ramp(z_m, Tbot, z_c, d=300.0, T0=T0, S0=S0, Sbot=SBOT_WARM):
    """Ensemble "ramp" profile (declared deviation from AD16, which has no thermocline):

    T = T0 + (Tbot - T0) r,  S = S0 + (Sbot - S0) r,  r = clip((-z - (z_c - d/2)) / d, 0, 1)

    with thermocline centre depth z_c [m, positive] and thickness d [m]. With z_c = 360 m,
    d = 720 m, r = z / zb,deep and the profile is AD16 Eq. 21-22 (Ocean0)."""
    r = np.clip((-np.asarray(z_m, dtype=float) - (z_c - d / 2.0)) / d, 0.0, 1.0)
    return T0 + (Tbot - T0) * r, S0 + (Sbot - S0) * r


def rho_linear(T, S):
    """AD16 Eq. 23: rho = rho_ref [1 - alpha (T - Tref) + beta (S - Sref)], kg m^-3."""
    return RHO_REF * (1.0 - ALPHA * (np.asarray(T) - T_REF) + BETA * (np.asarray(S) - S_REF))


def load_anomaly(T, S, dz, rho_const: float = RHO_REF, g: float = GRAVITY) -> np.ndarray:
    """Hydrostatic pressure anomaly p'(z_f) = g * integral_{z_f}^{0} (rho - rho_const) dz [Pa]
    at the nz + 1 level faces (p'[0] = 0 at the surface), for cell-centre profiles T, S (nz)
    and thicknesses dz (nz).

    The density is constant within each cell, as in MITgcm's hydrostatic integration for
    the initial load (verification/isomip/input/gendata.m:69-78). This is the quantity that
    pkg/shelfice adds to the surface pressure, shelficeLoadAnomaly = g (M_ice + rhoConst
    Ro_surf) (pkg/shelfice/shelfice_init_fixed.F:226-227), for an ice shelf floating in this
    profile. For profiles linear in z (Ocean0, COLD) the cell value is the exact cell mean."""
    rho = rho_linear(T, S)
    dz = np.asarray(dz, dtype=float)
    return np.concatenate([[0.0], np.cumsum(g * (rho - rho_const) * dz)])


def load_anomaly_at_draft(p_faces: np.ndarray, dz: np.ndarray, draft: np.ndarray) -> np.ndarray:
    """p' [Pa] at the ice draft (2-D, m, <= 0): linear in depth inside each cell, which
    is exact for the constant cell density of ``load_anomaly``."""
    depth_f = np.concatenate([[0.0], np.cumsum(np.asarray(dz, dtype=float))])
    return np.interp(-np.asarray(draft, dtype=float), depth_f, p_faces)
