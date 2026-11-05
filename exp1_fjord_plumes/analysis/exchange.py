r"""Fjord-shelf exchange diagnostics on a cross-fjord section (config (c)).

Exchange flux::

    Q_ex = 1/2 \int |u| dA                                         [m^3 s^-1]

over the wet part of a cross-section (u normal to it). Jackson et al. (2018)
measure the exchange with the inflowing/outflowing volume flux Q_in = -Q_out
(their Sec. 3c(3), and Q_out in the turnover time, Eq. 26). Since
\int |u| dA = Q_in - Q_out, Q_ex = (Q_in - Q_out)/2, which equals Q_in = -Q_out
when the net transport through the section is zero; with the plume's net
down-fjord transport Q0 it is the mean of the two directional fluxes.
[verificar resolved: the 1/2 \int |u| dA form is not written in Jackson et
al. 2018 Sec. 2 (Methods); it is our compact form of their Q_in/Q_out.]

Heat toward the front (ice face at x = 0, so u < 0 flows toward it)::

    H_front = rho0 c_p \int (-u) (T - T_ref) dA                     [W]

T_ref is a required argument: with a non-zero net volume transport the
value depends on it (e.g. the freezing point or the outflow temperature).
rho0 = 1027 kg m^-3 and c_p = 3974 J kg^-1 K^-1 match the runs' rhoConst and
HeatCapacity_Cp (config (a) ``data``).
"""
from __future__ import annotations

import numpy as np


def exchange_flux(u_section, dA) -> float:
    """Q_ex = 1/2 sum |u| dA [m^3 s^-1]; NaN velocities (land) are ignored.

    u_section [m s^-1] and dA [m^2] (cell face areas, hFac included) have the
    same shape.
    """
    u = np.asarray(u_section, dtype=float)
    dA = np.asarray(dA, dtype=float)
    if u.shape != dA.shape:
        raise ValueError(f"u {u.shape} and dA {dA.shape} differ")
    return float(0.5 * np.nansum(np.abs(u) * dA))


def heat_toward_front(u_section, T_section, dA, T_ref, rho0=1027.0,
                      cp=3974.0) -> float:
    """rho0 c_p sum (-u)(T - T_ref) dA [W]; > 0 means heat toward the front.

    The front is at x = 0 (``front_i = 0``), so u < 0 points toward it.
    NaN (land) entries are ignored.
    """
    u = np.asarray(u_section, dtype=float)
    T = np.asarray(T_section, dtype=float)
    dA = np.asarray(dA, dtype=float)
    if not (u.shape == T.shape == dA.shape):
        raise ValueError("u, T and dA must have the same shape")
    return float(rho0 * cp * np.nansum(-u * (T - T_ref) * dA))
