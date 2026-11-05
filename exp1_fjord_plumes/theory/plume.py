r"""Steady buoyant plume theory against a vertical ice face (top-hat profiles).

Line plume (per unit width along the face; Jenkins 2011 Eqs. 1-4 with
sin(alpha) = 1 and X = z)::

    d(b u)/dz     = E u + m_w                                       [m s^-1]
    d(b u^2)/dz   = b g' - Cd u^2                                   [m^2 s^-2]
    d(b u T)/dz   = E u Ta + m_w Tb - gT (T - Tb)                   [degC m s^-1]
    d(b u S)/dz   = E u Sa + m_w Sb - gS (S - Sb)                   [psu m s^-1]

Half-conical plume (Slater et al. 2016 Eqs. 1a-1d; flat side, of width 2b,
against the ice)::

    d(pi/2 b^2 u)/dz   = pi E b u + 2 b m_w                         [m^2 s^-1]
    d(pi/2 b^2 u^2)/dz = pi/2 b^2 g' - 2 Cd b u^2                   [m^3 s^-2]
    d(pi/2 b^2 u T)/dz = pi E b u Ta + 2 b [m_w Tb - gT (T - Tb)]
    d(pi/2 b^2 u S)/dz = pi E b u Sa + 2 b [m_w Sb - gS (S - Sb)]

Full cone (free plume, no wall: no drag, no melt; Morton et al. 1956 Eqs. 1
in top-hat form)::

    d(pi b^2 u)/dz = 2 pi E b u,  d(pi b^2 u^2)/dz = pi b^2 g',  heat/salt alike.

Here E = alpha is the entrainment coefficient, g' = g (rho_a - rho)/rho0 =
g [sBeta (Sa - S) - tAlpha (Ta - T)] (Jenkins 2011 Eqs. 5, 18), m_w the
meltwater volume flux (m s^-1, Boussinesq; ``three_eq_melt`` returns the
ice-thickness rate m = m_w rho0/rho_i), Tb, Sb, gT, gS the interface
values and transfer velocities of ``melt.three_eq_melt`` (Si = 0), and z the
elevation (m, negative below sea level), integrated upward from the source.

Numerics: the state is (Q, M^2, Q T, Q S) with Q the volume flux and M the
momentum flux. Writing the momentum equation for M^2 (as Morton et al. 1956
do with V^4 = (b u)^4) removes the 1/M singularity of the cone equations as
u -> 0, so the top of the plume (u = 0) is a regular zero crossing of M^2,
located with a terminal ``solve_ivp`` event. ``scipy.integrate.solve_ivp``,
RK45, rtol = 1e-8.

Top-hat vs Gaussian: in flux form the top-hat system equals the Gaussian
system of Morton et al. (1956) when alpha_tophat = sqrt(2) alpha_gauss.
Morton et al.'s alpha = 0.093 is a Gaussian value.

Linear equation of state (Jenkins 2011 Eq. 5; MITgcm ``eosType='LINEAR'``)::

    rho = rho0 [1 - tAlpha (T - Tref) + sBeta (S - Sref)]

tAlpha = 3.87e-5 K^-1 and sBeta = 7.86e-4 psu^-1 are Jenkins (2011,
Table 1). They are NOT the values of ``verification/isomip``, which uses
``eosType='JMD95Z'`` and ``rhoConst=1030`` (``isomip/input/data:18-20``).
rho0 = 1027 kg m^-3 [verificar]: no local source; the MITgcm runs must set
rhoNil = rhoConst = 1027, tAlpha, sBeta to these values for consistency.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Literal, Union

import numpy as np
from scipy.integrate import solve_ivp

from .melt import MeltParams, Mode, three_eq_melt, transfer_velocities

Ambient = Union[float, np.ndarray, Callable[[float], float]]

_RTOL = 1e-8


@dataclass(frozen=True)
class LinearEOS:
    """Linear equation of state (Jenkins 2011 Eq. 5, Table 1)."""

    rho0: float = 1027.0  # kg m^-3 [verificar]
    tAlpha: float = 3.87e-5  # K^-1    Jenkins (2011) Table 1
    sBeta: float = 7.86e-4  # psu^-1  Jenkins (2011) Table 1
    Tref: float = 0.0  # degC (only differences matter for g')
    Sref: float = 0.0  # psu
    g: float = 9.81  # m s^-2  Jenkins (2011) Table 1

    def rho(self, T, S):
        """rho = rho0 [1 - tAlpha (T - Tref) + sBeta (S - Sref)]  [kg m^-3]."""
        return self.rho0 * (1.0 - self.tAlpha * (np.asarray(T) - self.Tref)
                            + self.sBeta * (np.asarray(S) - self.Sref))

    def gprime(self, T, S, Ta, Sa):
        """Reduced gravity g' = g [sBeta (Sa - S) - tAlpha (Ta - T)]  [m s^-2]
        (Jenkins 2011 Eq. 18 times g)."""
        return self.g * (self.sBeta * (Sa - S) - self.tAlpha * (Ta - T))


@dataclass
class PlumeSolution:
    """Plume profiles on the requested grid (NaN above ``z_detach``).

    z [m], b [m] (thickness for the line plume, radius for cones),
    u [m s^-1], T [degC], S [psu], rho [kg m^-3], melt [m s^-1 of ice,
    positive = melting; zero for the full cone], Q [m^2 s^-1 line, m^3 s^-1
    cones], gprime [m s^-2]. ``z_detach``: elevation where u -> 0 (the
    surface z_grid[-1] if the plume gets there, ``reached_top`` = True);
    ``z_neutral``: first level of zero buoyancy (NaN if none).
    """

    z: np.ndarray
    b: np.ndarray
    u: np.ndarray
    T: np.ndarray
    S: np.ndarray
    rho: np.ndarray
    melt: np.ndarray
    Q: np.ndarray
    gprime: np.ndarray
    z_detach: float
    z_neutral: float
    reached_top: bool
    geometry: str
    mode: str


# geometry -> (k in Q = k b^n u, n, entrainment perimeter / b (cones),
#              wall contact width / b (cones))
_GEOM = {
    "line": (1.0, 1, None, None),
    "half_cone": (math.pi / 2.0, 2, math.pi, 2.0),
    "full_cone": (math.pi, 2, 2.0 * math.pi, 0.0),
}


def _as_callable(f: Ambient, z_grid: np.ndarray) -> Callable[[float], float]:
    if callable(f):
        return lambda z: float(f(z))
    arr = np.asarray(f, dtype=float)
    if arr.ndim == 0:
        val = float(arr)
        return lambda z: val
    if arr.shape != z_grid.shape:
        raise ValueError("ambient profile must be scalar, callable or match z_grid")
    return lambda z: float(np.interp(z, z_grid, arr))


def _solve(z_grid, Ta, Sa, Q0, T0, S0, p: MeltParams, mode: Mode, alpha,
           geometry, eos: LinearEOS, u0) -> PlumeSolution:
    z_grid = np.asarray(z_grid, dtype=float)
    if z_grid.ndim != 1 or z_grid.size < 2 or np.any(np.diff(z_grid) <= 0):
        raise ValueError("z_grid must be 1-D and strictly increasing (base first)")
    if Q0 <= 0:
        raise ValueError("source volume flux must be positive")
    if geometry not in _GEOM:
        raise ValueError(f"unknown geometry {geometry!r}")
    # validates mode early
    transfer_velocities(0.0, p, mode)
    k, n, perim, wall = _GEOM[geometry]
    Ta_f, Sa_f = _as_callable(Ta, z_grid), _as_callable(Sa, z_grid)
    zb, zt = z_grid[0], z_grid[-1]
    g0 = eos.gprime(T0, S0, Ta_f(zb), Sa_f(zb))
    if g0 <= 0:
        raise ValueError(f"source is not buoyant (g'0 = {g0:.3g} m s^-2)")

    # ---- source conditions
    if u0 is None:
        if geometry == "line":
            # Balance of buoyancy, entrainment and drag: U^3 = g'0 q0/(E + Cd)
            # (Jenkins 2011 Eq. 21 and Table 3 with sin(alpha)=1; u_x0 = 1).
            u0 = (g0 * Q0 / (alpha + p.Cd)) ** (1.0 / 3.0)
        else:
            # Pure plume G0 = 5 b g'/(8 alpha u^2) = 1 with Q0 = k b0^2 u0
            # -> b0^5 = 8 alpha Q0^2 / (5 k^2 g'0)
            # (Slater et al. 2016 Eq. 5 for k = pi/2; same algebra for k = pi).
            b0 = (8.0 * alpha * Q0**2 / (5.0 * k**2 * g0)) ** 0.2
            u0 = Q0 / (k * b0**2)
    M0 = Q0 * u0
    y0 = np.array([Q0, M0**2, Q0 * T0, Q0 * S0])
    atol = 1e-12 * np.array([Q0, M0**2, Q0 * 10.0, Q0 * 40.0])
    has_wall = geometry != "full_cone"
    use_melt = has_wall and not p.melt_off
    drag = 0.0 if geometry == "full_cone" else p.Cd
    ratio_w = p.rho_i / p.rho0  # ice-thickness rate -> meltwater volume flux

    def kinematics(Q, M):
        u = M / Q
        if geometry == "line":
            b = Q / u if u > 0 else math.inf
            entr = alpha * u  # Jenkins (2011) Eq. 6
            width = 1.0
            dragM = drag * u * u * M  # M * Cd u^2
        else:
            Ms = max(M, 1e-300)
            b = Q / math.sqrt(k * Ms)
            entr = alpha * perim * math.sqrt(M / k)  # perim*alpha*b*u
            width = wall * b
            dragM = 2.0 * drag * b * u * u * M if geometry == "half_cone" else 0.0
        return u, b, entr, width, dragM

    def rhs(z, y):
        Q, M2, QT, QS = y
        M = math.sqrt(M2) if M2 > 0 else 0.0
        T, S = QT / Q, QS / Q
        ta, sa = Ta_f(z), Sa_f(z)
        u, b, entr, width, dragM = kinematics(Q, M)
        gp = eos.gprime(T, S, ta, sa)
        dQ = entr
        dQT = entr * ta
        dQS = entr * sa
        if use_melt:
            m, Tb, Sb = three_eq_melt(T, S, z, u, p, mode)
            gT, gS = transfer_velocities(u, p, mode)
            m_w = m * ratio_w
            dQ += width * m_w
            dQT += width * (m_w * Tb - gT * (T - Tb))
            dQS += width * (m_w * Sb - gS * (S - Sb))
        # M dM/dz = Q^2 g' - M * drag  ->  d(M^2)/dz = 2 (Q^2 g' - drag M)
        dM2 = 2.0 * (Q * Q * gp - dragM)
        return [dQ, dM2, dQT, dQS]

    def ev_detach(z, y):
        return y[1]

    ev_detach.terminal = True
    ev_detach.direction = -1

    def ev_neutral(z, y):
        Q, _, QT, QS = y
        return eos.gprime(QT / Q, QS / Q, Ta_f(z), Sa_f(z))

    ev_neutral.terminal = False
    ev_neutral.direction = -1

    sol = solve_ivp(rhs, (zb, zt), y0, method="RK45", rtol=_RTOL, atol=atol,
                    t_eval=z_grid, events=[ev_detach, ev_neutral])
    if sol.status == -1:
        raise RuntimeError(f"plume integration failed: {sol.message}")

    reached_top = sol.status == 0
    z_detach = float(zt) if reached_top else float(sol.t_events[0][0])
    z_neutral = float(sol.t_events[1][0]) if sol.t_events[1].size else math.nan

    N = z_grid.size
    out = {name: np.full(N, np.nan) for name in
           ("b", "u", "T", "S", "rho", "melt", "Q", "gprime")}
    for j, (z, y) in enumerate(zip(sol.t, sol.y.T)):
        Q, M2, QT, QS = y
        M = math.sqrt(M2) if M2 > 0 else 0.0
        T, S = QT / Q, QS / Q
        u, b, _, _, _ = kinematics(Q, M)
        out["b"][j], out["u"][j], out["T"][j], out["S"][j] = b, u, T, S
        out["Q"][j] = Q
        out["gprime"][j] = eos.gprime(T, S, Ta_f(z), Sa_f(z))
        out["melt"][j] = three_eq_melt(T, S, z, u, p, mode)[0] if use_melt else 0.0
    out["rho"] = np.where(np.isfinite(out["T"]), eos.rho(out["T"], out["S"]), np.nan)
    return PlumeSolution(z=z_grid, z_detach=z_detach, z_neutral=z_neutral,
                         reached_top=reached_top, geometry=geometry, mode=mode,
                         **out)


def line_plume(z_grid, Ta: Ambient, Sa: Ambient, q0: float, T0: float,
               S0: float, p: MeltParams, mode: Mode, alpha: float = 0.1, *,
               eos: LinearEOS = LinearEOS(), u0: float | None = None
               ) -> PlumeSolution:
    """Line plume against a vertical ice face (Jenkins 2011 Eqs. 1-4).

    q0 : source volume flux per unit width [m^2 s^-1]; T0 [degC], S0 [psu].
    Ta, Sa : ambient profiles (scalar, array on z_grid or callable of z).
    u0 : source velocity [m s^-1]; default is the buoyancy/entrainment/drag
    balance of Jenkins (2011, Eq. 21), b0 = q0/u0.
    """
    return _solve(z_grid, Ta, Sa, q0, T0, S0, p, mode, alpha, "line", eos, u0)


def point_plume(z_grid, Ta: Ambient, Sa: Ambient, Q0: float, T0: float,
                S0: float, p: MeltParams, mode: Mode, alpha: float = 0.1,
                geometry: Literal["half_cone", "full_cone"] = "half_cone", *,
                eos: LinearEOS = LinearEOS(), u0: float | None = None
                ) -> PlumeSolution:
    """Point-source plume: half cone against the ice (Slater et al. 2016
    Eqs. 1a-1d) or free full cone (Morton et al. 1956, no wall terms).

    Q0 : source volume flux [m^3 s^-1]. Default source radius/velocity:
    pure plume G0 = 1 (Slater et al. 2016 Eq. 5).
    """
    if geometry not in ("half_cone", "full_cone"):
        raise ValueError(f"unknown point-plume geometry {geometry!r}")
    return _solve(z_grid, Ta, Sa, Q0, T0, S0, p, mode, alpha, geometry, eos, u0)
