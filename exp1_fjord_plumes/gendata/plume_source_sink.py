r"""Subglacial plume as MITgcm ``addMass`` sinks and sources (config (c)).

At 200 m horizontal resolution the hydrostatic fjord model cannot resolve the
plume, so it is parameterised with buoyant-plume theory, following Cowton et
al. (2015, Sec. 2.1.3, and their Fig. 1): "Water, heat, and salt are then
removed from MITgcm cells at the rate predicted by the plume model in those
vertical layers where the plume is entraining, and put into the cell at the
depth at which the plume is predicted to terminate", the plume "terminating
upon reaching neutral buoyancy or the fjord surface".

Per-level volume flux (``q``, m^3 s^-1, MITgcm order: k = 0 at the surface)::

    entraining cells (below k_n):  q_k = -[Q(z_top,k) - Q(z_bot,k)]   (= -dQ/dz dz)
    outflow cell k_n:              q_k = Q(z_lo)
    cells above k_n:               q_k = 0

with Q(z) the plume volume flux (``theory.plume.point_plume``, half cone,
Slater et al. 2016 Eqs. 1a-1d), z_e = z_neutral (first level of zero
buoyancy) or the surface if the plume is still buoyant there, k_n the cell
holding z_e and z_lo its bottom face. The sum telescopes,
``sum(q) = Q(z_lo) - [Q(z_lo) - Q(z_base)] = Q0`` (criterion C5), exactly up to
round-off.

Choices:

* **Outflow at the neutral level, one cell (as Cowton et al. 2015).** Above
  z_neutral the plume rises on its momentum only (their Fig. 3, dashed) and
  is not counted as entraining. In the (c) Sermilik profile with
  Q0 = 100 m^3 s^-1 the theory puts z_neutral at ~181 m depth but the
  momentum top (z_detach) at ~3 m (melt off), i.e. 18 cells higher;
  spreading the outflow over that overshoot would inject AW-derived plume
  water into the polar-water layer up to the surface, where it is denser
  than the ambient and would sink back. We therefore follow Cowton et al.
  and put it at the neutral level (z_detach is not used).
* **Cells are not split.** A cell either entrains or receives the outflow;
  the entrainment between z_lo and z_neutral (< one dz) is attributed to the
  outflow, whose T, S are those of the plume at z_lo. Every cell is then a
  pure sink (removes local water) or a pure source (adds plume water),
  which is what PD3 needs (``notes/pd3_addmass_sink.md``): a net cell value
  would hide the gross outflow that sets the outflow heat/salt tendency.
* **Melt.** Q(z) also contains the meltwater flux of the half cone in contact
  with the ice (``mode="mitgcm"``, same transfer law as configs (a)/(b)), so
  the sinks remove ``entrainment + melt`` from the ambient. The melt volume
  is a small fraction of the entrained volume (Task 7 report); it is kept
  so that ``sum(q) = Q0`` holds exactly.
* **Outflow levels (``n_out_levels``, default 1 = Cowton).** Optionally the
  outflow is split equally between the z_neutral cell and the
  ``n_out_levels - 1`` cells above it, i.e. into the overshoot, where the
  intruding plume water lies (z_neutral < z < z_detach). Config (c) uses
  2 for numerical reasons (outflow-cell turnover time, see
  ``notes/pd1_c_fjord3d.md``); the sinks and the sum are unchanged.
* **Integration top.** ``point_plume`` ends the half cone with melt at
  u = u_stop (ruling R15, ``theory/plume.py``), so one integration from the
  base to the surface suffices.

Outflow T, S: the plume's flux-weighted values at z_lo (top-hat profiles,
T = QT/Q, S = QS/Q), returned in the outflow cells and NaN elsewhere.
"""
from __future__ import annotations

import numpy as np

from exp1_fjord_plumes.theory.melt import MeltParams
from exp1_fjord_plumes.theory.plume import point_plume


def plume_source_sink(z, Ta, Sa, Q0, dz, *, T0=0.0, S0=0.0, alpha=0.1,
                      geometry="half_cone", mode="mitgcm",
                      p: MeltParams | None = None, n_out_levels=1, n_fine=2001,
                      return_solution=False):
    """Per-level plume sinks/sources for ``addMass`` (see module docstring).

    z : cell-centre elevations [m], negative, surface first (MITgcm ``rC``).
    Ta, Sa : ambient T [degC], S [psu] on ``z``. Q0 : discharge [m^3 s^-1]
    entering at the base of the deepest cell, with T0 [degC], S0 [psu].
    dz : cell thickness [m], scalar or per level. p : melt constants,
    default ``MeltParams(ice_depth=H)`` with H the water depth (as (a)/(b)).
    n_out_levels : number of outflow cells (z_neutral cell and those above).

    Returns ``(q, T_out, S_out)``: q [m^3 s^-1] per level (< 0 entrainment,
    > 0 outflow, sum = Q0); T_out [degC], S_out [psu] of the outflow in the
    outflow cells, NaN elsewhere. With ``return_solution`` the
    ``PlumeSolution`` on the fine grid is appended.
    """
    z = np.asarray(z, dtype=float)
    dz = np.broadcast_to(np.asarray(dz, dtype=float), z.shape).copy()
    Ta = np.broadcast_to(np.asarray(Ta, dtype=float), z.shape)
    Sa = np.broadcast_to(np.asarray(Sa, dtype=float), z.shape)
    if z.ndim != 1 or z.size < 2 or np.any(np.diff(z) >= 0):
        raise ValueError("z must be 1-D cell centres ordered surface first (decreasing)")
    if np.any(dz <= 0):
        raise ValueError("dz must be positive")
    top, bot = z + dz / 2.0, z - dz / 2.0
    if not np.allclose(bot[:-1], top[1:]):
        raise ValueError("z and dz are inconsistent (cells must be contiguous)")
    if Q0 <= 0:
        raise ValueError("Q0 must be positive")
    if int(n_out_levels) != n_out_levels or n_out_levels < 1:
        raise ValueError("n_out_levels must be a positive integer")
    z_base, z_surf = float(bot[-1]), float(top[0])
    if p is None:
        p = MeltParams(ice_depth=-z_base)

    zi, Tai, Sai = z[::-1], Ta[::-1], Sa[::-1]           # increasing for interp
    zg = np.unique(np.concatenate([np.linspace(z_base, z_surf, n_fine), top, bot]))
    sol = point_plume(zg, np.interp(zg, zi, Tai), np.interp(zg, zi, Sai), Q0, T0,
                      S0, p, mode, alpha=alpha, geometry=geometry)
    ok = np.isfinite(sol.Q)

    def at(field, zz):
        return float(np.interp(zz, sol.z[ok], field[ok]))

    z_e = float(sol.z_neutral) if np.isfinite(sol.z_neutral) else min(float(sol.z_detach), z_surf)
    k_n = int(np.flatnonzero((bot < z_e) & (z_e <= top))[0])
    z_lo = float(bot[k_n])

    q = np.zeros(z.size)
    for k in range(k_n + 1, z.size):          # entrainment sinks, whole cells
        q[k] = -(at(sol.Q, top[k]) - at(sol.Q, bot[k]))
    ks = np.arange(max(k_n - n_out_levels + 1, 0), k_n + 1)   # outflow cells
    q[ks] = at(sol.Q, z_lo) / ks.size         # outflow at the neutral level
    T_out = np.full(z.size, np.nan)
    S_out = np.full(z.size, np.nan)
    T_out[ks], S_out[ks] = at(sol.T, z_lo), at(sol.S, z_lo)
    if return_solution:
        return q, T_out, S_out, sol
    return q, T_out, S_out
