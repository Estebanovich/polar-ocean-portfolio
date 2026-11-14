"""Wind-stress pulse of the Exp. 4 fjord runs (Task 4).

The temporal shape is the thesis pulse (MSc thesis, input/wind_patch.py,
``make_meridional_wind_patch``): a half Gaussian of the surface stress,

    tau(t) = tau_max exp(-(t - t_c)^2 / (2 sigma^2))  for t <= t_c,   0 for t > t_c,

with tau_max = 0.1 N m^-2, sigma = 1 h and t_c = 2.5 h (spec section 3.4). In the fjord the
stress is spatially uniform and points along ``angle_deg``, measured in degrees from the UTM
grid east, counter-clockwise (the direction the stress points to, oceanographic convention):
180 deg is a stress towards the west (Review Focus 2).

Time sampling. MITgcm periodic external forcing (``periodicExternalForcing``) treats record n
(1-based) of the file as valid at (n - 1/2) * externForcingPeriod and interpolates linearly in
time (eesupp/src/get_periodic_interval.F:118-128; model/src/external_fields_load.F:69-72).
``record_times`` gives those instants and the pulse is evaluated there, so the forcing the
model applies is the nominal pulse sampled exactly (the thesis sampled linspace(0, 60 h, 180)
instead, 1206.7 s apart; notes/pd1_puente.md, Task 3). Before t = period/2 and after
(n - 1/2) period MITgcm interpolates between the last record and record 1 (the cycle wraps);
``cycle_records`` adds one zero record so the forced stage ends before the wrap.
"""
from __future__ import annotations

import numpy as np


def pulse(t_s, tau_max: float = 0.1, sigma_h: float = 1.0, center_h: float = 2.5) -> np.ndarray:
    """Half-Gaussian stress magnitude [N m^-2] at times ``t_s`` [s] (thesis shape)."""
    t_h = np.asarray(t_s, dtype=float) / 3600.0
    g = tau_max * np.exp(-(t_h - center_h) ** 2 / (2 * sigma_h ** 2))
    g[t_h > center_h] = 0.0
    return g


def record_times(n: int, period_s: float) -> np.ndarray:
    """Model times [s] of the n records of a periodic forcing file, (k + 1/2) * period."""
    return (np.arange(int(n)) + 0.5) * float(period_s)


def cycle_records(forced_h: float, period_s: float) -> int:
    """Number of records for a forced stage of ``forced_h`` hours: forced/period + 1, so the
    stage end (forced_h) is before (n - 1/2) * period, where the periodic interpolation would
    start to bring back record 1 (eesupp/src/get_periodic_interval.F:118-123)."""
    n = float(forced_h) * 3600.0 / float(period_s)
    if abs(n - round(n)) > 1e-9:
        raise ValueError("forced_h is not a whole number of forcing periods")
    return int(round(n)) + 1


def _cos_sin(angle_deg: float) -> tuple[float, float]:
    """cos and sin of an angle in degrees, exact (0, +-1) at multiples of 90 deg."""
    a = float(angle_deg) % 360.0
    exact = {0.0: (1.0, 0.0), 90.0: (0.0, 1.0), 180.0: (-1.0, 0.0), 270.0: (0.0, -1.0)}
    if a in exact:
        return exact[a]
    r = np.radians(a)
    return float(np.cos(r)), float(np.sin(r))


def uniform_stress(t_s, angle_deg: float, shape: tuple, **pulse_kw) -> tuple[np.ndarray, np.ndarray]:
    """Spatially uniform stress (taux, tauy), each (nt, *shape) [N m^-2], of magnitude
    ``pulse(t_s, **pulse_kw)`` pointing ``angle_deg`` from UTM east, counter-clockwise.
    MITgcm reads taux at U points and tauy at V points (model/src/external_fields_load.F:114-125);
    a uniform field needs no staggering."""
    g = pulse(t_s, **pulse_kw)
    c, s = _cos_sin(angle_deg)
    ones = np.ones(tuple(shape))
    return np.multiply.outer(g * c, ones), np.multiply.outer(g * s, ones)
