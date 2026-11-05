"""Initial T/S profiles for Exp. 1 (two-layer and continuous Sermilik-like).

Convention: ``z`` is elevation in metres, negative below sea level (MITgcm
``rC``); depth = -z. Temperature in degC (potential), salinity in psu.

Values (verified against ``literatura/*.md``):

* Sciascia et al. (2013), Table 2, winter Sermilik survey (two-layer runs
  WIN/TEMP/...): polar water (PW) T1 = -1.5 degC, S1 = 32.9 psu; Atlantic water
  (AW) T2 = 4 degC, S2 = 34.6 psu (AW thickness H2 = 450 m in their 600 m
  domain). Sciascia's text (Sec. 2.4) describes Sermilik as "a cold and fresh
  layer of polar waters (PW) overlying warm, salty Atlantic waters (AW)".
* The task brief listed polar water ~ -1.5 degC / 33.0 and Atlantic water
  ~ 3.5 degC / 34.8 as [verificar]. Sciascia (2013) gives 32.9 / 4 / 34.6, which
  we adopt (differences are small: 0.1 psu, 0.5 K, 0.2 psu). Neither
  Straneo & Cenedese (2015) nor Jackson et al. (2014) tabulate layer values.
* Pycnocline depth: Jackson et al. (2018, Sec. 4) centre the Sermilik pycnocline
  at 160 m (Delta rho = 1.1 kg m^-3). The brief's 200 m / 75 m tanh centre and
  width are kept as defaults [verificar: not tabulated in any local paper]; they
  are keyword arguments so campaigns can change them.
* R11: in config (b) (300 m domain, H2 = 200 m) make_case centres the
  continuous profile at H - H2 = 100 m depth (same tanh shape, width and layer
  values), so its lower layer is 200 m thick like the two-layer case. (a) and
  (c) keep the 200 m centre.
"""
from __future__ import annotations

import numpy as np

# Sciascia et al. (2013) Table 2 (winter Sermilik layer averages)
T_PW, S_PW = -1.5, 32.9  # degC, psu
T_AW, S_AW = 4.0, 34.6  # degC, psu


def two_layer(z, z_int, T1, S1, T2, S2, dz_int):
    """Two-layer profile with a tanh interface.

    z : elevation [m], negative down. z_int : interface depth [m, positive].
    (T1, S1) above, (T2, S2) below. dz_int : tanh scale of the interface [m],
    profile = X1 + (X2 - X1) (1 + tanh((depth - z_int)/dz_int))/2. Returns
    (T(z), S(z)).
    """
    depth = -np.asarray(z, dtype=float)
    w = 0.5 * (1.0 + np.tanh((depth - z_int) / dz_int))
    return T1 + (T2 - T1) * w, S1 + (S2 - S1) * w


def continuous_sermilik(z, z_c=200.0, width=75.0, T_pw=T_PW, S_pw=S_PW,
                        T_aw=T_AW, S_aw=S_AW):
    """Polar water over Atlantic water with a tanh transition centred at z_c.

    ``width`` [m] is the full width between the tanh = -1 and +1 points, i.e.
    profile weight = (1 + tanh(2 (depth - z_c)/width))/2. Defaults: see module
    docstring (Sciascia 2013 layer values; 200 m / 75 m from the brief).
    """
    return two_layer(z, z_c, T_pw, S_pw, T_aw, S_aw, width / 2.0)
