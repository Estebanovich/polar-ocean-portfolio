"""Vertical DIC gradients over the analysis region R, after Jarnikova et al. (2025, Methods,
eqs. 1 and 4): the gradient over the top 200 m and the deep-minus-surface difference at
1000 m. Surface = model level 1. mol m-3 -> umol kg-1 with rho0 = rhoConst = 1035 kg m-3."""
from __future__ import annotations

import numpy as np

RHO0 = 1035.0


def _profile(field, z, region_mask) -> np.ndarray:
    f = np.asarray(field, float)
    m = np.asarray(region_mask, bool)
    return np.array([f[k][m].mean() for k in range(f.shape[0])])


def _umolkg(x):
    return x * 1e6 / RHO0


def g200(dic, z, region_mask) -> float:
    """[DIC(200 m) - DIC(level 1)] / 200 m, umol kg-1 m-1 (Jarnikova et al. 2025, eq. 1)."""
    p = _profile(dic, z, region_mask)
    return float(_umolkg(np.interp(200.0, z, p) - p[0]) / 200.0)


def ddic1000(dic, z, region_mask) -> float:
    """DIC(1000 m) - DIC(level 1), umol kg-1 (Jarnikova et al. 2025, eq. 4)."""
    p = _profile(dic, z, region_mask)
    return float(_umolkg(np.interp(1000.0, z, p) - p[0]))
