"""Air-sea CO2 flux and wind sensitivity (spec Sec. 4-5).

Sign: positive into the ocean throughout. DICCFLX is FluxCO2 of pkg/dic (mol m-2 s-1),
Kwexch (ff pCO2atm - k0 fugf pCO2ocean)/permil, positive when the atmosphere is richer
(pkg/dic/dic_surfforcing.F:261-272; filled at pkg/dic/dic_biotic_forcing.F:413)."""
from __future__ import annotations

import numpy as np

YEAR_S = 31536000.0


def flux_into_ocean(diccflx) -> np.ndarray:
    """DICCFLX (mol m-2 s-1) -> mol C m-2 yr-1, positive into the ocean."""
    return np.asarray(diccflx, float) * YEAR_S


def channel_mean(f2d, area, wet) -> float:
    w = np.asarray(area, float) * np.asarray(wet, bool)
    return float((np.asarray(f2d, float) * w).sum() / w.sum())


def sensitivity(f_branch, f_parent, window=(6, 15)) -> float:
    """Mean of branch - parent over branch years window[0]..window[1] (1-based, inclusive)."""
    b, p = np.asarray(f_branch, float), np.asarray(f_parent, float)
    a, e = window
    return float((b[a - 1:e] - p[a - 1:e]).mean())
