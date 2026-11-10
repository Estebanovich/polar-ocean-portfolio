"""Physical melt parameterisations M1-M3 with one fitted coefficient each.

Constants (Favier et al. 2019, Table 2): rho_sw = 1028 kg m-3, c_p = 3974 J kg-1 K-1,
rho_i = 917 kg m-3, L_i = 3.34e5 J kg-1. c = rho_sw c_p / (rho_i L_i) [K-1].

Units of the coefficient. The literature prefactor gamma * c gives an ice-equivalent
melt in m s-1 (gamma in m s-1). The dataset melt is m a-1 of FRESHWATER
(exp2_isomip_cavity/notes/exp3_dataset.md). Because the coefficient is fitted by
least squares, the factors s->a and ice->freshwater (rho_i/rho_fw) are absorbed into
it: gamma / K are "effective" coefficients, with units m a-1 (M1: per unit c*TF, TF in K)
and m a-1 K-1 (M2, M3, times c^2 TF^2), so they equal the literature gamma times
(seconds per year) * rho_i/rho_fw, and are comparable in magnitude to it only up to
that documented factor. Signs follow R1: none uses TF+, refreezing (TF < 0) is allowed.

Fit: minimise sum over all training runs' shelf cells of (coef * basis - melt)^2
(uniform grid: area weighting is a plain pooled sum), coef = sum(basis melt)/sum(basis^2).
"""
from __future__ import annotations

import numpy as np

from exp3_melt_emulator.data import Run, cell_features

RHO_SW = 1028.0   # Favier 2019 Table 2
CP_O = 3974.0     # Favier 2019 Table 2
RHO_I = 917.0     # Favier 2019 Table 2
L_I = 3.34e5      # Favier 2019 Table 2
C = RHO_SW * CP_O / (RHO_I * L_I)


class _OneCoef:
    name = ""

    def __init__(self):
        self.coef: float | None = None

    def _basis(self, run: Run) -> np.ndarray:
        raise NotImplementedError

    def fit(self, runs: list[Run]) -> None:
        if not len(runs):
            raise ValueError("fit requires at least one run")
        num = den = 0.0
        for r in runs:
            b = self._basis(r)[r.shelf_mask]
            num += float(np.sum(b * r.melt[r.shelf_mask]))
            den += float(np.sum(b * b))
        if den == 0.0:
            raise ValueError("degenerate basis: all-zero thermal forcing")
        self.coef = num / den

    def predict(self, run: Run) -> np.ndarray:
        if self.coef is None:
            raise RuntimeError(f"{self.name}: predict called before fit")
        return np.where(run.shelf_mask, self.coef * self._basis(run), np.nan)


class LinearLocal(_OneCoef):
    """M1: m = gamma * c * TF (Favier 2019 Eq. 2, linear local). Signed."""
    name = "M1"

    def _basis(self, run):
        return C * cell_features(run)["TF"]


class QuadLocal(_OneCoef):
    """M2: m = K * c^2 * |TF| * TF (Burgard 2022 Eq. 14 with the constant slope and
    far-field factors folded into K; same form as Favier 2019 Eq. 4). Signed."""
    name = "M2"

    def _basis(self, run):
        TF = cell_features(run)["TF"]
        return C**2 * np.abs(TF) * TF


class QuadNonLocal(_OneCoef):
    """M3: m = K * c^2 * TF * |<TF>| (Jourdain 2020 Eq. 1 with delta T = 0;
    Favier 2019 Eq. 5, M+). <TF> is the shelf mean of the run (TF_mean)."""
    name = "M3"

    def _basis(self, run):
        f = cell_features(run)
        return C**2 * f["TF"] * np.abs(f["TF_mean"])
