"""Evaluation metrics for Exp. 3 (numpy/scipy only).

All map metrics take a prediction on the full grid and the evaluated run; only that
run's own ``shelf_mask`` cells are used (predictions elsewhere are ignored, may be
NaN). Targets are never clipped: refreezing (negative melt) cells count as they are.
The grid is uniform (4 km), so area weighting is the plain mean over mask cells.
A NaN prediction inside the mask raises ValueError.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

from exp3_melt_emulator.data import Run, cell_features


def _masked(pred: np.ndarray, run: Run, what: str = "prediction") -> np.ndarray:
    """Values of ``pred`` on the run's shelf mask; raise if any is not finite."""
    v = np.asarray(pred, float)[run.shelf_mask]
    if not np.isfinite(v).all():
        raise ValueError(f"non-finite {what} inside run {run.case!r} shelf mask")
    return v


def rmse_map(pred: np.ndarray, run: Run) -> float:
    """Area-weighted RMSE over the shelf mask (m/a); Burgard 2022 Eq. 32 on a uniform grid."""
    e = _masked(pred, run) - run.melt[run.shelf_mask]
    return float(np.sqrt(np.mean(e**2)))


def integrated_rel_error(pred: np.ndarray, run: Run) -> float:
    """(sum pred - sum melt) / sum melt over the mask.

    Single-shelf, single-run relative form of the integrated melt of Burgard 2022
    Eq. 34 (M = sum m_j a_j) that enters RMSE_int, Eq. 33; the uniform cell area cancels.
    """
    p = _masked(pred, run).sum()
    m = run.melt[run.shelf_mask].sum()
    return float((p - m) / m)


def _ordered_cells(run: Run, key: np.ndarray, descending: bool = False) -> np.ndarray:
    """Flat indices of shelf cells ordered by ``key``; ties broken by flat index."""
    idx = np.flatnonzero(run.shelf_mask)
    k = key.ravel()[idx]
    order = np.argsort(-k if descending else k, kind="stable")   # stable -> ties by flat index
    return idx[order]


def _top_fraction(run: Run, key: np.ndarray, frac: float, descending: bool = False) -> np.ndarray:
    n = int(round(frac * run.shelf_mask.sum()))
    out = np.zeros(run.shelf_mask.shape, bool)
    out.flat[_ordered_cells(run, key, descending)[:n]] = True
    return out


def gl_region(run: Run, frac: float = 0.10) -> np.ndarray:
    """Boolean map of the ``frac`` of shelf cells closest to the grounding line.

    Plays the role of the first box of the five-box setup (about 10 % of the shelf
    area) in Burgard 2022 Sect. 3.2.2.
    """
    return _top_fraction(run, cell_features(run)["dist_gl"], frac)


def rmse_gl(pred: np.ndarray, run: Run, frac: float = 0.10) -> float:
    """Burgard 2022 Eq. 35 for one shelf and one simulation.

    Eq. 35 is the RMS over shelves k and simulations n of the difference between the
    parameterised and reference GL-region mean melt m_GL[k, n]. With N_isf = N_simu = 1
    the root of the single squared difference is its absolute value.
    """
    reg = gl_region(run, frac)
    p = np.asarray(pred, float)[reg]
    if not np.isfinite(p).all():
        raise ValueError(f"non-finite prediction inside GL region of run {run.case!r}")
    return float(abs(p.mean() - run.melt[reg].mean()))


def skill(err: float, err_ref: float) -> float:
    """1 - err / err_ref."""
    return float(1.0 - err / err_ref)


def coverage(mu: np.ndarray, sigma: np.ndarray, run: Run, level: float) -> float:
    """Fraction of mask cells with |melt - mu| <= z * sigma, z = Phi^-1((1 + level) / 2)."""
    z = norm.ppf((1.0 + level) / 2.0)
    m, s = _masked(mu, run, "mu"), _masked(sigma, run, "sigma")
    return float(np.mean(np.abs(run.melt[run.shelf_mask] - m) <= z * s))


def crps_gaussian(mu: np.ndarray, sigma: np.ndarray, run: Run) -> float:
    """Area-mean CRPS of N(mu, sigma^2) against the melt.

    Closed form, Gneiting & Raftery 2007 (Eq. 21): CRPS = sigma [z (2 Phi(z) - 1)
    + 2 phi(z) - 1/sqrt(pi)], z = (y - mu) / sigma.
    """
    m, s = _masked(mu, run, "mu"), _masked(sigma, run, "sigma")
    y = run.melt[run.shelf_mask]
    z = (y - m) / s
    crps = s * (z * (2.0 * norm.cdf(z) - 1.0) + 2.0 * norm.pdf(z) - 1.0 / np.sqrt(np.pi))
    return float(np.mean(crps))


def lambert_sensitivities(melt_by_tbot: dict[float, float]) -> dict:
    """Melt sensitivity and its nonlinearity, Lambert & Burgard 2025 Eq. 1-2.

    ``melt_by_tbot`` maps the temperature forcing T (REF, +1, +2 degC warming) to the
    region-mean melt M. Exactly three entries are required.
      Eq. 1: dM/dT |(REF,+1) = (M_+1 - M_REF) / (T_+1 - T_REF)
      Eq. 2: d2M/dT2 = [dM/dT|(+1,+2) - dM/dT|(REF,+1)] / (T_+1.5 - T_+0.5)
    where T_+1.5 and T_+0.5 are the mid-point temperatures of the (+1,+2) and
    (REF,+1) intervals.
    """
    if len(melt_by_tbot) != 3:
        raise ValueError("need exactly three temperatures (REF, +1, +2)")
    t0, t1, t2 = sorted(melt_by_tbot)
    m0, m1, m2 = (melt_by_tbot[t] for t in (t0, t1, t2))
    s01 = (m1 - m0) / (t1 - t0)
    s12 = (m2 - m1) / (t2 - t1)
    t_hi, t_lo = 0.5 * (t1 + t2), 0.5 * (t0 + t1)
    return {"dMdT_ref_p1": s01, "dMdT_p1_p2": s12, "d2MdT2": (s12 - s01) / (t_hi - t_lo)}


def deep_amplification(pred_a, pred_b, run_a: Run, run_b: Run, frac: float = 0.10) -> float:
    """Deep amplification DA, Lambert & Burgard 2025 Eq. 3.

    DA = (dM10/dT)|(REF,+1) / (dM/dT)|(REF,+1), M10 being the mean melt over the deepest
    ``frac`` of the draft. a = REF run, b = +1 degC run; the common temperature
    difference T_+1 - T_REF cancels in the ratio.
    """
    def means(pred, run):
        p = _masked(pred, run)
        deep = _top_fraction(run, np.abs(run.draft), frac, descending=True)
        q = np.asarray(pred, float)[deep]
        return p.mean(), q.mean()

    ma, m10a = means(pred_a, run_a)
    mb, m10b = means(pred_b, run_b)
    return float((m10b - m10a) / (mb - ma))
