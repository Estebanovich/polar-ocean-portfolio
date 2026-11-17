"""Adjoint sensitivities from GA/GB (plan Task 5): adxx_* files written by pkg/ctrl, the tau_x
sensitivity per monthly record, lead times and scaling by a typical variability.

adxx_<name> holds dJ/du for control u (J in TW, as in cost_accumulate_mean.F); for xx_fu the
records are the gentim2d periods (xx_gentim2d_period, 30 d), (nrec, ny, nx)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from pfcommon.mitgcm_io import read_mds

J_WINDOW_DAYS = 30.0          # J is the mean over the last 30 days (PD0 amendment)


def load_adxx(run_dir: Path, name: str) -> np.ndarray:
    """dJ/du from adxx_<name>.<it>.data/.meta in run_dir (the first iteration found)."""
    hits = sorted(Path(run_dir).glob(f"adxx_{name}.[0-9]*.data"))
    if not hits:
        raise FileNotFoundError(f"no adxx_{name}.*.data in {run_dir}")
    it = int(hits[0].name.split(".")[1])
    a, _ = read_mds(Path(run_dir) / f"adxx_{name}", it)
    return np.asarray(a, dtype=float)


def tau_x_sensitivity(adxx_fu, grid) -> np.ndarray:
    """Sum over x of dJ/dtau_x per record, (nrec, ny): TW per (N m-2) applied on a whole row."""
    a = np.asarray(adxx_fu, dtype=float)
    if a.ndim == 2:
        a = a[None]
    return a.sum(axis=-1)


def lead_time_days(m: int, window_days: float, period_days: float = 30.0) -> float:
    """Lead time of record m (1-based) of a window_days run: centre of the J averaging window
    (the last J_WINDOW_DAYS) minus the centre of the record [(m-1)P, mP]."""
    return (window_days - 0.5 * J_WINDOW_DAYS) - (m - 0.5) * period_days


def scaled(sens, sigma: float) -> np.ndarray:
    """Sensitivity times a typical variability of the control (Mazloff 2012): dJ per sigma."""
    return np.asarray(sens, dtype=float) * sigma
