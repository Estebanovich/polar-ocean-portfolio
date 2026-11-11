"""Exp. 3 dataset loader and per-cell features.

The Exp. 2 export (exp2_isomip_cavity/notes/exp3_dataset.md) is read with
``allow_pickle=False`` after its sha256 is checked. Depths are negative downward
(``z``, ``draft``); ``z_c`` and ``d`` are positive depth and thickness (m).
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt

from exp2_isomip_cavity.gendata.profiles import GRAVITY
from exp2_isomip_cavity.theory.melt import (
    LAMBDA1, LAMBDA2, LAMBDA3, PA_TO_DBAR, RHO_SW,
)

DATASET_PATH = Path("/Volumes/LaCie/portfolio_runs/exp2/exp3_dataset.npz")
DATASET_SHA256 = "3d2eceb18be2d21c785b175e5a569dd2edbadd0b5437d9aba16ac71ae8d6200f"

FEATURES = ("draft", "T", "S", "TF", "TF_mean", "draft_minus_zc",
            "dist_gl", "dist_front", "slope", "s")


@dataclass(frozen=True)
class Run:
    case: str
    group: str
    split: str
    Tbot: float
    z_c: float
    d: float
    s: float
    n1_pass: bool
    melt: np.ndarray          # (ny, nx) m/a, NaN on dry cells
    draft: np.ndarray         # (ny, nx) m (negative), NaN on dry cells
    shelf_mask: np.ndarray    # (ny, nx) bool
    T_profile: np.ndarray     # (nz,)
    S_profile: np.ndarray     # (nz,)
    z: np.ndarray             # (nz,) level centres, negative downward
    x: np.ndarray             # (nx,)
    y: np.ndarray             # (ny,)


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_runs(path=DATASET_PATH, sha256=DATASET_SHA256) -> list[Run]:
    """Load all runs; ValueError if the file's sha256 differs from ``sha256``."""
    got = _sha256(path)
    if got != sha256:
        raise ValueError(f"sha256 mismatch for {path}: {got} != {sha256}")
    with np.load(path, allow_pickle=False) as d:
        a = {k: d[k] for k in d.files}
    return [
        Run(case=str(a["case"][i]), group=str(a["group"][i]), split=str(a["split"][i]),
            Tbot=float(a["Tbot"][i]), z_c=float(a["z_c"][i]), d=float(a["d"][i]),
            s=float(a["s"][i]), n1_pass=bool(a["N1_pass"][i]),
            melt=a["melt"][i], draft=a["draft"][i], shelf_mask=a["shelf_mask"][i],
            T_profile=a["T_profile"][i], S_profile=a["S_profile"][i],
            z=a["z"], x=a["x"], y=a["y"])
        for i in range(len(a["case"]))
    ]


def union_mask(runs) -> np.ndarray:
    """Cells that belong to the shelf of at least one run."""
    return np.logical_or.reduce([r.shelf_mask for r in runs])


def pressure_dbar(z):
    """Pressure (dbar) at depth z (m, negative down): RHO_SW g |z| SItodBar.

    Same rule as MITgcm shelfice (pLoc = mass g SItodBar) as used in theory/melt.py.
    """
    return RHO_SW * GRAVITY * np.abs(np.asarray(z, dtype=float)) * PA_TO_DBAR


def freezing_temperature(S, z):
    """AD16 Table 4 liquidus T_f = LAMBDA1 S + LAMBDA2 + b p, with p in dbar.

    LAMBDA3 is in degC/Pa; as in theory/melt.py (``b = LAMBDA3 / PA_TO_DBAR``) it is
    converted to degC/dbar before multiplying the pressure in dbar.
    """
    return LAMBDA1 * np.asarray(S) + LAMBDA2 + (LAMBDA3 / PA_TO_DBAR) * pressure_dbar(z)


def _centred_gradient(f, h, axis):
    """d f/dh with NaN-aware differences: centred where both neighbours are finite,
    else one-sided from the finite neighbour, else 0."""
    f = np.moveaxis(f, axis, 0)
    pad = np.full((1,) + f.shape[1:], np.nan)
    up = np.concatenate([f[1:], pad])      # f[i+1]
    dn = np.concatenate([pad, f[:-1]])     # f[i-1]
    fu, fd = np.isfinite(up), np.isfinite(dn)
    out = np.zeros_like(f)
    both = fu & fd
    out[both] = (up[both] - dn[both]) / (2 * h)
    m = fu & ~fd
    out[m] = (up[m] - f[m]) / h
    m = fd & ~fu
    out[m] = (f[m] - dn[m]) / h
    return np.moveaxis(out, 0, axis)


def cell_features(run: Run) -> dict[str, np.ndarray]:
    """Per-cell predictors, shape (ny, nx), NaN outside ``run.shelf_mask``.

    - T, S: profile linearly interpolated at the draft depth (clamped at the ends).
    - TF = T - T_f(S, p(draft)); TF_mean is the shelf mean of TF, broadcast.
    - draft_minus_zc = (-draft) - z_c: depth of the ice base below the thermocline
      centre (m, positive = deeper).
    - dist_gl: Euclidean distance (m) to the nearest grounded cell, via
      scipy.ndimage.distance_transform_edt(sampling=(dy, dx)). Grounded = dry cell
      (non-finite draft) other than the side-wall land row, which is the last row
      (ny-1) of the 4 km grid and is not grounded ice.
    - dist_front: the same distance to the nearest open-ocean wet cell
      (finite draft, not in shelf_mask).
    - slope = |grad draft| from centred differences over shelf cells only (one-sided
      where a neighbour is not a shelf cell, so the open-ocean draft of 0 and dry
      cells do not contaminate it).
    """
    sh = run.shelf_mask
    dx = float(run.x[1] - run.x[0])
    dy = float(run.y[1] - run.y[0])
    draft = np.where(sh, run.draft, np.nan)

    zi = run.z[::-1]                       # np.interp needs increasing abscissa
    T = np.interp(draft, zi, run.T_profile[::-1])
    S = np.interp(draft, zi, run.S_profile[::-1])
    TF = T - freezing_temperature(S, draft)
    TF_mean = np.full(sh.shape, np.nanmean(TF[sh]))

    wet = np.isfinite(run.draft)
    wall = np.zeros_like(sh)
    wall[-1, :] = True
    grounded = ~wet & ~wall
    openocean = wet & ~sh
    dist_gl = distance_transform_edt(~grounded, sampling=(dy, dx))
    dist_front = distance_transform_edt(~openocean, sampling=(dy, dx))

    slope = np.hypot(_centred_gradient(draft, dy, 0), _centred_gradient(draft, dx, 1))

    feats = {
        "draft": draft, "T": T, "S": S, "TF": TF, "TF_mean": TF_mean,
        "draft_minus_zc": -draft - run.z_c,
        "dist_gl": dist_gl, "dist_front": dist_front, "slope": slope,
        "s": np.full(sh.shape, run.s),
    }
    return {k: np.where(sh, v, np.nan) for k, v in feats.items()}
