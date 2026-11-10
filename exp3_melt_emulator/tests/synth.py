"""Small synthetic runs (ny=5, nx=12) with known masks and profiles for tests.

Layout per run: columns 0-2 grounded (dry), 3..(3+w-1) shelf, remaining columns open
ocean (draft 0); last row (y index 4) is the dry side-wall land row.
"""
from __future__ import annotations

import numpy as np

from exp3_melt_emulator.data import Run

NY, NX, DX = 5, 12, 4000.0


def _run(i, rng, width, negative, n):
    z = -(np.arange(8) * 100.0 + 50.0)                  # -50 ... -750
    Tbot, z_c, d, s = -0.5 + 0.4 * i, 300.0 + 20.0 * i, 200.0, 1.0
    r = np.clip((-z - (z_c - d / 2)) / d, 0, 1)
    T = -1.9 + (Tbot + 1.9) * r
    S = 34.0 + 0.7 * r
    wet = np.zeros((NY, NX), bool)
    wet[:4, 3:] = True
    sh = np.zeros((NY, NX), bool)
    sh[:4, 3:3 + width] = True
    draft = np.full((NY, NX), np.nan)
    draft[wet] = 0.0
    col = np.arange(NX)[None, :] - 3
    row = np.arange(NY)[:, None]
    draft[sh] = (-(100.0 + 60.0 * col) * s - 5.0 * row)[sh] * np.ones((NY, NX))[sh]
    melt = np.full((NY, NX), np.nan)
    melt[wet] = 0.0
    melt[sh] = 2.0 + 0.01 * np.abs(draft[sh]) + 0.1 * rng.standard_normal(sh.sum())
    if negative:
        melt[tuple(np.argwhere(sh)[:3].T)] = -0.5
    return Run(case=f"synth_{i}", group="factorial", split="test" if i == n - 1 else "train",
               Tbot=Tbot, z_c=z_c, d=d, s=s, n1_pass=(i != 2),
               melt=melt, draft=draft, shelf_mask=sh, T_profile=T, S_profile=S, z=z,
               x=320000.0 + DX * np.arange(NX), y=2000.0 + DX * np.arange(NY))


def synthetic_runs(n=6, seed=0) -> list[Run]:
    """n runs; run 1 has a smaller shelf mask (width 4 vs 5) and run 3 has negative melt cells."""
    rng = np.random.default_rng(seed)
    return [_run(i, rng, 4 if i == 1 else 5, negative=(i == 3), n=n) for i in range(n)]


def write_npz(runs, path):
    """Write runs in the Exp. 2 export layout (the arrays load_runs reads)."""
    st = lambda k: np.stack([getattr(r, k) for r in runs])
    np.savez(path,
             case=np.array([r.case for r in runs]), group=np.array([r.group for r in runs]),
             split=np.array([r.split for r in runs]),
             Tbot=st("Tbot"), z_c=st("z_c"), d=st("d"), s=st("s"),
             N1_pass=np.array([r.n1_pass for r in runs]),
             z=runs[0].z, x=runs[0].x, y=runs[0].y,
             T_profile=st("T_profile"), S_profile=st("S_profile"),
             melt=st("melt"), draft=st("draft"), shelf_mask=st("shelf_mask"))
