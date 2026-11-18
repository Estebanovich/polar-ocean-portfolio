"""Coarse (g50) -> fine (f10) initial tracer fields (plan Task 9).

Land filling and bilinear regridding reuse exp5_adjoint_slope.gendata.interp (_fill_land,
_regrid; same Grid attributes xc, yc, nx, dx_m). Then each tracer is clipped at >= 0 on wet
fine cells and rescaled so that its total inventory sum(C V) equals the coarse one."""
from __future__ import annotations

import numpy as np

from exp5_adjoint_slope.gendata.interp import _fill_land, _regrid


def tracers_to_fine(tr_c, gc, gf, vol_c, vol_f) -> np.ndarray:
    tr_c, vol_c, vol_f = (np.asarray(a, float) for a in (tr_c, vol_c, vol_f))
    wet_c, wet_f = vol_c > 0, vol_f > 0
    out = np.zeros((tr_c.shape[0], gf.nr, gf.ny, gf.nx))
    for k in range(wet_c.shape[0]):          # _fill_land fills along y only: it would never end
        cols = ~wet_c[k].any(axis=0)
        if cols.any() and wet_c[k].any():
            raise ValueError(f"level {k}: x-columns {np.nonzero(cols)[0].tolist()} have no wet cell in y")
    for i in range(tr_c.shape[0]):
        f = _regrid(_fill_land(tr_c[i], wet_c), gc, gf)
        f = np.where(wet_f, np.maximum(f, 0.0), 0.0)
        inv_c, inv_f = (tr_c[i] * vol_c).sum(), (f * vol_f).sum()
        out[i] = f * (inv_c / inv_f) if inv_f > 0 else f
    return out
