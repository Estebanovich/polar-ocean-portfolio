"""Coarse (g50) -> fine (f10) initial fields (plan Task 9).

Land filling and bilinear regridding reuse exp5_adjoint_slope.gendata.interp (_fill_land,
_regrid; same Grid attributes xc, yc, nx, dx_m). Each tracer is then clipped at >= 0 on wet
fine cells and rescaled so that its total inventory sum(C V) equals the coarse one. Theta is
regridded the same way, then convectively adjusted on wet cells (stabilize_wet; no rescaling; salinity is uniform, so a convex horizontal combination
keeps every column statically stable) and gets the Exp. 5 member noise (<= 1e-3 degC, spec
Sec. 4) in wet cells."""
from __future__ import annotations

import numpy as np

from exp5_adjoint_slope.gendata.interp import _fill_land, _regrid, perturb_member

MEMBER_AMP, MEMBER_SEED = 1e-3, 2026          # as Exp. 5 (spec Sec. 4)


def wet_volume(g, bathy, hfac_min: float = 0.1, hfac_min_dr: float = 5.0) -> np.ndarray:
    """Cell volumes (nr, ny, nx) with MITgcm partial cells (model/src/ini_masks_etc.F): the wet
    fraction of each level, rounded below hFacMnSz = max(hFacMin, min(hFacMinDr/drF, 1)) to 0 or
    hFacMnSz. Uniform dx = dy (Cartesian channel); bathy negative down, 0 = land."""
    depth = -np.asarray(bathy, float)
    top, dr = g.rf[:-1, None, None], g.dr[:, None, None]
    h = np.clip((np.minimum(depth[None], top + dr) - top) / dr, 0.0, 1.0)
    mn = np.maximum(hfac_min, np.minimum(hfac_min_dr / dr, 1.0)) * np.ones_like(h)
    small = h < mn
    h = np.where(small, np.where(h < 0.5 * mn, 0.0, mn), h)
    return h * dr * g.dx_m ** 2


def _check_fillable(wet_c) -> None:
    for k in range(wet_c.shape[0]):          # _fill_land fills along y only: it would never end
        cols = ~wet_c[k].any(axis=0)
        if cols.any() and wet_c[k].any():
            raise ValueError(f"level {k}: x-columns {np.nonzero(cols)[0].tolist()} have no wet cell in y")


def tracers_to_fine(tr_c, gc, gf, vol_c, vol_f) -> np.ndarray:
    tr_c, vol_c, vol_f = (np.asarray(a, float) for a in (tr_c, vol_c, vol_f))
    wet_c, wet_f = vol_c > 0, vol_f > 0
    _check_fillable(wet_c)
    out = np.zeros((tr_c.shape[0], gf.nr, gf.ny, gf.nx))
    for i in range(tr_c.shape[0]):
        f = _regrid(_fill_land(tr_c[i], wet_c), gc, gf)
        f = np.where(wet_f, np.maximum(f, 0.0), 0.0)
        inv_c, inv_f = (tr_c[i] * vol_c).sum(), (f * vol_f).sum()
        out[i] = f * (inv_c / inv_f) if inv_f > 0 else f
    return out


def _adjust_column(T: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Merge adjacent blocks until T is non-increasing downward, each block at its w-weighted
    mean (pool-adjacent-violators): the convective adjustment of a column, conserving sum(T w)."""
    blocks: list[list[float]] = []                          # [mean, weight, n]
    for x, wi in zip(T, w):
        blocks.append([float(x), float(wi), 1])
        while len(blocks) > 1 and blocks[-1][0] > blocks[-2][0]:
            m2, w2, n2 = blocks.pop()
            m1, w1, n1 = blocks.pop()
            blocks.append([(m1 * w1 + m2 * w2) / (w1 + w2), w1 + w2, n1 + n2])
    return np.concatenate([np.full(n, m) for m, _, n in blocks])


def stabilize_wet(theta, vol) -> np.ndarray:
    """Convective adjustment of theta on the wet cells of each column (the linear EOS has
    sBeta = 0, so a column is stable iff theta does not increase downward). Heat (theta x volume)
    is conserved; land cells are left as they are."""
    T, v = np.array(theta, float), np.asarray(vol, float)
    pair = (v[1:] > 0) & (v[:-1] > 0)
    for j, i in zip(*np.nonzero((pair & (T[1:] > T[:-1])).any(axis=0))):
        k = np.nonzero(v[:, j, i] > 0)[0]
        if k.size and (np.diff(k) == 1).all():
            T[k, j, i] = _adjust_column(T[k, j, i], v[k, j, i])
        else:
            raise ValueError(f"column ({j}, {i}): wet cells not contiguous in depth")
    return T


def fine_init(theta_c, tr_c, gc, gf, vol_c, vol_f, member: int | None = 0) -> dict:
    """make_case `init` for a fine run: {'theta': (nr, ny_f, nx_f), 1..6: tracers}. member None
    leaves theta unperturbed."""
    wet_c, wet_f = np.asarray(vol_c) > 0, np.asarray(vol_f) > 0
    _check_fillable(wet_c)
    th = stabilize_wet(_regrid(_fill_land(np.asarray(theta_c, float), wet_c), gc, gf), vol_f)
    if member is not None:
        th = np.where(wet_f, perturb_member(th, int(member), MEMBER_AMP, MEMBER_SEED), th)
    out: dict = {"theta": th}
    out.update({i + 1: a for i, a in enumerate(tracers_to_fine(tr_c, gc, gf, vol_c, vol_f))})
    return out
