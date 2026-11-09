"""ISOMIP+ cavity geometry: bed and ice draft on the 2 km / 4 km ocean grid.

Model geometry (ruling R7): both the bed and the ice draft come from the official ISOMIP+
geometry file ``Ocean1_input_geom_v1.01.nc`` (Cornford and Asay-Davis 2016,
doi:10.5880/PIK.2016.002; Asay-Davis et al. 2016, hereafter AD16, sect. 3.1.1), a 1 km
grid, block-averaged n x n with n = dx_km (2 x 2 for 2 km, 4 x 4 for 4 km). Bed and draft
use the same averaging, so the grounding line stays consistent.

``bathymetry`` is the analytic MISMIP+ bed (AD16 Eq. 1-3, Table 1). It is kept as a
documented reference: the file's ``bedrockTopography`` is the BISICLES-interpolated
version of it, flat at zb,deep = -720 m beyond x = 640 km.

Coordinates are in metres: x (along flow, "north") 320-800 km, y 0-80 km (AD16 sect. 3.1.1).
Depths are negative downwards.
"""
from __future__ import annotations

from collections import deque

import numpy as np

# AD16 Table 1 (MISMIP+ bedrock), m
B0, B2, B4, B6 = -150.0, -728.8, 343.91, -50.57
XBAR = 300e3            # characteristic along-flow length scale, m
FC = 4e3                # side-wall width, m
DC = 500.0              # trough depth relative to the side walls, m
WC = 24e3               # trough half-width, m
LY = 80e3               # domain width, m
ZB_DEEP = -720.0        # maximum bed depth, m
# AD16 sect. 3.1.1: ocean domain
X_MIN, X_MAX, Y_MIN, Y_MAX = 320e3, 800e3, 0.0, 80e3
# AD16 sect. 3.1.2: ice thinner than Hcalve has calved (draft set to 0)
H_CALVE = 100.0         # m
COORD_TOL = 1.0         # m, tolerance on the domain edges of the file


def bathymetry(x_m, y_m):
    """Analytic MISMIP+ bed zb(x, y) [m], AD16 Eq. 1-3 with Table 1 values.

    zb = max(Bx(x) + By(y), zb,deep)                                       (AD16 Eq. 1)
    Bx = B0 + B2 xt^2 + B4 xt^4 + B6 xt^6, xt = x / xbar                    (AD16 Eq. 2)
    By = dc / (1 + exp(-2 (y - Ly/2 - wc) / fc))
       + dc / (1 + exp( 2 (y - Ly/2 + wc) / fc))                            (AD16 Eq. 3)
    """
    x = np.asarray(x_m, dtype=float)
    y = np.asarray(y_m, dtype=float)
    xt = x / XBAR
    bx = B0 + B2 * xt**2 + B4 * xt**4 + B6 * xt**6
    by = (DC / (1.0 + np.exp(-2.0 * (y - LY / 2 - WC) / FC))
          + DC / (1.0 + np.exp(2.0 * (y - LY / 2 + WC) / FC)))
    return np.maximum(bx + by, ZB_DEEP)


# ---- reading the 1 km file --------------------------------------------------------------

def _check_axis(v: np.ndarray, lo: float, hi: float, name: str, path: str) -> float:
    if v.ndim != 1 or v.size < 2:
        raise ValueError(f"{path}: coordinate {name} is not a 1-D axis")
    d = np.diff(v)
    if not np.allclose(d, d[0], rtol=0, atol=1e-6 * abs(d[0])) or d[0] <= 0:
        raise ValueError(f"{path}: coordinate {name} is not uniform and increasing")
    lo_f, hi_f = v[0] - d[0] / 2, v[-1] + d[0] / 2
    if abs(lo_f - lo) > COORD_TOL or abs(hi_f - hi) > COORD_TOL:
        raise ValueError(
            f"{path}: {name} cells span {lo_f:g}-{hi_f:g}; expected metres covering "
            f"x 320-800 km, y 0-80 km ({name}: {lo:g}-{hi:g} m)")
    return float(d[0])


def _read_field(path: str, var: str):
    """(x, y, field) of the 1 km file, with dimensions and coordinates checked
    (Review Focus 1: (y, x) order, metres, domain 320-800 x 0-80 km)."""
    import netCDF4
    with netCDF4.Dataset(path) as nc:
        if var not in nc.variables:
            raise ValueError(f"{path}: no variable {var!r}")
        dims = nc.variables[var].dimensions
        if dims != ("y", "x"):
            raise ValueError(f"{path}: {var} has dimensions {dims}, expected ('y', 'x')")
        x = np.asarray(nc.variables["x"][:], dtype=float)
        y = np.asarray(nc.variables["y"][:], dtype=float)
        f = np.ma.filled(nc.variables[var][:].astype(float), np.nan)
    _check_axis(x, X_MIN, X_MAX, "x", path)
    _check_axis(y, Y_MIN, Y_MAX, "y", path)
    if f.shape != (y.size, x.size):
        raise ValueError(f"{path}: {var} shape {f.shape} != (ny, nx) = {(y.size, x.size)}")
    if not np.all(np.isfinite(f)):
        raise ValueError(f"{path}: {var} has missing or non-finite values")
    return x, y, f


def _factor(src: np.ndarray, dst: np.ndarray, name: str) -> int:
    ds, dd = src[1] - src[0], dst[1] - dst[0]
    r = dd / ds
    n = int(round(r))
    if n < 1 or abs(r - n) > 1e-9 or src.size != n * dst.size:
        raise ValueError(f"target {name} spacing {dd:g} m is not an integer multiple of the "
                         f"file's {ds:g} m covering the same domain")
    centres = src.reshape(dst.size, n).mean(axis=1)
    if not np.allclose(centres, dst, rtol=0, atol=1e-6):
        raise ValueError(f"target {name} centres are not the centres of {n}-cell blocks of the file grid")
    return n


def block_average(field: np.ndarray, x1: np.ndarray, y1: np.ndarray,
                  xc_m: np.ndarray, yc_m: np.ndarray) -> np.ndarray:
    """Mean of each n x n block of the 1 km ``field`` (y, x) onto centres (yc_m, xc_m)."""
    xc_m, yc_m = np.asarray(xc_m, float), np.asarray(yc_m, float)
    nx = _factor(x1, xc_m, "x")
    ny = _factor(y1, yc_m, "y")
    if nx != ny:
        raise ValueError(f"block factors differ in x ({nx}) and y ({ny})")
    return field.reshape(yc_m.size, ny, xc_m.size, nx).mean(axis=(1, 3))


def load_draft(path: str, xc_m: np.ndarray, yc_m: np.ndarray) -> np.ndarray:
    """Ice draft [m, <= 0] on the target grid: ``lowerSurface`` with calving applied at
    1 km (thickness upperSurface - lowerSurface < Hcalve = 100 m -> draft 0, AD16
    sect. 3.1.2), then block-averaged. Grounded ice keeps its draft (= bed), so the
    minimum-column rule of ``apply_scale_and_mask`` makes it land."""
    x, y, lower = _read_field(path, "lowerSurface")
    _, _, upper = _read_field(path, "upperSurface")
    draft = np.where(upper - lower < H_CALVE, 0.0, lower)
    return block_average(draft, x, y, xc_m, yc_m)


def load_bed(path: str, xc_m: np.ndarray, yc_m: np.ndarray) -> np.ndarray:
    """Bed [m] on the target grid: the file's ``bedrockTopography``, block-averaged
    exactly like ``load_draft`` (ruling R7)."""
    x, y, bed = _read_field(path, "bedrockTopography")
    return block_average(bed, x, y, xc_m, yc_m)


def load_floating(path: str, xc_m: np.ndarray, yc_m: np.ndarray) -> np.ndarray:
    """Floating-ice columns on the target grid (ruling R9): a coarse cell floats if a strict
    majority (> 1/2) of its 1 km sub-cells have ``floatingMask`` = 1. A tie (2 of 4 at
    2 km, 8 of 16 at 4 km) counts as not floating, so it keeps its ice unscaled. Calved
    sub-cells keep floatingMask = 1 in the file, but their draft is 0, so scaling has no
    effect on them."""
    x, y, fm = _read_field(path, "floatingMask")
    if not np.all((fm == 0) | (fm == 1)):
        raise ValueError(f"{path}: floatingMask is not 0/1")
    return block_average(fm, x, y, xc_m, yc_m) > 0.5


# ---- minimum water column and connectivity ----------------------------------------------

def cell_fractions(bathy: np.ndarray, draft: np.ndarray, dz: np.ndarray) -> np.ndarray:
    """Open fraction (nz, ny, nx) of each cell between the bed and the ice draft,
    with partial cells at top and bottom; level faces at 0, -dz[0], -dz[0]-dz[1], ..."""
    dz = np.asarray(dz, float)
    zf = np.concatenate([[0.0], -np.cumsum(dz)])
    top = np.minimum(draft[None], zf[:-1, None, None])
    bot = np.maximum(bathy[None], zf[1:, None, None])
    return np.clip(top - bot, 0.0, None) / dz[:, None, None]


def mitgcm_hfac(bathy: np.ndarray, draft: np.ndarray, dz: np.ndarray, hfac_min: float,
                hfac_min_dr: float = 0.0):
    """MITgcm's partial-cell factors and effective boundaries for bed ``bathy`` and ice
    draft ``draft`` (2-D, m, <= 0), line by line as model/src/ini_masks_etc.F:105-195:

    1. against the bed: h = clip((rF(k) - R_low) / drF(k), 0, 1); h = 0 if
       h < hFacMnSz / 2 or R_low >= Ro_surf, else max(h, hFacMnSz), with
       hFacMnSz = max(hFacMin, min(hFacMinDr / drF(k), 1)) (lines 105-122);
       R_low = rF(1) - sum(drF h) (lines 124-142);
    2. against the ice draft: h = max(h - max((rF(k) - Ro_surf) / drF(k), 0), 0), then the
       same minimum rule (lines 144-164); Ro_surf = R_low + sum(drF h) (lines 166-195).

    Returns (hFacC (nz, ny, nx), R_low, Ro_surf); R_low and Ro_surf are the bed and ice
    base the model actually uses (pkg/shelfice reads R_shelfIce from the draft and
    ini_masks_etc.F rewrites Ro_surf)."""
    dz = np.asarray(dz, float)
    r_low0 = np.asarray(bathy, float)
    ro_surf0 = np.asarray(draft, float)
    rF = np.concatenate([[0.0], -np.cumsum(dz)])
    recip = 1.0 / dz
    nz = dz.size
    h = np.zeros((nz,) + r_low0.shape)
    for k in range(nz):
        mn = max(hfac_min, min(hfac_min_dr * recip[k], 1.0))
        hl = np.minimum(np.maximum((rF[k] - r_low0) * recip[k], 0.0), 1.0)
        h[k] = np.where((hl < mn * 0.5) | (r_low0 >= ro_surf0), 0.0, np.maximum(hl, mn))
    r_low = rF[0] - np.sum(dz[:, None, None] * h, axis=0)
    for k in range(nz):
        mn = max(hfac_min, min(hfac_min_dr * recip[k], 1.0))
        hl = h[k] - np.maximum((rF[k] - ro_surf0) * recip[k], 0.0)
        hl = np.maximum(hl, 0.0)
        h[k] = np.where(hl < mn * 0.5, 0.0, np.maximum(hl, mn))
    ro_surf = r_low + np.sum(dz[:, None, None] * h, axis=0)
    return h, r_low, ro_surf


def _connected(frac: np.ndarray, wet: np.ndarray, seeds: np.ndarray) -> np.ndarray:
    """Wet columns reachable from ``seeds`` through 4-neighbours that share at least one
    level with open water in both columns. With useMin4hFacEdges, MITgcm opens a face at
    level k when hFacC > 0 on both sides: hFacW = min of the two hFacC
    (model/src/ini_masks_etc.F:251-252). Method 2 differs only when ice draft and bed of
    adjacent columns fall in the same level (ini_masks_etc.F:240-241)."""
    ny, nx = wet.shape
    open_k = frac > 0
    reach = np.zeros_like(wet)
    q = deque()
    for j, i in zip(*np.nonzero(seeds & wet)):
        reach[j, i] = True
        q.append((j, i))
    while q:
        j, i = q.popleft()
        for jj, ii in ((j - 1, i), (j + 1, i), (j, i - 1), (j, i + 1)):
            if (0 <= jj < ny and 0 <= ii < nx and wet[jj, ii] and not reach[jj, ii]
                    and np.any(open_k[:, j, i] & open_k[:, jj, ii])):
                reach[jj, ii] = True
                q.append((jj, ii))
    return reach


def apply_scale_and_mask(bathy: np.ndarray, draft: np.ndarray, s: float, dz: np.ndarray,
                         min_cells: int = 2, xc_m: np.ndarray | None = None,
                         x_connect_m: float = 790e3, floating: np.ndarray | None = None,
                         hfac_min: float | None = None):
    """Scale the draft of floating ice, enforce the minimum water column and remove wet
    regions not connected to the open ocean (Review Focus 2).

    - Ruling R9: draft' = s draft only where ``floating`` (``load_floating``: floating
      at s = 1); grounded columns keep their draft. Scaling all ice would open water under
      the grounded ice sheet for s < 1. ``floating`` is required when s != 1.
    - A column is wet if the open water between bed and draft' fills at least
      ``min_cells`` cells (sum of partial-cell fractions; 2 x 20 m = 40 m, AD16 sect. 3.1.4).
      Thinner columns become land, as grounded ice.
    - Only wet columns connected to a wet column with x >= ``x_connect_m`` (the
      restoring zone, AD16 Eq. 20) are kept. Without ``xc_m`` the last column is used.
    - With ``hfac_min`` (the namelist hFacMin), the open fractions are MITgcm's own
      partial cells (``mitgcm_hfac``, model/src/ini_masks_etc.F:105-195), not the exact
      ones: a sliver below hFacMin/2 is closed, as in the model, so a column connected
      only through such a sliver is not counted as connected (it would be a lake in the
      model). The returned bed and draft of wet columns are then snapped to the model's
      effective R_low / Ro_surf, so MITgcm's rounding leaves them unchanged and the load
      anomaly is evaluated at the ice base the model uses. Without ``hfac_min`` the exact
      fractions are used (reference behaviour).
    Dry columns get bed = draft = 0. Returns (bathy, draft, wet).
    """
    bathy = np.asarray(bathy, float)
    draft = np.asarray(draft, float)
    if floating is None:
        if s != 1.0:
            raise ValueError("apply_scale_and_mask: s != 1 needs the floating-ice mask")
        floating = np.zeros(draft.shape, dtype=bool)
    draft_s = np.where(floating, s * draft, draft)
    if hfac_min is None:
        frac = cell_fractions(bathy, draft_s, dz)
    else:
        frac, bathy, draft_s = mitgcm_hfac(bathy, draft_s, dz, hfac_min)
    wet = frac.sum(axis=0) >= min_cells - 1e-9
    seeds = np.zeros_like(wet)
    if xc_m is None:
        seeds[:, -1] = True
    else:
        seeds[:, np.asarray(xc_m) >= x_connect_m] = True
    wet = _connected(frac, wet, seeds)
    return np.where(wet, bathy, 0.0), np.where(wet, draft_s, 0.0), wet


def grounding_line(draft: np.ndarray, wet: np.ndarray, xc_m: np.ndarray, yc_m: np.ndarray) -> dict:
    """Grounding line on the trough centre line: in the two rows on either side of
    y = 40 km, the first wet column (smallest x) under ice (draft < 0). Returns the mean
    x [m] and depth (-draft) [m] of those cells, or None values if a row has none."""
    yc_m = np.asarray(yc_m)
    rows = np.argsort(np.abs(yc_m - LY / 2))[:2]
    xs, ds = [], []
    for j in rows:
        i = np.flatnonzero(wet[j] & (draft[j] < 0))
        if i.size == 0:
            return {"grounding_line_x_m": None, "grounding_line_depth_m": None}
        xs.append(float(xc_m[i[0]]))
        ds.append(float(-draft[j, i[0]]))
    return {"grounding_line_x_m": float(np.mean(xs)), "grounding_line_depth_m": float(np.mean(ds))}
