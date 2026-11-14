"""Model-ready bathymetry of Van Mijenfjorden for MITgcm (Exp. 4, Task 1b).

Pipeline (``build``): IBCAO v5.2 100 m grid (EPSG:3996) -> area-average onto the model grid
(UTM 33N, EPSG:32633) -> land/sea mask from the Norwegian Polar Institute S100 coastline ->
``clean`` (puddles, lakes, 1-cell inlets, minimum depth, partial cells) -> ``open_boundaries``
(zero normal depth gradient in the sponge) -> optional ``smooth_if_needed`` (off by default).
``report`` gives the controls written to figures/bathy_report.json.

Conventions: inside this module depth is in metres, positive downward, and 0 on land. Model arrays
are (ny, nx) with j = 0 the southern row and i = 0 the western column (MITgcm order); only
``write_bathy`` turns depth into the negative ``bathyFile`` field.

Data and sources:
- IBCAO v5.2 (Jakobsson et al. 2024, https://doi.org/10.1038/s41597-024-04278-w; grid DOI
  10.17043/ibcao-5.2): elevation (m, negative below mean sea level), SID and TID rasters, polar
  stereographic EPSG:3996, cell-centred values on a 100 m grid (GeoTIFF PixelIsArea).
- TID codes follow GEBCO (https://www.gebco.net/gebco-tid-grid): 10-17 direct measurements,
  40-46 indirect, 70-72 unknown. IBCAO also uses 47 and 100 (100 marks land in this window: 99.5 %
  of the TID = 100 cells are above sea level) [verificar en la documentacion de IBCAO 5].
- NPI S100 Kartdata (Norwegian Polar Institute, https://doi.org/10.21334/npolar.2014.645336c7,
  CC BY 4.0), land polygons ``S100_Land_f`` in ETRS89 / UTM 33N (EPSG:25833). ETRS89 and WGS 84
  differ by < 1 m in Svalbard [verificar], far below the 200 m grid, so the polygons are used as
  EPSG:32633.
- Partial cells follow MITgcm ``model/src/ini_masks_etc.F:105-121``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import shapefile
import shapely
from PIL import Image
from pyproj import Transformer
from scipy import ndimage
from shapely.geometry import Point, box, shape
from shapely.ops import nearest_points, unary_union

from . import vgrid

CRS_IBCAO = "EPSG:3996"
CRS_MODEL = "EPSG:32633"
TID_CLASSES = {"measured": range(10, 18), "indirect": range(40, 48), "unknown": range(70, 73), "land": (0, 100)}
IBCAO_FILES = {k: f"ibcao_v5_2_2026_{k}_100m_svalbard.tif" for k in ("depth", "sid", "tid")}
SIDES = ("W", "E", "S", "N")


# ----------------------------------------------------------------------------- reading IBCAO
def _read_geotiff(path: Path):
    """Read a single-band north-up GeoTIFF with PIL; return (array, x centres, y centres, epsg)."""
    Image.MAX_IMAGE_PIXELS = None
    im = Image.open(path)
    tags = dict(im.tag_v2)
    sx, sy, _ = tags[33550]                                   # ModelPixelScale
    _, _, _, X0, Y0, _ = tags[33922]                          # ModelTiepoint: pixel (0,0) corner -> (X0, Y0)
    keys = tags[34735]
    geo = {keys[4 + 4 * n]: keys[4 + 4 * n + 3] for n in range(keys[3])}
    if geo.get(1025, 1) != 1:
        raise ValueError(f"{path}: raster is not PixelIsArea")
    arr = np.array(im)
    ny, nx = arr.shape
    x = X0 + sx * (np.arange(nx) + 0.5)
    y = Y0 - sy * (np.arange(ny) + 0.5)
    return arr, x, y, geo.get(3072)


def bbox_utm_to_lonlat(bbox_utm, margin_m: float = 2000.0):
    """Lon/lat box that contains the UTM box (x0, x1, y0, y1) plus a margin."""
    x0, x1, y0, y1 = bbox_utm
    xs = np.concatenate([np.linspace(x0, x1, 50), np.full(50, x1), np.linspace(x1, x0, 50), np.full(50, x0)])
    ys = np.concatenate([np.full(50, y0), np.linspace(y0, y1, 50), np.full(50, y1), np.linspace(y1, y0, 50)])
    t = Transformer.from_crs(CRS_MODEL, "EPSG:4326", always_xy=True)
    lon, lat = t.transform(np.concatenate([xs - margin_m * np.sign(xs - (x0 + x1) / 2), xs]),
                           np.concatenate([ys - margin_m * np.sign(ys - (y0 + y1) / 2), ys]))
    return (float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max()))


def load_ibcao(path: Path, bbox_lonlat) -> dict:
    """Crop the IBCAO depth, SID and TID rasters to a lon/lat box (lon0, lon1, lat0, lat1).

    ``path`` is the directory with the three GeoTIFFs. Returns ``z`` (elevation, m, negative in the
    sea), ``sid``, ``tid`` as (ny, nx) arrays with y ascending, and the 1-D cell-centre coordinates
    ``x``, ``y`` in EPSG:3996 (m).
    """
    path = Path(path)
    lon0, lon1, lat0, lat1 = bbox_lonlat
    lons = np.concatenate([np.linspace(lon0, lon1, 60), np.full(60, lon1), np.linspace(lon1, lon0, 60), np.full(60, lon0)])
    lats = np.concatenate([np.full(60, lat0), np.linspace(lat0, lat1, 60), np.full(60, lat1), np.linspace(lat1, lat0, 60)])
    px, py = Transformer.from_crs("EPSG:4326", CRS_IBCAO, always_xy=True).transform(lons, lats)
    out = {}
    for key, fname in IBCAO_FILES.items():
        arr, x, y, epsg = _read_geotiff(path / fname)
        if epsg != 3996:
            raise ValueError(f"{fname}: EPSG {epsg}, expected 3996")
        ii = np.where((x >= px.min() - 100) & (x <= px.max() + 100))[0]
        jj = np.where((y >= py.min() - 100) & (y <= py.max() + 100))[0]
        if ii.size == 0 or jj.size == 0:
            raise ValueError("bbox outside the IBCAO window")
        sub = arr[jj[0]:jj[-1] + 1, ii[0]:ii[-1] + 1][::-1]        # flip to y ascending
        out[key] = sub
        out["x"], out["y"] = x[ii[0]:ii[-1] + 1], y[jj[0]:jj[-1] + 1][::-1]
    out["z"] = out.pop("depth").astype(np.float64)
    return out


# ----------------------------------------------------------------------------- regridding
def model_axes(bbox_utm, dx_m: float):
    """Cell-centre coordinates (xc, yc) of the model grid covering the UTM box (x0, x1, y0, y1)."""
    x0, x1, y0, y1 = bbox_utm
    nx, ny = round((x1 - x0) / dx_m), round((y1 - y0) / dx_m)
    if not np.isclose(nx * dx_m, x1 - x0) or not np.isclose(ny * dx_m, y1 - y0):
        raise ValueError("bbox is not a whole number of cells")
    return x0 + dx_m * (np.arange(nx) + 0.5), y0 + dx_m * (np.arange(ny) + 0.5)


def to_model_grid(src: dict, bbox_utm, dx_m: float) -> dict:
    """Area average of the IBCAO points onto the UTM model grid.

    Each source point (cell centre of the 100 m polar-stereographic grid) is transformed to
    EPSG:32633 with pyproj and assigned to the model cell that contains it; the model value is the
    mean over those points (not a sample at the cell centre). Land points count as depth 0.
    The sea volume is conserved only in aggregate: a model cell holds 3-5 source points of slightly
    different UTM area (the projection change aliases the 100 m points onto the 200 m cells), so
    ``depth_mean`` * dx^2 is not the exact sea volume of one cell, but its sum over a region
    matches the source volume (sum of depth x ``src_area``) to well under 1 % (tested).

    Returns depth_mean (m, >= 0), depth_wet_mean (mean over sea points only, NaN if none),
    wet_fraction (fraction of points below sea level), tid_mode, sid_mode, n_points, xc, yc,
    src_cell (flat model-cell index of every source point, -1 outside), src_area (area in UTM of
    each source cell, m^2, from the Jacobian of the projection change).
    """
    xc, yc = model_axes(bbox_utm, dx_m)
    nx, ny = xc.size, yc.size
    x0, _, y0, _ = bbox_utm
    XX, YY = np.meshgrid(src["x"], src["y"])
    t = Transformer.from_crs(CRS_IBCAO, CRS_MODEL, always_xy=True)
    U, V = t.transform(XX, YY)
    # source-cell area in UTM: |det J| * dx_src * dy_src, J by one-sided differences (h = 1 m)
    h = 1.0
    U1, V1 = t.transform(XX + h, YY)
    U2, V2 = t.transform(XX, YY + h)
    detJ = np.abs((U1 - U) * (V2 - V) - (V1 - V) * (U2 - U)) / h**2
    dxs = abs(src["x"][1] - src["x"][0]) if src["x"].size > 1 else 100.0
    dys = abs(src["y"][1] - src["y"][0]) if src["y"].size > 1 else 100.0
    area = detJ * dxs * dys

    i = np.floor((U - x0) / dx_m).astype(np.int64)
    j = np.floor((V - y0) / dx_m).astype(np.int64)
    inside = (i >= 0) & (i < nx) & (j >= 0) & (j < ny)
    cell = np.where(inside, j * nx + i, -1)
    c = cell[inside]
    ncell = nx * ny
    n = np.bincount(c, minlength=ncell).astype(float)
    if np.any(n == 0):
        raise ValueError(f"{int((n == 0).sum())} model cells contain no source point")
    z = src["z"][inside]
    d = np.maximum(-z, 0.0)
    sea = z < 0
    nsea = np.bincount(c, weights=sea, minlength=ncell)
    with np.errstate(invalid="ignore", divide="ignore"):
        dwet = np.bincount(c, weights=d * sea, minlength=ncell) / nsea

    def mode(v, nmax):
        counts = np.bincount(c * nmax + v.astype(np.int64), minlength=ncell * nmax).reshape(ncell, nmax)
        return counts.argmax(axis=1)

    return {
        "depth_mean": (np.bincount(c, weights=d, minlength=ncell) / n).reshape(ny, nx),
        "depth_wet_mean": dwet.reshape(ny, nx),
        "wet_fraction": (nsea / n).reshape(ny, nx),
        "tid_mode": mode(src["tid"][inside], 256).reshape(ny, nx),
        "sid_mode": mode(src["sid"][inside], int(src["sid"].max()) + 1).reshape(ny, nx),
        "n_points": n.reshape(ny, nx).astype(int),
        "xc": xc, "yc": yc, "dx_m": float(dx_m), "bbox_utm": tuple(float(v) for v in bbox_utm),
        "src_cell": cell, "src_area": area, "src_U": U, "src_V": V,
    }


# ----------------------------------------------------------------------------- coastline
def read_land(coastline_path: Path, bbox_utm, margin_m: float = 5000.0):
    """Union of the NPI land polygons that intersect the UTM box (+ margin), clipped to it."""
    x0, x1, y0, y1 = bbox_utm
    clip = box(x0 - margin_m, y0 - margin_m, x1 + margin_m, y1 + margin_m)
    geoms = []
    with shapefile.Reader(str(coastline_path)) as r:
        for sh in r.iterShapes():
            b = sh.bbox
            if b[2] < clip.bounds[0] or b[0] > clip.bounds[2] or b[3] < clip.bounds[1] or b[1] > clip.bounds[3]:
                continue
            g = shapely.make_valid(shape(sh.__geo_interface__))
            geoms.append(g.intersection(clip))
    land = unary_union([g for g in geoms if not g.is_empty])
    return land


def land_fraction(xc, yc, land, dx_m: float, tile: int = 25) -> np.ndarray:
    """Exact fraction of each model cell covered by the land geometry (tiles of tile x tile cells)."""
    nx, ny = len(xc), len(yc)
    frac = np.zeros((ny, nx))
    h = dx_m / 2
    shapely.prepare(land)
    for j0 in range(0, ny, tile):
        for i0 in range(0, nx, tile):
            js, is_ = slice(j0, min(j0 + tile, ny)), slice(i0, min(i0 + tile, nx))
            tb = box(xc[is_][0] - h, yc[js][0] - h, xc[is_][-1] + h, yc[js][-1] + h)
            if not land.intersects(tb):
                continue
            lt = land.intersection(tb)
            if abs(lt.area - tb.area) < 1e-6 * tb.area:
                frac[js, is_] = 1.0
                continue
            XC, YC = np.meshgrid(xc[is_], yc[js])
            cells = shapely.box(XC - h, YC - h, XC + h, YC + h)
            frac[js, is_] = shapely.area(shapely.intersection(cells, lt)) / dx_m**2
    return frac


def coast_mask(xc, yc, coastline_path, dx_m: float) -> np.ndarray:
    """True where more than half of the cell area is land in the NPI S100 polygons."""
    xc, yc = np.asarray(xc, float), np.asarray(yc, float)
    bb = (xc[0] - dx_m / 2, xc[-1] + dx_m / 2, yc[0] - dx_m / 2, yc[-1] + dx_m / 2)
    return land_fraction(xc, yc, read_land(coastline_path, bb), dx_m) > 0.5 + 1e-9


def mask_disagreements(wet_fraction, coast) -> dict:
    """Cells where IBCAO (sea if > half of its points are below sea level) and the NPI mask differ."""
    ibcao_sea = np.asarray(wet_fraction) > 0.5
    coast = np.asarray(coast, bool)
    a = int((ibcao_sea & coast).sum())
    b = int((~ibcao_sea & ~coast).sum())
    return {"ibcao_sea_npi_land": a, "ibcao_land_npi_sea": b, "total": a + b}


def make_wet(depth_mean, wet_fraction, coast, rule: str = "npi") -> np.ndarray:
    """Final sea mask. rule 'npi': the NPI coastline decides wherever IBCAO and NPI differ;
    'ibcao': more than half of the IBCAO points below sea level; 'both': sea in both."""
    ibcao_sea = np.asarray(wet_fraction) > 0.5
    npi_sea = ~np.asarray(coast, bool)
    if np.shape(depth_mean) != ibcao_sea.shape or ibcao_sea.shape != npi_sea.shape:
        raise ValueError("shape mismatch")
    if rule == "npi":
        return npi_sea.copy()
    if rule == "ibcao":
        return ibcao_sea
    if rule == "both":
        return ibcao_sea & npi_sea
    raise ValueError(f"unknown rule {rule!r}")


# ----------------------------------------------------------------------------- partial cells
def hfac_c(depth, faces, hfac_min: float, hfac_min_dr: float) -> np.ndarray:
    """MITgcm open fraction hFacC of each level, shape (nz,) + depth.shape.

    model/src/ini_masks_etc.F:108-121: hFacMnSz = max(hFacMin, min(hFacMinDr/drF(k), 1));
    hFac = clip((depth - rF_k)/drF_k, 0, 1); hFac = 0 if hFac < hFacMnSz/2, else max(hFac, hFacMnSz).
    Columns with depth <= 0 are land (hFac = 0).
    """
    f = np.asarray(faces, float)
    drf = np.diff(f)
    d = np.asarray(depth, float)
    shp = (-1,) + (1,) * d.ndim
    loc = np.clip((d[None] - f[:-1].reshape(shp)) / drf.reshape(shp), 0.0, 1.0)
    mn = np.maximum(hfac_min, np.minimum(hfac_min_dr / drf, 1.0)).reshape(shp)
    h = np.where(loc < 0.5 * mn, 0.0, np.maximum(loc, mn))
    return np.where(d[None] > 0, h, 0.0)


def partial_cells(depth, faces, hfac_min: float, hfac_min_dr: float) -> np.ndarray:
    """Depth that MITgcm actually represents, sum_k hFacC drF (m); idempotent."""
    return np.tensordot(np.diff(np.asarray(faces, float)), hfac_c(depth, faces, hfac_min, hfac_min_dr), axes=1)


# ----------------------------------------------------------------------------- topology
def _edge_mask(shape_, sides) -> np.ndarray:
    m = np.zeros(shape_, bool)
    if "W" in sides:
        m[:, 0] = True
    if "E" in sides:
        m[:, -1] = True
    if "S" in sides:
        m[0, :] = True
    if "N" in sides:
        m[-1, :] = True
    return m


def connected_to_ob(wet, ob_sides=("W", "N", "S")) -> np.ndarray:
    """Wet cells 4-connected (C-grid faces) to a wet cell on an open-boundary edge."""
    lab, _ = ndimage.label(wet)                    # default 2-D structure = 4-connectivity
    keep = np.unique(lab[_edge_mask(wet.shape, ob_sides) & wet])
    return np.isin(lab, keep[keep > 0])


def inlet_cells(wet, ob_sides=("W", "N", "S")) -> np.ndarray:
    """Wet cells with >= 3 land faces (only one open C-grid face). Outside an open edge is sea."""
    p = np.zeros((wet.shape[0] + 2, wet.shape[1] + 2), bool)
    p[1:-1, 1:-1] = wet
    if "W" in ob_sides:
        p[1:-1, 0] = True
    if "E" in ob_sides:
        p[1:-1, -1] = True
    if "S" in ob_sides:
        p[0, 1:-1] = True
    if "N" in ob_sides:
        p[-1, 1:-1] = True
    nwet = p[:-2, 1:-1].astype(int) + p[2:, 1:-1] + p[1:-1, :-2] + p[1:-1, 2:]
    return wet & (nwet <= 1)


def clean(depth, wet, dz_faces, min_depth_m: float, min_cells: int = 2, ob_sides=("W", "N", "S"), *,
          hfac_min: float, hfac_min_dr: float) -> dict:
    """Clean a wet mask and depth for MITgcm; every change is counted.

    ``dz_faces`` are the interface depths of the vertical grid (vgrid.faces(delR), m). ``hfac_min``
    and ``hfac_min_dr`` are hFacMin and hFacMinDr of the ``data`` namelist (from cases.yaml).
    1. puddles and lakes: wet cells not 4-connected to an open boundary become land;
    2. 1-cell inlets (wet cells with three land faces) become land; 1 and 2 repeat until stable;
    3. minimum depth: wet cells shallower than max(min_depth_m, depth of the bottom of level
       ``min_cells``) are deepened to it;
    4. partial cells: depth is replaced by the depth MITgcm represents (``partial_cells``).
    """
    faces = np.asarray(dz_faces, float)
    wet = np.asarray(wet, bool).copy()
    depth = np.where(wet, np.asarray(depth, float), 0.0)
    counts = {"removed_disconnected_cells": 0, "removed_components": 0, "filled_inlet_cells": 0}
    removed = []
    while True:
        conn = connected_to_ob(wet, ob_sides)
        lost = wet & ~conn
        if lost.any():
            lab, n = ndimage.label(lost)
            for k in range(1, n + 1):
                jj, ii = np.nonzero(lab == k)
                removed.append({"n_cells": int(jj.size), "j_mean": float(jj.mean()), "i_mean": float(ii.mean()),
                                "depth_max_m": float(depth[jj, ii].max())})
            counts["removed_disconnected_cells"] += int(lost.sum())
            counts["removed_components"] += int(n)
            wet &= conn
        inl = inlet_cells(wet, ob_sides)
        if not inl.any() and not lost.any():
            break
        counts["filled_inlet_cells"] += int(inl.sum())
        wet &= ~inl
    depth[~wet] = 0.0
    dmin = max(float(min_depth_m), float(faces[min_cells]))
    shallow = wet & (depth < dmin)
    counts["deepened_cells"] = int(shallow.sum())
    counts["deepened_volume_m_per_cell_mean"] = float((dmin - depth[shallow]).mean()) if shallow.any() else 0.0
    depth[shallow] = dmin
    deep = wet & (depth > faces[-1])
    counts["clipped_to_bottom_cells"] = int(deep.sum())
    counts["clipped_to_bottom_max_m"] = float(depth[deep].max() - faces[-1]) if deep.any() else 0.0
    depth = np.minimum(depth, faces[-1])
    new = partial_cells(depth, faces, hfac_min, hfac_min_dr)
    ch = wet & (np.abs(new - depth) > 1e-9)
    counts["partial_cell_adjusted_cells"] = int(ch.sum())
    counts["partial_cell_max_change_m"] = float(np.abs(new - depth)[wet].max()) if wet.any() else 0.0
    counts["partial_cell_mean_change_m"] = float((new - depth)[ch].mean()) if ch.any() else 0.0
    if np.any(new[wet] <= 0):
        raise RuntimeError("partial cells emptied a wet column")
    return {"depth": np.where(wet, new, 0.0), "wet": wet, "counts": counts, "removed": removed, "dmin_m": dmin}


def open_boundaries(depth, wet, n_sponge: int = 10, ob_sides=("W", "N", "S")) -> np.ndarray:
    """Depth (0 on land) with zero normal gradient in the sponge band of every open boundary.

    In the band of ``n_sponge`` cells along each open edge the depth (and so the mask) is copied
    from the first interior row/column, so d h / d n = 0 exactly there (W/E first, then S/N, so the
    corner blocks are uniform). Raises ValueError if land falls in an open corner (where two open
    boundaries meet): no island or coast may touch an open corner.
    """
    d = np.where(wet, np.asarray(depth, float), 0.0).copy()
    n = int(n_sponge)
    if "W" in ob_sides:
        d[:, :n] = d[:, [n]]
    if "E" in ob_sides:
        d[:, -n:] = d[:, [-n - 1]]
    if "S" in ob_sides:
        d[:n, :] = d[[n], :]
    if "N" in ob_sides:
        d[-n:, :] = d[[-n - 1], :]
    corners = {("S", "W"): (slice(0, n), slice(0, n)), ("N", "W"): (slice(-n, None), slice(0, n)),
               ("S", "E"): (slice(0, n), slice(-n, None)), ("N", "E"): (slice(-n, None), slice(-n, None))}
    for (a, b), sl in corners.items():
        if a in ob_sides and b in ob_sides and not np.all(d[sl] > 0):
            raise ValueError(f"land in the open {a}{b} corner")
    return d


def smooth_if_needed(depth, wet, max_step_cells: int | None, dz_faces=None) -> tuple[np.ndarray, dict]:
    """Optional local limiter of bottom steps; off by default (``max_step_cells=None``).

    If enabled, a wet column whose number of (partially) open levels exceeds that of a wet
    4-neighbour by more than ``max_step_cells`` is made shallower, to the face that satisfies the
    limit; repeated until no step exceeds it. Only the changed cells move; the change is reported.
    """
    depth = np.asarray(depth, float)
    if max_step_cells is None:
        return depth.copy(), {"applied": False, "changed_cells": 0, "max_change_m": 0.0, "mean_change_m": 0.0}
    faces = np.asarray(dz_faces, float)
    wet = np.asarray(wet, bool)
    d = depth.copy()
    big = 10**6
    for _ in range(10_000):
        k = np.where(wet, np.searchsorted(faces, d - 1e-9), big)
        p = np.pad(k, 1, constant_values=big)
        knb = np.minimum.reduce([p[:-2, 1:-1], p[2:, 1:-1], p[1:-1, :-2], p[1:-1, 2:]])
        over = wet & (knb < big) & (k > knb + max_step_cells)
        if not over.any():
            break
        d[over] = faces[np.minimum(knb[over] + max_step_cells, faces.size - 1)]
    ch = np.abs(d - depth)
    m = ch > 0
    return d, {"applied": True, "changed_cells": int(m.sum()), "max_change_m": float(ch.max()),
               "mean_change_m": float(ch[m].mean()) if m.any() else 0.0}


# ----------------------------------------------------------------------------- sounds
def _component(geom, pt):
    for g in getattr(geom, "geoms", [geom]):
        if g.contains(Point(pt)):
            return g
    raise ValueError(f"no land polygon contains {pt}")


def sound_sections(land, sounds: dict) -> dict:
    """Section line of each sound: the shortest segment between the island and the opposite shore.

    ``sounds`` maps name -> {"island_point": [x, y], "shore_point": [x, y], "search_box": [x0, x1,
    y0, y1]} (UTM; each point on land inside its polygon). Both polygons are clipped to the search
    box before the shortest segment is taken (the mainland wraps around the fjord). Returns
    name -> {p0 (island end), p1 (shore end), length_m, width_npi_m}; width_npi_m is the length of
    the segment in NPI water (smaller islets on the segment are subtracted).
    """
    out = {}
    for name, s in sounds.items():
        isl = _component(land, s["island_point"])
        sho = _component(land, s["shore_point"])
        if isl.equals(sho):
            raise ValueError(f"{name}: island and shore are the same polygon")
        if "search_box" in s:
            x0, x1, y0, y1 = s["search_box"]
            sb = box(x0, y0, x1, y1)
            isl, sho = isl.intersection(sb), sho.intersection(sb)
        a, b = nearest_points(isl, sho)
        seg = shapely.LineString([a, b])
        out[name] = {"p0": [a.x, a.y], "p1": [b.x, b.y], "length_m": float(seg.length),
                     "width_npi_m": float(seg.length - seg.intersection(land).length)}
        if s.get("islet_points"):
            # passages of a sound split by islets: shortest gaps along the chain island-islets-shore
            chain = [isl] + [_component(land, q) for q in s["islet_points"]] + [sho]
            out[name]["passages_npi_m"] = [float(g1.distance(g2)) for g1, g2 in zip(chain[:-1], chain[1:])]
    return out


def _samples(p0, p1, ds):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    L = float(np.hypot(*(p1 - p0)))
    n = max(int(np.ceil(L / ds)), 1)
    s = (np.arange(n) + 0.5) / n
    return p0[None] + s[:, None] * (p1 - p0)[None], L / n


def section_stats(depth, wet, sections: dict, dx_m: float, xc, yc, src: dict | None = None) -> dict:
    """Width, maximum depth and cross-sectional area along each section line.

    The line is sampled every dx/20; each sample takes the depth of the model cell that contains it
    (and, if ``src`` is given, of the nearest IBCAO point). Area = sum(depth * ds) (m^2), normal to
    the line. wet_cells_on_section counts the distinct wet model cells the line crosses.
    """
    x0, y0 = xc[0] - dx_m / 2, yc[0] - dx_m / 2
    out = {}
    for name, sec in sections.items():
        P, ds = _samples(sec["p0"], sec["p1"], dx_m / 20)
        i = np.clip(((P[:, 0] - x0) // dx_m).astype(int), 0, len(xc) - 1)
        j = np.clip(((P[:, 1] - y0) // dx_m).astype(int), 0, len(yc) - 1)
        dd = np.where(wet[j, i], depth[j, i], 0.0)
        cells = {(a, b) for a, b in zip(j, i) if wet[a, b]}
        r = {"length_m": float(len(P) * ds),
             "grid": {"wet_cells_on_section": len(cells), "wet_length_m": float((dd > 0).sum() * ds),
                      "max_depth_m": float(dd.max()), "area_m2": float(dd.sum() * ds)}}
        if src is not None:
            to_ps = Transformer.from_crs(CRS_MODEL, CRS_IBCAO, always_xy=True)
            px, py = to_ps.transform(P[:, 0], P[:, 1])
            ii = np.abs(src["x"][None] - px[:, None]).argmin(axis=1)
            jj = np.abs(src["y"][None] - py[:, None]).argmin(axis=1)
            zs = np.maximum(-src["z"][jj, ii], 0.0)
            r["ibcao"] = {"wet_length_m": float((zs > 0).sum() * ds), "max_depth_m": float(zs.max()),
                          "area_m2": float(zs.sum() * ds)}
        out[name] = r
    return out


def _cell_of(pt, xc, yc, dx_m):
    return (int((pt[1] - (yc[0] - dx_m / 2)) // dx_m), int((pt[0] - (xc[0] - dx_m / 2)) // dx_m))


def closure_disks(xc, yc, sections: dict, margin_m: float) -> dict:
    """Per sound, the model cells within (half the section length + margin) of the section midpoint."""
    XC, YC = np.meshgrid(xc, yc)
    out = {}
    for name, s in sections.items():
        m = (np.asarray(s["p0"]) + np.asarray(s["p1"])) / 2
        r = s["length_m"] / 2 + margin_m
        out[name] = (XC - m[0]) ** 2 + (YC - m[1]) ** 2 <= r**2
    return out


def sill_depth(depth, wet, a, b, blocked=None) -> float:
    """Bottleneck depth between cells a and b: the largest h such that cells with depth >= h
    (4-connected, outside ``blocked``) join a and b. 0 if they are not connected at all."""
    ok = wet & (depth > 0) if blocked is None else wet & (depth > 0) & ~blocked
    lab, _ = ndimage.label(ok)
    if lab[a] == 0 or lab[a] != lab[b]:
        return 0.0
    levels = np.unique(depth[ok])
    lo, hi = 0, levels.size - 1                       # levels[lo] always connects
    while lo < hi:
        mid = (lo + hi + 1) // 2
        lab, _ = ndimage.label(ok & (depth >= levels[mid]))
        if lab[a] > 0 and lab[a] == lab[b]:
            lo = mid
        else:
            hi = mid - 1
    return float(levels[lo])


def open_width_test(wet, a, b, blocked=None, min_width: int = 3) -> bool:
    """True if a and b stay 4-connected after eroding the sea by (min_width-1)/2 cells with a cross
    structure, i.e. the passage is at least ``min_width`` cells wide everywhere (odd min_width)."""
    ok = wet if blocked is None else wet & ~blocked
    r = (min_width - 1) // 2
    er = ndimage.binary_erosion(ok, structure=ndimage.generate_binary_structure(2, 1), iterations=r,
                                border_value=1) if r > 0 else ok
    lab, _ = ndimage.label(er)
    return bool(lab[a] > 0 and lab[a] == lab[b])


# ----------------------------------------------------------------------------- output
def write_bathy(path, depth, wet) -> None:
    """MITgcm bathyFile: big-endian float64, -depth in the sea and 0 on land, shape (ny, nx) with
    x varying fastest and j = 0 (south) first, as MITgcm reads it."""
    d = np.where(np.asarray(wet, bool), -np.abs(np.asarray(depth, float)), 0.0)
    np.ascontiguousarray(d, dtype=">f8").tofile(Path(path))


# ----------------------------------------------------------------------------- PD1 mask edits
def region_cells(xc, yc, region: dict, fill_mask=None) -> np.ndarray:
    """Model cells of a force_land / force_sea region of cases.yaml.

    ``region`` has ``polygon`` (list of [x, y] vertices, EPSG:32633; cells whose centre is inside)
    and/or ``seed`` ([x, y]). With a seed, only the cells of ``fill_mask`` 4-connected to the seed
    cell (inside the polygon, if one is given) are selected.
    """
    XC, YC = np.meshgrid(xc, yc)
    if "polygon" in region:
        sel = shapely.contains_xy(shapely.Polygon(region["polygon"]), XC, YC)
    else:
        sel = np.ones(XC.shape, bool)
    if "seed" in region:
        if fill_mask is None:
            raise ValueError("a seeded region needs a fill mask")
        dx = xc[1] - xc[0]
        js, is_ = _cell_of(region["seed"], xc, yc, dx)
        lab, _ = ndimage.label(fill_mask & sel)
        if lab[js, is_] == 0:
            raise ValueError(f"seed {region['seed']} is not in the fill mask")
        sel = lab == lab[js, is_]
    return sel


def apply_force_regions(wet, g: dict, faces, force_land: dict | None, force_sea: dict | None):
    """PD1 edits of the sea mask, all from cases.yaml (``bathy: force_land / force_sea``).

    force_sea regions (seed fill on the IBCAO sea, wet_fraction > 0.5) turn land cells into sea,
    optionally only where ``require_measured_tid`` (GEBCO TID 10-17) and the IBCAO area-mean depth
    exceeds ``min_depth_cells`` model levels. force_land regions (seed fill on the current sea)
    then turn sea into land. Returns (wet, counts per region: cells changed and km^2).
    """
    wet = np.asarray(wet, bool).copy()
    dx = g["xc"][1] - g["xc"][0]
    counts = {}
    for name, r in (force_sea or {}).items():
        sel = region_cells(g["xc"], g["yc"], r, fill_mask=(g["wet_fraction"] > 0.5) | wet) & ~wet
        if r.get("require_measured_tid"):
            sel &= np.isin(g["tid_mode"], list(TID_CLASSES["measured"]))
        if r.get("min_depth_cells") is not None:
            sel &= g["depth_mean"] > np.asarray(faces)[int(r["min_depth_cells"])]
        wet |= sel
        counts[f"force_sea:{name}"] = {"cells": int(sel.sum()), "area_km2": float(sel.sum() * dx**2 / 1e6)}
    for name, r in (force_land or {}).items():
        sel = region_cells(g["xc"], g["yc"], r, fill_mask=wet) & wet
        wet &= ~sel
        counts[f"force_land:{name}"] = {"cells": int(sel.sum()), "area_km2": float(sel.sum() * dx**2 / 1e6)}
    return wet, counts


def check_sponge(depth, n_sponge: int, ob_sides) -> None:
    """Raise if the sponge band of an open boundary does not have d h / d n = 0 exactly."""
    n = int(n_sponge)
    d = np.asarray(depth)
    ok = True
    if "W" in ob_sides:
        ok &= bool(np.all(d[:, :n + 1] == d[:, [n]]))
    if "E" in ob_sides:
        ok &= bool(np.all(d[:, -n - 1:] == d[:, [-n - 1]]))
    if "S" in ob_sides:
        ok &= bool(np.all(d[:n + 1, :] == d[[n], :]))
    if "N" in ob_sides:
        ok &= bool(np.all(d[-n - 1:, :] == d[[-n - 1], :]))
    if not ok:
        raise RuntimeError("sponge band lost d h / d n = 0")


# ----------------------------------------------------------------------------- pipeline
def build(cfg: dict, data_dir: Path) -> tuple[dict, dict]:
    """Run the whole pipeline for a config of cases.yaml (configs: fjord_vmf:)."""
    data_dir = Path(data_dir)
    b = cfg["bathy"]
    bbox, dx = tuple(cfg["bbox_utm"]), float(cfg["dx_m"])
    faces = vgrid.faces(cfg["delR"])
    ob = tuple(cfg["ob_sides"])
    src = load_ibcao(data_dir / b["ibcao_dir"], bbox_utm_to_lonlat(bbox))
    g = to_model_grid(src, bbox, dx)
    land = read_land(data_dir / b["coastline"], bbox)
    frac = land_fraction(g["xc"], g["yc"], land, dx)
    coast = frac > 0.5 + 1e-9
    wet0 = make_wet(g["depth_mean"], g["wet_fraction"], coast, rule=b["mask_rule"])
    wet1, forced = apply_force_regions(wet0, g, faces, b.get("force_land"), b.get("force_sea"))
    kw = dict(hfac_min=cfg["hFacMin"], hfac_min_dr=cfg["hFacMinDr"])
    c1 = clean(np.where(wet1, g["depth_mean"], 0.0), wet1, faces, cfg["min_depth_m"], 2, ob, **kw)
    d_ob = open_boundaries(c1["depth"], c1["wet"], cfg["n_sponge"], ob)
    sponge = {"changed_cells": int((np.abs(d_ob - c1["depth"]) > 1e-9).sum()),
              "land_to_sea": int(((d_ob > 0) & ~c1["wet"]).sum()),
              "sea_to_land": int(((d_ob == 0) & c1["wet"]).sum())}
    c2 = clean(d_ob, d_ob > 0, faces, cfg["min_depth_m"], 2, ob, **kw)
    check_sponge(c2["depth"], cfg["n_sponge"], ob)
    depth, info = smooth_if_needed(c2["depth"], c2["wet"], cfg.get("max_step_cells"), faces)
    sections = sound_sections(land, b["sounds"])
    g.update({"land_fraction": frac, "coast": coast, "wet_npi": wet0, "wet_forced": wet1, "forced": forced, "clean1": c1, "clean2": c2,
              "sponge": sponge, "smooth": info, "depth": depth, "wet": depth > 0, "faces": faces,
              "sections": sections, "land": land, "cfg": cfg})
    return src, g


def _hyps(depths, areas, edges):
    h, _ = np.histogram(depths, bins=edges, weights=areas)
    return h


def inner_fjord_mask(grid: dict):
    """Inner fjord (model sea east of the two sounds): the wet cells 4-connected to the east seed
    with both sounds closed by their ``closure_disks``. Returns (inner, outer, disks), where outer
    is the sea connected to the west seed (Bellsund and the shelf) and disks the closure disks."""
    b = grid["cfg"]["bathy"]
    xc, yc, dx, wet = grid["xc"], grid["yc"], grid["dx_m"], grid["wet"]
    west = _cell_of(b["seeds"]["west"], xc, yc, dx)
    east = _cell_of(b["seeds"]["east"], xc, yc, dx)
    disks = closure_disks(xc, yc, grid["sections"], b["closure_margin_m"])
    both = np.logical_or.reduce(list(disks.values()))
    lab, _ = ndimage.label(wet & ~both)
    if lab[east] == 0 or lab[west] == 0:
        raise RuntimeError("a seed is not in the sea")
    if lab[east] == lab[west]:
        raise RuntimeError("inner fjord and Bellsund still connected with both sounds closed")
    return lab == lab[east], lab == lab[west], disks


def report(src: dict, grid: dict) -> dict:
    """Controls of the model bathymetry (figures/bathy_report.json)."""
    cfg, b = grid["cfg"], grid["cfg"]["bathy"]
    xc, yc, dx = grid["xc"], grid["yc"], grid["dx_m"]
    depth, wet, faces = grid["depth"], grid["wet"], grid["faces"]
    secs = grid["sections"]
    west = _cell_of(b["seeds"]["west"], xc, yc, dx)
    east = _cell_of(b["seeds"]["east"], xc, yc, dx)
    gin, outer, disks = inner_fjord_mask(grid)
    both = np.logical_or.reduce(list(disks.values()))

    # For the IBCAO side the inner region is extended to the land cells within
    # ``inner_buffer_cells`` whose nearest model sea cell is in the inner fjord (so IBCAO sea points
    # under NPI land along the fjord coast count, but none across Akseloya into Bellsund).
    dist, (jn, in_) = ndimage.distance_transform_edt(~(wet & ~both), return_indices=True)
    inner = gin[jn, in_] & (dist <= b["inner_buffer_cells"]) & ~both
    cell_inner = np.concatenate([inner.ravel(), [False]])
    pin = cell_inner[grid["src_cell"]]                    # src_cell = -1 -> last element (False)
    z, area = src["z"], grid["src_area"]
    d_src = np.maximum(-z[pin], 0.0)
    a_src = area[pin]
    v_ibcao = float((d_src * a_src).sum())
    v_grid = float(depth[gin].sum() * dx**2)
    v_npi_mask = float(np.where(grid["wet_npi"] & inner, grid["depth_mean"], 0.0).sum() * dx**2)

    # hypsometry and area/volume below each model level, inner fjord
    edges = np.arange(0.0, faces[-1] + 10.0, 10.0)
    hyps = {"bin_edges_m": edges.tolist(),
            "area_ibcao_m2": _hyps(d_src[d_src > 0], a_src[d_src > 0], edges).tolist(),
            "area_grid_m2": _hyps(depth[gin], np.full(gin.sum(), dx**2), edges).tolist()}
    below = {"level_face_m": faces.tolist(),
             "area_ibcao_m2": [float(a_src[d_src > f].sum()) for f in faces],
             "area_grid_m2": [float((depth[gin] > f).sum() * dx**2) for f in faces],
             "volume_ibcao_m3": [float((np.maximum(d_src - f, 0) * a_src).sum()) for f in faces],
             "volume_grid_m3": [float(np.maximum(depth[gin] - f, 0).sum() * dx**2) for f in faces]}

    # sounds: section statistics, sill (bottleneck) depth, 3-cell open-width test
    stats = section_stats(depth, wet, secs, dx, xc, yc, src)
    # IBCAO-native sill depth: the same bottleneck search on the 100 m grid, with the source points
    # inside the closure disks removed and sea = below sea level in IBCAO
    U, V = grid["src_U"], grid["src_V"]
    ps_sea = z < 0
    ps_west = np.unravel_index(np.argmin((U - b["seeds"]["west"][0]) ** 2 + (V - b["seeds"]["west"][1]) ** 2), z.shape)
    ps_east = np.unravel_index(np.argmin((U - b["seeds"]["east"][0]) ** 2 + (V - b["seeds"]["east"][1]) ** 2), z.shape)
    sounds = {}
    for name in secs:
        other = np.logical_or.reduce([m for k, m in disks.items() if k != name] or [np.zeros_like(wet)])
        s = dict(secs[name])
        s.update(stats[name])
        s["sill_depth_grid_m"] = sill_depth(depth, wet, west, east, blocked=other)
        s["open_with_min_width_3"] = open_width_test(wet, west, east, blocked=other, min_width=3)
        s["open_with_min_width_1"] = open_width_test(wet, west, east, blocked=other, min_width=1)
        o = secs[[k for k in secs if k != name][0]] if len(secs) > 1 else None
        if o is not None:
            m = (np.asarray(o["p0"]) + np.asarray(o["p1"])) / 2
            r = o["length_m"] / 2 + b["closure_margin_m"]
            blocked_ps = (U - m[0]) ** 2 + (V - m[1]) ** 2 <= r**2
        else:
            blocked_ps = np.zeros_like(ps_sea)
        s["sill_depth_ibcao_m"] = sill_depth(np.maximum(-z, 0.0), ps_sea, ps_west, ps_east, blocked=blocked_ps)
        sounds[name] = s

    # data quality: TID of the inner-fjord sea cells and points
    def tid_frac(t):
        t = np.asarray(t).ravel()
        return {k: float(np.isin(t, list(v)).mean()) if t.size else 0.0 for k, v in TID_CLASSES.items()}

    tid = {"cells_inner": tid_frac(grid["tid_mode"][gin]),
           "points_inner_sea": tid_frac(src["tid"][pin][z[pin] < 0]),
           "cells_all_sea": tid_frac(grid["tid_mode"][wet])}

    # mask disagreements (whole box, and inside the inner fjord)
    dis = mask_disagreements(grid["wet_fraction"], grid["coast"])
    dis_inner = mask_disagreements(grid["wet_fraction"][inner], grid["coast"][inner])

    # Vallunden
    vp = _cell_of(b["vallunden_point"], xc, yc, dx)
    vall = {"cell": [int(vp[0]), int(vp[1])], "npi_sea": bool(grid["wet_npi"][vp]),
            "ibcao_wet_fraction": float(grid["wet_fraction"][vp]), "final_wet": bool(wet[vp]),
            "forced_land": bool(grid["wet_npi"][vp] and not grid["wet_forced"][vp])}
    if vall["npi_sea"]:
        lab0, _ = ndimage.label(grid["wet_npi"])
        comp = lab0 == lab0[vp]
        vall["npi_component_cells"] = int(comp.sum())
        vall["npi_component_connected_to_ob"] = bool(connected_to_ob(grid["wet_npi"], tuple(cfg["ob_sides"]))[vp])
    if vall["final_wet"]:
        # how Vallunden is joined to the fjord in the final grid: passage width and sill (bottleneck) depth
        vall["open_to_fjord_min_width_3"] = open_width_test(wet, vp, east, min_width=3)
        vall["open_to_fjord_min_width_1"] = open_width_test(wet, vp, east, min_width=1)
        vall["sill_depth_to_fjord_grid_m"] = sill_depth(depth, wet, vp, east)
        vall["depth_grid_at_point_m"] = float(depth[vp])

    # maximum depths
    maxd = {"inner_grid_m": float(depth[gin].max()), "inner_ibcao_m": float(d_src.max()),
            "box_grid_m": float(depth[wet].max()), "box_depth_mean_npi_sea_m": float(grid["depth_mean"][grid["wet_npi"]].max()),
            "outer_bellsund_grid_m": float(depth[outer].max()) if outer.any() else 0.0}

    rep = {
        "config": {k: cfg[k] for k in ("bbox_utm", "dx_m", "hFacMin", "hFacMinDr", "min_depth_m", "n_sponge", "ob_sides")},
        "grid_shape_ny_nx": list(depth.shape), "nz": len(faces) - 1, "max_face_m": float(faces[-1]),
        "wet_cells": int(wet.sum()),
        "inner_fjord": {"cells": int(gin.sum()), "volume_ibcao_m3": v_ibcao, "volume_grid_m3": v_grid,
                        "volume_after_npi_mask_m3": v_npi_mask,
                        "volume_ratio_grid_to_ibcao": v_grid / v_ibcao, "area_grid_m2": float(gin.sum() * dx**2),
                        "area_ibcao_m2": float(a_src[d_src > 0].sum())},
        "hypsometry_inner": hyps, "below_levels_inner": below, "max_depth": maxd, "sounds": sounds,
        "tid_fraction": tid, "clean_pass1": grid["clean1"]["counts"], "clean_pass2": grid["clean2"]["counts"],
        "removed_components_pass1": sorted(grid["clean1"]["removed"], key=lambda r: -r["n_cells"])[:15],
        "force_regions": grid["forced"], "sponge": grid["sponge"], "smoothing": grid["smooth"],
        "mask_disagreements": dis, "mask_disagreements_inner": dis_inner, "vallunden": vall,
    }
    for r in rep["removed_components_pass1"]:
        r["x_utm"] = float(xc[0] + r["i_mean"] * dx)
        r["y_utm"] = float(yc[0] + r["j_mean"] * dx)
    rep["controls"] = controls(rep, b.get("controls", {}), b.get("classifications", {}) or {})
    return rep


def _get(rep, dotted):
    v = rep
    for k in dotted.split("."):
        v = v[k]
    return v


def controls(rep: dict, spec: dict, classifications: dict) -> list:
    """Evaluate the pre-declared controls of cases.yaml (bathy: controls:) on the report.

    Each control: value (dotted key into the report), reference (number or dotted key), rel_tol or
    abs_tol, source. pass = |value - reference| <= tolerance. A failed control keeps pass = False
    and carries its classification (a)-(d) from cases.yaml if one has been written.
    """
    out = []
    for name, c in spec.items():
        val = float(_get(rep, c["value"]))
        ref = c["reference"]
        ref = float(_get(rep, ref)) if isinstance(ref, str) else float(ref)
        tol = c["abs_tol"] if "abs_tol" in c else c["rel_tol"] * abs(ref)
        out.append({"name": name, "value": val, "reference": ref, "tolerance": float(tol),
                    "diff": val - ref, "pass": bool(abs(val - ref) <= tol), "source": c["source"],
                    "classification": classifications.get(name)})
    return out
