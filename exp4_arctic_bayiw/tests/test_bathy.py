"""Tests for exp4_arctic_bayiw.gendata.bathy (Task 1b: model-ready Van Mijenfjorden bathymetry).

Synthetic tests need no external data. The tests marked "real IBCAO" run the full pipeline on the
IBCAO v5.2 window and the NPI S100 coastline in data_external/ and are skipped if those files are missing.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest
import shapefile
import yaml
from pyproj import Transformer
from shapely.geometry import box

from exp4_arctic_bayiw.gendata import bathy, vgrid
from pfcommon.mitgcm_io import read_bin

PKG = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((PKG / "cases.yaml").read_text())["configs"]["fjord_vmf"]

# a small UTM box in Van Mijenfjorden (inside zone 33X), used by the synthetic tests
X0, Y0 = 510_000.0, 8_630_000.0


def synthetic_src(x0, x1, y0, y1, zfun, step=100.0, margin=1500.0):
    """IBCAO-like source: 100 m grid in EPSG:3996 covering the UTM box [x0,x1]x[y0,y1] plus a margin.

    zfun(X_utm, Y_utm) gives the elevation (m, negative in the sea) at each source point.
    """
    to_ps = Transformer.from_crs("EPSG:32633", "EPSG:3996", always_xy=True)
    cx, cy = to_ps.transform([x0 - margin, x1 + margin, x0 - margin, x1 + margin],
                             [y0 - margin, y0 - margin, y1 + margin, y1 + margin])
    xs = np.arange(np.floor(min(cx) / step) * step + step / 2, max(cx), step)
    ys = np.arange(np.floor(min(cy) / step) * step + step / 2, max(cy), step)
    XX, YY = np.meshgrid(xs, ys)
    to_utm = Transformer.from_crs("EPSG:3996", "EPSG:32633", always_xy=True)
    U, V = to_utm.transform(XX, YY)
    z = zfun(U, V).astype(np.float64)
    return {"z": z, "sid": np.ones_like(z, dtype=np.uint16), "tid": np.full(z.shape, 11, dtype=np.uint8),
            "x": xs, "y": ys}


# --------------------------------------------------------------------------- to_model_grid
def test_to_model_grid_ramp_is_exact_point_mean():
    bb = (X0, X0 + 4000.0, Y0, Y0 + 3000.0)
    ramp = lambda U, V: -(60.0 + 0.004 * (U - X0) + 0.009 * (V - Y0))  # noqa: E731
    src = synthetic_src(*bb, ramp)
    g = bathy.to_model_grid(src, bb, 200.0)
    assert g["depth_mean"].shape == (15, 20)            # (ny, nx), j = 0 is the southern row
    # independent brute-force mean of the source points whose UTM position falls in each cell
    to_utm = Transformer.from_crs("EPSG:3996", "EPSG:32633", always_xy=True)
    XX, YY = np.meshgrid(src["x"], src["y"])
    U, V = to_utm.transform(XX, YY)
    d = -src["z"]
    err = 0.0
    for j in range(15):
        for i in range(20):
            m = ((U >= X0 + 200 * i) & (U < X0 + 200 * (i + 1))
                 & (V >= Y0 + 200 * j) & (V < Y0 + 200 * (j + 1)))
            err = max(err, abs(d[m].mean() - g["depth_mean"][j, i]))
    assert err < 1e-6
    # and the point mean of a linear ramp is within half a cell of the ramp value at the cell centre
    XC, YC = np.meshgrid(g["xc"], g["yc"])
    centre = 60.0 + 0.004 * (XC - X0) + 0.009 * (YC - Y0)
    assert np.max(np.abs(g["depth_mean"] - centre)) < (0.004 + 0.009) * 100.0
    assert np.all(g["wet_fraction"] == 1.0)
    assert np.all(g["tid_mode"] == 11)


def test_to_model_grid_conserves_basin_volume():
    bb = (X0, X0 + 12000.0, Y0, Y0 + 10000.0)
    xm, ym = X0 + 6000.0, Y0 + 5000.0
    basin = lambda U, V: -100.0 * np.exp(-((U - xm) ** 2 + (V - ym) ** 2) / (2 * 2000.0 ** 2))  # noqa: E731
    src = synthetic_src(*bb, basin)
    g = bathy.to_model_grid(src, bb, 200.0)
    v_grid = g["depth_mean"].sum() * 200.0 ** 2
    inside = g["src_cell"] >= 0
    v_src = np.sum(np.maximum(-src["z"], 0.0)[inside] * g["src_area"][inside])
    assert v_grid == pytest.approx(v_src, rel=0.01)
    # the source point area in UTM is the 100 m x 100 m polar-stereographic cell times the scale ratio
    assert 0.9e4 < g["src_area"].mean() < 1.1e4


def test_to_model_grid_axes_and_hemisphere():
    lon, lat = 15.5, 77.8
    xu, yu = Transformer.from_crs("EPSG:4326", "EPSG:32633", always_xy=True).transform(lon, lat)
    # independent check of the projection: zone 33 central meridian is 15 E, so 15.5 E is ~11.8 km east;
    # 77.8 N is ~8636 km north of the equator (scale 0.9996)
    assert xu == pytest.approx(511_760.0, abs=1500.0)
    assert yu == pytest.approx(8_636_000.0, abs=3000.0)
    bb = (xu - 2000.0, xu + 3000.0, yu - 1000.0, yu + 2000.0)   # asymmetric so a flipped axis fails
    xp, yp = Transformer.from_crs("EPSG:4326", "EPSG:3996", always_xy=True).transform(lon, lat)

    def one_wet_point(U, V):
        z = np.full(U.shape, 50.0)
        return z

    src = synthetic_src(*bb, one_wet_point)
    jj = np.argmin(np.abs(src["y"] - yp))
    ii = np.argmin(np.abs(src["x"] - xp))
    src["z"][jj, ii] = -10.0
    xs_u, ys_u = Transformer.from_crs("EPSG:3996", "EPSG:32633", always_xy=True).transform(
        src["x"][ii], src["y"][jj])
    g = bathy.to_model_grid(src, bb, 200.0)
    j, i = np.argwhere(g["wet_fraction"] > 0)[0]
    assert np.count_nonzero(g["wet_fraction"]) == 1
    assert i == int((xs_u - bb[0]) // 200.0) and j == int((ys_u - bb[2]) // 200.0)
    assert abs(g["xc"][i] - xu) < 200.0 and abs(g["yc"][j] - yu) < 200.0


# --------------------------------------------------------------------------- coast and mask
def write_land_shapefile(path: Path, polygons):
    w = shapefile.Writer(str(path), shapeType=shapefile.POLYGON)
    w.field("OBJTYPE", "C", size=32)
    for poly in polygons:
        ext = list(poly.exterior.coords)[::-1] if poly.exterior.is_ccw else list(poly.exterior.coords)
        w.poly([ext])  # shapefile outer rings are clockwise
        w.record("Landareal")
    w.close()
    return path


def test_make_wet_npi_wins_and_counts_disagreement():
    wet_fraction = np.ones((3, 3))
    coast = np.zeros((3, 3), dtype=bool)
    coast[1, 1] = True                                    # NPI says land, IBCAO says sea
    depth = np.full((3, 3), 20.0)
    wet = bathy.make_wet(depth, wet_fraction, coast, rule="npi")
    assert not wet[1, 1] and wet.sum() == 8
    dis = bathy.mask_disagreements(wet_fraction, coast)
    assert dis["ibcao_sea_npi_land"] == 1 and dis["ibcao_land_npi_sea"] == 0 and dis["total"] == 1


def test_coast_mask_half_area_rule(tmp_path):
    # land covers x < X0 + 500: cells 0-1 are all land, cell 2 is 50 % land (not land), 3+ water
    shp = write_land_shapefile(tmp_path / "land", [box(X0 - 5000, Y0 - 5000, X0 + 500, Y0 + 5000)])
    xc = X0 + 100.0 + 200.0 * np.arange(6)
    yc = Y0 + 100.0 + 200.0 * np.arange(4)
    land = bathy.coast_mask(xc, yc, shp, 200.0)
    assert land.shape == (4, 6)
    assert land[:, :2].all() and not land[:, 2:].any()
    shp2 = write_land_shapefile(tmp_path / "land2", [box(X0 - 5000, Y0 - 5000, X0 + 501, Y0 + 5000)])
    assert bathy.coast_mask(xc, yc, shp2, 200.0)[:, 2].all()


# --------------------------------------------------------------------------- clean
FACES = vgrid.faces(np.array([1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 20.0, 30.0, 40.0]))


def wet_levels(depth, faces, hfac_min=0.2, hfac_min_dr=2.0):
    h = bathy.hfac_c(depth, faces, hfac_min, hfac_min_dr)
    return (h > 0).sum(axis=0)


def clean_case():
    ny, nx = 30, 40
    depth = np.full((ny, nx), 50.0)
    wet = np.ones((ny, nx), dtype=bool)
    wet[:, 10:15] = False                       # barrier with two 3-cell channels
    wet[5:8, 10:15] = True
    wet[20:23, 10:15] = True
    wet[:, -1] = False                          # closed east
    wet[25:30, 30:39] = False                   # land block in the SE ...
    wet[27, 33] = True                          # ... with an isolated puddle inside it
    wet[0:5, 30:39] = False
    wet[2, 32] = True                           # puddle touching nothing
    # 2x2 lake connected to the basin only through the diagonal of cell (12, 30)
    wet[12:18, 30:39] = False
    wet[13:15, 31:33] = True
    wet[12, 30] = True                          # basin notch (2 land sides), diagonal to lake cell (13, 31)
    # 1-cell inlet into the northern land block: W, E and N neighbours are land
    wet[25, 35] = True
    depth[~wet] = 0.0
    depth[wet & (np.arange(nx)[None, :] > 30)] = 1.0     # shallow cells must be deepened
    return depth, wet


def test_clean_removes_puddles_lakes_inlets_and_sets_min_depth():
    depth, wet = clean_case()
    out = bathy.clean(depth, wet, FACES, min_depth_m=2.0, min_cells=2, ob_sides=("W", "N", "S"),
                      hfac_min=0.2, hfac_min_dr=2.0)
    w, d, c = out["wet"], out["depth"], out["counts"]
    assert not w[27, 33] and not w[2, 32]            # isolated puddles removed
    assert not w[13:15, 31:33].any()                 # diagonal-only lake removed
    assert not w[25, 35]                             # 1-cell inlet filled
    assert c["removed_disconnected_cells"] >= 6 and c["filled_inlet_cells"] >= 1
    assert c["deepened_cells"] > 0
    # every wet column has at least min_cells wet levels and every land column is zero
    assert wet_levels(d[w], FACES).min() >= 2
    assert np.all(d[~w] == 0.0) and np.all(d[w] >= FACES[2] - 1e-12)
    # all wet cells are 4-connected to an open boundary
    assert bathy.connected_to_ob(w, ("W", "N", "S")).sum() == w.sum()
    # partial cells are already MITgcm-consistent (idempotent)
    assert np.allclose(bathy.partial_cells(d, FACES, 0.2, 2.0), d)
    # no wet cell is left with three land sides
    assert bathy.inlet_cells(w, ("W", "N", "S")).sum() == 0


def test_two_600m_channels_stay_open(tmp_path):
    # land wall 2 km thick at 200 m resolution with two 600 m gaps; one aligned with the grid,
    # one shifted by half a cell
    wall = box(X0 + 2000, Y0 - 1000, X0 + 4000, Y0 + 7000)
    gaps = [box(X0 + 1900, Y0 + 1000, X0 + 4100, Y0 + 1600), box(X0 + 1900, Y0 + 4100, X0 + 4100, Y0 + 4700)]
    land = wall.difference(gaps[0]).difference(gaps[1])
    shp = write_land_shapefile(tmp_path / "wall", list(land.geoms))
    xc = X0 + 100.0 + 200.0 * np.arange(40)
    yc = Y0 + 100.0 + 200.0 * np.arange(30)
    coast = bathy.coast_mask(xc, yc, shp, 200.0)
    coast[:, -1] = True
    depth = np.full(coast.shape, 30.0)
    wet = bathy.make_wet(depth, np.ones(coast.shape), coast)
    out = bathy.clean(np.where(wet, depth, 0.0), wet, FACES, min_depth_m=2.0, hfac_min=0.2, hfac_min_dr=2.0)
    w = out["wet"]
    i_wall = slice(10, 20)
    rows_open = np.where(w[:, i_wall].all(axis=1))[0]
    assert len(rows_open) >= 6
    for lo, hi in ((0, 15), (15, 30)):
        assert ((rows_open >= lo) & (rows_open < hi)).sum() >= 3
    assert w[:, 25:39].all()                           # east basin kept: connected through the gaps


# --------------------------------------------------------------------------- PD1 force regions
def force_grid():
    """10 x 12 grid at 200 m: sea in the west (i < 4), a lagoon (i >= 7, j >= 5) joined to the sea
    by a channel on row 7, and a second unconnected pond (row 1, i = 9) inside the lagoon polygon."""
    xc = X0 + 100.0 + 200.0 * np.arange(12)
    yc = Y0 + 100.0 + 200.0 * np.arange(10)
    wet = np.zeros((10, 12), bool)
    wet[:, :4] = True
    wet[7, 4:7] = True                       # channel
    wet[5:9, 7:11] = True                    # lagoon
    wet[1, 9] = True                         # separate pond
    g = {"xc": xc, "yc": yc, "wet_fraction": np.where(wet, 1.0, 0.0),
         "depth_mean": np.where(wet, 10.0, 0.0), "tid_mode": np.full(wet.shape, 11)}
    return wet, g


def test_force_land_seed_fills_only_connected_region_in_polygon():
    wet, g = force_grid()
    poly = [[X0 + 1400, Y0], [X0 + 2400, Y0], [X0 + 2400, Y0 + 2000], [X0 + 1400, Y0 + 2000]]   # i >= 7
    region = {"Lagoon": {"seed": [X0 + 1700, Y0 + 1300], "polygon": poly}}
    out, counts = bathy.apply_force_regions(wet, g, FACES, region, None)
    assert not out[5:9, 7:11].any()                      # lagoon closed
    assert out[1, 9]                                     # pond in the polygon, not connected: kept
    assert out[7, 4:7].all() and out[:, :4].all()        # channel and sea untouched
    assert counts["force_land:Lagoon"] == {"cells": 16, "area_km2": pytest.approx(16 * 0.04)}


def test_force_sea_with_measured_tid_and_depth_condition():
    wet, g = force_grid()
    npi = wet.copy()
    npi[:, 2:4] = False                                  # NPI land over IBCAO sea in columns 2-3
    g["tid_mode"][:, 2] = 41                             # column 2 interpolated
    g["depth_mean"][0:3, 3] = 2.5                        # rows 0-2 of column 3 shallower than 2 levels
    region = {"Head": {"seed": [X0 + 100, Y0 + 100], "require_measured_tid": True, "min_depth_cells": 2}}
    out, counts = bathy.apply_force_regions(npi, g, FACES, None, region)
    assert not out[:, 2].any()                           # no measured data: NPI keeps winning
    assert not out[0:3, 3].any() and out[3:, 3].all()    # measured and deeper than faces[2] only
    assert counts["force_sea:Head"]["cells"] == 7
    out2, _ = bathy.apply_force_regions(npi, g, FACES, None, {"Head": {"seed": [X0 + 100, Y0 + 100]}})
    assert out2[:, 2:4].all()                            # unconditional: the whole IBCAO-sea fill


# --------------------------------------------------------------------------- open boundaries
def test_open_boundaries_zero_normal_gradient_and_wet_corners():
    rng = np.random.default_rng(0)
    ny, nx, n = 40, 50, 10
    depth = 30.0 + 20.0 * rng.random((ny, nx))
    depth[:, -3:] = 0.0                                 # closed east
    out = bathy.open_boundaries(depth, depth > 0, n_sponge=n, ob_sides=("W", "N", "S"))
    # W: in the sponge, every row is constant along x (d h / d x = 0) and equals column n
    assert np.all(out[:, :n + 1] == out[:, [n]])
    # S and N: every column is constant along y in the sponge
    assert np.all(out[:n + 1, :] == out[[n], :])
    assert np.all(out[-n - 1:, :] == out[[-n - 1], :])
    # interior unchanged
    assert np.array_equal(out[n + 1:-n - 1, n + 1:], depth[n + 1:-n - 1, n + 1:])
    # the open corners (SW, NW) are water
    assert (out[:n, :n] > 0).all() and (out[-n:, :n] > 0).all()


def test_open_boundaries_rejects_land_in_open_corner():
    depth = np.full((30, 30), 40.0)
    depth[10:13, 10:13] = 0.0                           # island at the inner edge of the SW corner
    with pytest.raises(ValueError, match="corner"):
        bathy.open_boundaries(depth, depth > 0, n_sponge=10, ob_sides=("W", "N", "S"))


def test_smooth_if_needed_default_is_identity():
    depth = np.array([[10.0, 80.0], [20.0, 0.0]])
    out, info = bathy.smooth_if_needed(depth, depth > 0, None)
    assert np.array_equal(out, depth) and info["applied"] is False


def test_smooth_if_needed_limits_step():
    depth = np.array([[2.0, 2.0, 2.0, 2.0], [2.0, 2.0, 100.0, 2.0]]) * 1.0 + 1.0
    out, info = bathy.smooth_if_needed(depth, depth > 0, 2, dz_faces=FACES)
    k = np.searchsorted(FACES, out - 1e-9)
    assert info["applied"] is True and info["max_change_m"] > 0
    assert np.abs(np.diff(k, axis=1)).max() <= 2 and np.abs(np.diff(k, axis=0)).max() <= 2


# --------------------------------------------------------------------------- partial cells
def test_partial_cells_follow_mitgcm_rule():
    faces = vgrid.faces(np.array([2.0, 4.0, 10.0]))          # faces 0, 2, 6, 16
    # hFacMnSz (ini_masks_etc.F:108) = max(0.2, min(2/drF, 1)) -> 1.0, 0.5, 0.2 for the three levels
    d = np.array([2.5, 3.5, 6.5, 7.5, 16.0, 30.0])
    out = bathy.partial_cells(d, faces, 0.2, 2.0)
    # 2.5: level 2 hFac 0.125 < 0.5/2 -> 0 -> depth 2.0; 3.5: hFac 0.375 >= 0.25 -> max(0.375, 0.5) -> 4.0
    # 6.5: level 3 hFac 0.05 < 0.1 -> 6.0; 7.5: hFac 0.15 >= 0.1 -> max(0.15, 0.2) = 0.2 -> 8.0
    # 16: full column; 30: clipped to the bottom face
    assert np.allclose(out, [2.0, 4.0, 6.0, 8.0, 16.0, 16.0])


# --------------------------------------------------------------------------- vertical grid
def test_vertical_grid_in_cases_yaml():
    delR = np.asarray(CFG["delR"], dtype=float)
    assert delR[0] == pytest.approx(1.5)
    assert 36 <= delR.size <= 44
    assert np.all(np.diff(delR) >= -1e-12)                     # dz monotonic non-decreasing
    assert delR.sum() == pytest.approx(CFG["vgrid"]["H_m"], abs=1e-6)
    # the stored list is what the generator produces from the stored recipe
    v = CFG["vgrid"]
    regen = vgrid.stretched_delR(v["nz"], v["dz_top"], v["H_m"], v["p"])
    assert np.allclose(regen, delR, atol=1e-6)
    # fine and nearly uniform through the pycnocline: >= 20 levels and dz <= 6 m in the top 60 m,
    # and smooth stretching below (ratio of neighbouring thicknesses <= 1.15)
    f = vgrid.faces(delR)
    assert (f[1:] <= 60.0).sum() >= 20 and delR[f[1:] <= 60.0].max() <= 6.0
    assert (delR[1:] / delR[:-1]).max() <= 1.15


# --------------------------------------------------------------------------- binary
def test_write_bathy_roundtrip(tmp_path):
    depth = np.array([[0.0, 12.5, 30.0], [5.0, 0.0, 100.25]])
    p = tmp_path / "bathy.bin"
    bathy.write_bathy(p, depth, depth > 0)
    back = read_bin(p, depth.shape)
    assert np.array_equal(back, -depth)
    assert np.all(back[depth > 0] < 0) and np.all(back[depth == 0] == 0)
    assert p.stat().st_size == depth.size * 8
    assert np.fromfile(p, dtype=">f8")[1] == -12.5            # big-endian, x fastest, j = 0 first


# --------------------------------------------------------------------------- real IBCAO
DATA = PKG / "data_external"
HAVE_DATA = (DATA / "ibcao" / "ibcao_v5_2_2026_depth_100m_svalbard.tif").exists() and \
    (DATA / "npi_coast" / "extracted" / "NP_S100_SHP" / "S100_Land_f.shp").exists()
real = pytest.mark.skipif(not HAVE_DATA, reason="IBCAO v5.2 window or NPI S100 coastline not in data_external/")


@pytest.fixture(scope="module")
def pipeline():
    if not HAVE_DATA:
        pytest.skip("external data missing")
    src, grid = bathy.build(CFG, DATA)
    return src, grid, bathy.report(src, grid)


@real
def test_real_both_sounds_open_three_cells(pipeline):
    _, _, rep = pipeline
    for name in ("Akselsundet", "Mariasundet"):
        s = rep["sounds"][name]
        assert s["open_with_min_width_3"], name
        assert s["grid"]["wet_cells_on_section"] >= 3, name


@real
def test_real_all_wet_connected_and_min_cells(pipeline):
    _, grid, _ = pipeline
    w = grid["wet"]
    assert bathy.connected_to_ob(w, tuple(CFG["ob_sides"])).sum() == w.sum()
    faces = vgrid.faces(np.asarray(CFG["delR"]))
    assert wet_levels(grid["depth"][w], faces, CFG["hFacMin"], CFG["hFacMinDr"]).min() >= 2
    assert grid["depth"][w].max() <= faces[-1] + 1e-9


@real
def test_real_controls_vs_stoylen_fer(pipeline):
    _, _, rep = pipeline
    assert rep["controls"], "no controls"
    for c in rep["controls"]:
        if not c["pass"]:
            cl = c["classification"] or ""
            assert re.match(r"^\([a-d]\)", cl) and "post hoc" in cl, c


@real
def test_real_vallunden_closed_and_force_regions_counted(pipeline):
    _, grid, rep = pipeline
    xc, yc = grid["xc"], grid["yc"]
    j, i = bathy._cell_of(CFG["bathy"]["vallunden_point"], xc, yc, CFG["dx_m"])
    assert not grid["wet"][j, i]                          # ruling R6
    assert rep["vallunden"]["forced_land"]
    assert rep["force_regions"]["force_land:Vallunden"]["cells"] > 0
    fs = rep["force_regions"]["force_sea:VanKeulen_head"]
    assert fs["cells"] >= 0
    forced = grid["wet_forced"] & ~grid["wet_npi"]
    assert np.isin(grid["tid_mode"][forced], list(range(10, 18))).all()


@real
def test_real_sponge_zero_normal_gradient(pipeline):
    _, grid, _ = pipeline
    bathy.check_sponge(grid["depth"], CFG["n_sponge"], tuple(CFG["ob_sides"]))


@real
def test_real_inner_fjord_volume_within_5pct(pipeline):
    _, _, rep = pipeline
    v = rep["inner_fjord"]
    assert abs(v["volume_grid_m3"] / v["volume_ibcao_m3"] - 1.0) < 0.05


@real
def test_real_box_max_depth_bin_matches_vgrid(pipeline):
    _, grid, _ = pipeline
    H = CFG["vgrid"]["H_m"]
    dmax = grid["depth"][grid["wet"]].max()               # final model sea, after the sponge copy
    assert dmax <= H < dmax + CFG["vgrid"]["bin_m"]


@real
def test_real_report_json_serialisable(pipeline):
    _, _, rep = pipeline
    json.dumps(rep)
