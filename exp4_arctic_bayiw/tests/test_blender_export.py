"""Tests of the Blender scene export (analysis/blender_export.py, Task 10)."""
from __future__ import annotations

import json

import numpy as np
import pytest

from exp4_arctic_bayiw.analysis import blender_export as be


def _column(nz=60, dz=2.0):
    RC = -(np.arange(nz) + 0.5) * dz
    return RC, np.full(nz, dz)


# ---- isopycnal depth -----------------------------------------------------------------------------

def test_isopycnal_depth_is_exact_on_a_linear_profile():
    RC, drF = _column()
    rho = (1025.0 - 0.01 * RC)[:, None, None] * np.ones((1, 2, 3))      # 0.01 kg m^-3 per m
    hfac = np.ones_like(rho)
    for zt in (-7.3, -31.0, -100.0):
        sigma = 1025.0 - 0.01 * zt
        z = be.isopycnal_depth(rho, RC, hfac, sigma)
        assert z.shape == (2, 3)
        assert np.allclose(z, zt, atol=1e-9)


def test_isopycnal_depth_interpolates_between_the_level_centres_of_a_step():
    # two-layer column: 1025 above level 10, 1026 from level 10 down; sigma half-way sits half-way
    # between the two centres
    RC, drF = _column()
    rho = np.where(np.arange(60) < 10, 1025.0, 1026.0)[:, None, None]
    z = be.isopycnal_depth(rho, RC, np.ones_like(rho), 1025.5)
    assert z[0, 0] == pytest.approx(0.5 * (RC[9] + RC[10]))
    z = be.isopycnal_depth(rho, RC, np.ones_like(rho), 1025.25)
    assert z[0, 0] == pytest.approx(RC[9] + 0.25 * (RC[10] - RC[9]))


def test_isopycnal_depth_is_nan_where_outcropped_too_shallow_or_land():
    RC, drF = _column(nz=10, dz=10.0)
    rho = np.repeat((1025.0 - 0.01 * RC)[:, None, None], 4, axis=2)    # 1025.05 .. 1025.95
    hfac = np.ones_like(rho)
    hfac[3:, 0, 1] = 0.0                                                # column 1: 30 m deep
    hfac[:, 0, 2] = 0.0                                                 # column 2: land
    rho[:, 0, 3] += 1.0                                                 # column 3: denser everywhere
    z = be.isopycnal_depth(rho, RC, hfac, 1025.5)                       # sigma at z = -50 m
    assert z[0, 0] == pytest.approx(-50.0)
    assert np.isnan(z[0, 1])          # bottom lighter than sigma: the isopycnal meets the bed
    assert np.isnan(z[0, 2])          # land
    assert np.isnan(z[0, 3])          # outcropped: the top cell already denser than sigma


def test_isopycnal_depth_ignores_dry_levels_under_a_partial_bottom():
    RC, drF = _column(nz=10, dz=10.0)
    rho = (1025.0 - 0.01 * RC)[:, None, None].copy()
    hfac = np.ones_like(rho)
    hfac[6:] = 0.0
    rho[6:] = 0.0                      # MITgcm writes 0 (+ rhoConst offset) below the bed
    assert np.isnan(be.isopycnal_depth(rho, RC, hfac, 1025.8)[0, 0])


# ---- displacement sign ---------------------------------------------------------------------------

def test_displacement_is_positive_when_the_isopycnals_are_lifted():
    RC, drF = _column()
    d = 3.0
    rho0 = (1025.0 - 0.01 * RC)[:, None, None]
    rho1 = (1025.0 - 0.01 * (RC - d))[:, None, None]       # rho1(z) = rho0(z - d): lifted by d
    sigma = 1025.0 + 0.01 * 40.0
    z0 = be.isopycnal_depth(rho0, RC, np.ones_like(rho0), sigma)
    z1 = be.isopycnal_depth(rho1, RC, np.ones_like(rho1), sigma)
    zeta = be.displacement(z1, z0)
    assert zeta[0, 0] == pytest.approx(+d)
    assert be.displacement(z0, z1)[0, 0] == pytest.approx(-d)


# ---- the pycnocline rule -------------------------------------------------------------------------

def test_max_n2_sigma_below_the_mixed_layer():
    # mixed top 4 m, a strong step at ~11 m and a weaker one at ~31 m; with min depth 5 m the
    # step at 11 m is chosen and sigma is the mean of the two densities around it
    RC, drF = _column(nz=30, dz=2.0)
    rho_b = np.full(30, 1025.0)
    rho_b[5:] += 0.5                   # step between levels 4 (z = -9) and 5 (z = -11)
    rho_b[15:] += 0.2                  # step between levels 14 and 15
    rho_b[0] -= 1.0                    # a surface step stronger than both, inside the top 5 m
    s = be.max_n2_sigma(rho_b, RC, min_depth_m=5.0, rho0=1025.0)
    assert s["k_upper"] == 4 and s["k_lower"] == 5
    assert s["sigma"] == pytest.approx(0.5 * (rho_b[4] + rho_b[5]))
    assert s["z_m"] == pytest.approx(-10.0)
    s0 = be.max_n2_sigma(rho_b, RC, min_depth_m=0.0, rho0=1025.0)
    assert s0["k_upper"] == 0                                     # without the cut: the surface step


def test_max_n2_sigma_respects_a_maximum_depth():
    RC, drF = _column(nz=30, dz=2.0)
    rho_b = 1025.0 + np.where(np.arange(30) >= 20, 1.0, 0.0) + np.where(np.arange(30) >= 8, 0.1, 0.0)
    s = be.max_n2_sigma(rho_b, RC, min_depth_m=2.0, rho0=1025.0, max_depth_m=30.0)
    assert (s["k_upper"], s["k_lower"]) == (7, 8)


# ---- the inertial filter and its edge drop ------------------------------------------------------

def test_filtered_colour_drops_the_frames_within_half_an_inertial_period_of_the_ends():
    f = 1.4266e-4
    Tf = 2 * np.pi / f
    t = np.arange(0, 60 * 3600 + 1, 1200.0)
    x = np.sin(2 * np.pi * t / (8 * 3600.0))[:, None, None] * np.ones((1, 2, 2)) + 5.0 + 1e-5 * t[:, None, None]
    tk, xf = be.filtered_colour(x.astype(np.float32), t, f)
    assert tk[0] >= t[0] + Tf / 2 - 1e-6 and tk[0] - 1200.0 < t[0] + Tf / 2
    assert tk[-1] <= t[-1] - Tf / 2 + 1e-6 and tk[-1] + 1200.0 > t[-1] - Tf / 2
    assert xf.shape == (tk.size, 2, 2)
    assert abs(float(np.mean(xf))) < 0.1                     # constant and trend removed
    # the same as the Task 9 filter
    from exp4_arctic_bayiw.analysis import animate as an
    tk9, x9 = an.inertial_filter(x.astype(np.float32), t, Tf)
    assert np.array_equal(tk9, tk) and np.allclose(x9, xf)


def test_select_frames_picks_the_nearest_kept_time_and_reports_it():
    t = np.arange(22800.0, 200000.0, 1200.0)
    sel = be.select_frames(t, onset_s=0.0, hours=[6.0, 24.0])
    assert sel[0]["t_s"] == t[0] and sel[0]["requested_h"] == 6.0
    assert sel[0]["shifted"] is True                     # 6 h is inside the dropped edge
    assert sel[1]["t_s"] == pytest.approx(86400.0) and sel[1]["shifted"] is False


# ---- labels --------------------------------------------------------------------------------------

def test_time_labels_in_hours_and_inertial_periods():
    f = 1.4266e-4
    lab = be.time_label(24 * 3600.0, 0.0, f)
    assert lab == "t = 24.0 h = 1.96 inertial periods"
    assert be.time_label(6.3333 * 3600.0, 0.0, f) == "t = 6.3 h = 0.52 inertial periods"
    assert be.time_label(48 * 3600.0 + 1800.0, 1800.0, f).startswith("t = 48.0 h = 3.92")


# ---- colour scale and colour map -----------------------------------------------------------------

def test_symmetric_limit_is_the_nice_ceiling_of_the_percentile_over_all_frames():
    a = {"F3": np.linspace(-2.7, 2.7, 1001), "F2": np.array([0.1, np.nan, -0.2])}
    L = be.symmetric_limit(a, pctl=100.0)
    assert L == 3.0


def test_lut_colour_maps_the_ends_and_the_centre_and_greys_nan():
    lut = np.array([[0, 0, 1], [1, 1, 1], [1, 0, 0]], float)
    c = be.lut_colour(np.array([-5.0, 0.0, 5.0, 99.0, np.nan]), 5.0, lut, nan_rgb=(0.5, 0.5, 0.5))
    assert np.allclose(c[0], [0, 0, 1]) and np.allclose(c[1], [1, 1, 1]) and np.allclose(c[2], [1, 0, 0])
    assert np.allclose(c[3], [1, 0, 0])                          # clipped, not wrapped
    assert np.allclose(c[4], [0.5, 0.5, 0.5])


# ---- heightfield mesh ----------------------------------------------------------------------------

def test_quad_faces_only_where_all_four_corners_are_valid():
    valid = np.array([[1, 1, 1],
                      [1, 1, 0]], bool)
    f = be.quad_faces(valid)
    # vertex index = j * nx + i; the quad (0,0)-(1,1) is valid, (0,1)-(1,2) is not
    assert f.shape == (1, 4)
    assert sorted(f[0].tolist()) == [0, 1, 3, 4]


def test_export_metadata_round_trip(tmp_path):
    meta = be.Meta(shot="A", runs=["F1"], sigma={"F1": {"sigma": 1024.8, "rule": "x"}}, ve=10.0, ve_land=2.0,
                   colour_limit_m=3.0, colour_field="filt", onset_s=0.0, f={"F1": 1.4266e-4},
                   frames=[{"t_s": 86400.0, "label": "t = 24.0 h"}], crop={"j0": 10, "j1": 290,
                                                                          "i0": 10, "i1": 500})
    p = tmp_path / "meta.json"
    be.write_meta(meta, p)
    d = json.loads(p.read_text())
    assert d["ve"] == 10.0 and d["sigma"]["F1"]["sigma"] == 1024.8 and d["colour_field"] == "filt"
    assert "lut_diverging" in d and len(d["lut_diverging"]) == 256


# ---- render helper (numpy only) ------------------------------------------------------------------

def test_terrain_height_uses_the_land_and_sea_exaggerations_separately():
    from exp4_arctic_bayiw.analysis import blender_render as br
    ter = {"land": np.array([[True, False]]), "elev_m": np.array([[500.0, 0.0]]),
           "depth_m": np.array([[0.0, 100.0]])}
    z, rgb = br.terrain_height_rgb(ter, ve=20.0, ve_land=2.0)
    assert z[0, 0] == pytest.approx(1.0) and z[0, 1] == pytest.approx(-2.0)      # km
    assert rgb.shape == (1, 2, 3)


# ---- full render: open-boundary sponge mask and render-time saturation ---------------------------

DATA_OBCS = """ &OBCS_PARM01
 OB_Iwest=6*1,
 OB_Jnorth=8*-1,
 OB_Jsouth=8*1,
 useOBCSsponge=.TRUE.,
 &
 &OBCS_PARM03
 spongeThickness=2,
 &
"""


def test_sponge_mask_from_data_obcs_covers_the_boundary_and_sponge_cells():
    m = be.sponge_mask(DATA_OBCS, (8, 8))
    exp = np.zeros((8, 8), bool)
    # the OB row plus spongeThickness = 2 cells inward (pkg/obcs/obcs_sponge.F:68-69, jsl = 1..T)
    exp[:, :3] = True                     # west
    exp[:3, :] = True                     # south
    exp[-3:, :] = True                    # north
    assert np.array_equal(m, exp)
    assert not m[3:5, 3:].any()           # no east boundary in this namelist


def test_sponge_mask_is_empty_without_the_sponge():
    txt = DATA_OBCS.replace("useOBCSsponge=.TRUE.", "useOBCSsponge=.FALSE.")
    m = be.sponge_mask(txt, (6, 8))
    exp = np.zeros((6, 8), bool)
    exp[:, 0] = exp[0, :] = exp[-1, :] = True        # the open-boundary rows themselves only
    assert np.array_equal(m, exp)


def test_saturation_stats_mean_and_worst_frame_per_run():
    a = np.zeros((3, 2, 2))
    a[1, 0, 0] = 2.0                       # frame 1: 1 of 4 pixels saturated
    a[2, :, 0] = -1.5                      # frame 2: 2 of 4 (|x| >= L counts)
    a[2, 1, 1] = np.nan                    # NaN pixels are not counted
    s = be.saturation_stats({"F3": a}, 1.5)
    assert s["F3"]["worst"] == pytest.approx(2 / 3)
    assert s["F3"]["mean"] == pytest.approx((0 + 0.25 + 2 / 3) / 3)
    assert s["F3"]["worst_index"] == 2
