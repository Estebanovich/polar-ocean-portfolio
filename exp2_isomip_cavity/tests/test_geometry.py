"""Tests for exp2_isomip_cavity/gendata/geometry.py (Task 1; Review Focus 1-2)."""
from __future__ import annotations

import numpy as np
import pytest
from scipy import ndimage

from exp2_isomip_cavity.gendata import geometry as geo

from exp2_testdata import REAL_GEOM, synthetic_fields, write_geom

DZ = np.full(36, 20.0)


def grid(dx_km):
    dx = dx_km * 1e3
    xc = 320e3 + dx * (np.arange(int(480 / dx_km)) + 0.5)
    yc = dx * (np.arange(int(80 / dx_km)) + 0.5)
    return xc, yc


# ---- analytic reference bathymetry (AD16 Eq. 1-3) -------------------------------

def test_bathymetry_reference_point():
    # B0 + B2 + B4 + B6 = -585.46 m at x = xbar; By(40 km) ~ 2 dc / (1 + e^12) ~ 0.006 m
    assert geo.bathymetry(300e3, 40e3) == pytest.approx(-585.46, abs=0.5)


def test_bathymetry_never_below_zb_deep():
    X, Y = np.meshgrid(np.linspace(0, 800e3, 801), np.linspace(0, 80e3, 161))
    zb = geo.bathymetry(X, Y)
    assert zb.min() >= -720.0
    assert zb.shape == X.shape


def test_bathymetry_symmetric_in_y():
    X, Y = np.meshgrid(np.linspace(320e3, 800e3, 97), np.linspace(0, 40e3, 41))
    np.testing.assert_allclose(geo.bathymetry(X, Y), geo.bathymetry(X, 80e3 - Y), atol=1e-9)


@pytest.mark.skipif(not REAL_GEOM.exists(), reason="Ocean1_input_geom_v1.01.nc not present")
def test_file_bed_consistent_with_analytic():
    """Sanity check that the right field is read with the right orientation (ruling R8);
    not a physics criterion. The model uses the file's bed, not the analytic one (R7).

    Asserted: over the trough interior (|y - 40 km| < 16 km, x < 640 km) median
    |file - analytic| < 0.05 m, p99 < 5 m, fraction(|diff| < 0.5 m) >= 0.90; over the whole
    domain Pearson r(file, analytic) > 0.99 (orientation guard).

    Measured (2026-11-08, 1 km centres):
    - trough interior (10240 cells): median 0.0074 m, p99 1.79 m, fraction 0.921, max 3.36 m;
    - x < 640 km (25600 cells): median 0.2346 m, p99 15.38 m, fraction 0.642, max 90.08 m;
      the side walls differ by up to ~90 m (the file is BISICLES' interpolated bed);
    - whole domain (38400 cells): median 0.0071 m, p99 17.05 m, max 103.92 m (the file is
      flat at -720 m beyond x = 640 km, the analytic walls are not until ~650 km);
    - Pearson r = 0.99966 (with x reversed: -0.047).
    """
    import netCDF4
    with netCDF4.Dataset(REAL_GEOM) as nc:
        assert nc.variables["bedrockTopography"].dimensions == ("y", "x")
        x, y = np.asarray(nc["x"][:]), np.asarray(nc["y"][:])
        b = np.asarray(nc["bedrockTopography"][:], dtype=float)
    X, Y = np.meshgrid(x, y)
    a = geo.bathymetry(X, Y)
    d = np.abs(a - b)[(np.abs(Y - 40e3) < 16e3) & (X < 640e3)]
    assert np.median(d) < 0.05
    assert np.percentile(d, 99) < 5.0
    assert (d < 0.5).mean() >= 0.90
    assert np.corrcoef(a.ravel(), b.ravel())[0, 1] > 0.99


# ---- file reading: block averages, calving, coordinate checks --------------------

@pytest.mark.parametrize("dx_km", [2, 4])
def test_load_draft_block_average_and_calving(synth_geom, dx_km):
    f = synthetic_fields()
    n = dx_km
    lower = np.where(f["upperSurface"] - f["lowerSurface"] < 100.0, 0.0, f["lowerSurface"])
    want = lower.reshape(80 // n, n, 480 // n, n).mean(axis=(1, 3))
    xc, yc = grid(dx_km)
    got = geo.load_draft(str(synth_geom), xc, yc)
    np.testing.assert_allclose(got, want, rtol=0, atol=1e-10)
    # the 80 m thick strip 630-640 km has calved (draft 0); just upstream it has not
    strip = (xc > 630e3) & (xc < 640e3)
    assert np.all(got[:, strip] == 0.0)
    assert np.all(got[:, (xc > 600e3) & (xc < 628e3)] < 0.0)


@pytest.mark.parametrize("dx_km", [2, 4])
def test_load_bed_same_block_average(synth_geom, dx_km):
    f = synthetic_fields()
    n = dx_km
    want = f["bedrockTopography"].reshape(80 // n, n, 480 // n, n).mean(axis=(1, 3))
    xc, yc = grid(dx_km)
    np.testing.assert_allclose(geo.load_bed(str(synth_geom), xc, yc), want, rtol=0, atol=1e-10)


def test_load_draft_rejects_transposed_axes(tmp_path):
    p = write_geom(tmp_path / "t.nc", synthetic_fields(), dims=("x", "y"))
    with pytest.raises(ValueError, match="dimensions"):
        geo.load_draft(str(p), *grid(2))


def test_load_draft_rejects_km_coordinates(tmp_path):
    p = write_geom(tmp_path / "km.nc", synthetic_fields(), xscale=1e-3)
    with pytest.raises(ValueError, match="320-800 km"):
        geo.load_draft(str(p), *grid(2))


def test_load_draft_rejects_short_domain(tmp_path):
    f = synthetic_fields()
    g = {k: (v[:, :400] if np.ndim(v) == 2 else v) for k, v in f.items()}
    g["x"] = f["x"][:400]
    p = write_geom(tmp_path / "short.nc", g)
    with pytest.raises(ValueError, match="320-800 km"):
        geo.load_draft(str(p), *grid(2))


def test_load_draft_rejects_misaligned_target(synth_geom):
    xc, yc = grid(2)
    with pytest.raises(ValueError):
        geo.load_draft(str(synth_geom), xc + 500.0, yc)
    with pytest.raises(ValueError):                       # 3 km does not divide 80 km
        geo.load_draft(str(synth_geom), xc[:160], 1.5e3 + 3e3 * np.arange(26))


# ---- minimum water column and connectivity ---------------------------------------

def _geometry(dx_km, s, path):
    xc, yc = grid(dx_km)
    bed, draft = geo.load_bed(path, xc, yc), geo.load_draft(path, xc, yc)
    fl = geo.load_floating(path, xc, yc)
    return xc, yc, bed, draft, geo.apply_scale_and_mask(bed, draft, s, DZ, xc_m=xc, floating=fl)


def _paths(synth_geom):
    return [str(synth_geom)] + ([str(REAL_GEOM)] if REAL_GEOM.exists() else [])


@pytest.mark.parametrize("dx_km", [2, 4])
@pytest.mark.parametrize("s", [0.8, 1.0, 1.2])
def test_wet_columns_have_two_cells_and_reach_open_ocean(synth_geom, dx_km, s):
    for path in _paths(synth_geom):
        xc, yc, bed, draft, (b2, d2, wet) = _geometry(dx_km, s, path)
        assert wet.any()
        thick = d2 - b2
        assert np.all(thick[wet] >= 2 * 20.0 - 1e-9)           # >= 2 cells of 20 m
        assert np.all((b2[~wet] == 0) & (d2[~wet] == 0))        # dry columns are land
        lab, n = ndimage.label(wet)                             # 4-neighbour components
        for k in range(1, n + 1):
            assert np.any((lab == k) & (xc[None, :] >= 790e3)), f"isolated wet region {k}"


@pytest.mark.parametrize("dx_km", [2, 4])
def test_deeper_draft_never_adds_wet_cells(synth_geom, dx_km):
    for path in _paths(synth_geom):
        n1 = _geometry(dx_km, 1.0, path)[-1][2].sum()
        n12 = _geometry(dx_km, 1.2, path)[-1][2].sum()
        assert n12 <= n1


def test_min_column_threshold_exact():
    bed = np.full((1, 4), -100.0)
    draft = np.array([[-60.0 - 1e-6, -60.0, 0.0, 0.0]])  # columns: 39.999999 m, 40 m, 100 m, 100 m
    b2, d2, wet = geo.apply_scale_and_mask(bed, draft, 1.0, DZ)
    assert wet.tolist() == [[False, True, True, True]]


def test_isolated_lake_is_dried():
    bed = np.full((3, 6), -300.0)
    draft = np.zeros((3, 6))
    draft[:, 1:3] = -300.0                    # grounded barrier at i = 1, 2 ...
    draft[1, 0] = -100.0                      # ... cutting off a wet pocket at i = 0
    b2, d2, wet = geo.apply_scale_and_mask(bed, draft, 1.0, DZ)
    assert not wet[:, :3].any() and wet[:, 3:].all()


def test_columns_without_a_shared_level_are_not_connected():
    # column 0 is wet only from -700 to -640 m, column 1 only from -100 to 0 m: no shared level
    bed = np.array([[-700.0, -100.0, -100.0]])
    draft = np.array([[-640.0, 0.0, 0.0]])
    b2, d2, wet = geo.apply_scale_and_mask(bed, draft, 1.0, DZ)
    assert wet.tolist() == [[False, True, True]]


def test_scale_multiplies_floating_draft_only():
    bed = np.full((1, 4), -600.0)
    draft = np.array([[-300.0, -200.0, -100.0, 0.0]])
    fl = np.array([[False, True, True, False]])
    b2, d2, wet = geo.apply_scale_and_mask(bed, draft, 1.2, DZ, floating=fl)
    np.testing.assert_allclose(d2, [[-300.0, -240.0, -120.0, 0.0]])
    np.testing.assert_allclose(b2, bed)


def test_scale_without_floating_mask_refused():
    with pytest.raises(ValueError, match="floating"):
        geo.apply_scale_and_mask(np.full((1, 2), -600.0), np.full((1, 2), -100.0), 0.8, DZ)


@pytest.mark.parametrize("dx_km", [2, 4])
def test_load_floating_strict_majority(synth_geom, dx_km):
    f = synthetic_fields()
    n = dx_km
    frac = f["floatingMask"].reshape(80 // n, n, 480 // n, n).mean(axis=(1, 3))
    got = geo.load_floating(str(synth_geom), *grid(dx_km))
    np.testing.assert_array_equal(got, frac > 0.5)
    assert got.dtype == bool and got.any() and not got.all()


def test_floating_tie_is_not_floating(tmp_path):
    f = synthetic_fields()
    f["floatingMask"] = np.zeros_like(f["floatingMask"])
    f["floatingMask"][0, 0] = f["floatingMask"][1, 0] = 1.0      # 2 of 4 sub-cells
    f["floatingMask"][0, 2] = f["floatingMask"][1, 2] = f["floatingMask"][0, 3] = 1.0  # 3 of 4
    p = write_geom(tmp_path / "tie.nc", f)
    got = geo.load_floating(str(p), *grid(2))
    assert not got[0, 0] and got[0, 1] and got.sum() == 1


@pytest.mark.parametrize("dx_km", [2, 4])
def test_thinner_draft_does_not_move_grounding_line(synth_geom, dx_km):
    """R9 (a): with s = 0.8 the grounded columns (ice that is not floating) are exactly
    those of s = 1: same set, same draft, same wet/dry state. Columns that become wet
    at s = 0.8 are all floating ones, whose water column was thinner than the minimum at
    s = 1 (real file: the centre-line wet edge moves 4-6 km, one cell, for that reason)."""
    for path in _paths(synth_geom):
        xc, yc = grid(dx_km)
        fl = geo.load_floating(path, xc, yc)
        r1, r08 = _geometry(dx_km, 1.0, path)[-1], _geometry(dx_km, 0.8, path)[-1]
        draft = geo.load_draft(path, xc, yc)
        grounded = ~fl & (draft < 0)
        assert grounded.any()
        np.testing.assert_array_equal(r1[2][grounded], r08[2][grounded])
        np.testing.assert_array_equal(r1[1][grounded], r08[1][grounded])
        dry_ice_1 = ~r1[2] & ~fl & (draft < 0)
        dry_ice_08 = ~r08[2] & ~fl & (draft < 0)
        np.testing.assert_array_equal(dry_ice_1, dry_ice_08)
        newly_wet = r08[2] & ~r1[2]
        assert np.all(fl[newly_wet])
