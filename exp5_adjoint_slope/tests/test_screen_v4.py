"""analysis.screen_v4: the registered V4 screening criterion (notes/v4_revision_posthoc.md)."""
import numpy as np
import pytest

from exp5_adjoint_slope.analysis import screen_v4 as sv


def test_shelf_mean_is_volume_weighted_and_masked():
    nr, ny, nx = 2, 4, 3
    th = np.zeros((nr, ny, nx))
    th[0] = 1.0
    hfc = np.ones((nr, ny, nx))
    hfc[1, 0] = 0.0                       # land under the first row
    drf = np.array([10.0, 30.0])
    rac = np.full((ny, nx), 2.0)
    yc = np.array([10e3, 30e3, 90e3, 110e3])
    got = sv.shelf_mean(th, hfc, drf, rac, yc, y_max=80e3)
    vol_top = 2 * 3 * 10.0 * 2.0          # two shelf rows, top level
    vol_bot = 1 * 3 * 30.0 * 2.0          # only the second row is wet at depth
    assert got == pytest.approx(vol_top / (vol_top + vol_bot))


def test_verdict_needs_both_conditions():
    assert sv.passes(shelf_theta=-1.2, J=0.01)
    assert not sv.passes(shelf_theta=-0.9, J=0.5)
    assert not sv.passes(shelf_theta=-1.5, J=-0.01)


def test_bottom_theta_skips_land_columns():
    nr, ny, nx = 3, 3, 2
    th = np.full((nr, ny, nx), 5.0)
    th[1, 1:] = -1.0                      # bottom wet level of rows 1-2
    hfc = np.ones((nr, ny, nx))
    hfc[:, 0] = 0.0                       # row 0 is land (south wall)
    hfc[2, 1:] = 0.0
    yc = np.array([5e3, 15e3, 25e3])
    assert sv.shelf_bottom_mean(th, hfc, yc) == pytest.approx(-1.0)
