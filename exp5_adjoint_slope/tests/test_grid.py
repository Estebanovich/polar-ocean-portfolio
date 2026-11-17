"""gendata.grid and the J-line index (Exp. 5 Task 1)."""
import numpy as np
import pytest

from exp5_adjoint_slope.gendata import fields, grid


def test_vertical_spacing_sums_and_grows():
    dz = grid.vertical_spacing(35, 10.0, 3000.0)
    assert dz.shape == (35,) and dz[0] == pytest.approx(10.0)
    assert dz.sum() == pytest.approx(3000.0, rel=1e-12)
    assert np.all(np.diff(dz) > 0)


@pytest.mark.parametrize("nx,ny,dx", [(40, 45, 10e3), (160, 180, 2.5e3)])
def test_j_qb_on_v_face_at_170km(nx, ny, dx):
    # notes/pd1b_config_dettling.md: J line at 170 km, ~920 m isobath, a v face on both grids.
    g = grid.Grid(nx, ny, dx, grid.vertical_spacing())
    j = fields.j_qb(g, 170e3)
    assert g.yg[j] == pytest.approx(170e3)
    assert 900.0 <= fields.depth(np.array([g.yg[j]]))[0] <= 940.0


def test_j_qb_rejects_line_off_grid():
    g = grid.Grid(40, 45, 10e3, grid.vertical_spacing())
    with pytest.raises(ValueError):
        fields.j_qb(g, 95e3)
