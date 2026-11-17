"""analysis.perturb: adjoint-guided perturbation patterns for V2/V3 (PD0) and their xx_* files."""
import numpy as np
import pytest

from exp5_adjoint_slope.analysis import perturb
from pfcommon.mitgcm_io import read_mds


def test_pattern_normalised_to_amplitude():
    s = np.array([[0.5, -2.0], [1.0, 0.0]])
    p = perturb.pattern(s, 0.01)
    assert np.abs(p).max() == pytest.approx(0.01)
    assert p[0, 1] == pytest.approx(-0.01) and p[0, 0] == pytest.approx(0.0025)


def test_predicted_dj_is_inner_product():
    s = np.array([1.0, 2.0, -1.0])
    du = np.array([0.1, 0.0, 0.2])
    assert perturb.predicted_dj(s, du) == pytest.approx(0.1 - 0.2)


def test_write_xx_round_trip(tmp_path):
    a = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)       # (nr, ny, nx) for a 3-D control
    perturb.write_xx(tmp_path, "theta", a)
    b, meta = read_mds(tmp_path / "xx_theta", 0)
    assert np.array_equal(np.asarray(b)[0], a)              # read_mds keeps the record axis
    t = (tmp_path / "xx_theta.0000000000.meta").read_text()
    assert "nDims = [   3 ]" in t and "nrecords = [          1 ]" in t


def test_write_xx_2d_constant_in_time(tmp_path):
    a = np.ones((3, 4))
    perturb.write_xx(tmp_path, "fu", a)
    t = (tmp_path / "xx_fu.0000000000.meta").read_text()
    assert "nDims = [   2 ]" in t and "nrecords = [          1 ]" in t
