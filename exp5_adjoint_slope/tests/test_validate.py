"""analysis.validate: V1-V3 criteria as frozen in notes/pd0_preregistro.md."""
from exp5_adjoint_slope.analysis import validate


def test_v1_fails_if_one_point_over_one_percent():
    pts = [{"rel_err": 1e-4}] * 9 + [{"rel_err": 0.012}]
    r = validate.v1(pts)
    assert r["pass"] is False and r["value"] == 0.012
    assert validate.v1([{"rel_err": 0.009}] * 10)["pass"] is True


def test_v1_needs_ten_points():
    assert validate.v1([{"rel_err": 1e-6}] * 9)["pass"] is False


def test_v2():
    r = validate.v2(1.0, -0.95)
    assert r["pass"] is True and abs(r["value"] - 0.05 / 1.95) < 1e-12
    assert validate.v2(1.0, -0.5)["pass"] is False


def test_v3():
    assert validate.v3(1.25, 1.0)["pass"] is True
    r = validate.v3(1.4, 1.0)
    assert r["pass"] is False and abs(r["value"] - 0.4) < 1e-12


def test_top_points_largest_and_separated():
    import numpy as np
    s = np.zeros((3, 10, 10))
    s[0, 5, 5], s[0, 5, 6], s[2, 1, 1], s[1, 8, 2] = 10.0, 9.0, -8.0, 1e-12
    pts = validate.top_points(s, 2, min_sep=2)
    assert pts == [(0, 5, 5), (2, 1, 1)]                    # (0,5,6) too close to (0,5,5)


def test_top_points_skip_points_without_signal():
    import numpy as np
    s = np.zeros((1, 4, 4))
    s[0, 0, 0], s[0, 3, 3] = 1.0, 1e-11                     # < 1e-10 * max: no signal
    assert validate.top_points(s, 2, min_sep=1) == [(0, 0, 0)]
