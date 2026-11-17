"""analysis.hypotheses: H5a / H5c interface of plan Task 6 with the PD1b scales and PD2 inputs."""
import numpy as np
import pytest

from exp5_adjoint_slope.analysis import hypotheses as hy

Y = (np.arange(45) + 0.5) * 10e3


def test_normalize_rms_is_one_over_shelf_and_slope():
    for p in (hy.pattern_uniform(Y), hy.pattern_dipole(Y, 37.5e3)):
        q = hy.normalize_rms(p, Y)
        assert np.sqrt((q[Y < hy.Y_MAX] ** 2).mean()) == pytest.approx(1.0)


def test_dipole_integrates_to_zero():
    # Odd about the break: integral < 1e-6 of its norm on a symmetric domain. On the channel the
    # southern wall truncates it 3.3 half-widths from y0 (tail ~exp(-11)): residual < 1e-4.
    w = 37.5e3
    ys = np.linspace(hy.Y0 - 10 * w, hy.Y0 + 10 * w, 8001)
    d = hy.pattern_dipole(ys, w)
    assert abs(d.sum()) < 1e-6 * np.abs(d).sum()
    yf = np.linspace(0, 450e3, 9001)
    d = hy.pattern_dipole(yf, w)
    assert abs(d.sum()) < 1e-4 * np.abs(d).sum()


def test_h5a_dipole_sensitivity_passes_uniform_does_not():
    assert hy.h5a(hy.pattern_dipole(Y, 37.5e3), Y)["pass"]
    assert not hy.h5a(hy.pattern_uniform(Y), Y)["pass"]


def test_h5c_exponential_half_life():
    lead = np.arange(0.0, 181.0, 10.0)
    r = hy.h5c(np.exp(-lead / 100.0), lead)
    assert r["t_half_days"] == pytest.approx(69.3, abs=1.0) and r["pass"]
