import numpy as np

from exp5_adjoint_slope.analysis import animate as an


def test_frame_times_are_centres_since_t0():
    t0 = 100 * 86400.0
    t = an.frame_times([105 * 86400.0, 110 * 86400.0], t0)
    assert np.allclose(t, [2.5, 7.5])


def test_eke_months_restart_at_the_extension():
    t, e = an.eke_months([1.0, 2.0], [3.0], f1_days=365.0)
    assert np.allclose(t, [15.0, 45.0, 380.0]) and np.allclose(e, [1.0, 2.0, 3.0])
