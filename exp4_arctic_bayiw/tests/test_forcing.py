"""Tests for exp4_arctic_bayiw.gendata.forcing (Task 4): the wind pulse of the fjord runs.

The pulse must have the thesis shape (thesis input/wind_patch.py, make_meridional_wind_patch,
read from the frozen copy in thesis_copy/), the stress direction must follow the UTM
convention (Review Focus 2: angle in degrees from east, counter-clockwise) and the peak must
be 0.1 N m^-2.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from exp4_arctic_bayiw.gendata import forcing
from exp4_arctic_bayiw.tests.thesis_helpers import thesis_wind_patch_module

PKG = Path(__file__).resolve().parents[1]
SPEC = yaml.safe_load((PKG / "cases.yaml").read_text())


def _thesis_module():
    return thesis_wind_patch_module()


def test_pulse_matches_thesis_make_meridional_wind_patch():
    """Same half Gaussian in time as the thesis function, to < 1e-15 N m^-2, on the thesis
    time axis (linspace(0, 60 h, 180)), at the patch centre (spatial factor 1)."""
    mod = _thesis_module()
    out = mod.make_meridional_wind_patch(nt=180, end_time=60, max_tau=0.10, sigma_t=1.0,
                                         center_t=2.5, nx=9, ny=9, cx=4, cy=4, radius=3,
                                         write_bin=False, make_plots=False)
    tau = out[0] if isinstance(out, tuple) else out
    thesis_t = np.asarray(tau)[4, 4, :]                     # (nx, ny, nt) at the centre
    t_h = np.linspace(0, 60, 180)
    mine = forcing.pulse(t_h * 3600.0, tau_max=0.1, sigma_h=1.0, center_h=2.5)
    assert mine.shape == thesis_t.shape
    assert np.max(np.abs(mine - thesis_t)) < 1e-15
    assert np.all(mine[t_h > 2.5] == 0.0)


def test_pulse_defaults_and_peak():
    t = np.array([0.0, 9000.0, 9000.0 + 1e-6, 20000.0])
    p = forcing.pulse(t)
    assert p[1] == 0.1                                      # exp(0) at the centre, 2.5 h
    assert p[2] == 0.0 and p[3] == 0.0                      # zero after the centre
    assert p[0] == pytest.approx(0.1 * np.exp(-2.5 ** 2 / 2), rel=1e-15)


def test_record_times_follow_mitgcm_periodic_forcing():
    """Record n (1-based) is valid at (n - 1/2) * period (eesupp/src/get_periodic_interval.F:
    locTime = t - period/2, tRec1 = 1 + INT(locTime/period))."""
    t = forcing.record_times(4, 1200.0)
    assert np.array_equal(t, [600.0, 1800.0, 3000.0, 4200.0])


def test_cycle_records_keep_the_wrap_out_of_the_forced_stage():
    """With n = forced/period + 1 records, the forced stage ends at (n - 1) * period, before
    the time (n - 1/2) * period at which MITgcm starts to interpolate back to record 1."""
    n = forcing.cycle_records(60.0, 1200.0)
    assert n == 181
    tr = forcing.record_times(n, 1200.0)
    p = forcing.pulse(tr)
    assert p[-1] == 0.0 and p[-2] == 0.0
    assert 60.0 * 3600 < tr[-1]


def test_uniform_stress_west_wind_sign():
    """Review Focus 2: angle 180 deg is a stress towards the west: tau_x < 0, tau_y == 0."""
    t = forcing.record_times(20, 1200.0)
    tx, ty = forcing.uniform_stress(t, 180.0, (3, 4))
    assert tx.shape == ty.shape == (20, 3, 4)
    k = int(np.argmax(np.abs(tx[:, 0, 0])))
    assert tx[k, 0, 0] < 0
    assert np.all(ty == 0.0)
    assert np.all(tx <= 0.0)


def test_uniform_stress_cardinal_and_spatially_uniform():
    t = forcing.record_times(20, 1200.0)
    for ang, sx, sy in ((0.0, 1, 0), (90.0, 0, 1), (270.0, 0, -1), (-90.0, 0, -1)):
        tx, ty = forcing.uniform_stress(t, ang, (2, 2))
        k = 7                                                  # record at 9000 s, the peak
        assert np.sign(tx[k, 0, 0]) == sx and np.sign(ty[k, 0, 0]) == sy
        assert np.ptp(tx, axis=(1, 2)).max() == 0.0 and np.ptp(ty, axis=(1, 2)).max() == 0.0


def test_uniform_stress_peak_is_0p1():
    t = forcing.record_times(181, 1200.0)
    for ang in (180.0, 192.412, 282.412):
        tx, ty = forcing.uniform_stress(t, ang, (2, 3))
        mag = np.hypot(tx, ty)
        assert mag.max() == pytest.approx(0.1, abs=1e-15)
        assert int(np.argmax(mag[:, 0, 0])) == 7               # 9000 s = 2.5 h


def test_down_fjord_angle_matches_pd0_vector():
    """PD0 section 7: stress unit vector (-0.977, -0.215), the down-fjord axis of theory.json."""
    th = json.loads((PKG / "figures" / "theory.json").read_text())["wind"]
    ang = SPEC["configs"]["fjord_vmf"]["wind"]["angle_deg"]
    u = (np.cos(np.radians(ang)), np.sin(np.radians(ang)))
    assert np.allclose(u, th["down_fjord_stress_unit_vector_xy"], atol=1e-5)
    assert np.allclose(u, (-0.977, -0.215), atol=5e-4)
    f5 = SPEC["cases"]["F5"]["wind_angle_deg"]
    assert (f5 - ang) % 360 == pytest.approx(90.0)              # 90 deg counter-clockwise
