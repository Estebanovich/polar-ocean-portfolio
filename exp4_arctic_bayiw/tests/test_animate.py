"""Tests of the descriptive animations (analysis/animate.py, Task 9)."""
from __future__ import annotations

import numpy as np
import pytest

from exp4_arctic_bayiw.analysis import animate as an


# ---- mode-1 displacement amplitude on a synthetic stratification ----------------------------------

def _tanh_setup(nz=30):
    drF = np.full(nz, 4.0)
    RF = -np.r_[0.0, np.cumsum(drF)]
    RC = 0.5 * (RF[:-1] + RF[1:])
    rho_b = 1025.0 - 1.5 * np.tanh((RC + 25.0) / 6.0)
    return drF, RC, rho_b


def test_mode1_amplitude_recovers_a_mode1_displacement_in_columns_of_different_depth():
    drF, RC, rho_b = _tanh_setup()
    rho0, g = 1025.0, 9.81
    nz, ny, nx = len(drF), 2, 3
    kb = np.array([[30, 30, 20], [12, 2, 25]])                # wet levels per column (2: too shallow)
    hfac = (np.arange(nz)[:, None, None] < kb[None]).astype(float)
    A = np.array([[1.5, -0.7, 2.0], [0.3, 9.0, -1.1]])        # mode-1 amplitude [m]
    rho = np.tile(rho_b[:, None, None], (1, ny, nx))
    for j in range(ny):
        for i in range(nx):
            sh = an.mode1_shape(rho_b, RC, drF, int(kb[j, i]), rho0, g)
            if sh is None:
                continue
            w, n2 = sh
            # N_b^2 zeta = (g / rho0) rho', zeta = A W_1
            rho[:kb[j, i], j, i] += rho0 / g * n2 * A[j, i] * w
    a = an.mode1_amplitude(rho, rho_b, RC, drF, hfac, rho0, g)
    assert np.isnan(a[1, 1])
    m = ~np.isnan(a)
    np.testing.assert_allclose(a[m], A[m], rtol=1e-10)


def test_mode1_amplitude_of_a_small_uniform_lift_in_linear_stratification():
    # linear rho_b: zeta = d at every level; projection on W_1 = sin(pi z / H) gives
    # d int W dz / int W^2 dz = d (2H/pi) / (H/2) = 4 d / pi
    nz = 200
    drF = np.full(nz, 0.5)
    RC = -(np.arange(nz) + 0.5) * 0.5
    rho_b = 1020.0 - 0.01 * RC
    d = 0.2
    rho = (1020.0 - 0.01 * (RC - d))[:, None, None]
    a = an.mode1_amplitude(rho, rho_b, RC, drF, np.ones((nz, 1, 1)), 1020.0)
    assert a[0, 0] == pytest.approx(4 * d / np.pi, rel=2e-3)


# ---- clock --------------------------------------------------------------------------------------

def test_bay_crop_is_the_wind_patch_up_to_the_bay_head():
    assert an.bay_crop({"cx": 280, "cy": 310, "radius": 60}, 352, 560) == (250, 352, 220, 341)
    assert an.bay_crop({"cx": 10, "cy": 5, "radius": 60}, 352, 560, margin=1.5) == (0, 352, 0, 101)


def test_pulse_onset_is_centre_minus_two_and_a_half_sigma_clipped_at_zero():
    assert an.pulse_onset_s(center_h=2.5, sigma_h=1.0) == 0.0
    assert an.pulse_onset_s(center_h=10.0, sigma_h=2.0) == 5.0 * 3600.0
    assert an.pulse_onset_s(center_h=1.0, sigma_h=1.0) == 0.0


def test_clock_hours_and_inertial_periods():
    f = 1.4266e-4
    ti = 2 * np.pi / f                                      # 44 043 s = 12.23 h
    h, n = an.clock(t_s=ti + 7200.0, onset_s=7200.0, f=f)
    assert h == pytest.approx(ti / 3600.0)
    assert n == pytest.approx(1.0)
    h, n = an.clock(t_s=36 * 3600.0, onset_s=0.0, f=-6.97e-5)   # sign of f does not matter
    assert h == pytest.approx(36.0)
    assert n == pytest.approx(36 * 3600.0 * 6.97e-5 / (2 * np.pi))
    assert an.clock_text(36 * 3600.0, 0.0, f).startswith("t = 36.0 h")


def test_common_frame_times_subsample_both_records():
    ta = np.arange(0, 121 * 3600 + 1, 900.0)                # bay: every 900 s
    tb = np.arange(0, 280 * 3600 + 1, 1200.0)               # fjord: every 1200 s
    t = an.common_times([ta, tb], step_s=3600.0, t0_s=0.0, t1_s=96 * 3600.0)
    np.testing.assert_array_equal(t, np.arange(0, 96 * 3600 + 1, 3600.0))
    with pytest.raises(ValueError):
        an.common_times([ta, tb], step_s=1000.0, t0_s=0.0, t1_s=3600.0)   # 1000 s not in both


# ---- scales -------------------------------------------------------------------------------------

def test_nice_ceiling():
    assert an.nice_ceil(0.37) == pytest.approx(0.4)
    assert an.nice_ceil(1.2) == 1.5
    assert an.nice_ceil(2.0) == 2.0
    assert an.nice_ceil(2.91) == 3.0
    assert an.nice_ceil(7.1) == 8.0
    assert an.nice_ceil(9.5) == 10.0


def test_limits_shared_or_separate():
    rng = np.random.default_rng(0)
    a = {"A": rng.normal(0, 1.0, (5, 10, 10)), "B": rng.normal(0, 4.0, (5, 10, 10))}
    a["A"][0, 0, 0] = np.nan
    sh = an.colour_limits(a, shared=True)
    assert sh["A"] == sh["B"]
    assert sh["A"][0] == -sh["A"][1] > -np.inf
    sep = an.colour_limits(a, shared=False)
    assert sep["A"][1] < sep["B"][1]
    for v in sep.values():
        assert v[0] == -v[1]


def test_colour_limits_fixed_across_frames_and_shared_across_panels():
    ny, nx, nt = 6, 8, 5
    rng = np.random.default_rng(1)
    frames = {"P": rng.normal(0, 1, (nt, ny, nx)), "Q": rng.normal(0, 3, (nt, ny, nx))}
    lims = an.colour_limits(frames, shared=True)
    panels = [an.Panel(name=r, title=r, x_km=np.arange(nx + 1.0), y_km=np.arange(ny + 1.0),
                       land=np.zeros((ny, nx), bool), shallow=np.zeros((ny, nx), bool),
                       eta=frames[r], lim=lims[r], f=1e-4) for r in ("P", "Q")]
    fig, draw = an.build_figure(panels, t_s=np.arange(nt) * 3600.0, onset_s=0.0,
                                suptitle="test", footer="footer")
    clims = []
    for n in range(nt):
        draw(n)
        clims.append(tuple(p.mesh.get_clim() for p in panels))
    assert len(set(clims)) == 1
    assert clims[0][0] == clims[0][1] == lims["P"]
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_saturated_fraction():
    eta = np.array([[0.0, 1.0, -2.0, np.nan]])
    assert an.saturated_fraction(eta, 1.0) == pytest.approx(2 / 3)


# ---- GIF budget ---------------------------------------------------------------------------------

def test_stride_doubles_until_the_gif_fits_the_budget(tmp_path):
    ny, nx, nt = 30, 40, 24
    rng = np.random.default_rng(2)
    eta = rng.normal(0, 1, (nt, ny, nx))                    # noise: poorly compressible frames
    panels = [an.Panel(name="N", title="noise", x_km=np.arange(nx + 1.0), y_km=np.arange(ny + 1.0),
                       land=np.zeros((ny, nx), bool), shallow=np.zeros((ny, nx), bool),
                       eta=eta, lim=(-2.0, 2.0), f=1e-4)]
    fig, draw = an.build_figure(panels, t_s=np.arange(nt) * 1200.0, onset_s=0.0,
                                suptitle="noise", footer="")
    out = tmp_path / "a.gif"
    full = an.write_gif(fig, draw, nt, out, max_mb=1e3, dpi=40)
    assert full["stride"] == 1 and full["frames"] == nt
    budget = full["mb"] / 3.0
    res = an.write_gif(fig, draw, nt, out, max_mb=budget, dpi=40)
    assert res["mb"] <= budget
    assert res["stride"] in (4, 8, 16) and not res["over_budget"]
    assert out.stat().st_size / 1e6 == pytest.approx(res["mb"])
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_mode1_amplitude_is_undefined_in_columns_shallower_than_the_minimum_depth():
    nz = 20
    drF = np.full(nz, 2.0)
    RC = -(np.arange(nz) + 0.5) * 2.0
    rho_b = 1020.0 - 0.01 * RC
    kb = np.array([[4, 5, 20]])                               # 8 m, 10 m, 40 m deep
    hfac = (np.arange(nz)[:, None, None] < kb[None]).astype(float)
    rho = np.tile(rho_b[:, None, None], (1, 1, 3)) + 0.001
    a = an.mode1_amplitude(rho, rho_b, RC, drF, hfac, 1020.0)
    assert an.MIN_DEPTH_M == 10.0
    assert np.isnan(a[0, 0]) and np.isfinite(a[0, 1]) and np.isfinite(a[0, 2])
