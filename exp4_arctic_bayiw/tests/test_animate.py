"""Tests of the descriptive animations (analysis/animate.py, Task 9)."""
from __future__ import annotations

import numpy as np
import pytest

from exp4_arctic_bayiw.analysis import animate as an


# ---- mode-1 displacement amplitude on a synthetic stratification ----------------------------------

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

def _linear_column(nz, dz):
    drF = np.full(nz, dz)
    RC = -(np.arange(nz) + 0.5) * dz
    return drF, RC, 1020.0 - 0.01 * RC                      # rho_b linear in z (constant N^2)


def test_mode1_amplitude_uniform_lift_gives_4d_over_pi_at_several_depths():
    # independent of mode1_shape: linear N^2 -> W_1 = sin(pi z / H) for every H, and a uniform lift
    # zeta = d projects to d int W dz / int W^2 dz = 4 d / pi (positive: lifted isopycnals, + up)
    nz = 200
    drF, RC, rho_b = _linear_column(nz, 0.5)
    kb = np.array([[200, 120, 60]])                          # H = 100, 60, 30 m
    hfac = (np.arange(nz)[:, None, None] < kb[None]).astype(float)
    d = 0.2
    rho = np.tile((rho_b + 0.01 * d)[:, None, None], (1, 1, 3))   # rho_b(z - d) = rho_b + 0.01 d
    a = an.mode1_amplitude(rho, rho_b, RC, drF, hfac, 1020.0)
    np.testing.assert_allclose(a[0], 4 * d / np.pi, rtol=3e-3)
    a_down = an.mode1_amplitude(np.tile((rho_b - 0.01 * d)[:, None, None], (1, 1, 3)), rho_b, RC,
                                drF, hfac, 1020.0)
    np.testing.assert_allclose(a_down[0], -4 * d / np.pi, rtol=3e-3)


def test_mode1_amplitude_with_a_partial_bottom_cell_uses_the_wet_thickness():
    # column of 12 cells of 2.5 m with hFac = 0.2 in the last one: wet depth H' = 27.5 + 0.5 = 28 m.
    # zeta = d sin(pi z / H') is the exact mode 1 of the wet column (linear N^2), so a_1 = d.
    nz = 40
    drF, RC, rho_b = _linear_column(nz, 2.5)
    hfac = np.zeros((nz, 1, 1))
    hfac[:11] = 1.0
    hfac[11] = 0.2
    dz = drF * hfac[:, 0, 0]
    zc = -(np.cumsum(dz) - dz / 2)                           # centres of the wet thickness
    Hw = dz.sum()
    assert Hw == pytest.approx(28.0)
    d = 0.5
    zeta = np.where(dz > 0, d * np.sin(np.pi * -zc / Hw), 0.0)
    rho = (rho_b + 0.01 * zeta)[:, None, None]
    a = an.mode1_amplitude(rho, rho_b, RC, drF, hfac, 1020.0)
    assert a[0, 0] == pytest.approx(d, rel=1e-2)


# ---- inertial-period filter ----------------------------------------------------------------------

def test_inertial_filter_removes_trend_and_constant_and_keeps_the_oscillation():
    f = 1.4266e-4
    T = 2 * np.pi / f                                        # 12.23 h, not a multiple of the step
    t = np.arange(0, 96 * 3600 + 1, 1200.0)
    A = 1.3
    osc = A * np.sin(2 * np.pi * t / T)
    x = (osc + 2.0e-5 * t - 4.0)[:, None, None] * np.ones((1, 2, 3))
    tk, xf = an.inertial_filter(x, t, T)
    keep = (tk[:, None] == t[None, :]).any(axis=0)
    assert xf.shape == (len(tk), 2, 3)
    np.testing.assert_allclose(xf, np.broadcast_to(osc[keep][:, None, None], xf.shape), atol=0.01 * A)


def test_inertial_filter_drops_the_frames_without_a_full_window():
    T = 10 * 3600.0
    t = np.arange(0, 48 * 3600 + 1, 3600.0)
    tk, xf = an.inertial_filter(np.zeros((len(t), 1, 1)), t, T)
    np.testing.assert_array_equal(tk, np.arange(5, 44) * 3600.0)      # 5 h .. 43 h
    assert xf.shape[0] == len(tk)
    tk2, _ = an.inertial_filter(np.zeros((len(t), 1, 1)), t, T + 1.0)   # 5 h no longer complete
    assert tk2[0] == 6 * 3600.0 and tk2[-1] == 42 * 3600.0
    with pytest.raises(ValueError):
        an.inertial_filter(np.zeros((5, 1, 1)), t[:5], T)                  # record shorter than T


# ---- views ---------------------------------------------------------------------------------------

def test_bay_crop_is_the_bay_plus_one_bay_width_and_one_bay_length():
    ny, nx = 60, 80
    land = np.zeros((ny, nx), bool)
    land[40:] = True                                         # coast row 40, land to the north
    land[40:55, 30:36] = False                               # bay: 6 cells wide, 15 long
    bay = an.bay_mask(land)
    assert bay.sum() == 6 * 15 and bay[40:55, 30:36].all()
    assert an.bay_crop(land) == (40 - 15, 60, 30 - 6, 36 + 6)


def test_bay_runs_are_told_apart_by_one_rule():
    assert an.is_bay_run({"params": {"nx": 560}})
    assert not an.is_bay_run({"params": {"bbox_utm": [0, 1, 2, 3]}})


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
    drF, RC, rho_b = _linear_column(nz, 2.0)
    kb = np.array([[4, 5, 20]])                               # 8 m, 10 m, 40 m deep
    hfac = (np.arange(nz)[:, None, None] < kb[None]).astype(float)
    rho = np.tile(rho_b[:, None, None], (1, 1, 3)) + 0.001
    a = an.mode1_amplitude(rho, rho_b, RC, drF, hfac, 1020.0)
    assert an.MIN_DEPTH_M == 10.0
    assert np.isnan(a[0, 0]) and np.isfinite(a[0, 1]) and np.isfinite(a[0, 2])


# ---- fix round 2: footer notes computed at render time --------------------------------------------

def test_filter_note_names_the_sub_inertial_attenuation():
    s = an.FILTER_NOTE
    for w in ("sub-inertial", "Kelvin", "seiches"):
        assert w in s


def test_saturation_note_only_where_the_shallow_columns_are_over_represented():
    key = f"share_of_saturated_in_columns_shallower_than_{an.SHALLOW_REPORT_M:g}m"
    area = f"fraction_of_shown_pixels_in_columns_shallower_than_{an.SHALLOW_REPORT_M:g}m"
    scales = {"R1": {key: 0.0, area: 0.1}, "F1": {key: 0.4564, area: 0.15}, "shared": False,
              "step_s": 3600.0}
    note = an.saturation_note(scales)
    assert "F1" in note and "R1" not in note
    assert f"shallower than {an.SHALLOW_REPORT_M:g} m" in note and "46 %" in note and "15 %" in note
    # 4.5 % of the saturated pixels in 7 % of the area: not concentrated
    assert an.saturation_note({"F3": {key: 0.045, area: 0.07}, "F2": {key: None, area: 0.07}}) == ""


# ---- Task 10 fix round 1: a missing sample blanks only the windows that contain it --------------

def test_inertial_filter_one_missing_sample_blanks_only_its_windows():
    f = 1.4266e-4
    T = 2 * np.pi / f
    dt = 1200.0
    t = np.arange(0.0, 100 * 3600.0 + 1, dt)
    x = (np.sin(2 * np.pi * t / (9 * 3600.0)) + 1e-5 * t)[:, None]
    ref_t, ref = an.inertial_filter(x, t, T)
    xn = x.copy()
    s = int(round(10 * 3600.0 / dt))                    # the sample at 10 h
    xn[s] = np.nan
    tk, xf = an.inertial_filter(xn, t, T)
    assert np.array_equal(tk, ref_t)
    d = np.abs(tk - t[s])
    # windows that contain 10 h are blank; the piecewise-linear integral also uses the sample over
    # one step on either side, so windows reaching within one step of it are blank too
    assert np.all(np.isnan(xf[d <= T / 2, 0]))
    far = d >= T / 2 + dt
    assert far.sum() > 200                               # almost every frame stays valid
    assert np.all(np.isfinite(xf[far, 0]))
    assert np.allclose(xf[far, 0], ref[far, 0])          # and unchanged
