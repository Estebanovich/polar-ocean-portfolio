"""Post hoc V2 diagnostics (analysis/pd2_v2_diagnostics.py) on synthetic series."""
from __future__ import annotations

import numpy as np
import pytest

from exp4_arctic_bayiw.analysis import pd2_v2_diagnostics as dg

F, DT = 1.4266e-4, 1200.0


def test_peak_map_finds_each_column_frequency():
    t = np.arange(840) * DT
    w = np.array([1.0, 1.1, 1.15]) * F
    U = np.cos(t[:, None] * w[None]); V = -np.sin(t[:, None] * w[None])
    r, df = dg.peak_map(U, V, t, F)
    np.testing.assert_allclose(r, w / F, atol=df)


def test_lag_speed_of_a_travelling_signal():
    t = np.arange(840) * DT
    s = np.array([0.0, 10e3, 20e3, 30e3])
    c = 0.4
    z = np.stack([np.sin(2 * np.pi * (t - x / c) / (60 * 3600)) for x in s], axis=1)
    cl, lags = dg.lag_speed(z, t, s)
    assert cl == pytest.approx(c, rel=0.05)


def test_v4_section_summary_net_gross_and_bands():
    from exp4_arctic_bayiw.analysis import pd2_v4_diagnostics as d4, validate as va
    from exp4_arctic_bayiw.tests.test_validate import _channel_ds, BANDS_HZ
    ds, mask, groups, rho_b = _channel_ds(nt=240)
    fp = va.face_profiles(ds, mask, groups, rho_b=rho_b)
    s = d4.summarise(fp, BANDS_HZ)["mouth"]
    assert s["gross_J"] >= abs(s["net_out_J"]) > 0
    assert sum(s["band_out_J"].values()) == pytest.approx(s["mean_out_W"] * 240 * DT, rel=1e-9)


def test_v4_descriptive_shares_and_signs():
    from exp4_arctic_bayiw.analysis import pd2_v4_diagnostics as d4
    band = lambda a, b, c, e, low: {"kelvin_subinertial": a, "near_inertial": b, "superinertial": c,
                                    "seiche": e, "low": low}
    diag = {"F1": {"S": {"gross_J": 10.0, "band_out_J": band(2.0, -1.0, 0.5, 0.5, -5.0)}},
            "V4": {"S": {"gross_J": 4.0, "band_out_J": band(4.0, -1.0, -0.5, 0.5, 1.0)}}}
    r = d4.descriptive(diag)["S"]
    assert r["gross_J"]["V4_over_F1"] == 0.4
    assert r["bands"]["kelvin_subinertial"]["abs_share_F1"] == 0.5
    assert r["bands"]["kelvin_subinertial"]["abs_share_V4"] == pytest.approx(4 / 6)
    assert r["bands"]["superinertial"]["same_sign"] is False and r["bands"]["near_inertial"]["same_sign"]
    assert r["max_abs_share_diff"] == pytest.approx(4 / 6 - 0.5)


def test_free_wave_speed_ignores_the_forced_period_and_tracks_the_first_extremum():
    """PD2 V2b descriptive (post hoc): a uniform forced bump before 6 h plus a free wave of 60 h
    period travelling at 0.4 m/s; the speed comes from the first crest or trough after 6 h."""
    t = np.arange(840) * DT
    s = np.arange(12) * 5e3
    c = 0.4
    forced = 2.0 * np.exp(-0.5 * ((t - 2.5 * 3600) / 3600) ** 2)
    free = np.stack([np.where(t > x / c + 8 * 3600,
                              np.sin(2 * np.pi * (t - x / c - 8 * 3600) / (60 * 3600)), 0.0) for x in s], axis=1)
    z = forced[:, None] + free
    band = (2 / (280 * 3600), 0.85 * F / (2 * np.pi))
    r = dg.free_wave_speed(z, t, s, band, t_min=6 * 3600.0)
    assert r["c_m_s"] == pytest.approx(c, rel=0.1)
    assert r["c_stderr_m_s"] > 0 and r["n_stations"] == 12 and r["extremum"] in ("crest", "trough")


def _fw_case(forced_amp=2.0):
    """Uniform forced bump of amplitude forced_amp before 6 h plus a unit free wave at 0.4 m/s."""
    t = np.arange(840) * DT
    s = np.arange(12) * 5e3
    c = 0.4
    forced = forced_amp * np.exp(-0.5 * ((t - 2.5 * 3600) / 3600) ** 2)
    free = np.stack([np.where(t > x / c + 8 * 3600,
                              np.sin(2 * np.pi * (t - x / c - 8 * 3600) / (60 * 3600)), 0.0) for x in s], axis=1)
    return t, s, c, forced[:, None] + free, (2 / (280 * 3600), 0.85 * F / (2 * np.pi))


def test_free_wave_needs_the_cut_before_filtering():
    """A forced pulse 6 times the free wave (F1: |zeta| 1.4-1.9 m during the pulse): filtering the
    whole record rings past 6 h above the 1/2-max floor; cutting first recovers c."""
    t, s, c, z, band = _fw_case(forced_amp=6.0)
    assert dg.free_wave_speed(z, t, s, band)["c_m_s"] == pytest.approx(c, rel=0.1)
    r = dg.free_wave_speed(z, t, s, band, cut_before_filter=False)
    assert not (r["c_m_s"] == pytest.approx(c, rel=0.25))


def test_free_wave_needs_the_amplitude_floor():
    """Without the 1/2-max floor the first ripple after the cut counts as the crest."""
    t, s, c, z, band = _fw_case()
    assert dg.free_wave_speed(z, t, s, band)["c_m_s"] == pytest.approx(c, rel=0.1)
    r = dg.free_wave_speed(z, t, s, band, frac=0.0)
    assert not (r["c_m_s"] == pytest.approx(c, rel=0.25))
