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
