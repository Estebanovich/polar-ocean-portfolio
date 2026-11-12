"""Tests for exp4_arctic_bayiw.analysis.spectra (synthetic series)."""
from __future__ import annotations

import numpy as np
import pytest

from exp4_arctic_bayiw.analysis import spectra

DT = 1200.0                                   # s, fjord diagnostics interval
T5D = 5 * 86400.0                             # 5 days


def _t():
    return np.arange(int(T5D / DT)) * DT


def test_psd_peak_of_inertial_sinusoid_within_one_bin():
    f0 = 1.0 / (12.23 * 3600.0)               # Hz
    f, P, df = spectra.psd(np.cos(2 * np.pi * f0 * _t()), DT)
    assert df == pytest.approx(1.0 / T5D)
    assert abs(f[np.argmax(P)] - f0) <= df


def test_psd_parseval_with_window():
    rng = np.random.default_rng(0)
    x = rng.standard_normal(500)
    f, P, df = spectra.psd(x, DT)
    w = np.hanning(x.size)
    xd = x - x.mean()
    assert P.sum() * df == pytest.approx(np.sum((xd * w) ** 2) / np.sum(w ** 2), rel=1e-12)


def test_band_energy_rejects_overlap_and_narrow_bands():
    f, P, df = spectra.psd(np.cos(2 * np.pi * _t() / 43200.0), DT)
    with pytest.raises(ValueError, match="overlap"):
        spectra.band_energy(f, P, {"a": (0.0, 10 * df), "b": (9 * df, 20 * df)}, df)
    with pytest.raises(ValueError, match="narrower"):
        spectra.band_energy(f, P, {"a": (5 * df, 6.5 * df)}, df)
    e = spectra.band_energy(f, P, {"semi": (1 / 43200.0 - 2 * df, 1 / 43200.0 + 2 * df)}, df)
    assert e["semi"] == pytest.approx(0.5, rel=0.05)          # variance of a unit cosine


def test_rotary_inertial_motion_is_clockwise_in_the_north():
    f_cor = 1.42e-4                                           # s^-1, ~78 N
    t = _t()
    u, v = np.cos(f_cor * t), -np.sin(f_cor * t)
    f, Pcw, Pccw, df = spectra.rotary(u, v, DT)
    assert Pcw.sum() > 100 * Pccw.sum()
    assert abs(f[np.argmax(Pcw)] - f_cor / (2 * np.pi)) <= df
    _, Pcw2, Pccw2, _ = spectra.rotary(u, -v, DT)             # reversed rotation
    assert Pccw2.sum() > 100 * Pcw2.sum()
    w = np.hanning(t.size)
    tot = (np.sum(((u - u.mean()) * w) ** 2) + np.sum(((v - v.mean()) * w) ** 2)) / np.sum(w ** 2)
    assert (Pcw + Pccw).sum() * df == pytest.approx(tot, rel=1e-12)
