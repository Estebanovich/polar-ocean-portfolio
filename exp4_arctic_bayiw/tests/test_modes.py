"""Tests for exp4_arctic_bayiw.theory.modes (analytic limits only)."""
from __future__ import annotations

import json

import numpy as np
import pytest

from exp4_arctic_bayiw.theory.modes import (
    burger,
    coriolis,
    kelvin_circuit_period,
    merian_period,
    modal_inner_product,
    rossby_radius,
    theory_summary,
    vertical_modes,
)


def test_coriolis_78N():
    assert coriolis(78.0) == pytest.approx(1.4266e-4, abs=1e-8)


def test_constant_N_analytic():
    H, N = 164.0, 0.01
    z = np.linspace(0.0, -H, 200)
    c, _ = vertical_modes(z, np.full_like(z, N**2), H, n_modes=3)
    for n in (1, 2, 3):
        exact = N * H / (n * np.pi)
        assert abs(c[n - 1] - exact) / exact < 0.005


def test_two_layer_limit():
    H, h1, gp = 100.0, 20.0, 0.02
    z = np.linspace(0.0, -H, 4001)
    N2 = 1e-8 + gp * np.exp(-0.5 * ((z + h1) / 0.15) ** 2) / (0.15 * np.sqrt(2 * np.pi))
    c, _ = vertical_modes(z, N2, H, n_modes=2)
    exact = np.sqrt(gp * h1 * (H - h1) / H)
    assert abs(c[0] - exact) / exact < 0.02


def test_orthogonality_and_order():
    H = 150.0
    z = np.linspace(0.0, -H, 300)
    N2 = 1e-4 * (1 + 3 * np.exp(z / 30.0))
    c, phi = vertical_modes(z, N2, H, n_modes=4)
    assert np.all(np.diff(c) < 0)
    for m in range(4):
        for n in range(m):
            assert abs(modal_inner_product(phi[m], phi[n], z, N2)) < 1e-6
        assert modal_inner_product(phi[m], phi[m], z, N2) == pytest.approx(1.0)
        assert abs(phi[m][0]) < 1e-12 and abs(phi[m][-1]) < 1e-12  # w = 0 at both ends


def test_merian_thesis_number():
    T = merian_period(119e3, 164.0, open_end=True)
    assert T == pytest.approx(4 * 119e3 / np.sqrt(9.81 * 164.0))
    assert T / 3600 == pytest.approx(3.29, abs=0.01)
    assert merian_period(119e3, 164.0, open_end=False) == pytest.approx(T / 2)


def test_kelvin_period():
    assert kelvin_circuit_period(5e5, 0.5) == pytest.approx(1e6)


def test_rossby_and_burger():
    f = coriolis(78.0)
    assert rossby_radius(1.0, f) == pytest.approx(1.0 / f)
    assert burger(1.0, f, 1e4) == pytest.approx(1.0 / f / 1e4)


def test_summary_json_serialisable():
    H = 164.0
    z = np.linspace(0.0, -H, 100)
    s = theory_summary(z, np.full_like(z, 1e-4), H, 78.0, 20e3, 80e3)
    json.dumps(s)
    assert len(s["c_n"]) == len(s["rossby_radius_m"]) == len(s["burger"])
