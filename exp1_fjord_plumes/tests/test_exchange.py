"""Tests for the fjord-shelf exchange diagnostics of config (c) (Task 7)."""
from __future__ import annotations

import numpy as np

from exp1_fjord_plumes.analysis.exchange import exchange_flux, heat_toward_front


def _section(nz=20, ny=8, dz=10.0, dy=200.0):
    return np.full((nz, ny), dz * dy)


def test_exchange_flux_two_layer_synthetic():
    """u = +U in the upper half, -U in the lower half of a section of area A:
    Q_ex = U A / 2, the transport of one layer."""
    dA = _section()
    U = 0.05
    u = np.where(np.arange(dA.shape[0])[:, None] < dA.shape[0] // 2, U, -U)
    u = np.broadcast_to(u, dA.shape)
    A = dA.sum()
    assert np.isclose(exchange_flux(u, dA), 0.5 * U * A, rtol=1e-12)
    assert np.isclose(exchange_flux(u, dA), (u * dA)[u > 0].sum(), rtol=1e-12)


def test_exchange_flux_zero_for_rest():
    dA = _section()
    assert exchange_flux(np.zeros_like(dA), dA) == 0.0


def test_exchange_flux_ignores_land_nan():
    dA = _section()
    u = np.full(dA.shape, 0.1)
    u[-1, :] = np.nan
    dA[-1, :] = 0.0
    assert np.isclose(exchange_flux(u, dA), 0.5 * 0.1 * dA.sum())


def test_heat_toward_front_sign_and_value():
    """Front at x = 0: u < 0 carries heat toward it. Uniform u = -U, T - Tref
    = dT -> rho0 cp U dT A (W), positive."""
    dA = _section()
    U, dT = 0.02, 3.0
    H = heat_toward_front(np.full(dA.shape, -U), np.full(dA.shape, 1.0 + dT),
                          dA, T_ref=1.0)
    assert np.isclose(H, 1027.0 * 3974.0 * U * dT * dA.sum(), rtol=1e-12)
    assert heat_toward_front(np.full(dA.shape, U), np.full(dA.shape, 1.0 + dT),
                             dA, T_ref=1.0) < 0
