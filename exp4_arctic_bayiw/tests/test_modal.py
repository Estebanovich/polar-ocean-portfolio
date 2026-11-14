"""Tests for exp4_arctic_bayiw.analysis.modal (ruling R4: u-modes = d phi_n / dz)."""
from __future__ import annotations

import numpy as np
import pytest

from exp4_arctic_bayiw.analysis import modal

H = 120.0


def _grid(nz=40):
    zf = -np.linspace(0, H, nz + 1) ** 1.2 / H ** 0.2          # stretched faces, 0 ... -H
    N2 = 1e-4 * np.exp(zf / 20.0) + 1e-6
    return zf, N2


def test_u_modes_orthonormal_and_depth_mean_free():
    zf, N2 = _grid()
    c, U, dz = modal.u_modes(zf, N2, n_modes=6)
    assert np.all(np.diff(c) < 0)
    G = (U * dz) @ U.T / H                                      # (1/H) int U_m U_n dz
    np.testing.assert_allclose(G, np.eye(U.shape[0]), atol=1e-10)
    np.testing.assert_allclose(U[0], 1.0)                       # barotropic mode 0
    np.testing.assert_allclose((U[1:] * dz).sum(1), 0.0, atol=1e-10)


def test_project_reconstructs_modes_1_and_2():
    zf, N2 = _grid()
    _, U, dz = modal.u_modes(zf, N2, n_modes=6)
    u = 0.02 * U[0] + 0.7 * U[1] - 0.3 * U[2]
    a = modal.project(u, U, dz)
    np.testing.assert_allclose(a[:3], [0.02, 0.7, -0.3], atol=1e-12)
    assert np.max(np.abs(a @ U - u)) < 1e-6


def test_q_local_alford_fraction():
    zf, N2 = _grid()
    _, U, dz = modal.u_modes(zf, N2, n_modes=20)
    assert modal.q_local(modal.project(U[1], U, dz)) == pytest.approx(0.0, abs=1e-10)
    assert modal.q_local(modal.project(U[5], U, dz)) == pytest.approx(1.0, abs=1e-10)
    zc = 0.5 * (zf[1:] + zf[:-1])
    q = modal.q_local(modal.project(modal.slab(zc, 15.0), U, dz))
    assert 0.0 < q < 1.0
