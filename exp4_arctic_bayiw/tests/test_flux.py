"""Tests for exp4_arctic_bayiw.analysis.flux (Nash et al. 2005 as in the thesis script)."""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

from exp4_arctic_bayiw.analysis import flux

PKG = Path(__file__).resolve().parents[1]
G, RHO0 = 9.81, 999.8
H, NZ = 600.0, 200
DZ = np.full(NZ, H / NZ)
ZC = -(np.cumsum(DZ) - 0.5 * DZ)
OMEGA = 2 * np.pi / 43200.0
T = np.arange(5 * 48) * 900.0                       # 5 whole periods, endpoint excluded
DL = np.array([100.0, 200.0])
RHO_BG = 1025.0 - 0.01 * ZC                         # kg m^-3, stable


def _mode1(U0=0.1, P0=50.0, sign=1.0):
    """u' = U0 cos(mz) cos(wt), p' = P0 cos(mz) cos(wt) from rho' = (P0 m/g) sin(mz) cos(wt)
    (hydrostatic, Nash et al. 2005 eq. 3 with eta = 0); <F> = P0 U0 H / 4 [W m^-1]."""
    m = np.pi / H
    ct = np.cos(OMEGA * T)[:, None, None]
    u = sign * U0 * np.cos(m * ZC)[None, :, None] * ct + 0.05 + 0.01 * ct     # + mean and bt
    rho = RHO_BG[None, :, None] + (P0 * m / G) * np.sin(m * ZC)[None, :, None] * ct
    u = np.repeat(u, 2, axis=2)
    rho = np.repeat(rho, 2, axis=2)
    eta = np.zeros((T.size, 2))
    return u, rho, eta, P0 * U0 * H / 4


def test_mode1_flux_within_2_percent_and_sign_flips():
    u, rho, eta, F_exact = _mode1()
    r = flux.section_flux(u, rho, eta, DZ, DL, T, rho0=RHO0, g=G)
    s = r["stages"][0]
    np.testing.assert_allclose(s["F_bc_mean"], F_exact, rtol=0.02)
    assert s["F_bc_total_mean"] == pytest.approx(F_exact * DL.sum(), rel=0.02)
    u2, rho2, eta2, _ = _mode1(sign=-1.0)
    r2 = flux.section_flux(u2, rho2, eta2, DZ, DL, T, rho0=RHO0, g=G)
    assert r2["stages"][0]["F_bc_total_mean"] == pytest.approx(-F_exact * DL.sum(), rel=0.02)


def test_pure_barotropic_has_no_baroclinic_flux():
    U, eta0 = 0.2, 0.3
    ct = np.cos(OMEGA * T)
    u = np.broadcast_to(U * ct[:, None, None], (T.size, NZ, 2)).copy()
    rho = np.broadcast_to(RHO_BG[None, :, None], (T.size, NZ, 2)).copy()
    eta = np.repeat((eta0 * ct)[:, None], 2, axis=1)
    s = flux.section_flux(u, rho, eta, DZ, DL, T, rho0=RHO0, g=G)["stages"][0]
    assert np.max(np.abs(s["F_bc"])) < 1e-9 * RHO0 * G * eta0 * U * H
    # thesis F_bt = p_bt' u_bt' H (stage perturbations): <F_bt> = rho0 g eta0 U H / 2
    np.testing.assert_allclose(s["F_bt_mean"], RHO0 * G * eta0 * U * H / 2, rtol=1e-12)


def test_stages_share_the_boundary_instant():
    u, rho, eta, _ = _mode1()
    t_split = T[T.size // 2]
    r = flux.section_flux(u, rho, eta, DZ, DL, T, rho0=RHO0, g=G,
                          stages=[(T[0], t_split), (t_split, T[-1])])
    a, b = r["stages"]
    assert a["t"][-1] == b["t"][0] == t_split
    assert a["t"].size + b["t"].size == T.size + 1


def _thesis_nash():
    sys.dont_write_bytecode = True
    p = str(PKG / "thesis_copy" / "Output_IWBay")
    if p not in sys.path:
        sys.path.insert(0, p)
    pytest.importorskip("netCDF4")
    return importlib.import_module("energy_flux_IWBay")


def test_equals_thesis_functions_on_random_section():
    """Same numbers as energy_flux_IWBay.py (presion_baroclínica, flujo_energia_seccion and the
    stage de-meaning of procesar_etapa) on a random one-row tile, to round-off."""
    th = _thesis_nash()
    rng = np.random.default_rng(1)
    nt, nz, nx = 30, 8, 5
    dz = rng.uniform(5, 20, nz)
    Traw = 10 + rng.standard_normal((nt, nz, 1, nx))
    Sraw = 33 + 0.1 * rng.standard_normal((nt, nz, 1, nx))
    Uraw = 0.1 * rng.standard_normal((nt, nz, 1, nx + 1))
    Vraw = 0.1 * rng.standard_normal((nt, nz, 2, nx))
    Eta = 0.1 * rng.standard_normal((nt, 1, nx))
    eos = dict(tRef=16.582, sRef=33.6, tAlpha=2e-4, sBeta=7.4e-4, rhoNil=RHO0)
    rho = th.eos(Traw, Sraw, **eos)
    # thesis: means over the whole run (calcular_medias_run), then per time
    rho_ref = rho[:, :, 0, :].mean(0).mean(1)
    UT, VT = th._interpolar_UV_Tpoint(Uraw, Vraw)
    u_mean, v_mean = UT[:, :, 0].mean(0), VT[:, :, 0].mean(0)
    Fy, pbt, vbt = [], [], []
    for n in range(nt):
        p_bc, p_bar = th.presion_baroclínica(rho[n], rho_ref, Eta[n], dz, RHO0)
        _, fy, pb, _, vb = th.flujo_energia_seccion(Uraw[n], Vraw[n], p_bc, p_bar, dz,
                                                     u_mean, v_mean, 0, slice(0, nx))
        Fy.append(fy); pbt.append(pb); vbt.append(vb)
    Fy, pbt, vbt = np.array(Fy), np.array(pbt), np.array(vbt)
    Fbt = (pbt - pbt.mean(0)) * (vbt - vbt.mean(0)) * dz.sum()
    t = np.arange(nt) * 900.0
    r = flux.section_flux(VT[:, :, 0], rho[:, :, 0], Eta[:, 0], dz, np.ones(nx), t, rho0=RHO0, g=G)
    s = r["stages"][0]
    np.testing.assert_allclose(r["rho_ref"], rho_ref, rtol=1e-14)
    np.testing.assert_allclose(s["F_bc"], Fy, rtol=1e-9, atol=1e-12 * np.abs(Fy).max())
    np.testing.assert_allclose(s["F_bt"], Fbt, rtol=1e-9, atol=1e-12 * np.abs(Fbt).max())


def test_section_T_point_velocity_like_thesis():
    """Normal velocity at T points of row j: 0.5 (V[j] + V[j+1]) (energy_flux_IWBay
    _interpolar_UV_Tpoint); x-normal: 0.5 (U[i] + U[i+1])."""
    V = np.arange(2 * 4 * 6, dtype=float).reshape(2, 4, 6)
    np.testing.assert_allclose(flux.t_point_normal(V, "y", 1, slice(2, 5)),
                               0.5 * (V[:, 1, 2:5] + V[:, 2, 2:5]))
    np.testing.assert_allclose(flux.t_point_normal(V, "x", 1, slice(2, 4)),
                               0.5 * (V[:, 2:4, 1] + V[:, 2:4, 2]))
