"""Tests for exp4_arctic_bayiw.analysis.budget (Kang & Fringer 2012 as in the thesis script)."""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from exp4_arctic_bayiw.analysis import budget

PKG = Path(__file__).resolve().parents[1]
G, RHO0 = 9.81, 999.8


def _thesis(name):
    sys.dont_write_bytecode = True
    p = str(PKG / "thesis_copy" / "Output_IWBay")
    if p not in sys.path:
        sys.path.insert(0, p)
    return importlib.import_module(name)


# ---- pure pieces against the thesis functions ----------------------------------------------

def test_pieces_equal_thesis_functions():
    kf = _thesis("balance_kang_fringer")
    bt = _thesis("conversion_bt_2D")
    rng = np.random.default_rng(2)
    nz, ny, nx = 6, 4, 5
    dz = rng.uniform(2, 10, (nz, ny, nx))
    rp = 0.1 * rng.standard_normal((nz, ny, nx))
    eta = 0.1 * rng.standard_normal((ny, nx))
    np.testing.assert_allclose(budget.pressure_tilde(rp, eta, dz, RHO0, G),
                               kf.presion_tilde(rp, eta, dz, RHO0), rtol=1e-14)
    v, p, d2 = rng.standard_normal((nz, nx)), rng.standard_normal((nz, nx)), dz[:, 0]
    for a, b in zip(budget.face_pressure_work(v, p, d2), kf.flujo_cara(v, p, d2)):
        np.testing.assert_allclose(a, b, rtol=1e-14)
    u, ep = rng.standard_normal((nz, nx)), rng.random((nz, nx))
    for a, b in zip(budget.face_advection(u, v, d2, ep, RHO0),
                    kf.flujo_advectivo_cara(u, v, d2, ep, RHO0)):
        np.testing.assert_allclose(a, b, rtol=1e-14)
    uc, vc, vol = rng.standard_normal((3, nz, ny, nx))
    vol = np.abs(vol)
    np.testing.assert_allclose(budget.ke_bt_bc(uc, vc, vol, RHO0), kf.ke_bt_bc(uc, vc, vol, RHO0),
                               rtol=1e-14)
    zb = -(np.cumsum(np.full(nz, 10.0)) - 5.0)
    rb = 1025 + 0.01 * np.arange(nz) ** 1.5
    rt = rb[:, None, None] + 0.05 * rng.standard_normal((nz, ny, nx))
    z3 = zb[:, None, None] * np.ones((nz, ny, nx))
    np.testing.assert_allclose(budget.ape_exact(rt, z3, rb, zb, G), kf.ape_exacta(rt, z3, rb, zb),
                               rtol=1e-12, atol=1e-12)
    hW, hS = rng.random((nz, ny, nx + 1)), rng.random((nz, ny + 1, nx))
    dyG, dxG, rA = rng.uniform(50, 100, (ny, nx + 1)), rng.uniform(50, 100, (ny + 1, nx)), \
        rng.uniform(1e3, 1e4, (ny, nx))
    drF = dz[:, 0, 0]
    AW, AS = budget.areas_below_face(drF, hW, hS, dyG, dxG)
    AW2, AS2 = bt.areas_bajo_cara(drF, hW, hS, dyG, dxG)
    np.testing.assert_allclose(AW, AW2); np.testing.assert_allclose(AS, AS2)
    Ub, Vb = rng.standard_normal((ny, nx + 1)), rng.standard_normal((ny + 1, nx))
    Wc = budget.w_bt_center(Ub, Vb, AW, AS, rA)
    for k in range(nz):
        np.testing.assert_allclose(Wc[k], bt.w_bt_centro(Ub[None], Vb[None], AW, AS, rA, k)[0],
                                   rtol=1e-13, atol=1e-18)


def test_ape_exact_linear_background_equals_quadratic_form():
    """Kang & Fringer (2012) eq. 12 with linear rho_b gives exactly g^2 rho'^2 / (2 rho0 N^2)."""
    zb = -np.arange(0.5, 100, 1.0)
    rb = 1020 + 0.02 * (-zb)                                # d rho_b / dz = -0.02 kg m^-4
    N2 = G / RHO0 * 0.02
    rp = np.array([-0.3, -0.05, 0.0, 0.07, 0.4])
    z = np.full(5, -50.5)
    ep = budget.ape_exact(rb[50] + rp, z, rb, zb, G)
    np.testing.assert_allclose(ep, G ** 2 * rp ** 2 / (2 * RHO0 * N2), rtol=1e-10)


# ---- box geometry ------------------------------------------------------------------------------

def test_box_faces_staircase_signs():
    m = np.zeros((5, 6), bool)
    m[1, 1:3] = m[2, 1:4] = m[3, 2:5] = True
    f = budget.box_faces(m)
    # every open face separates a box cell from a non-box cell, with sign +1 = +x/+y into the box
    for j, i, s in f["x"]:
        assert m[j, i] != m[j, i - 1] and s == (1 if m[j, i] else -1)
    for j, i, s in f["y"]:
        assert m[j, i] != m[j - 1, i] and s == (1 if m[j, i] else -1)
    assert len(f["x"]) == 6 and len(f["y"]) == 8


def test_box_from_section_closes_a_channel_with_a_staircase():
    ny, nx = 20, 20
    xc, yc = np.arange(nx) + 0.5, np.arange(ny) + 0.5
    wet = np.zeros((ny, nx), bool)
    wet[5:15] = True
    sec = {"S": ((9.3, 4.0), (10.7, 16.0))}
    mask, groups = budget.box_from_sections(wet, xc, yc, sec, seed_in=(17.5, 10.5),
                                            seed_out=(2.5, 10.5), extend_m=0.0)
    assert mask[10, 17] and not mask[10, 2] and not (mask & ~wet).any()
    f = budget.box_faces(mask, wet=wet)
    assert len(f["x"]) + len(f["y"]) == groups["S"]["x"].sum() + groups["S"]["y"].sum() > 0
    for j, i, _ in f["x"]:
        assert groups["S"]["x"][j, i]
    for j, i, _ in f["y"]:
        assert groups["S"]["y"][j, i]
    with pytest.raises(RuntimeError, match="not closed"):
        budget.box_from_sections(wet, xc, yc, {"S": ((9.3, 6.0), (10.7, 12.0))},
                                 seed_in=(17.5, 10.5), seed_out=(2.5, 10.5), extend_m=0.0)


# ---- closed synthetic budget --------------------------------------------------------------------

def _filling_basin(nt=13, Vin=1e-4, eta0=0.0):
    """A 3 x 3 basin (box) fed through one open face by a uniform inflow Vin, eta = a t uniform
    with a A_box = Vin H dx (continuity), rho = rho_b: dE_bt/dt = F_bt,in (+ the U^3 advection)."""
    nz, ny, nx, d = 3, 6, 5, 100.0
    drF = np.full(nz, 10.0)
    H = drF.sum()
    hC = np.zeros((nz, ny, nx))
    hC[:, 2:5, 1:4] = 1.0
    hC[:, 1, 2] = 1.0                                          # the feeding cell
    wet = hC[0] > 0
    hW = np.zeros((nz, ny, nx)); hS = np.zeros((nz, ny, nx))
    hW[:, :, 1:] = (wet[:, 1:] & wet[:, :-1])[None]
    hS[:, 1:, :] = (wet[1:] & wet[:-1])[None]
    mask = np.zeros((ny, nx), bool); mask[2:5, 1:4] = True
    a = Vin * H * d / (mask.sum() * d * d)                     # m s^-1
    t = np.arange(nt) * 900.0
    V = np.zeros((nt, nz, ny, nx)); V[:, :, 2, 2] = Vin
    U = np.zeros_like(V)
    rho_b = 1025.0 + 0.01 * np.arange(nz)
    rho = np.broadcast_to(rho_b[None, :, None, None], V.shape).copy()
    eta = np.where(wet, 1.0, 0.0)[None] * (eta0 + a * t)[:, None, None]
    one = np.ones((ny, nx))
    ds = xr.Dataset(
        {"U": (("time", "k", "j", "i"), U), "V": (("time", "k", "j", "i"), V),
         "RHO": (("time", "k", "j", "i"), rho), "ETA": (("time", "j", "i"), eta),
         "hFacC": (("k", "j", "i"), hC), "hFacW": (("k", "j", "i"), hW),
         "hFacS": (("k", "j", "i"), hS), "drF": (("k",), drF),
         "RC": (("k",), -(np.cumsum(drF) - 0.5 * drF)),
         "rA": (("j", "i"), d * d * one), "rAs": (("j", "i"), d * d * one),
         "rAw": (("j", "i"), d * d * one), "dxG": (("j", "i"), d * one),
         "dyG": (("j", "i"), d * one)},
        coords={"time": t}, attrs={"rho0": RHO0, "g": G})
    return ds, mask, a


def test_closed_synthetic_budget_has_zero_residual():
    ds, mask, a = _filling_basin()
    r = budget.kf_budget(ds, mask)
    c = budget.closure(r)
    F = c["bt"]["F"][-1]
    assert F > 0                                                # inflow work enters the box
    assert c["bt"]["dE"][-1] == pytest.approx(F, rel=1e-6)
    assert abs(c["bt"]["D"][-1]) < 1e-6 * F
    assert np.abs(r["C"]).max() < 1e-9 * np.abs(r["F_bt"]).max() and np.all(r["P_bt"] == 0)
    # a draining basin above its rest level loses energy through the face: F < 0, D still 0
    ds2, mask2, _ = _filling_basin(Vin=-1e-4, eta0=1.0)
    c2 = budget.closure(budget.kf_budget(ds2, mask2))["bt"]
    assert c2["F"][-1] < 0 and abs(c2["D"][-1]) < 1e-6 * abs(c2["F"][-1])


F1_RUN = Path("/Volumes/LaCie/portfolio_runs/exp4/F1_pd1final_48h")


@pytest.mark.skipif(not (F1_RUN / "hFacC.meta").exists(), reason="fjord run not mounted")
def test_fjord_box_closes_at_the_two_sounds():
    from exp4_arctic_bayiw.analysis import load
    ds = load.run(F1_RUN)
    mask, groups = budget.fjord_box(ds)
    assert 13000 < mask.sum() < 13500                  # inner fjord: 13169 cells without the disks
    f = budget.box_faces(mask, wet=ds["hFacC"].values[0] > 0)
    n_open = len(f["x"]) + len(f["y"])
    assert n_open == sum(int(g["x"].sum() + g["y"].sum()) for g in groups.values())
    assert all(int(g["x"].sum() + g["y"].sum()) > 0 for g in groups.values())
    with pytest.raises(RuntimeError, match="not closed"):
        budget.fjord_box(ds, extend_m=0.0)


def test_laplacian_dissipation_equals_thesis_function():
    th = _thesis("disipacion_viscosa_bahia")
    rng = np.random.default_rng(3)
    nz, ny, nx = 5, 4, 3
    drF = rng.uniform(1, 5, nz)
    g = {"drF": drF, "drC": rng.uniform(1, 5, nz + 1),
         "hC": (rng.random((nz, ny, nx)) > 0.2).astype(float),
         "hW": (rng.random((nz, ny, nx + 1)) > 0.2).astype(float),
         "hS": (rng.random((nz, ny + 1, nx)) > 0.2).astype(float),
         "rA": rng.uniform(1e3, 2e3, (ny, nx)), "dxF": rng.uniform(30, 60, (ny, nx)),
         "dyF": rng.uniform(30, 60, (ny, nx)), "rAw": rng.uniform(1e3, 2e3, (ny, nx + 1)),
         "rAs": rng.uniform(1e3, 2e3, (ny + 1, nx)), "dyU": rng.uniform(30, 60, (ny - 1, nx + 1)),
         "dxV": rng.uniform(30, 60, (ny + 1, nx - 1)), "rAzU": rng.uniform(1e3, 2e3, (ny - 1, nx + 1)),
         "rAzV": rng.uniform(1e3, 2e3, (ny + 1, nx - 1))}
    dz = drF[:, None, None]
    g["HC"], g["HW"], g["HS"] = (g["hC"] * dz).sum(0), (g["hW"] * dz).sum(0), (g["hS"] * dz).sum(0)
    U, V = rng.standard_normal((nz, ny, nx + 1)), rng.standard_normal((nz, ny + 1, nx))
    np.testing.assert_allclose(budget.laplacian_dissipation(U, V, g, RHO0, th.NU_H, th.NU_V),
                               th.disipacion(U, V, g, RHO0), rtol=1e-13)


def test_conversion_mask_outside_the_region_is_refused_and_inside_is_computed():
    ds, mask, _ = _filling_basin(nt=3)
    feed = np.zeros_like(mask); feed[1, 2] = True
    r = budget.kf_budget(ds, mask, conversion_masks={"feed": feed})
    assert r["C_feed"].shape == r["C"].shape
    big = np.zeros((6, 5), bool); big[3, 2] = True          # box region: rows 2-4, cols 1-3
    corner = np.zeros((6, 5), bool); corner[0, 0] = True
    with pytest.raises(ValueError, match="leaves"):
        budget.kf_budget(ds, big, conversion_masks={"corner": corner})
