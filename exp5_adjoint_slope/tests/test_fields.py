"""gendata.fields: bathymetry, wind, polynya salt, initial T/S, surface and sponge relaxation,
GM map. Values from Dettling et al. (2024, Table 1) and Stewart & Thompson (2015); see
notes/pd1b_config_dettling.md."""
import numpy as np
import pytest

from exp5_adjoint_slope.gendata import fields, grid

G10 = grid.Grid(40, 45, 10e3, grid.vertical_spacing())
G2P5 = grid.Grid(160, 180, 2.5e3, grid.vertical_spacing())
ALPHA, BETA = 3.0e-5, 7.4e-4


def test_depth_profile():
    # Table 1: 500 m shelf, 3000 m deep, slope centred at 200 km, 150 km wide (shelf 125 km).
    y = np.linspace(0, 450e3, 901)
    d = fields.depth(y)
    assert d[0] == pytest.approx(500.0, abs=15.0)
    assert d[-1] == pytest.approx(3000.0, abs=15.0)
    assert fields.depth(np.array([200e3]))[0] == pytest.approx(1750.0, abs=1e-6)
    assert np.all(np.diff(d) >= 0)
    # tanh scale W/4 = 37.5 km [verificar]: the shelf edge (125 km) is within 50 m of 500 m and
    # the mean slope over the 150 km slope is 0.016, near Dettling's ~0.02 (p. 1676).
    assert fields.depth(np.array([125e3]))[0] < 550.0
    s = fields.topo_slope(y)
    assert s.max() == pytest.approx(1250.0 / 37.5e3, rel=1e-3)
    mean = (fields.depth(np.array([275e3])) - fields.depth(np.array([125e3])))[0] / 150e3
    assert 0.015 < mean < 0.02


def test_bathy_uniform_in_x():
    b = fields.bathymetry(G10)
    assert b.shape == (45, 40) and np.all(b[1:-1] < 0)
    assert np.all(b == b[:, :1])


@pytest.mark.parametrize("g", [G10, G2P5])
def test_bathy_has_walls_south_and_north(g):
    # MITgcm is periodic in y unless land closes it: without these walls the 3000 m northern
    # ocean connects to the shelf across the boundary (first S0: 0.17 Sv net through every row).
    b = fields.bathymetry(g)
    assert np.all(b[0] == 0) and np.all(b[-1] == 0)


def test_wind_easterly_piecewise_sinusoid():
    # S&T p. 434 / Table 1: easterly only, 0 at y = 50 and 450 km, max -0.075 N m-2 at 225 km.
    t = fields.wind_tau_x(np.array([0.0, 30e3, 50e3, 225e3, 450e3]))
    assert t[0] == 0.0 and t[1] == 0.0 and t[2] == pytest.approx(0.0, abs=1e-12)
    assert t[3] == pytest.approx(-0.075) and t[4] == pytest.approx(0.0, abs=1e-12)
    y = np.linspace(0, 450e3, 4501)
    tau = fields.wind_tau_x(y)
    assert np.all(tau <= 0) and tau.min() == pytest.approx(-0.075, abs=1e-6)
    # Ekman convergence over the slope: d(tau)/dy < 0 shoreward of the maximum (downwelling
    # for f < 0 and easterlies, Thompson et al. 2018 p. 750).
    on_slope = (y > 125e3) & (y < 225e3)
    assert np.all(np.diff(tau)[on_slope[1:]] < 0)


def test_polynya_salt_flux():
    # Table 1: 2.5 mg m-2 s-1 over the shoreward 50 km; MITgcm saltFlux is upward (FFIELDS.h).
    f = fields.polynya_salt_flux(G10)
    assert f.shape == (45, 40)
    assert np.all(f[G10.yc < 50e3] == pytest.approx(-2.5e-3))
    assert np.all(f[G10.yc > 50e3] == 0)


def test_initial_ts_uniform_in_y_and_stable():
    T, S = fields.initial_ts(G10)
    assert T.shape == S.shape == (35, 45, 40)
    assert np.all(T == T[:, :1, :1]) and np.all(S == S[:, :1, :1])
    rho = -ALPHA * T + BETA * S
    assert np.all(np.diff(rho, axis=0) >= -1e-15)
    assert T[0, 0, 0] == pytest.approx(-1.8, abs=0.01) and T[-1, 0, 0] == pytest.approx(0.8, abs=0.01)


def test_surface_relaxation_ice_like():
    # notes/pd1b: SST relaxed to theta_f with tau = dz1/gamma (gamma = 6e-5 m s-1, Schmidt 2004);
    # SSS relaxed to 34.2 in 30 days except in the polynya (y < 50 km), where lambda = 0.
    sst, sss, lam_t, lam_s = fields.surface_relaxation(G10)
    assert np.all(sst == -1.9) and np.all(sss == pytest.approx(34.2))
    assert np.all(lam_t == pytest.approx(6e-5 / G10.dz_m[0]))
    assert np.all(lam_s[G10.yc < 50e3] == 0)
    assert np.all(lam_s[G10.yc > 50e3] == pytest.approx(1.0 / (30 * 86400.0)))


def test_sponge_mask_50km_linear():
    # Table 1 / S&T p. 434: rate linear from 0 at y = 400 km to 1/tau at the northern wall.
    m = fields.sponge_mask(G10)
    assert m.shape == (35, 45, 40)
    assert np.all(m[:, G10.yc < 400e3] == 0) and m.max() <= 1
    assert m[0, 43, 0] == pytest.approx((G10.yc[43] - 400e3) / 50e3)


def test_gm_factor_slope_aware():
    # Dettling Eq. (24): kappa_GM = kappa_high / (1 + 800 s_topo); ~1 on shelf and offshore,
    # ~1/14 at the slope centre.
    f = fields.gm_factor(np.array([20e3, 200e3, 430e3]))
    assert f[1] == pytest.approx(1.0 / (1.0 + 800.0 * 1250.0 / 37.5e3), rel=1e-6)
    assert f[0] > 0.9 and f[2] > 0.9
