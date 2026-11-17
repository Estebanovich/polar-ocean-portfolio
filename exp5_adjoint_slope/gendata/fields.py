"""Input fields of the Exp. 5 channel: the configuration of Dettling et al. (2024, Table 1) and
Stewart & Thompson (2015, p. 434), adapted as listed in notes/pd1b_config_dettling.md (linear EOS,
ice-like surface relaxation instead of a sea-ice model, idealized profiles [verificar])."""
from __future__ import annotations

import numpy as np

from .grid import Grid

# Bathymetry (Dettling Table 1): 500 m shelf, 3000 m deep, slope centred at 200 km, 150 km wide
# (shelf width 125 km). "tanh-like" (S&T p. 434) with scale W/4 [verificar]: the shelf edge at
# 125 km is then 545 m deep and kappa_GM ~ kappa_high up to ~75 km, as in Dettling Fig. 5.
H_SHELF, H_DEEP = 500.0, 3000.0
Y_SLOPE, W_SLOPE = 200e3, 150e3
# Wind (S&T p. 434; Table 1): easterly only, piecewise sinusoid, 0 at Y_0 and Y_N, max at Y_W.
TAU_MAX, Y_0, Y_W, Y_N = -0.075, 50e3, 225e3, 450e3
# Polynya salt input (Table 1): 2.5 mg m-2 s-1 over the shoreward 50 km; saltFlux is upward
# (model/inc/FFIELDS.h), so the input is negative.
SALT_POLYNYA, Y_POLYNYA = -2.5e-3, 50e3
# Sponge (Table 1): 50 km at the northern wall (rbcs; time scales in data.rbcs).
W_SPONGE = 50e3
# Surface relaxation (notes/pd1b): SST to the freezing point with the conductance of the
# two-equation ice model (gamma = 6e-5 m s-1, Schmidt et al. 2004); SSS to S_SURF in 30 days
# outside the polynya.
THETA_FREEZE, GAMMA_T, TAU_SSS = -1.9, 6e-5, 30 * 86400.0
# Slope-aware GM (Dettling Eq. 24): kappa_GM = kappa_high / (1 + EPS_C s_topo).
EPS_C = 800.0
# Profiles [verificar]: degC and g/kg (psu), depths in m; uniform in y.
T_SURF, T_CDW = -1.8, 0.8
S_SURF, S_CDW = 34.2, 34.7
Z_CLINE, W_CLINE = 300.0, 100.0


def depth(y_m) -> np.ndarray:
    """Water depth (m, positive) as a function of y only."""
    y = np.asarray(y_m, dtype=float)
    return H_SHELF + 0.5 * (H_DEEP - H_SHELF) * (1.0 + np.tanh((y - Y_SLOPE) / (0.25 * W_SLOPE)))


def topo_slope(y_m) -> np.ndarray:
    """|dH/dy| (dimensionless)."""
    y = np.asarray(y_m, dtype=float)
    L = 0.25 * W_SLOPE
    return 0.5 * (H_DEEP - H_SHELF) / L / np.cosh((y - Y_SLOPE) / L) ** 2


def bathymetry(g: Grid) -> np.ndarray:
    """MITgcm bathyFile (ny, nx): negative depths, identical in every column x, with land
    (0) in the first and last rows: MITgcm is periodic in y unless land closes the domain, so
    these rows are the coast (south) and the northern wall. The channel stays periodic in x."""
    b = np.repeat(-depth(g.yc)[:, None], g.nx, axis=1)
    b[0] = 0.0
    b[-1] = 0.0
    return b


def j_qb(g: Grid, y_qb_m: float = 90e3) -> int:
    """Index (0-based, Python) of the v face at y = y_qb_m; the Fortran row is j_qb + 1."""
    j = int(round(y_qb_m / g.dx_m))
    if not np.isclose(g.yg[j], y_qb_m, rtol=0, atol=1e-6):
        raise ValueError(f"y = {y_qb_m} m is not a v face of a {g.dx_m} m grid")
    return j


def wind_tau_x(y_m) -> np.ndarray:
    """Zonal wind stress (N m-2): 0 for y < Y_0, quarter sines rising to TAU_MAX at Y_W and
    falling to 0 at Y_N [verificar the exact shape]."""
    y = np.asarray(y_m, dtype=float)
    up = np.sin(0.5 * np.pi * np.clip((y - Y_0) / (Y_W - Y_0), 0.0, 1.0))
    down = np.sin(0.5 * np.pi * np.clip((Y_N - y) / (Y_N - Y_W), 0.0, 1.0))
    return TAU_MAX * np.where(y <= Y_W, up, down)


def polynya_salt_flux(g: Grid) -> np.ndarray:
    """saltFluxFile (ny, nx), g m-2 s-1, upward positive."""
    f = np.where(g.yc < Y_POLYNYA, SALT_POLYNYA, 0.0)
    return np.repeat(f[:, None], g.nx, axis=1)


def initial_ts(g: Grid) -> tuple[np.ndarray, np.ndarray]:
    """Initial (and sponge) theta (degC) and salinity, (nr, ny, nx), uniform in y and x:
    a tanh thermocline from the cold, fresh surface to CDW; statically stable (linear EOS)."""
    z = g.zc[:, None, None]
    c = 0.5 * (1.0 + np.tanh((z - Z_CLINE) / W_CLINE))
    T = T_SURF + (T_CDW - T_SURF) * c
    S = S_SURF + (S_CDW - S_SURF) * c
    shape = (g.nr, g.ny, g.nx)
    return np.broadcast_to(T, shape).copy(), np.broadcast_to(S, shape).copy()


def surface_relaxation(g: Grid) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """thetaClimFile, saltClimFile, lambdaThetaFile, lambdaSaltFile (ny, nx); lambdas in s-1."""
    shape = (g.ny, g.nx)
    sst = np.full(shape, THETA_FREEZE)
    sss = np.full(shape, S_SURF)
    lam_t = np.full(shape, GAMMA_T / float(g.dz_m[0]))
    lam_s = np.repeat(np.where(g.yc < Y_POLYNYA, 0.0, 1.0 / TAU_SSS)[:, None], g.nx, axis=1)
    return sst, sss, lam_t, lam_s


def sponge_mask(g: Grid, width_m: float = W_SPONGE) -> np.ndarray:
    """rbcs mask (nr, ny, nx): 0 south of Ly - width_m, rising linearly to 1 at the wall."""
    y0 = g.ny * g.dx_m - width_m
    m = np.clip((g.yc - y0) / width_m, 0.0, 1.0)
    return np.broadcast_to(m[None, :, None], (g.nr, g.ny, g.nx)).copy()


def gm_factor(y_m) -> np.ndarray:
    """Slope-aware GM factor (Dettling Eq. 24): 1 / (1 + EPS_C s_topo); multiplies
    GM_background_K through GM_bol2dFile. The Redi coefficient stays uniform."""
    return 1.0 / (1.0 + EPS_C * topo_slope(y_m))
