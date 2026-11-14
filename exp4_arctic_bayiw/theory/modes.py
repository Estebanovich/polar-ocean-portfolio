"""Linear internal-wave theory for Exp. 4: vertical modes and basic scales.

Conventions
-----------
* ``z_m`` is the vertical coordinate in metres, positive up, so the water
  column spans ``z = 0`` (surface) to ``z = -H`` (bottom); consistent with
  ``exp2_isomip_cavity`` profiles.  Either ordering of the array is accepted.
* ``N2`` is the squared buoyancy frequency in s^-2 at the same points.
* Hydrostatic, rigid lid, flat bottom: w = 0 at z = 0 and z = -H.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

OMEGA = 7.2921e-5  # s^-1, Earth rotation rate (sidereal) [verificar]
G = 9.81  # m s^-2
N2_FLOOR = 1e-10  # s^-2, floor to keep the eigenproblem positive definite


def coriolis(lat_deg: float) -> float:
    """f = 2 Omega sin(lat). Gill 1982 Sec. 4.5 [verificar]. Units: s^-1."""
    return float(2.0 * OMEGA * np.sin(np.deg2rad(lat_deg)))


def vertical_modes(z_m, N2, H: float, n_modes: int = 4):
    """Hydrostatic rigid-lid vertical modes.

    Solves  W'' + (N^2/c^2) W = 0,  W(0) = W(-H) = 0   (Gill 1982 Sec. 6.11;
    Cushman-Roisin & Beckers 2011 Sec. 11 [verificar]) for the vertical-velocity
    structure W_n(z).  Discretised on the (possibly non-uniform) grid with
    linear finite elements and a lumped mass matrix:
        A W = lambda M W,  lambda = 1/c^2 [s^2 m^-2],
        A = stiffness (dz^-1 second difference), M = diag(N2 * cell weight).
    ``H`` is the total depth; the grid is rescaled to span exactly H if it
    differs by less than 1 % (otherwise ValueError).

    Returns
    -------
    c : (n_modes,) phase speeds in m s^-1, descending (mode 1 = fastest).
    phi : (n_modes, nz) modes on the input grid, normalised so that
        <phi, phi> = 1 with the N2-weighted inner product
        ``modal_inner_product`` (<a,b> = int N2 a b dz / int N2 dz); sign
        fixed so the largest-magnitude lobe is positive.  Orthogonality of
        W modes holds under this weight (Sturm-Liouville weight N2).
    """
    z = np.asarray(z_m, dtype=float)
    N2 = np.asarray(N2, dtype=float)
    if z.shape != N2.shape or z.ndim != 1 or z.size < 4:
        raise ValueError("z_m and N2 must be 1-D arrays of equal length >= 4")
    order = np.argsort(z)  # ascending: bottom first
    zs, n2s = z[order], np.maximum(N2[order], N2_FLOOR)
    span = zs[-1] - zs[0]
    if abs(span - H) / H > 0.01:
        raise ValueError(f"grid spans {span:.3f} m but H = {H:.3f} m")
    dz = np.diff(zs)
    if np.any(dz <= 0):
        raise ValueError("z_m must be strictly monotonic")
    w = np.zeros_like(zs)  # cell weights
    w[:-1] += 0.5 * dz
    w[1:] += 0.5 * dz
    n = zs.size
    # Stiffness matrix on interior nodes (Dirichlet W = 0 at both ends).
    Afull = np.zeros((n, n))
    for i, d in enumerate(dz):
        Afull[i, i] += 1 / d
        Afull[i + 1, i + 1] += 1 / d
        Afull[i, i + 1] -= 1 / d
        Afull[i + 1, i] -= 1 / d
    A = Afull[1:-1, 1:-1]
    M = np.diag(n2s[1:-1] * w[1:-1])
    lam, vec = eigh(A, M, subset_by_index=[0, n_modes - 1])
    c = 1.0 / np.sqrt(lam)  # m s^-1
    phi_s = np.zeros((n_modes, n))
    norm_den = np.sum(n2s * w)
    for k in range(n_modes):
        p = np.zeros(n)
        p[1:-1] = vec[:, k]
        p /= np.sqrt(np.sum(n2s * w * p * p) / norm_den)
        if p[np.argmax(np.abs(p))] < 0:
            p = -p
        phi_s[k] = p
    phi = np.empty_like(phi_s)
    phi[:, order] = phi_s
    return c, phi


def modal_inner_product(a, b, z_m, N2) -> float:
    """<a,b> = int N2 a b dz / int N2 dz (trapezoid), the W-mode weight."""
    z = np.asarray(z_m, float)
    n2 = np.maximum(np.asarray(N2, float), N2_FLOOR)
    o = np.argsort(z)
    zs, n2s = z[o], n2[o]
    a, b = np.asarray(a, float)[o], np.asarray(b, float)[o]
    integ = getattr(np, "trapezoid", None) or np.trapz
    return float(integ(n2s * a * b, zs) / integ(n2s, zs))


def rossby_radius(c: float, f: float) -> float:
    """Internal Rossby radius R = c/|f|. Gill 1982 Sec. 7.3 [verificar]. Units: m."""
    return float(c / abs(f))


def burger(c: float, f: float, L: float) -> float:
    """Burger number Bu = (c/f)/L = R/L (here defined as the ratio, not squared).
    Cushman-Roisin & Beckers 2011 [verificar]. Dimensionless."""
    return float(rossby_radius(c, f) / L)


def kelvin_circuit_period(perimeter_m: float, c: float) -> float:
    """Period of a Kelvin wave circuit of a closed basin, T = P/c.
    Gill 1982 ch. 10 (coastal Kelvin waves) [verificar]. Units: s."""
    return float(perimeter_m / c)


def merian_period(L: float, H: float, n: int = 1, open_end: bool = True) -> float:
    """Merian seiche period with c = sqrt(g H) (long-wave speed, Gill 1982 ch. 5 (shallow-water waves) [verificar]).
    open_end (closed-open basin, quarter wave): T_n = 4 L / ((2n-1) c);
    closed-closed (half wave): T_n = 2 L / (n c).  Units: s."""
    c = np.sqrt(G * H)  # m s^-1
    if open_end:
        return float(4.0 * L / ((2 * n - 1) * c))
    return float(2.0 * L / (n * c))


def theory_summary(z, N2, H, lat_deg, L_basin, perimeter, n_modes: int = 4) -> dict:
    """All scales as plain floats/lists (JSON-serialisable) for figures/theory.json."""
    f = coriolis(lat_deg)
    c, _ = vertical_modes(z, N2, H, n_modes=n_modes)
    cl = [float(x) for x in c]
    return {
        "lat_deg": float(lat_deg),
        "H_m": float(H),
        "L_basin_m": float(L_basin),
        "perimeter_m": float(perimeter),
        "f_s-1": f,
        "inertial_period_h": float(2 * np.pi / abs(f) / 3600.0),
        "c_n": cl,
        "rossby_radius_m": [rossby_radius(x, f) for x in cl],
        "burger": [burger(x, f, L_basin) for x in cl],
        "kelvin_circuit_period_h": [kelvin_circuit_period(perimeter, x) / 3600.0 for x in cl],
        "c_barotropic_m_s-1": float(np.sqrt(G * H)),
        "merian_period_open_h": [merian_period(L_basin, H, n, True) / 3600.0 for n in (1, 2, 3)],
        "merian_period_closed_h": [merian_period(L_basin, H, n, False) / 3600.0 for n in (1, 2, 3)],
    }
