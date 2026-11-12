"""Modal projection of horizontal velocity profiles and the local fraction q (Alford 2020).

Ruling R4: the horizontal-velocity modes are U_n = d phi_n / dz, with phi_n the vertical-velocity
(displacement) modes of ``theory.modes.vertical_modes`` (Gill 1982 sec. 6.11; Alford 2020
section 2: "the horizontal modes F(z) for velocity are calculated by continuity as the vertical
derivative"). On the linear finite-element grid of vertical_modes, phi_n lives on the cell faces
and d phi_n / dz is piecewise constant, i.e. one value per cell (the u points). Because
(d phi/dz)^T diag(dz) (d phi/dz) is the stiffness matrix of that eigenproblem, the discrete U_n
are orthogonal with the plain inner product (1/H) sum U_m U_n dz, to round-off. Mode 0 is the
barotropic mode U_0 = 1; sum U_n dz = phi_n(0) - phi_n(-H) = 0 for n >= 1.

Normalisation: (1/H) sum U_n^2 dz = 1 (dimensionless modes; coefficients in the units of u).

Reference: Alford, M. H. (2020), Global calculations of local and remote near-inertial-wave
dissipation, J. Phys. Oceanogr. 50 (November 2020), p. 3157 ff., doi:10.1175/JPO-D-20-0106.1
(thesis library file "Alford_2020b"); q is the fraction of the flux in modes 4 and higher.
"""
from __future__ import annotations

import numpy as np

from exp4_arctic_bayiw.theory.modes import vertical_modes


def u_modes(z_faces, N2_faces, n_modes: int = 10):
    """Horizontal-velocity modes on the cells between ``z_faces`` (0 ... -H, descending).

    Returns (c [m s^-1], U (n_modes + 1, nz) with U[0] = 1 barotropic, dz (nz,) [m])."""
    zf = np.asarray(z_faces, float)
    if not np.all(np.diff(zf) < 0):
        raise ValueError("z_faces must descend from the surface (0) to the bottom (-H)")
    H = zf[0] - zf[-1]
    c, phi = vertical_modes(zf, np.asarray(N2_faces, float), H, n_modes=n_modes)
    dz = -np.diff(zf)
    U = np.empty((n_modes + 1, dz.size))
    U[0] = 1.0
    U[1:] = (phi[:, :-1] - phi[:, 1:]) / dz[None]        # d phi / dz at the cell centres
    U[1:] /= np.sqrt((U[1:] ** 2 * dz).sum(1) / H)[:, None]
    for n in range(1, n_modes + 1):                      # sign: positive at the surface
        if U[n, 0] < 0:
            U[n] = -U[n]
    return c, U, dz


def project(u_profile, modes, dz) -> np.ndarray:
    """Modal coefficients a_n of u(z) = sum a_n U_n(z): a_n = sum u U_n dz / sum U_n^2 dz
    (least squares on orthogonal modes). Units of u."""
    U, dz, u = np.asarray(modes, float), np.asarray(dz, float), np.asarray(u_profile, float)
    return (U * dz) @ u / (U ** 2 * dz).sum(1)


def slab(z_c, H_ml: float) -> np.ndarray:
    """Slab mixed-layer step S(z) = 1 above z = -H_ml, 0 below (Alford 2020 eq. 3; Gill 1984)."""
    return (np.asarray(z_c, float) > -H_ml).astype(float)


def q_local(coeffs, n_far: int = 3) -> float:
    """Local fraction q = sum_{n > n_far} |a_n| / sum_{n >= 1} |a_n|: the share of the
    near-inertial flux in the modes that do not propagate far (Alford 2020 eqs. 4-5, flux into a
    mode proportional to |a_n|; modes 1-3 far, after Vic et al. 2019). The barotropic
    coefficient a_0 is excluded. Dimensionless."""
    a = np.abs(np.asarray(coeffs, float)[1:])
    tot = a.sum()
    if tot == 0:
        raise ValueError("profile has no baroclinic projection")
    return float(a[n_far:].sum() / tot)
