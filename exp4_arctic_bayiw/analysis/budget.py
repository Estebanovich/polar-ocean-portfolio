"""Barotropic and baroclinic energy budget of a box, Kang & Fringer (2012), as the MSc thesis.

PD0 Enmienda 1 (2026-11-12) and ruling R17: every flux of H4a-H4d and the budget use the
hydrostatic Kang & Fringer (2012, eqs. 13-17; q = 0) decomposition of
``thesis_copy/Output_IWBay/balance_kang_fringer.py``:

    dE_bt/dt + div(F_bt + A_bt) = P_bt - C - D_bt
    dE_bc/dt + div(F_bc + A_bc) = P_bc + C - D_bc

with u = u_bar + u' (u_bar the instantaneous depth mean of the column) and rho' = rho - rho_b(z),
rho_b the initial rest profile (fondo_kang_fringer.py), the same in every term:

    E_bt = 1/2 rho0 sum H |u_bar|^2 dA + 1/2 rho0 g sum eta^2 dA                       [J]
    E_bc = 1/2 rho0 sum |u'|^2 dV + sum E'_p dV, E'_p the exact local APE (KF12 eq. 12)  [J]
           (and, for comparison only, the linear APE 1/2 g^2 rho'^2 / (rho0 N_b^2))
    p~   = rho0 g eta + g int_z^0 rho' dz',  p_bt = depth mean of p~,  p' = p~ - p_bt     [Pa]
    F    = pressure work on the box faces: F_bt = H u_n,bar p_bt, F_bc = sum u_n' p' dz  [W]
    A    = energy advection on the faces (KF12 eqs. 15-16): A_bt = U_n H E_k0,
           A_bc = sum u_n (E'_k + E'_k0 + E'_p) dz                                        [W]
    C    = sum rho' g W_bt dV, W_bt = -div((H + z) u_bar) in the discrete C-grid form of
           MITgcm's integrate_for_w with u replaced by the column mean (conversion_bt_2D.py) [W]
    P_bt = sum tau . u_bar dA,  P_bc = sum tau . (u_surf - u_bar) dA, on the U/V points of the
           box; a face on the box edge counts with weight 1/2                              [W]
    D    residual (dissipation, energy diffusion, numerical error); D includes bottom drag,
         Smagorinsky viscosity and KPP in the fjord runs.

Face values (as the thesis at the bay mouth): u_n on the face, layer thickness drF hFacW/S of the
face, p~ and E'_p as the mean of the two adjacent cells, the tangential velocity as the mean of
the four surrounding faces. Cell-centre velocities are the mean of the two faces (hydrostatic:
no w^2 in the KE).

Sign convention (project): F and A are positive INTO the box (bay / inner fjord), as the Nash
flux of flux.py. balance_kang_fringer.py writes + = out of the bay, so F_thesis = -F_here. D is
the same in both: D = P -/+ C - dE + F_in.
"""
from __future__ import annotations

from collections import deque

import numpy as np

G = 9.81


# ---- pure pieces (balance_kang_fringer.py / conversion_bt_2D.py) ---------------------------------

def pressure_tilde(rp, eta, dz, rho0, g=G):
    """p~ at cell centres [Pa]: rho0 g eta + g (sum_{m<k} rho' dz + 1/2 rho'_k dz_k)
    (Nash et al. 2005 eq. 3; thesis presion_tilde). rp, dz: (nz, ...)."""
    layer = rp * dz
    above = np.cumsum(layer, axis=0) - layer
    return rho0 * g * eta[None] + g * (above + 0.5 * layer)


def face_pressure_work(un, p, dz):
    """(H u_bar p_bar, sum u' p' dz) per face column [W m^-1]; un, p, dz: (nz, n)
    (KF12 eqs. 15-16 pressure work; thesis flujo_cara)."""
    H = dz.sum(0)
    Hs = np.where(H > 0, H, 1.0)
    ub, pb = (un * dz).sum(0) / Hs, (p * dz).sum(0) / Hs
    return H * ub * pb, ((un - ub) * (p - pb) * dz).sum(0)


def face_advection(ut, un, dz, ep, rho0):
    """Energy advection through a face [W m^-1] (KF12 eqs. 15-16; thesis flujo_advectivo_cara):
    (U_n H 1/2 rho0 (U_t^2 + U_n^2), sum u_n (1/2 rho0 |u'|^2 + rho0 U.u' + E'_p) dz).
    ut tangential, un normal velocity; ut, un, dz, ep: (nz, n)."""
    H = dz.sum(0)
    Hs = np.where(H > 0, H, 1.0)
    Ut, Un = (ut * dz).sum(0) / Hs, (un * dz).sum(0) / Hs
    utp, unp = ut - Ut, un - Un
    F0 = Un * H * 0.5 * rho0 * (Ut ** 2 + Un ** 2)
    E1 = 0.5 * rho0 * (utp ** 2 + unp ** 2) + rho0 * (Ut * utp + Un * unp) + ep
    return F0, (un * E1 * dz).sum(0)


def ke_bt_bc(u, v, vol, rho0):
    """(1/2 rho0 sum H |u_bar|^2 dA, 1/2 rho0 sum |u'|^2 dV) [J]; u, v, vol (nz, ...) at the
    cell centres (KF12 eq. 13-14 kinetic parts; thesis ke_bt_bc). H dA = sum_k vol."""
    H = vol.sum(0)
    Hs = np.where(H > 0, H, 1.0)
    ub, vb = (u * vol).sum(0) / Hs, (v * vol).sum(0) / Hs
    bt = 0.5 * rho0 * ((ub ** 2 + vb ** 2) * H).sum()
    bc = 0.5 * rho0 * (((u - ub) ** 2 + (v - vb) ** 2) * vol).sum()
    return bt, bc


def pe_bt(eta, dA, rho0, g=G):
    """1/2 rho0 g sum eta^2 dA [J]."""
    return 0.5 * rho0 * g * float((eta ** 2 * dA).sum())


def ape_linear(rp, n2, vol, rho0, g=G):
    """sum 1/2 g^2 rho'^2 / (rho0 N_b^2) dV [J] (linear APE, comparison only)."""
    n2 = np.asarray(n2).reshape(np.shape(n2) + (1,) * (np.ndim(rp) - np.ndim(n2)))
    return float((0.5 * g ** 2 * rp ** 2 / (rho0 * n2) * vol).sum())


def ape_exact(rho_t, z, rho_b, z_b, g=G):
    """Exact local APE [J m^-3] (Kang & Fringer 2012 eq. 12): g int_{z_r}^{z} [rho - rho_b(z')] dz'
    with z_r the rest level where rho_b(z_r) = rho. rho_b is piecewise linear between the centres
    z_b (negative, surface to bottom) and constant outside them (thesis ape_exacta, same
    arithmetic). With linear rho_b it equals g^2 rho'^2 / (2 rho0 N^2)."""
    d_b = -np.asarray(z_b, float)
    rb = np.maximum.accumulate(np.asarray(rho_b, float))       # non-decreasing with depth
    slope = np.diff(rb) / np.diff(d_b)
    Rn = np.concatenate(([0.0], np.cumsum(0.5 * (rb[1:] + rb[:-1]) * np.diff(d_b))))

    def R(d):                                                  # int_{d_b[0]}^{d} rho_b, exact
        k = np.clip(np.searchsorted(d_b, d, side="right") - 1, 0, d_b.size - 2)
        x = np.clip(d, d_b[0], d_b[-1]) - d_b[k]
        inside = Rn[k] + rb[k] * x + 0.5 * slope[k] * x ** 2
        return (inside + rb[0] * np.minimum(d - d_b[0], 0.0)
                + rb[-1] * np.maximum(d - d_b[-1], 0.0))

    rt = np.asarray(rho_t, float)
    d = -np.asarray(z, float) * np.ones_like(rt)
    k = np.clip(np.searchsorted(rb, rt, side="right") - 1, 0, rb.size - 2)
    rise = rb[k + 1] - rb[k]
    d_r = d_b[k] + np.where(rise > 0, (rt - rb[k]) / np.where(rise > 0, rise, 1.0), 0.0) \
        * (d_b[k + 1] - d_b[k])
    d_r = np.where(rt <= rb[0], d_b[0], np.where(rt >= rb[-1], d_b[-1], d_r))
    return g * (rt * (d_r - d) - (R(d_r) - R(d)))


def areas_below_face(drF, hFacW, hFacS, dyG, dxG):
    """A_W,k, A_S,k: area of each face below the top of level k [m^2] (conversion_bt_2D.py).
    hFacW (nz, ny, nx+1), hFacS (nz, ny+1, nx)."""
    dW = np.asarray(drF)[:, None, None] * hFacW
    dS = np.asarray(drF)[:, None, None] * hFacS
    return dW[::-1].cumsum(0)[::-1] * dyG[None], dS[::-1].cumsum(0)[::-1] * dxG[None]


def w_bt_center(Ub, Vb, AW, AS, rA):
    """Barotropic vertical velocity at the cell centres (nz, ny, nx) [m s^-1]:
    W_top,k = -[d_x(U_bar A_W,k) + d_y(V_bar A_S,k)] / rA (top face of level k, the discrete
    -div((H + z) u_bar)); centre = 1/2 (W_top,k + W_top,k+1), bottom face 0
    (conversion_bt_2D w_bt_caras, w_bt_centro). Ub (ny, nx+1), Vb (ny+1, nx)."""
    fx = Ub[None] * AW
    fy = Vb[None] * AS
    Wt = -((fx[..., 1:] - fx[..., :-1]) + (fy[:, 1:, :] - fy[:, :-1, :])) / rA[None]
    Wb = np.concatenate((Wt[1:], np.zeros_like(Wt[:1])), axis=0)
    return 0.5 * (Wt + Wb)


def wind_work(tau, u_sup, u_bar, dA):
    """(sum tau u_bar dA, sum tau (u_sup - u_bar) dA) [W] (thesis trabajo_viento)."""
    return float((tau * u_bar * dA).sum()), float((tau * (u_sup - u_bar) * dA).sum())


def background_profile(rho0_field, wet, tol=1e-6):
    """rho_b(z) [kg m^-3] = mean of the wet cells of each level at t = 0 (fondo_kang_fringer
    perfil_fondo); ValueError if a level varies by more than ``tol`` in the horizontal."""
    R = np.where(wet, np.asarray(rho0_field, float), np.nan).reshape(wet.shape[0], -1)
    rb = np.nanmean(R, axis=1)
    dev = np.nanmax(np.abs(R - rb[:, None]), axis=1)
    if np.nanmax(dev) > tol:
        raise ValueError(f"initial state not horizontally uniform: deviation {np.nanmax(dev):.3g}")
    return rb


def n2_background(rho_b, RC, rho0, g=G):
    """N_b^2 = -(g / rho0) d rho_b / dz at the centres RC [s^-2] (fondo_kang_fringer n2_fondo)."""
    return -g / rho0 * np.gradient(np.asarray(rho_b, float), np.asarray(RC, float))


# ---- box geometry -----------------------------------------------------------------------------

def box_faces(mask, wet=None, allow_edge: bool = False) -> dict:
    """Open faces of a box. {'x': [(j, i, s)] west faces of cell (j, i) (U points),
    'y': [(j, i, s)] south faces (V points)}; s = +1 if +x / +y points into the box. With
    ``wet``, only the faces whose outer cell is wet (water). A box touching the domain edge has
    faces on the model boundary that this list cannot hold (an open boundary would carry a flux):
    ValueError unless ``allow_edge=True``, in which case those faces are left out."""
    m = np.asarray(mask, bool)
    if not allow_edge and (m[0].any() or m[-1].any() or m[:, 0].any() or m[:, -1].any()):
        raise ValueError("the box touches the domain edge; its boundary faces would be dropped "
                         "(pass allow_edge=True to accept that)")
    w = np.ones_like(m) if wet is None else np.asarray(wet, bool)
    out = {"x": [], "y": []}
    jj, ii = np.nonzero(m[:, 1:] != m[:, :-1])
    for j, i in zip(jj, ii + 1):
        if w[j, i] and w[j, i - 1]:
            out["x"].append((int(j), int(i), 1 if m[j, i] else -1))
    jj, ii = np.nonzero(m[1:] != m[:-1])
    for j, i in zip(jj + 1, ii):
        if w[j, i] and w[j - 1, i]:
            out["y"].append((int(j), int(i), 1 if m[j, i] else -1))
    return out


def _orient(ax, ay, bx, by, cx, cy):
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def _crosses(c1x, c1y, c2x, c2y, p0, p1):
    """Does the segment between two cell centres cross the section p0-p1?"""
    o1 = _orient(p0[0], p0[1], p1[0], p1[1], c1x, c1y)
    o2 = _orient(p0[0], p0[1], p1[0], p1[1], c2x, c2y)
    o3 = _orient(c1x, c1y, c2x, c2y, p0[0], p0[1])
    o4 = _orient(c1x, c1y, c2x, c2y, p1[0], p1[1])
    return (o1 * o2 <= 0) & (o3 * o4 <= 0) & ~((o1 == 0) & (o2 == 0))


def box_from_sections(wet, xc, yc, sections: dict, seed_in, seed_out, extend_m: float = 0.0):
    """Box closed by straight sections on grid faces (a staircase), e.g. the inner fjord closed
    at Akselsundet and Mariasundet.

    A face is cut if the segment joining its two cell centres crosses a section line (extended
    by ``extend_m`` at each end, so a line that ends at the coast in the coastline polygon still
    reaches the model land). The box is the wet region 4-connected to ``seed_in`` without
    crossing a cut face. RuntimeError if it reaches ``seed_out`` (not closed) or if an open face
    of the box is not a cut face. Returns (mask (ny, nx), {name: {'x': bool (ny, nx) U faces,
    'y': bool (ny, nx) V faces}}) with the faces of each section on the box edge."""
    wet = np.asarray(wet, bool)
    ny, nx = wet.shape
    X, Y = np.meshgrid(np.asarray(xc, float), np.asarray(yc, float))
    cut_x = np.zeros((ny, nx), bool)          # west face of (j, i)
    cut_y = np.zeros((ny, nx), bool)          # south face of (j, i)
    per = {}
    for name, (p0, p1) in sections.items():
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
        e = (p1 - p0) / np.hypot(*(p1 - p0))
        q0, q1 = p0 - extend_m * e, p1 + extend_m * e
        cx = np.zeros((ny, nx), bool)
        cy = np.zeros((ny, nx), bool)
        cx[:, 1:] = _crosses(X[:, :-1], Y[:, :-1], X[:, 1:], Y[:, 1:], q0, q1)
        cy[1:] = _crosses(X[:-1], Y[:-1], X[1:], Y[1:], q0, q1)
        per[name] = (cx, cy)
        cut_x |= cx
        cut_y |= cy
    js, is_ = (int(np.abs(np.asarray(yc) - seed_in[1]).argmin()),
               int(np.abs(np.asarray(xc) - seed_in[0]).argmin()))
    jo, io = (int(np.abs(np.asarray(yc) - seed_out[1]).argmin()),
              int(np.abs(np.asarray(xc) - seed_out[0]).argmin()))
    if not wet[js, is_] or not wet[jo, io]:
        raise RuntimeError("a seed is on land")
    mask = np.zeros((ny, nx), bool)
    mask[js, is_] = True
    q = deque([(js, is_)])
    while q:
        j, i = q.popleft()
        for dj, di, open_ in ((0, 1, i + 1 < nx and not cut_x[j, min(i + 1, nx - 1)]),
                              (0, -1, i > 0 and not cut_x[j, i]),
                              (1, 0, j + 1 < ny and not cut_y[min(j + 1, ny - 1), i]),
                              (-1, 0, j > 0 and not cut_y[j, i])):
            a, b = j + dj, i + di
            if open_ and wet[a, b] and not mask[a, b]:
                mask[a, b] = True
                q.append((a, b))
    if mask[jo, io]:
        raise RuntimeError("the sections do not close the box: not closed (seed_out reached)")
    f = box_faces(mask, wet=wet)
    groups = {n: {"x": np.zeros((ny, nx), bool), "y": np.zeros((ny, nx), bool)} for n in sections}
    for kind, lst in (("x", f["x"]), ("y", f["y"])):
        for j, i, _ in lst:
            hit = [n for n, c in per.items() if c[0 if kind == "x" else 1][j, i]]
            if not hit:
                raise RuntimeError(f"open {kind} face ({j}, {i}) of the box is not on a section")
            groups[hit[0]][kind][j, i] = True
    return mask, groups


# ---- the budget on a run ------------------------------------------------------------------------

def _ext(a, j0, j1, i0, i1):
    """a[..., j0:j1, i0:i1] with zero padding where the range leaves the domain."""
    ny, nx = a.shape[-2:]
    out = np.zeros(a.shape[:-2] + (j1 - j0, i1 - i0))
    ja, jb, ia, ib = max(j0, 0), min(j1, ny), max(i0, 0), min(i1, nx)
    out[..., ja - j0:jb - j0, ia - i0:ib - i0] = a[..., ja:jb, ia:ib]
    return out


def kf_budget(ds, mask, tau=None, rho_b=None, face_groups=None, time=slice(None),
              bg_tol=1e-6, conversion_masks=None) -> dict:
    """Kang & Fringer (2012) budget terms of the box ``mask`` (ny, nx) for every record.

    ds: a ``load.run`` dataset (U, V, RHO, ETA and the grid; attrs rho0, g). tau: None (TAUX,
    TAUY of ds if present, else no wind), or a function t [s] -> (taux at U points, tauy at V
    points) (ny, nx) [N m^-2], e.g. ``lambda t: load.forcing_stress(ds, t)``. rho_b: background
    profile (default: background_profile of RHO at the first record). face_groups: {name: {'x',
    'y': bool (ny, nx)}} to split F and A by section. Returns time series (nt,): t [s], ke_bt,
    pe_bt, ke_bc, ape, ape_linear [J]; F_bt, F_bc, A_bt, A_bc (+ into the box), C, P_bt, P_bc [W];
    and F_bt_<name>, ... per group; C_<name> for each extra ``conversion_masks`` {name: mask}
    (cells within one cell of the box, e.g. the thesis C_escalon row in front of the bay mouth);
    plus rho_b, n2 and the face list."""
    rho0, g = float(ds.attrs["rho0"]), float(ds.attrs.get("g", G))
    m_full = np.asarray(mask, bool)
    ny, nx = m_full.shape
    jj, ii = np.nonzero(m_full)
    j0, j1, i0, i1 = jj.min() - 1, jj.max() + 2, ii.min() - 1, ii.max() + 2   # cells + halo
    drF = ds["drF"].values
    RC = ds["RC"].values
    nz = drF.size
    hC = _ext(ds["hFacC"].values, j0, j1, i0, i1)
    hW = _ext(ds["hFacW"].values, j0, j1, i0, i1 + 1)
    hS = _ext(ds["hFacS"].values, j0, j1 + 1, i0, i1)
    rA = _ext(ds["rA"].values, j0, j1, i0, i1)
    dyG = _ext(ds["dyG"].values, j0, j1, i0, i1 + 1)
    dxG = _ext(ds["dxG"].values, j0, j1 + 1, i0, i1)
    rAw = _ext(ds["rAw"].values, j0, j1, i0, i1 + 1)
    rAs = _ext(ds["rAs"].values, j0, j1 + 1, i0, i1)
    m = _ext(m_full.astype(float), j0, j1, i0, i1) > 0
    nyr, nxr = m.shape
    dzC = drF[:, None, None] * hC
    wetC = dzC > 0
    vol = dzC * rA[None] * m[None]
    dW = drF[:, None, None] * hW
    dS = drF[:, None, None] * hS
    HW, HS = dW.sum(0), dS.sum(0)
    AW, AS = areas_below_face(drF, hW, hS, dyG, dxG)
    # wind weights: 1 inside, 1/2 on a face of the box edge (thesis peso_v)
    mp = np.zeros((nyr + 2, nxr + 2))
    mp[1:-1, 1:-1] = m
    wV = rAs * hS[0] * 0.5 * (mp[:-1, 1:-1] + mp[1:, 1:-1])
    wU = rAw * hW[0] * 0.5 * (mp[1:-1, :-1] + mp[1:-1, 1:])
    # open faces (local indices)
    faces = box_faces(m_full, wet=None)
    fy = np.array([(j - j0, i - i0, s) for j, i, s in faces["y"]], int).reshape(-1, 3)
    fx = np.array([(j - j0, i - i0, s) for j, i, s in faces["x"]], int).reshape(-1, 3)
    gy = {n: np.array([gr["y"][j, i] for j, i, _ in faces["y"]], bool)
          for n, gr in (face_groups or {}).items()}
    gx = {n: np.array([gr["x"][j, i] for j, i, _ in faces["x"]], bool)
          for n, gr in (face_groups or {}).items()}
    if rho_b is None:
        rho_b = background_profile(ds["RHO"].isel(time=0).values, ds["hFacC"].values > 0, bg_tol)
    rho_b = np.asarray(rho_b, float)
    n2 = n2_background(rho_b, RC, rho0, g)
    RC3 = RC[:, None, None] * np.ones((1, nyr, nxr))
    lin_ok = bool(np.all(n2[(vol > 0).any(axis=(1, 2))] > 0))
    times = ds["time"].values[time]
    tidx = np.arange(ds.sizes["time"])[time]
    has_tau = tau is None and "TAUX" in ds and "TAUY" in ds
    keys = ("ke_bt", "pe_bt", "ke_bc", "ape", "ape_linear", "F_bt", "F_bc", "A_bt", "A_bc", "C",
            "P_bt", "P_bc")
    out = {k: np.zeros(times.size) for k in keys}
    for n in gy:
        for k in ("F_bt", "F_bc", "A_bt", "A_bc"):
            out[f"{k}_{n}"] = np.zeros(times.size)
    vol_c = {}
    for name, cm in (conversion_masks or {}).items():
        cm = np.asarray(cm, bool)
        region = np.zeros((ny, nx), bool)
        region[max(j0, 0):j1, max(i0, 0):i1] = True
        if (cm & ~region).any():
            raise ValueError(f"conversion mask {name} leaves the box region (box + 1 cell)")
        vol_c[name] = dzC * rA[None] * (_ext(cm.astype(float), j0, j1, i0, i1) > 0)[None]
        out[f"C_{name}"] = np.zeros(times.size)
    for q, (n, tn) in enumerate(zip(tidx, times)):
        U = _ext_read(ds["U"], n, j0, j1, i0, i1 + 1, ny, nx)
        V = _ext_read(ds["V"], n, j0, j1 + 1, i0, i1, ny, nx)
        R = _ext_read(ds["RHO"], n, j0, j1, i0, i1, ny, nx)
        E = _ext_read(ds["ETA"], n, j0, j1, i0, i1, ny, nx)
        U = U * (hW > 0)
        V = V * (hS > 0)
        rp = np.where(wetC, R - rho_b[:, None, None], 0.0)
        pt = pressure_tilde(rp, E, dzC, rho0, g)
        uc = 0.5 * (U[..., :-1] + U[..., 1:])
        vc = 0.5 * (V[:, :-1] + V[:, 1:])
        out["ke_bt"][q], out["ke_bc"][q] = ke_bt_bc(uc, vc, vol, rho0)
        out["pe_bt"][q] = pe_bt(E * m, rA, rho0, g)
        ep = np.where(wetC, ape_exact(np.where(wetC, R, rho_b[:, None, None]), RC3, rho_b, RC, g),
                      0.0)
        out["ape"][q] = float((ep * vol).sum())
        # linear APE (comparison only) needs N_b^2 > 0 at every level that the box holds
        out["ape_linear"][q] = ape_linear(rp, n2, vol, rho0, g) if lin_ok else np.nan
        Ub = (U * dW).sum(0) / np.where(HW > 0, HW, np.inf)
        Vb = (V * dS).sum(0) / np.where(HS > 0, HS, np.inf)
        Wc = w_bt_center(Ub, Vb, AW, AS, rA)
        out["C"][q] = g * float((rp * Wc * vol).sum())
        for name, vc_ in vol_c.items():
            out[f"C_{name}"][q] = g * float((rp * Wc * vc_).sum())
        # faces normal to y (V points): south neighbour j-1, north neighbour j
        for kind, fl, gg in (("y", fy, gy), ("x", fx, gx)):
            if not len(fl):
                continue
            j, i, s = fl[:, 0], fl[:, 1], fl[:, 2].astype(float)
            if kind == "y":
                un, dz, wdt = V[:, j, i], dS[:, j, i], dxG[j, i]
                p = 0.5 * (pt[:, j - 1, i] + pt[:, j, i])
                epf = 0.5 * (ep[:, j - 1, i] + ep[:, j, i])
                ut = 0.25 * (U[:, j - 1, i] + U[:, j - 1, i + 1] + U[:, j, i] + U[:, j, i + 1])
            else:
                un, dz, wdt = U[:, j, i], dW[:, j, i], dyG[j, i]
                p = 0.5 * (pt[:, j, i - 1] + pt[:, j, i])
                epf = 0.5 * (ep[:, j, i - 1] + ep[:, j, i])
                ut = 0.25 * (V[:, j, i - 1] + V[:, j + 1, i - 1] + V[:, j, i] + V[:, j + 1, i])
            fbt, fbc = face_pressure_work(un, p, dz)
            a0, a1 = face_advection(ut, un, dz, epf, rho0)
            cols = {"F_bt": s * fbt * wdt, "F_bc": s * fbc * wdt, "A_bt": s * a0 * wdt,
                    "A_bc": s * a1 * wdt}
            for k, c in cols.items():
                out[k][q] += float(c.sum())
                for name, sel in gg.items():
                    out[f"{k}_{name}"][q] += float(c[sel].sum())
        # wind work
        if tau is not None:
            tx, ty = tau(float(tn))
        elif has_tau:
            tx = ds["TAUX"].isel(time=n).values
            ty = ds["TAUY"].isel(time=n).values
        else:
            tx = ty = None
        if tx is not None:
            txe = _ext(np.asarray(tx), j0, j1, i0, i1 + 1)
            tye = _ext(np.asarray(ty), j0, j1 + 1, i0, i1)
            a, b = wind_work(tye, V[0], Vb, wV)
            c, d = wind_work(txe, U[0], Ub, wU)
            out["P_bt"][q], out["P_bc"][q] = a + c, b + d
    out["t"] = times.astype(float)
    out["rho_b"] = rho_b
    out["n2"] = n2
    # PD0 Enmienda 2: the H4d numerator is the outward baroclinic pressure work alone
    out["F_bc_pressure_out"] = -out["F_bc"]
    for n in gy:
        out[f"F_bc_pressure_out_{n}"] = -out[f"F_bc_{n}"]
    out["ape_linear_defined"] = lin_ok
    out["faces"] = faces
    return out


def _ext_read(var, n, j0, j1, i0, i1, ny, nx):
    """Record n of a lazy (time, [k,] j, i) variable on rows j0:j1, columns i0:i1 (zero pad)."""
    ja, jb, ia, ib = max(j0, 0), min(j1, ny), max(i0, 0), min(i1, nx)
    a = var.isel(time=n, j=slice(ja, jb), i=slice(ia, ib)).values
    out = np.zeros(a.shape[:-2] + (j1 - j0, i1 - i0))
    out[..., ja - j0:jb - j0, ia - i0:ib - i0] = a
    return out


def cumulative(t_s, x):
    """Trapezoidal running integral from the first record [J for W] (thesis _acumular)."""
    t_s, x = np.asarray(t_s, float), np.asarray(x, float)
    return np.concatenate(([0.0], np.cumsum(0.5 * (x[1:] + x[:-1]) * np.diff(t_s))))


def f_bc_pressure_out(r: dict, section: str | None = None):
    """H4d numerator [J] (PD0 Enmienda 2, 2026-11-12): the running integral of the outward
    Kang & Fringer baroclinic pressure work, -int F_bc dt with F_bc = sum u_n' p' dz dl over the
    open faces (+ into the box in kf_budget), so + = out of the box. Advection (A_bc) is NOT
    included; F_bc + A_bc is descriptive only. ``section``: one face group (e.g. 'Akselsundet').
    eta = f_bc_pressure_out[-1] / int (P_bt + P_bc) dt; eta' with P_bc alone is descriptive."""
    key = "F_bc" if section is None else f"F_bc_{section}"
    return -cumulative(r["t"], r[key])


def closure(r: dict, ape: str = "ape") -> dict:
    """Cumulative terms [J] of both budgets and the residual D (thesis
    figura_energia_balance_kf.terminos, with F_in = -F_out):
    bt: D = P - C - dE + F,  bc: D = P + C - dE + F, F = pressure work + advection into the box."""
    t = r["t"]
    E = {"bt": r["ke_bt"] + r["pe_bt"], "bc": r["ke_bc"] + r[ape]}
    C = cumulative(t, r["C"])
    out = {}
    for comp, sgn in (("bt", -1.0), ("bc", 1.0)):
        Fp = cumulative(t, r[f"F_{comp}"])
        A = cumulative(t, r[f"A_{comp}"])
        x = dict(P=cumulative(t, r[f"P_{comp}"]), Fp=Fp, A=A, F=Fp + A, C=C,
                 dE=E[comp] - E[comp][0])
        x["D"] = x["P"] + sgn * C - x["dE"] + x["F"]
        out[comp] = x
    return out


def fjord_box(ds, report=None, cases=None, extend_m=None):
    """Inner fjord box of a fjord_vmf run: closed by the Akselsundet and Mariasundet section
    lines of the bathymetry (figures/bathy_report.json, ``sounds.*.p0/p1``, UTM 33N = the model
    x, y; gendata/bathy.sound_sections) with the inner/outer seeds of cases.yaml
    (configs.fjord_vmf.bathy.seeds east/west). The lines end on the NPI coastline, which the
    200 m land mask does not follow exactly: with no extension the box leaks past the line ends,
    so each line is extended by one grid cell (``extend_m`` default = dx; chosen in Task 5 before
    any campaign run, after the 0 m test leaked). Returns (mask, face groups)."""
    import json
    from pathlib import Path

    import yaml
    pkg = Path(__file__).resolve().parents[1]
    rep = json.loads(Path(report or pkg / "figures" / "bathy_report.json").read_text())
    spec = yaml.safe_load(Path(cases or pkg / "cases.yaml").read_text())
    seeds = spec["configs"]["fjord_vmf"]["bathy"]["seeds"]
    secs = {n: (v["p0"], v["p1"]) for n, v in rep["sounds"].items()}
    wet = ds["hFacC"].values[0] > 0
    xc, yc = ds["XC"].values[0], ds["YC"].values[:, 0]
    if extend_m is None:
        extend_m = float(np.median(np.diff(xc)))
    return box_from_sections(wet, xc, yc, secs, seeds["east"], seeds["west"], extend_m)


def fjord_budget(ds, mask=None, sections=None, seeds=None, extend_m=None, **kw) -> dict:
    """Budget of the inner fjord. The box is ``mask`` or, if None, the box closed by ``sections``
    ({name: (p0, p1)} in the model x, y [m]) with ``seeds`` = (inside, outside) points, or, if
    both are None, ``fjord_box(ds)`` (Akselsundet and Mariasundet). F and A are also given per
    section. Returns kf_budget's series plus 'closure' and 'mask'."""
    groups = kw.pop("face_groups", None)
    if mask is None and sections is None:
        mask, groups = fjord_box(ds, extend_m=extend_m)
    elif mask is None:
        wet = ds["hFacC"].values[0] > 0
        xc, yc = ds["XC"].values[0], ds["YC"].values[:, 0]
        mask, groups = box_from_sections(wet, xc, yc, sections, seeds[0], seeds[1],
                                         extend_m or 0.0)
    r = kf_budget(ds, mask, face_groups=groups, **kw)
    r["closure"] = closure(r)
    r["mask"] = mask
    return r


# ---- explicit Laplacian dissipation (descriptive; R1/R1k, ruling R17) -------------------------

def laplacian_dissipation(U, V, g, rho0, nu_h, nu_v):
    """Explicit dissipation of constant Laplacian viscosity in a rectangle [W], as the thesis
    disipacion_viscosa_bahia.disipacion (free-slip walls and bottom, no drag):
    eps = rho0 [nu_h (|grad_h u|^2 + |grad_h v|^2) + nu_v (|d_z u|^2 + |d_z v|^2)], C-grid
    differences (d_x u, d_y v at centres; d_y u, d_x v at interior corners between two wet faces;
    d_z at the interfaces between two wet faces). Split u = u_bar + u': returns (eps_bt,
    eps_bc_h, eps_bc_v, cross term, part of eps_bc in the top cell). Applies to constant
    viscAh/viscAr only (R1); with Smagorinsky or KPP the explicit dissipation is not this, and
    it stays in the residual D. U (nz, ny, nx+1), V (nz, ny+1, nx)."""
    mW, mS = g["hW"] > 0, g["hS"] > 0
    wW, wS = g["hW"] * g["drF"][:, None, None], g["hS"] * g["drF"][:, None, None]

    def mean_v(a, w):
        W = w.sum(0)
        return (a * w).sum(0) / np.where(W > 0, W, np.inf)

    Ub, Vb = mean_v(U, wW), mean_v(V, wS)
    Up, Vp = np.where(mW, U - Ub, 0.0), np.where(mS, V - Vb, 0.0)
    volC = g["hC"] * g["drF"][:, None, None] * g["rA"]

    def grad_h(u, v, vol, column=False):
        dux = (u[..., 1:] - u[..., :-1]) / g["dxF"]
        dvy = (v[..., 1:, :] - v[..., :-1, :]) / g["dyF"]
        s = ((dux ** 2 + dvy ** 2) * vol).sum()
        mu = (g["HW"][1:] > 0) & (g["HW"][:-1] > 0) if column else mW[..., 1:, :] & mW[..., :-1, :]
        mv = (g["HS"][:, 1:] > 0) & (g["HS"][:, :-1] > 0) if column else mS[..., :, 1:] & mS[..., :, :-1]
        duy = (u[..., 1:, :] - u[..., :-1, :]) / g["dyU"]
        dvx = (v[..., :, 1:] - v[..., :, :-1]) / g["dxV"]
        if column:
            s += (np.where(mu, duy ** 2, 0) * g["rAzU"] * np.minimum(g["HW"][1:], g["HW"][:-1])).sum()
            s += (np.where(mv, dvx ** 2, 0) * g["rAzV"] * np.minimum(g["HS"][:, 1:], g["HS"][:, :-1])).sum()
        else:
            dz = g["drF"][:, None, None]
            s += (np.where(mu, duy ** 2, 0) * g["rAzU"] * dz).sum()
            s += (np.where(mv, dvx ** 2, 0) * g["rAzV"] * dz).sum()
        return s

    bt = rho0 * nu_h * grad_h(Ub, Vb, g["HC"] * g["rA"], column=True)
    bc_h = rho0 * nu_h * grad_h(Up, Vp, volC)
    tot_h = rho0 * nu_h * grad_h(np.where(mW, U, 0.0), np.where(mS, V, 0.0), volC)

    def grad_v(u, m, rAf):
        ok = m[1:] & m[:-1]
        du = (u[1:] - u[:-1]) / g["drC"][1:-1, None, None]
        e = np.where(ok, du ** 2, 0.0) * g["drC"][1:-1, None, None] * rAf
        return e.sum(), e[0].sum()

    vu, vu0 = grad_v(Up, mW, g["rAw"])
    vv, vv0 = grad_v(Vp, mS, g["rAs"])
    bc_v = rho0 * nu_v * (vu + vv)
    top = np.arange(U.shape[0])[:, None, None] == 0
    sup = rho0 * nu_h * grad_h(np.where(top, Up, 0.0), np.where(top, Vp, 0.0), volC) \
        + rho0 * nu_v * (vu0 + vv0)
    return bt, bc_h, bc_v, tot_h - bt - bc_h, sup


def explicit_dissipation(ds, rows, cols, nu_h, nu_v, time=slice(None)) -> dict:
    """laplacian_dissipation of the rectangle rows [j0, j1) x cols [i0, i1) for every record.
    Series [W]: bt, bc_h, bc_v, cross, sup; t [s]."""
    (j0, j1), (i0, i1) = rows, cols
    G2 = lambda v, a, b, c, d: ds[v].isel(j=slice(a, b), i=slice(c, d)).values   # noqa: E731
    g = {"drF": ds["drF"].values, "drC": ds["drC"].values,
         "hC": G2("hFacC", j0, j1, i0, i1), "hW": G2("hFacW", j0, j1, i0, i1 + 1),
         "hS": G2("hFacS", j0, j1 + 1, i0, i1), "rA": G2("rA", j0, j1, i0, i1),
         "dxF": G2("dxF", j0, j1, i0, i1), "dyF": G2("dyF", j0, j1, i0, i1),
         "rAw": G2("rAw", j0, j1, i0, i1 + 1), "rAs": G2("rAs", j0, j1 + 1, i0, i1),
         "dyU": G2("dyU", j0 + 1, j1, i0, i1 + 1), "dxV": G2("dxV", j0, j1 + 1, i0 + 1, i1),
         "rAzU": G2("rAz", j0 + 1, j1, i0, i1 + 1), "rAzV": G2("rAz", j0, j1 + 1, i0 + 1, i1)}
    dz = g["drF"][:, None, None]
    g["HC"], g["HW"], g["HS"] = (g["hC"] * dz).sum(0), (g["hW"] * dz).sum(0), (g["hS"] * dz).sum(0)
    rho0 = float(ds.attrs["rho0"])
    tidx = np.arange(ds.sizes["time"])[time]
    out = {k: np.zeros(tidx.size) for k in ("bt", "bc_h", "bc_v", "cross", "sup")}
    for q, n in enumerate(tidx):
        U = ds["U"].isel(time=n, j=slice(j0, j1), i=slice(i0, i1 + 1)).values
        V = ds["V"].isel(time=n, j=slice(j0, j1 + 1), i=slice(i0, i1)).values
        for k, v in zip(("bt", "bc_h", "bc_v", "cross", "sup"),
                        laplacian_dissipation(U, V, g, rho0, nu_h, nu_v)):
            out[k][q] = v
    out["t"] = ds["time"].values[tidx].astype(float)
    return out
