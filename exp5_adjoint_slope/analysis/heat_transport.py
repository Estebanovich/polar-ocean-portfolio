"""J: onshore heat transport across the J line (spec Sec. 5; option A, notes/v4_revision_posthoc.md),
in TW.

    J = -(rho0 cp / N) sum_t [ sum_i sum_k (v hFacS drF + v* hFacS drF) (thetabar - theta_f) dxG
                               + sum_i sum_k F_gm ] * 1e-12,

over the N model steps of the cost window (the last `lastinterval` of the run). thetabar is the
mean of theta in the two cells that share the v face (row j, the face south of cell j). v* is the
GM bolus velocity, v* hFacS drF = PsiY(k+1) - PsiY(k) (pkg/gmredi/gmredi_residual_flow.F:76-81;
flipSign4LHCoord = +1). F_gm is the meridional theta flux of the GM/Redi tensor through the face,
degC m3 s-1 (pkg/gmredi/gmredi_ytransport.F:128-186); in the advective form it is the Redi part:

    F_gm = -yA [ Kvy (theta_j - theta_j-1) / dyC + Kvz dTdz ],   yA = dxG drF hFacS.

J > 0 is
southward (onshore). This is the same sum as config/chan10_tap/code_tap/cost_accumulate_mean.F;
units W * 1e-12 = TW.

Arrays are (nt, nr, ny, nx) for time series, (nr, ny, nx) for hFacS and maskC; j is the 0-based
row of the face."""
from __future__ import annotations

import numpy as np

TW = 1e-12


def bolus_transport(psi_y: np.ndarray) -> np.ndarray:
    """v* hFacS drF (m2 s-1) at every v point from GM_PsiY (same shape, (..., nr, ny, nx))."""
    psi = np.asarray(psi_y, dtype=float)
    below = np.zeros_like(psi)
    below[..., :-1, :, :] = psi[..., 1:, :, :]          # PsiY(k+1); PsiY(Nr+1) = 0
    return below - psi


def gm_yflux(theta, kvy, hfac_s, drf, dxg, dyc, j: int, kvz=None, drc=None,
             mask_c=None) -> np.ndarray:
    """GM/Redi meridional theta flux through the faces of row j, (nt, nr, nx), degC m3 s-1.

    Kvy, Kvz are the GM_Kvy, GM_Kvz diagnostics (m2 s-1). The Kvz term is used only when kvz
    is given (GM_ExtraDiag; it is zero when the Redi and GM coefficients are equal in the
    skew-flux form, pkg/gmredi/gmredi_readparms.F:252-261); it then needs drC (nr,) and maskC."""
    th = np.asarray(theta, dtype=float)
    hf = np.asarray(hfac_s, dtype=float)[:, j, :]                     # (nr, nx)
    ya = np.asarray(dxg, dtype=float) * np.asarray(drf, dtype=float)[:, None] * hf
    dthdy = (th[:, :, j, :] - th[:, :, j - 1, :]) / np.asarray(dyc, dtype=float)
    f = -ya * np.asarray(kvy, dtype=float)[:, :, j, :] * dthdy
    if kvz is not None:
        nr = th.shape[1]
        mc = np.asarray(mask_c, dtype=float)
        rdrc = 1.0 / np.asarray(drc, dtype=float)
        dtdz = np.zeros_like(f)
        for k in range(nr):
            km1, kp1 = max(k - 1, 0), min(k + 1, nr - 1)
            mp1 = 0.0 if k >= nr - 1 else 1.0
            up = 0.0
            dn = 0.0
            for jj in (j - 1, j):
                up = up + mc[km1, jj] * mc[k, jj] * (th[:, km1, jj] - th[:, k, jj])
                dn = dn + mc[k, jj] * mc[kp1, jj] * mp1 * (th[:, k, jj] - th[:, kp1, jj])
            dtdz[:, k] = 0.5 * (0.5 * rdrc[k] * up + 0.5 * rdrc[kp1] * dn)
        f = f - ya * np.asarray(kvz, dtype=float)[:, :, j, :] * dtdz
    return f * (hf > 0)


def J_series(v, theta, hfac_s, drf, dxg, j: int, gm_flux=None, psi_y=None,
             theta_f: float = -1.9, rho0: float = 1027.5, cp: float = 3994.0) -> np.ndarray:
    """J of every record (TW), shape (nt,). gm_flux: `gm_yflux` output, (nt, nr, nx);
    psi_y: GM_PsiY, (nt, nr, ny, nx), for the bolus term."""
    v = np.asarray(v, dtype=float)
    th = np.asarray(theta, dtype=float)
    hf = np.asarray(hfac_s, dtype=float)[:, j, :]           # (nr, nx)
    mask = hf > 0
    drf = np.asarray(drf, dtype=float)[:, None]
    tbar = 0.5 * (th[:, :, j, :] + th[:, :, j - 1, :])       # (nt, nr, nx)
    trans = v[:, :, j, :] * hf * drf                          # m2 s-1
    if psi_y is not None:
        trans = trans + bolus_transport(psi_y)[:, :, j, :]
    flux = trans * (tbar - theta_f) * mask * np.asarray(dxg, dtype=float)
    if gm_flux is not None:
        flux = flux + np.asarray(gm_flux, dtype=float) * mask
    return -rho0 * cp * flux.sum(axis=(1, 2)) * TW


def J(v, theta, hfac_s, drf, dxg, j: int, gm_flux=None, psi_y=None, theta_f: float = -1.9,
      rho0: float = 1027.5, cp: float = 3994.0) -> float:
    """Mean of `J_series` over the records (TW)."""
    return float(J_series(v, theta, hfac_s, drf, dxg, j, gm_flux, psi_y, theta_f, rho0,
                          cp).mean())


def J_from_means(vvelth, vvel, hfac_s, drf, dxg, j: int, dfy=None, vbt=None,
                 theta_f: float = -1.9, rho0: float = 1027.5, cp: float = 3994.0) -> float:
    """J (TW) from time means of the diagnostics VVELTH (v * face-mean theta), VVEL and, for
    the coarse model, DFyE_TH (Redi flux, degC m3 s-1; diffKhT = 0) and GM_vbT (bolus heat flux,
    degC m3 s-1; the bolus carries no net volume, so theta_f drops out), (nr, ny, nx) each.
    Equal to `J` of the per-step fields. The fine model has no GM: dfy = vbt = None."""
    hf = np.asarray(hfac_s, dtype=float)[:, j, :]
    mask = hf > 0
    drf = np.asarray(drf, dtype=float)[:, None]
    f = (np.asarray(vvelth)[:, j, :] - theta_f * np.asarray(vvel)[:, j, :]) * hf * drf * mask
    f = f * np.asarray(dxg, dtype=float)
    if dfy is not None:
        f = f + np.asarray(dfy, dtype=float)[:, j, :] * mask
    if vbt is not None:
        f = f + np.asarray(vbt, dtype=float)[:, j, :] * mask
    return float(-rho0 * cp * f.sum() * TW)
