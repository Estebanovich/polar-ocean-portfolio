"""Section energy flux of Nash et al. (2005), exactly as the MSc thesis computes it.

Replica of ``thesis_copy/Output_IWBay/energy_flux_IWBay.py`` (ruling R16). Given the time
series at one section (velocity normal to the section at the T points, in-situ density from
the EOS of T and S, free surface):

1. rho_ref(z) = time mean of rho over all records, then mean over the section columns
   (Nash et al. 2005 eq. 2; thesis calcular_medias_run);
2. u_mean(z, x) = time mean of u per column (Nash eq. 5);
3. rho' = rho - rho_ref(z);  p~ = rho0 g eta + g int_z^0 rho' dz'  (Nash eq. 3, cell centres:
   layers above plus half the own layer, thesis presion_baroclinica);
4. p' = p~ - H^-1 int p~ dz  (Nash eq. 4), p_bt = H^-1 int p~ dz;
5. u' = (u - u_mean) - H^-1 int (u - u_mean) dz  (Nash eq. 5);
6. F_bc(t, x) = int p' u' dz  [W m^-1]  (Nash eq. 10);
7. per stage: F_bt(t, x) = p_bt' u_bt' H with p_bt', u_bt' the deviations from the stage mean
   (thesis procesar_etapa; this is the thesis definition, not Nash's); <F> = plain time mean
   over the stage; F_total(t) = sum_x F dl [W] with the real cell width dl.

Records are unique times; a stage is a closed interval [t0, t1], so two stages share their
boundary instant as in the thesis (chunks 6300 and 7200 both hold t = 216000 s).

Density: the thesis takes rho from the EOS of T and S; ``section_from_dataset`` passes
rhoConst + RHOAnoma instead (see its docstring: equal to 2e-13 kg m^-3 under the linear EOS,
class (c)). Dry columns (H = 0) contribute exactly zero.

Sign: + along the normal the caller passes. The project convention is + into the bay / fjord,
the same as energy_flux_IWBay.py (+y into the bay).
"""
from __future__ import annotations

import numpy as np


def _col_dz(dz, nz, nx):
    dz = np.asarray(dz, float)
    return np.broadcast_to(dz[:, None], (nz, nx)) if dz.ndim == 1 else dz


def section_flux(u_n, rho, eta, dz, dl, t, rho0: float = 999.8, g: float = 9.81,
                 stages=None) -> dict:
    """Nash et al. (2005) flux through a section, as energy_flux_IWBay.py.

    u_n (nt, nz, nx) [m s^-1]  normal velocity at the T points (+ along the chosen normal);
    rho (nt, nz, nx) [kg m^-3] in-situ density;  eta (nt, nx) [m];
    dz (nz,) or (nz, nx) [m]  wet layer thickness (0 = dry); H = sum dz per column;
    dl (nx,) [m] cell width along the section;  t (nt,) [s];
    stages: list of closed intervals (t0, t1) [s]; default one stage with every record.

    Returns {rho_ref (nz,), u_mean (nz, nx), F_bc (nt, nx), stages: [{t, F_bc, F_bt
    (nt_s, nx) [W m^-1], F_bc_mean, F_bt_mean (nx,) [W m^-1], F_bc_total, F_bt_total (nt_s,) [W],
    F_bc_total_mean, F_bt_total_mean [W]}]}.
    """
    u = np.asarray(u_n, float)
    rho = np.asarray(rho, float)
    eta = np.asarray(eta, float)
    t = np.asarray(t, float)
    dl = np.asarray(dl, float)
    nt, nz, nx = u.shape
    if np.any(np.diff(t) <= 0):
        raise ValueError("t must be strictly increasing (unique records)")
    d = _col_dz(dz, nz, nx)
    wet = d > 0
    H = d.sum(0)                                              # local depth [m]
    Hs = np.where(H > 0, H, 1.0)              # dry columns (H = 0): every term below is 0 there
    # Nash eq. 2: reference profile, time mean then mean over the (wet) columns
    rmean = rho.mean(0)
    rho_ref = np.where(wet, rmean, 0.0).sum(1) / np.maximum(wet.sum(1), 1)
    u_mean = u.mean(0)                                        # Nash eq. 5, u_bar(z, x)
    rp = np.where(wet[None], rho - rho_ref[None, :, None], 0.0)
    # Nash eq. 3 at cell centres: rho0 g eta + g (sum of the layers above + half the own layer)
    above = np.cumsum(rp * d[None], axis=1) - rp * d[None]
    p_t = rho0 * g * eta[:, None, :] + g * (above + 0.5 * rp * d[None])        # Pa
    p_bt = (p_t * d[None]).sum(1) / Hs                                         # Nash eq. 4
    p_bc = p_t - p_bt[:, None, :]
    up = u - u_mean[None]
    u_bc = up - ((up * d[None]).sum(1) / Hs)[:, None, :]                       # Nash eq. 5
    F_bc = (p_bc * u_bc * d[None]).sum(1)                                      # eq. 10, W m^-1
    u_bt = (u * d[None]).sum(1) / Hs
    if stages is None:
        stages = [(t[0], t[-1])]
    out = {"rho_ref": rho_ref, "u_mean": u_mean, "F_bc": F_bc, "stages": []}
    for t0, t1 in stages:
        m = (t >= t0) & (t <= t1)
        if not m.any():
            raise ValueError(f"stage [{t0}, {t1}] s holds no record")
        Fb = F_bc[m]
        # thesis F_bt = p_bt' u_bt' H, deviations from the stage mean (procesar_etapa)
        Ft = (p_bt[m] - p_bt[m].mean(0)) * (u_bt[m] - u_bt[m].mean(0)) * H
        out["stages"].append({
            "t": t[m], "F_bc": Fb, "F_bt": Ft,
            "F_bc_mean": Fb.mean(0), "F_bt_mean": Ft.mean(0),
            "F_bc_total": Fb @ dl, "F_bt_total": Ft @ dl,
            "F_bc_total_mean": float((Fb @ dl).mean()), "F_bt_total_mean": float((Ft @ dl).mean()),
        })
    return out


def t_point_normal(vel, normal: str, idx: int, along: slice):
    """Velocity normal to a section at its T points, as the thesis interpolation
    (energy_flux_IWBay._interpolar_UV_Tpoint): for a row section ('y', row idx) the mean of the
    V faces south and north of the cell, 0.5 (V[j] + V[j+1]); for a column section ('x', column
    idx) 0.5 (U[i] + U[i+1]). vel: (..., ny, nx) with MITgcm MDS face indexing (V[j] = south
    face of cell j, U[i] = west face of cell i)."""
    v = np.asarray(vel)
    if normal == "y":
        return 0.5 * (v[..., idx, along] + v[..., idx + 1, along])
    if normal == "x":
        return 0.5 * (v[..., along, idx] + v[..., along, idx + 1])
    raise ValueError("normal must be 'x' or 'y'")


def section_from_dataset(ds, normal: str, idx: int, along: slice, sign: float = 1.0,
                         time=slice(None)) -> dict:
    """Inputs of section_flux for a straight grid section of a ``load.run`` dataset.

    normal 'y': row ``idx``, columns ``along``, positive +y times ``sign``; 'x': column ``idx``,
    rows ``along``, positive +x times ``sign``. dz = drF hFacC, dl = dxF (row) or dyF (column).

    Density: rho = rhoConst + RHOAnoma from the dataset (load.run), not EOS(T, S) as the thesis
    script (R16 point 2): a literal departure, needed because the Exp. 4 MDS output has no salt.
    Under the linear EOS with uniform tRef/sRef the two agree to 2e-13 kg m^-3 on R1 (checked
    on every wet cell at five records, Task 5); difference class (c), numerical choice."""
    if normal == "y":
        sel_c = dict(j=idx, i=along)
        vel = ds["V"].isel(time=time, j=slice(idx, idx + 2), i=along).values
        un = t_point_normal(vel, "y", 0, slice(None))
        dl = ds["dxF"].isel(**sel_c).values
    elif normal == "x":
        sel_c = dict(j=along, i=idx)
        vel = ds["U"].isel(time=time, j=along, i=slice(idx, idx + 2)).values
        un = t_point_normal(vel, "x", 0, slice(None))
        dl = ds["dyF"].isel(**sel_c).values
    else:
        raise ValueError("normal must be 'x' or 'y'")
    dz = ds["drF"].values[:, None] * ds["hFacC"].isel(**sel_c).values
    return {"u_n": sign * un, "rho": ds["RHO"].isel(time=time, **sel_c).values,
            "eta": ds["ETA"].isel(time=time, **sel_c).values, "dz": dz, "dl": dl,
            "t": ds["time"].isel(time=time).values.astype(float)}
