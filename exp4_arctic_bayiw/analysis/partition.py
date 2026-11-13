"""Task 7, Step 1: band partition of the strait fluxes, box energy, wind work and station spectra
of every Exp. 4 run (figures/partition.json, figures/budget.json).

Definitions (PD0 sections 3 and 5, Enmiendas 1-2; task-7 rulings R25, fixed before any H4
number was computed):

* **Fluxes** are Kang & Fringer (2012, eqs. 13-17, hydrostatic) as ``budget.kf_budget``: pressure
  work F (F_bt = H u_n,bar p_bt, F_bc = sum u_n' p' dz) plus energy advection A (A_bt =
  U_n H 1/2 rho0 (U_t^2 + U_n^2), A_bc = sum u_n (1/2 rho0 |u'|^2 + rho0 U.u' + E'_p) dz) on the
  open faces of the box, summed per strait [W]. ``face_series`` keeps them per face and level so
  that they can be split by frequency. H4a-H4d are decided by the pressure work F_bc (and F_bt
  in H4b(ii)): PD0 Enmienda 2 (H4d) and Enmienda 3 (2026-11-13; H4a, H4b(ii), H4c(ii)).
  F_bc + A_bc is computed and reported as descriptive.
* **Band split of a flux**: the exact cospectral (Parseval) split of the time mean of a product,
  ``validate.band_flux`` (Task 6), applied to (u_n', p') for F_bc and to (u_n, E1) for A_bc,
  where E1 is the advected energy density per layer. No window; the low band holds the product
  of the means. Band sums equal -dt sum F over the records (the band total is the dt sum, not the
  trapezoid; they differ by the end terms dt (F_0 + F_n)/2).
* **Box energy without the low band** (H4b(i)): each field is high-passed at 2 df per point
  before the energy is formed (u, v, rho'; FFT, no window, bins with f >= the upper edge of the
  low band kept), then KE_bc = 1/2 rho0 sum |u' - u'_bar|^2 dV and the exact APE (KF12 eq. 12)
  are averaged over the records. The KE is split by band with Parseval, using
  sum_k dz (u_k)^2 - H (u_bar)^2 per band (the depth mean is linear, so it commutes with the
  filter). The linear APE 1/2 g^2 rho'^2 / (rho0 N_b^2) per band is a comparison only.
* **Wind work in a band** (H4c(i)): W_b = int sum tau . u_b dA dt, u_b the surface-level (k = 1)
  velocity at the U/V points band-passed per point, tau the stress the model applies
  (``load.forcing_stress``), with the U/V points and weights of ``budget.kf_budget`` (1/2 on the
  box edge), so that W over every band is int (P_bt + P_bc) dt.
* **Station band energy**: one-sided Hann periodogram (``spectra.psd``) of the surface u and v
  at the T point; band sums of P df (half-open bands, the last one closed at Nyquist).

Every equation cites its source in the function that evaluates it. Sign convention: kf_budget
and the face series are + INTO the box; the *_out_J values of this module are + OUT of the box
(as validate.flux_summary).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp4_arctic_bayiw.analysis import budget, spectra, validate as va  # noqa: E402

# Arneborg & Liljebladh (2001, JPO 31:2549, Gullmar fjord, sec. 5; thesis library
# Arneborg_Liljebladh_2001 ... .md, l. 601-605): "bandpass filtering the time series in the band
# 0.01-0.07 cph before determining the energies and energy fluxes". POST HOC descriptive (R24).
AL_BAND_CPH = (0.01, 0.07)


def _progress(what: str, r: int, n: int, every: int = 120, on: bool = True) -> None:
    if on and (r % every == 0 or r == n - 1):
        print(f"[partition] {what}: {r + 1}/{n}", file=sys.stderr, flush=True)


# ---- spectral helpers ---------------------------------------------------------------------------

def _band_of_bins(fk, bands: dict) -> dict:
    """{band: bool mask of the bins} with half-open bands [lo, hi), the highest one closed at its
    upper edge (Nyquist), as validate.band_flux. ValueError if a bin is in no band."""
    items = sorted(bands.items(), key=lambda kv: kv[1][0])
    out, cover = {}, np.zeros(fk.shape, bool)
    for q, (name, (lo, hi)) in enumerate(items):
        last = q == len(items) - 1
        sel = (fk >= lo) & ((fk <= hi + 1e-15) if last else (fk < hi))
        out[name] = sel
        cover |= sel
    if not cover.all():
        raise ValueError("the bands do not cover every frequency bin")
    return out


def _rfft_weights(n: int) -> np.ndarray:
    """Parseval weights of the rfft bins: 1 for DC (and Nyquist if n is even), 2 otherwise."""
    w = np.full(n // 2 + 1, 2.0)
    w[0] = 1.0
    if n % 2 == 0:
        w[-1] = 1.0
    return w


def band_power(x, dt: float, bands: dict) -> dict:
    """{band: time mean of x^2 carried by the band} along axis 0 (Parseval, no window, mean
    included): mean(x^2) = sum_k |X_k|^2 / n^2 over all FFT bins (x^2 units)."""
    x = np.asarray(x, float)
    n = x.shape[0]
    X = np.fft.rfft(x, axis=0)
    p = np.abs(X) ** 2 / n ** 2 * _rfft_weights(n).reshape((-1,) + (1,) * (x.ndim - 1))
    sel = _band_of_bins(np.fft.rfftfreq(n, d=dt), bands)
    return {b: p[s].sum(axis=0) for b, s in sel.items()}


def _keep(x, dt: float, keep_fn):
    x = np.asarray(x, float)
    n = x.shape[0]
    X = np.fft.rfft(x, axis=0)
    k = keep_fn(np.fft.rfftfreq(n, d=dt)).reshape((-1,) + (1,) * (x.ndim - 1))
    return np.fft.irfft(X * k, n=n, axis=0)


def band_pass(x, dt: float, lo: float, hi: float):
    """x with only the FFT bins lo <= f < hi kept (axis 0; no window; hi at Nyquist is closed)."""
    nyq = 0.5 / dt
    return _keep(x, dt, lambda f: (f >= lo) & ((f <= hi + 1e-15) if hi >= nyq - 1e-15 else (f < hi)))


def high_pass(x, dt: float, cut: float):
    """x without the bins f < cut (axis 0): the low band [0, 2 df) removed (H4b(i), ruling)."""
    return _keep(x, dt, lambda f: f >= cut)


# ---- strait fluxes per face (Kang & Fringer 2012, as budget.kf_budget) --------------------------

def face_series(ds, mask, face_groups: dict, rho_b=None, time=slice(None), progress=True) -> dict:
    """Per face and level, the pieces of the KF fluxes on the open faces of ``mask`` that belong to
    ``face_groups`` ({name: {'x', 'y': bool (ny, nx)}}), computed as budget.kf_budget:

    * up = u_n - u_n,bar, pp = p~ - p_bt (KF12 eqs. 15-16; p~ = rho0 g eta + g int rho' dz,
      Nash et al. 2005 eq. 3) [m s^-1, Pa]: F_bc = sum up pp dz [W m^-1];
    * un and e1 = 1/2 rho0 (u_t'^2 + u_n'^2) + rho0 (U_t u_t' + U_n u_n') + E'_p [J m^-3]
      (KF12 eq. 16, thesis flujo_advectivo_cara): A_bc = sum un e1 dz [W m^-1];
    * fbt = H u_n,bar p_bt and abt = U_n H 1/2 rho0 (U_t^2 + U_n^2) [W m^-1] (barotropic).

    Returns up, pp, un, e1 (nt, nz, nf); fbt, abt (nt, nf); dz (nz, nf) [m]; w (nf,) signed face
    width [m] (+ when +x/+y points into the box); group (nf,); t (nt,) [s]. Only a small
    rectangle around each group, and only its wet levels, is read per record."""
    rho0, g = float(ds.attrs["rho0"]), float(ds.attrs.get("g", budget.G))
    drF, RC = ds["drF"].values, ds["RC"].values
    hC, hW, hS = ds["hFacC"].values, ds["hFacW"].values, ds["hFacS"].values
    dxG, dyG = ds["dxG"].values, ds["dyG"].values
    ny, nx = hC.shape[1:]
    if rho_b is None:
        rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, hC > 0)
    rho_b = np.asarray(rho_b, float)
    faces = budget.box_faces(np.asarray(mask, bool), wet=None)
    sel = []
    for kind in ("x", "y"):
        for j, i, s in faces[kind]:
            hit = [n for n, gr in face_groups.items() if gr[kind][j, i]]
            if hit:
                sel.append((kind, j, i, s, hit[0]))
    if not sel:
        raise ValueError("no face of the box is in face_groups")
    tidx = np.arange(ds.sizes["time"])[time]
    nt, nz, nf = tidx.size, drF.size, len(sel)
    up, pp, un_, e1 = (np.zeros((nt, nz, nf)) for _ in range(4))
    fbt, abt = np.zeros((nt, nf)), np.zeros((nt, nf))
    dz = np.zeros((nz, nf))
    w = np.zeros(nf)
    for q, (kind, j, i, s, _) in enumerate(sel):
        dz[:, q] = drF * (hS[:, j, i] if kind == "y" else hW[:, j, i])
        w[q] = s * (dxG[j, i] if kind == "y" else dyG[j, i])
    groups: dict = {}
    for q, x in enumerate(sel):
        groups.setdefault(x[4], []).append(q)
    for name, qs in groups.items():
        js = [sel[q][1] for q in qs]
        is_ = [sel[q][2] for q in qs]
        j0, j1 = max(min(js) - 1, 0), min(max(js) + 2, ny)
        i0, i1 = max(min(is_) - 1, 0), min(max(is_) + 2, nx)
        nk = int(np.flatnonzero((hC[:, j0:j1, i0:i1] > 0).any(axis=(1, 2))).max()) + 1
        if np.any(dz[nk:, qs] != 0):
            raise AssertionError("face thickness below the last wet level of the rectangle")
        dzC = drF[:nk, None, None] * hC[:nk, j0:j1, i0:i1]
        wetC = dzC > 0
        mW, mS = hW[:nk, j0:j1, i0:i1] > 0, hS[:nk, j0:j1, i0:i1] > 0
        z3 = RC[:nk, None, None] * np.ones((1, j1 - j0, i1 - i0))
        rb3 = rho_b[:nk, None, None]
        for r, n_rec in enumerate(tidx):
            win = dict(time=int(n_rec), k=slice(0, nk), j=slice(j0, j1), i=slice(i0, i1))
            U = ds["U"].isel(**win).values * mW            # kf_budget: U (hW > 0), V (hS > 0)
            V = ds["V"].isel(**win).values * mS
            R = ds["RHO"].isel(**win).values
            E = ds["ETA"].isel(time=int(n_rec), j=slice(j0, j1), i=slice(i0, i1)).values
            _progress(f"face_series {name}", r, nt, on=progress)
            rp = np.where(wetC, R - rb3, 0.0)
            pt = budget.pressure_tilde(rp, E, dzC, rho0, g)                    # Pa, Nash 2005 eq. 3
            ep = np.where(wetC, budget.ape_exact(np.where(wetC, R, rb3), z3, rho_b, RC, g), 0.0)
            for q in qs:                                                       # KF12 eq. 12, J m^-3
                kind, j, i, _, _ = sel[q]
                a, b = j - j0, i - i0
                if kind == "y":
                    un = V[:, a, b]
                    p = 0.5 * (pt[:, a - 1, b] + pt[:, a, b])
                    epf = 0.5 * (ep[:, a - 1, b] + ep[:, a, b])
                    ut = 0.25 * (U[:, a - 1, b] + U[:, a - 1, b + 1] + U[:, a, b] + U[:, a, b + 1])
                else:
                    un = U[:, a, b]
                    p = 0.5 * (pt[:, a, b - 1] + pt[:, a, b])
                    epf = 0.5 * (ep[:, a, b - 1] + ep[:, a, b])
                    ut = 0.25 * (V[:, a, b - 1] + V[:, a + 1, b - 1] + V[:, a, b] + V[:, a + 1, b])
                dzq = dz[:nk, q]
                Hq = dzq.sum()
                Hs = Hq if Hq > 0 else 1.0
                ub, pb = (un * dzq).sum() / Hs, (p * dzq).sum() / Hs
                Ut = (ut * dzq).sum() / Hs
                up[r, :nk, q], pp[r, :nk, q] = un - ub, p - pb                 # KF12: u', p'
                un_[r, :nk, q] = un
                utp, unp = ut - Ut, un - ub
                # KF12 eq. 16 (thesis flujo_advectivo_cara), J m^-3
                e1[r, :nk, q] = 0.5 * rho0 * (utp ** 2 + unp ** 2) + rho0 * (Ut * utp + ub * unp) + epf
                fbt[r, q] = Hq * ub * pb                                       # W m^-1
                abt[r, q] = ub * Hq * 0.5 * rho0 * (Ut ** 2 + ub ** 2)          # W m^-1
    return {"up": up, "pp": pp, "un": un_, "e1": e1, "fbt": fbt, "abt": abt, "dz": dz, "w": w,
            "group": np.array([x[4] for x in sel]), "t": ds["time"].values[time].astype(float)}


def section_series(fs: dict) -> dict:
    """{group: {F_bc, A_bc, F_bt, A_bt: (nt,) [W], + into the box}} from face_series."""
    Fbc = (fs["up"] * fs["pp"] * fs["dz"][None]).sum(1) * fs["w"][None]
    Abc = (fs["un"] * fs["e1"] * fs["dz"][None]).sum(1) * fs["w"][None]
    out = {}
    for n in np.unique(fs["group"]):
        s = fs["group"] == n
        out[str(n)] = {"F_bc": Fbc[:, s].sum(1), "A_bc": Abc[:, s].sum(1),
                       "F_bt": (fs["fbt"] * fs["w"][None])[:, s].sum(1),
                       "A_bt": (fs["abt"] * fs["w"][None])[:, s].sum(1)}
    return out


def band_flux_series(fs: dict, bands: dict) -> dict:
    """Band split [J, + OUT of the box] of F_bc, A_bc and F_bc + A_bc per group and for all the
    faces ('total'): band value = -(band share of the time mean, validate.band_flux) x n dt, so
    the bands add up to -dt sum over the records ('dt_sum_out_J'). The record must start at
    t = 0 (PD0: 0-T_fin; ruling: assert t[0] = 0)."""
    t = fs["t"]
    if t[0] != 0.0:
        raise ValueError(f"band_flux_series: t[0] = {t[0]} s, the record must start at 0")
    dt = float(t[1] - t[0])
    n = t.size
    ser = section_series(fs)
    out = {}
    names = [str(x) for x in np.unique(fs["group"])]
    for g in names + ["total"]:
        s = np.ones(fs["group"].size, bool) if g == "total" else fs["group"] == g
        F = va.band_flux(fs["up"][:, :, s], fs["pp"][:, :, s], fs["dz"][:, s], fs["w"][s], dt, bands)
        A = va.band_flux(fs["un"][:, :, s], fs["e1"][:, :, s], fs["dz"][:, s], fs["w"][s], dt, bands)
        Fo = {b: -v * n * dt for b, v in F.items()}
        Ao = {b: -v * n * dt for b, v in A.items()}
        sl = names if g == "total" else [g]
        dts = {k: float(-dt * sum(ser[x][k].sum() for x in sl)) for k in ("F_bc", "A_bc")}
        dts["FA_bc"] = dts["F_bc"] + dts["A_bc"]
        out[g] = {"F_bc_out_J": Fo, "A_bc_out_J": Ao,
                  "FA_bc_out_J": {b: Fo[b] + Ao[b] for b in Fo}, "dt_sum_out_J": dts}
    return out


# ---- box energy level by level (H4b(i), maps of Fig. 4, trapped fraction) -----------------------

def level_energy(ds, mask, rho_b, bands: dict, time=slice(None), progress=True) -> dict:
    """Time-mean baroclinic energy of the box per cell, by band, reading one level at a time.

    * KE_bc per band [J m^-2] = 1/2 rho0 [sum_k dz_k (P_b(u_k) + P_b(v_k)) - H (P_b(u_bar) +
      P_b(v_bar))], P_b the Parseval band power (band_power), u_k the cell-centre velocity (mean
      of the two faces, as kf_budget), u_bar = sum u_k dz_k / H (KF12 eq. 13-14 kinetic parts:
      sum dz (u - u_bar)^2 = sum dz u^2 - H u_bar^2).
    * exact APE [J m^-2] = sum_k dz_k mean_t E'_p(rho_b,k + rho'), KF12 eq. 12 (budget.ape_exact),
      for rho' unfiltered ('all') and high-passed at the upper edge of the 'low' band ('hp').
    * linear APE per band [J m^-2] = sum_k dz_k 1/2 g^2 P_b(rho') / (rho0 N_b,k^2) (comparison
      only; NaN where N_b^2 <= 0).

    Returns per-cell arrays (cells (jj, ii), area rA [m^2]) and box totals [J]: ke_bc_band_J,
    ke_bc_all_J, ke_bc_hp_J, ape_exact_all_J, ape_exact_hp_J, ape_linear_band_J, E_in_hp_J."""
    rho0, g = float(ds.attrs["rho0"]), float(ds.attrs.get("g", budget.G))
    m = np.asarray(mask, bool)
    jj, ii = np.nonzero(m)
    j0, j1, i0, i1 = jj.min(), jj.max() + 1, ii.min(), ii.max() + 1
    lj, li = jj - j0, ii - i0
    nc = jj.size
    drF, RC = ds["drF"].values, ds["RC"].values
    hC, hW, hS = ds["hFacC"].values, ds["hFacW"].values, ds["hFacS"].values
    rA = ds["rA"].values[jj, ii]
    rho_b = np.asarray(rho_b, float)
    n2 = budget.n2_background(rho_b, RC, rho0, g)
    names = [b for b, _ in sorted(bands.items(), key=lambda kv: kv[1][0])]
    low = names[0]
    if bands[low][0] != 0.0:
        raise ValueError("the lowest band must start at 0 (the unresolved band)")
    cut = bands[low][1]
    tidx = np.arange(ds.sizes["time"])[time]
    t = ds["time"].values[tidx].astype(float)
    if t[0] != 0.0:
        raise ValueError(f"level_energy: t[0] = {t[0]} s, the record must start at 0")
    dt = float(t[1] - t[0])
    nt = tidx.size
    acc = np.zeros((nc, len(names)))
    ub_n, vb_n = np.zeros((nt, nc)), np.zeros((nt, nc))
    Hc = np.zeros(nc)
    ape_all, ape_hp = np.zeros(nc), np.zeros(nc)
    apelin = np.zeros((nc, len(names)))
    nk = int(np.flatnonzero((hC[:, jj, ii] > 0).any(axis=1)).max()) + 1
    for k in range(nk):
        dz = drF[k] * hC[k, jj, ii]
        if not (dz > 0).any():
            continue
        wet = dz > 0
        u, v, r = np.zeros((nt, nc)), np.zeros((nt, nc)), np.zeros((nt, nc))
        mW = hW[k, j0:j1, i0:i1 + 1] > 0
        mS = hS[k, j0:j1 + 1, i0:i1] > 0
        for q, n in enumerate(tidx):
            U = ds["U"].isel(time=int(n), k=k, j=slice(j0, j1), i=slice(i0, i1 + 1)).values * mW
            V = ds["V"].isel(time=int(n), k=k, j=slice(j0, j1 + 1), i=slice(i0, i1)).values * mS
            R = ds["RHO"].isel(time=int(n), k=k, j=slice(j0, j1), i=slice(i0, i1)).values
            u[q] = 0.5 * (U[lj, li] + U[lj, li + 1])
            v[q] = 0.5 * (V[lj, li] + V[lj + 1, li])
            r[q] = np.where(wet, R[lj, li] - rho_b[k], 0.0)
            _progress(f"level_energy k={k}", q, nt, on=progress)
        Pu, Pv, Pr = band_power(u, dt, bands), band_power(v, dt, bands), band_power(r, dt, bands)
        for c, b in enumerate(names):
            acc[:, c] += dz * (Pu[b] + Pv[b])
            apelin[:, c] += (dz * 0.5 * g ** 2 * Pr[b] / (rho0 * n2[k])) if n2[k] > 0 else np.nan
        ub_n += u * dz
        vb_n += v * dz
        Hc += dz
        z = np.full(r.shape, RC[k])
        ape_all += dz * budget.ape_exact(rho_b[k] + r, z, rho_b, RC, g).mean(0)   # J m^-3 -> J m^-2
        ape_hp += dz * budget.ape_exact(rho_b[k] + high_pass(r, dt, cut), z, rho_b, RC, g).mean(0)
        del u, v, r
    Hs = np.where(Hc > 0, Hc, 1.0)
    Pub, Pvb = band_power(ub_n / Hs, dt, bands), band_power(vb_n / Hs, dt, bands)
    ke = np.stack([0.5 * rho0 * (acc[:, c] - Hc * (Pub[b] + Pvb[b])) for c, b in enumerate(names)], 1)
    tot_band = {b: float((ke[:, c] * rA).sum()) for c, b in enumerate(names)}
    ke_hp = float(sum(v for b, v in tot_band.items() if b != low))
    tot = {"ke_bc_band_J": tot_band, "ke_bc_all_J": float(sum(tot_band.values())),
           "ke_bc_hp_J": ke_hp, "ape_exact_all_J": float((ape_all * rA).sum()),
           "ape_exact_hp_J": float((ape_hp * rA).sum()),
           "ape_linear_band_J": {b: float((apelin[:, c] * rA).sum()) for c, b in enumerate(names)}}
    tot["E_in_hp_J"] = tot["ke_bc_hp_J"] + tot["ape_exact_hp_J"]
    tot["E_in_all_J"] = tot["ke_bc_all_J"] + tot["ape_exact_all_J"]
    return {"cells": (jj, ii), "area_m2": rA, "band_names": names, "ke_bc_band_J_m-2": ke,
            "ape_exact_all_J_m-2": ape_all, "ape_exact_hp_J_m-2": ape_hp,
            "ape_linear_band_J_m-2": apelin, "low_cut_hz": cut, "totals": tot, "n_records": nt}


# ---- surface series and wind work (H4c(i)) ------------------------------------------------------

def wind_points(ds, mask) -> dict:
    """U and V points of the box with the kf_budget wind weights [m^2]: rAw hFacW(k=1) and
    rAs hFacS(k=1) times 1 inside the box and 1/2 on a face of the box edge (thesis peso_v)."""
    m = np.asarray(mask, bool).astype(float)
    ny, nx = m.shape
    mp = np.zeros((ny + 2, nx + 2))
    mp[1:-1, 1:-1] = m
    hW0, hS0 = ds["hFacW"].values[0], ds["hFacS"].values[0]
    wU = ds["rAw"].values * hW0 * 0.5 * (mp[1:-1, :-2] + mp[1:-1, 1:-1])   # faces (j, i): cells i-1, i
    wV = ds["rAs"].values * hS0 * 0.5 * (mp[:-2, 1:-1] + mp[1:-1, 1:-1])   # faces (j, i): cells j-1, j
    # east/north edge faces of the last column/row (index nx / ny) do not exist in a box that
    # does not touch the domain edge (budget.box_faces refuses such a box)
    ju, iu = np.nonzero(wU > 0)
    jv, iv = np.nonzero(wV > 0)
    return {"u_ji": (ju, iu), "wU": wU[ju, iu], "v_ji": (jv, iv), "wV": wV[jv, iv]}


def surface_series(ds, pts: dict, time=slice(None), progress=True):
    """(u (nt, nU), v (nt, nV)) [m s^-1] of the surface level (k = 1) at the U points
    pts['u_ji'] and V points pts['v_ji'] (masked by hFacW/S > 0 as kf_budget), one read of the
    bounding rectangle per record."""
    tidx = np.arange(ds.sizes["time"])[time]
    out = []
    for var, key, hf in (("U", "u_ji", "hFacW"), ("V", "v_ji", "hFacS")):
        jj, ii = (np.asarray(a, int) for a in pts[key])
        a = np.zeros((tidx.size, jj.size))
        if jj.size:
            j0, j1, i0, i1 = jj.min(), jj.max() + 1, ii.min(), ii.max() + 1
            wetf = ds[hf].values[0, jj, ii] > 0
            for q, n in enumerate(tidx):
                blk = ds[var].isel(time=int(n), k=0, j=slice(j0, j1), i=slice(i0, i1)).values
                a[q] = blk[jj - j0, ii - i0] * wetf
                _progress(f"surface {var}", q, tidx.size, on=progress)
        out.append(a)
    return out[0], out[1]


def wind_work(t, tau, pts: dict, u, v, band=None) -> float:
    """W = int sum (tau_x u w_U + tau_y v w_V) dt [J] (trapezoid, budget.cumulative), with
    tau(t) -> (taux at U points, tauy at V points) (ny, nx) [N m^-2] (load.forcing_stress) and,
    with ``band`` = (lo, hi) [Hz], u and v band-passed per point first (H4c(i): W_NI)."""
    t = np.asarray(t, float)
    if band is not None:
        dt = float(t[1] - t[0])
        u, v = band_pass(u, dt, *band), band_pass(v, dt, *band)
    ju, iu = pts["u_ji"]
    jv, iv = pts["v_ji"]
    p = np.zeros(t.size)
    for q, tt in enumerate(t):
        tx, ty = tau(float(tt))
        p[q] = (np.asarray(tx)[ju, iu] * u[q] * pts["wU"]).sum() + \
            (np.asarray(ty)[jv, iv] * v[q] * pts["wV"]).sum()
    return float(budget.cumulative(t, p)[-1])


# ---- stations -----------------------------------------------------------------------------------

def station_band_energy(u, v, dt: float, bands: dict) -> dict:
    """Surface horizontal KE spectrum P = P_u + P_v (spectra.psd: one Hann window, mean removed,
    window-corrected; [m^2 s^-2 Hz^-1]) and its band sums P df [m^2 s^-2] (half-open bands, the
    highest closed at Nyquist). Bands narrower than 2 df are not refused here (R1 at 120 h holds
    a 1.7-bin near-inertial band; run_bands flags it)."""
    fr, Pu, df = spectra.psd(u, dt)
    _, Pv, _ = spectra.psd(v, dt)
    P = Pu + Pv
    sel = _band_of_bins(fr, bands)
    e = {b: float(P[s].sum() * df) for b, s in sel.items()}
    tot = float(P.sum() * df)
    return {"band_m2_s-2": e, "total_m2_s-2": tot,
            "fraction": {b: (x / tot if tot > 0 else None) for b, x in e.items()},
            "freq_hz": fr, "psd": P, "df_hz": df}


def trapped_fraction(e, dist_m, r1_m: float, names) -> dict:
    """PD0 3 (secondary, descriptive): fraction of the energy of each band at less than 2 R1 from
    the coast. e (ncells, nbands) [J]; dist_m (ncells,) distance to land [m]."""
    e = np.asarray(e, float)
    near = np.asarray(dist_m, float) < 2.0 * r1_m
    tot = e.sum(0)
    return {b: (float(e[near, c].sum() / tot[c]) if tot[c] != 0 else None) for c, b in enumerate(names)}


# ---- bands per run (PD0 section 3) --------------------------------------------------------------

FJORD_RUNS = ("F1", "F2", "F3", "F4", "F5", "F1c", "F1g", "V3", "V4")
BAY_RUNS = ("R1", "R1k", "R2", "R3")


def run_bands(name: str, th: dict) -> dict:
    """Bands [Hz] of a run (PD0 section 3 via theory_numbers.spectral_bands; ruling: T_fin = 280 h
    for the fjord and R2/R3, the thesis-bay bands with thesis f for R1/R1k).

    * fjord runs: theory.json fjord_78N_280h;
    * R2/R3: f at 78 N, 280 h, 900 s records, seiche edge 0.5 x the lowest Merian frequency of the
      thesis bay (PD0 3: 0.149 cph);
    * R1k: theory.json bay_thesis_f_240h;
    * R1 and R1k over 0-120 h ('R1k_0_120h'): the same formula at 120 h; its near-inertial band
      holds 1.7 bins and is NOT resolvable (PD0 3), which ``resolvable`` records."""
    from exp4_arctic_bayiw.analysis import theory_numbers as tn
    b = th["bands"]
    if name in FJORD_RUNS:
        src = b["fjord_78N_280h"]
    elif name in ("R2", "R3"):
        src = tn.spectral_bands(b["fjord_78N_280h"]["f_cph"], 280.0, 900.0,
                                b["bay_thesis_f_240h"]["merian_lowest_cph"])
    elif name == "R1k":
        src = b["bay_thesis_f_240h"]
    elif name in ("R1", "R1k_0_120h"):
        ref = b["bay_thesis_f_240h"]
        df = 1.0 / 120.0
        f = ref["f_cph"]
        edges = [0.0, 2 * df, ref["ni_factors"][0] * f, ref["ni_factors"][1] * f,
                 ref["seiche_factor"] * ref["merian_lowest_cph"], ref["nyquist_cph"]]
        nm = ["low", "kelvin_subinertial", "near_inertial", "superinertial", "seiche"]
        src = {"record_h": 120.0, "dt_out_s": ref["dt_out_s"], "df_cph": df, "f_cph": f,
               "bands": [{"name": n, "lo_cph": lo, "hi_cph": hi, "n_bins": (hi - lo) / df}
                         for n, lo, hi in zip(nm, edges[:-1], edges[1:])]}
    else:
        raise KeyError(name)
    bands = {x["name"]: (x["lo_cph"] / 3600.0, x["hi_cph"] / 3600.0) for x in src["bands"]}
    nb = {x["name"]: (x["hi_cph"] - x["lo_cph"]) / src["df_cph"] for x in src["bands"]}
    return {"bands_hz": bands, "n_bins": nb,
            "resolvable": {k: (k == "low" or v >= 2.0) for k, v in nb.items()},
            "record_h": float(src["record_h"]), "dt_out_s": float(src["dt_out_s"]),
            "f_cph": float(src["f_cph"])}


# ---- local store of the box: one record-major LaCie pass per run ---------------------------------

from xarray.backends import BackendArray as _BackendArray  # noqa: E402


class _StoreArray(_BackendArray):
    """(time, [k,] j, i) view of a stored field on the cropped grid: the stored cells hold their
    values, every other cell (and every level >= nk) is 0. Outer indexing as load._Records."""

    def __init__(self, mm, jl, il, nyc: int, nxc: int, nz: int | None):
        self.mm, self.jl, self.il, self.nz = mm, jl, il, nz
        self.nk = mm.shape[1] if nz is not None else 0
        self.shape = (mm.shape[0],) + ((nz,) if nz is not None else ()) + (nyc, nxc)
        self.dtype = np.dtype("float64")

    def __getitem__(self, key):
        from xarray.core import indexing
        return indexing.explicit_indexing_adapter(key, self.shape, indexing.IndexingSupport.OUTER,
                                                  self._get)

    def _plane(self, n, k=None):
        a = np.zeros(self.shape[-2:])
        if k is None:
            a[self.jl, self.il] = self.mm[n]
        elif k < self.nk:
            a[self.jl, self.il] = self.mm[n, k]
        return a

    def _get(self, key):
        tidx = np.arange(self.shape[0])[key[0]]
        rest = list(key[1:])
        out = []
        for n in np.atleast_1d(tidx):
            if self.nz is not None:
                kk = np.arange(self.nz)[rest[0]]
                a = (self._plane(int(n), int(kk)) if np.ndim(kk) == 0
                     else np.stack([self._plane(int(n), int(k)) for k in kk]))
                r2 = rest[1:]
            else:
                a, r2 = self._plane(int(n)), rest
            lead = a.ndim - 2
            jk, ik = r2
            a = a[(slice(None),) * lead + (jk if isinstance(jk, (slice, int, np.integer)) else slice(None),)]
            if not isinstance(jk, (slice, int, np.integer)):
                a = np.take(a, jk, axis=lead)
            lead_i = a.ndim - 1
            a = a[(slice(None),) * lead_i + (ik if isinstance(ik, (slice, int, np.integer)) else slice(None),)]
            if not isinstance(ik, (slice, int, np.integer)):
                a = np.take(a, ik, axis=lead_i)
            out.append(np.array(a, dtype=np.float64))
        return out[0] if np.ndim(tidx) == 0 else np.stack(out)


def build_store(ds, mask, store_dir: Path, halo: int = 2, progress: bool = True) -> dict:
    """Copy U, V, RHO (the levels wet anywhere in the set) and ETA of the cells within ``halo``
    cells of the box (Chebyshev distance) to float64 memory maps in ``store_dir``, reading each
    record of the run once (one rectangle per field and record). The kernels then run on
    ``store_dataset`` with identical results (tests) and without re-reading LaCie."""
    import json
    from scipy import ndimage
    store_dir = Path(store_dir)
    meta_f = store_dir / "store.json"
    if meta_f.exists():
        m = json.loads(meta_f.read_text())
        if m.get("complete"):
            return m
    store_dir.mkdir(parents=True, exist_ok=True)
    S = ndimage.binary_dilation(np.asarray(mask, bool), structure=np.ones((3, 3), bool), iterations=halo)
    jS, iS = np.nonzero(S)
    J0, J1, I0, I1 = int(jS.min()), int(jS.max()) + 1, int(iS.min()), int(iS.max()) + 1
    hC = ds["hFacC"].values
    nk = int(np.flatnonzero((hC[:, jS, iS] > 0).any(axis=1)).max()) + 1
    nt = ds.sizes["time"]
    jl, il = jS - J0, iS - I0
    meta = {"crop": [J0, J1, I0, I1], "nk": nk, "nt": nt, "nS": int(jS.size), "nz": int(hC.shape[0]),
            "jl": jl.tolist(), "il": il.tolist(), "dir": str(store_dir), "complete": False}
    mms = {v: np.lib.format.open_memmap(store_dir / f"{v}.npy", mode="w+", dtype="float64",
                                        shape=(nt, nk, jS.size)) for v in ("U", "V", "RHO")}
    mms["ETA"] = np.lib.format.open_memmap(store_dir / "ETA.npy", mode="w+", dtype="float64",
                                           shape=(nt, jS.size))
    win = dict(k=slice(0, nk), j=slice(J0, J1), i=slice(I0, I1))
    for n in range(nt):
        for v in ("U", "V", "RHO"):
            mms[v][n] = ds[v].isel(time=n, **win).values[:, jl, il]
        mms["ETA"][n] = ds["ETA"].isel(time=n, j=slice(J0, J1), i=slice(I0, I1)).values[jl, il]
        _progress("store", n, nt, on=progress)
    for m_ in mms.values():
        m_.flush()
    meta["complete"] = True
    meta_f.write_text(json.dumps(meta))
    return meta


def store_dataset(ds, st: dict, mask, groups: dict | None = None):
    """(local dataset, local mask, local groups, crop) on the cropped grid of a store: grid fields
    cropped, U/V/RHO/ETA from the store (0 outside the stored cells and below the wet levels),
    attrs kept (rho0, g, eos, stages)."""
    import xarray as xr
    from xarray.core import indexing
    J0, J1, I0, I1 = st["crop"]
    d = Path(st["dir"])
    jl, il = np.asarray(st["jl"]), np.asarray(st["il"])
    nyc, nxc = J1 - J0, I1 - I0
    keep = {}
    for v in ds.data_vars:
        if v in ("U", "V", "RHO", "ETA", "T", "TAUX", "TAUY"):
            continue
        x = ds[v]
        if "j" in x.dims and "i" in x.dims:
            x = x.isel(j=slice(J0, J1), i=slice(I0, I1))
        keep[v] = x.variable
    for v in ("U", "V", "RHO"):
        mm = np.load(d / f"{v}.npy", mmap_mode="r")
        keep[v] = xr.Variable(("time", "k", "j", "i"),
                              indexing.LazilyIndexedArray(_StoreArray(mm, jl, il, nyc, nxc, st["nz"])))
    mm = np.load(d / "ETA.npy", mmap_mode="r")
    keep["ETA"] = xr.Variable(("time", "j", "i"), indexing.LazilyIndexedArray(_StoreArray(mm, jl, il, nyc, nxc, None)))
    out = xr.Dataset(keep, coords={"time": ds["time"].values})
    out.attrs.update(ds.attrs)
    ml = np.asarray(mask, bool)[J0:J1, I0:I1]
    gl = None if groups is None else {n: {k: np.asarray(a)[J0:J1, I0:I1] for k, a in g.items()}
                                      for n, g in groups.items()}
    if ml.sum() != np.asarray(mask, bool).sum():
        raise ValueError("the box is not inside the stored rectangle")
    return out, ml, gl, (J0, J1, I0, I1)


def crop_tau(tau, crop):
    """tau(t) on the cropped grid."""
    if tau is None:
        return None
    J0, J1, I0, I1 = crop

    def f(t):
        tx, ty = tau(t)
        return np.asarray(tx)[J0:J1, I0:I1], np.asarray(ty)[J0:J1, I0:I1]
    return f


def remove_store(st: dict) -> None:
    import shutil
    shutil.rmtree(st["dir"], ignore_errors=True)


# ---- run driver ---------------------------------------------------------------------------------

ICE_RUNS = ("F2", "F4")
LEVEL_RUNS = ("F1", "F2", "F3")          # H4b(i) needs F2, F3; F1 for the trapped fraction/maps
SUMMER = ("F1", "F5", "F1c", "F1g", "V3", "V4")
WINTER = ("F2", "F3", "F4")
J_BAY, I_BAY = (324, 348), (274, 286)    # thesis_contrast.J_BAY/I_BAY (balance_kang_fringer)


def bay_box(ds, rows=J_BAY, cols=I_BAY):
    """The thesis bay box (rows x cols, closed at the mouth, as the Task 5 contrast) and its open
    faces as one group 'mouth' (faces with water on both sides). ValueError if the box holds no
    wet cell (R3: the bay is land)."""
    wet = ds["hFacC"].values[0] > 0
    ny, nx = wet.shape
    m = np.zeros((ny, nx), bool)
    m[rows[0]:rows[1], cols[0]:cols[1]] = True
    m &= wet
    if not m.any():
        raise ValueError("no wet cell in the bay box (no bay in this run)")
    f = budget.box_faces(m, wet=wet)
    g = {"mouth": {"x": np.zeros((ny, nx), bool), "y": np.zeros((ny, nx), bool)}}
    for kind in ("x", "y"):
        for j, i, _ in f[kind]:
            g["mouth"][kind][j, i] = True
    return m, g


def tau_from_diagnostics(ds):
    """tau(t) -> (oceTAUX, oceTAUY) of the record at time t [N m^-2]: the stress the ocean
    receives in an ice run (under landfast ice the wind file is not what the ocean feels;
    notes/pd1_seaice.md). KeyError if t is not a record time."""
    t = ds["time"].values.astype(float)
    idx = {float(x): n for n, x in enumerate(t)}

    def tau(tt):
        n = idx[float(tt)]
        return ds["TAUX"].isel(time=n).values, ds["TAUY"].isel(time=n).values
    return tau


def window(r: dict, t_end_s: float) -> dict:
    """kf_budget series cut to 0..t_end (t_end must be a record time): R1k over 0-120 h (H4d)."""
    t = np.asarray(r["t"], float)
    hit = np.flatnonzero(np.isclose(t, t_end_s, rtol=0, atol=1e-6))
    if not hit.size:
        raise ValueError(f"{t_end_s} s is not a record time")
    n = int(hit[0]) + 1
    return {k: (v[:n] if isinstance(v, np.ndarray) and v.ndim == 1 and v.size == t.size else v)
            for k, v in r.items()}


def budget_summary(r: dict, sections, every: int | None = None) -> dict:
    """Integrals over the record [J] of the kf_budget series (trapezoid, budget.cumulative; + INTO
    the box for F and A), the closure of both budgets (budget.closure: D the residual), the H4d
    numerator F_bc_pressure_out (Enmienda 2) and the wind work P = P_bt + P_bc. ``every``: also
    the series every n records for Fig. 3."""
    t = np.asarray(r["t"], float)
    if t[0] != 0.0:
        raise ValueError(f"budget_summary: t[0] = {t[0]} s, the record must start at 0")
    cum = lambda x: float(budget.cumulative(t, x)[-1])                     # noqa: E731
    c = budget.closure(r)
    out = {"window_h": [0.0, float(t[-1] / 3600.0)], "n_records": int(t.size),
           "closure_J": {comp: {k: float(v[-1]) for k, v in c[comp].items()} for comp in c},
           "totals_in_J": {k: cum(r[k]) for k in ("F_bt", "F_bc", "A_bt", "A_bc")},
           "per_section_in_J": {s: {k: cum(r[f"{k}_{s}"]) for k in ("F_bt", "F_bc", "A_bt", "A_bc")}
                                for s in sections},
           "F_bc_pressure_out_J": float(budget.f_bc_pressure_out(r)[-1]),
           "F_bc_pressure_out_per_section_J": {s: float(budget.f_bc_pressure_out(r, s)[-1])
                                               for s in sections},
           "P_bt_J": cum(r["P_bt"]), "P_bc_J": cum(r["P_bc"]),
           "P_J": cum(r["P_bt"] + r["P_bc"]), "C_J": cum(r["C"]),
           "energy_J": {k: {"start": float(r[k][0]), "end": float(r[k][-1]), "mean": float(np.mean(r[k]))}
                        for k in ("ke_bt", "pe_bt", "ke_bc", "ape", "ape_linear")}}
    tin = out["totals_in_J"]
    tin["FA_bc"], tin["FA_bt"] = tin["F_bc"] + tin["A_bc"], tin["F_bt"] + tin["A_bt"]
    for s in sections:
        p = out["per_section_in_J"][s]
        p["FA_bc"], p["FA_bt"] = p["F_bc"] + p["A_bc"], p["F_bt"] + p["A_bt"]
    bc = out["closure_J"]["bc"]
    gross = abs(bc["P"]) + abs(bc["C"]) + abs(bc["F"])
    out["D_bc_over_gross_sources"] = bc["D"] / gross if gross else None
    if every:
        sl = slice(None, None, every)
        cb = c["bc"]
        out["series"] = {"t_h": (t[sl] / 3600.0).tolist(),
                         "FA_bc_in_W": {s: (r[f"F_bc_{s}"] + r[f"A_bc_{s}"])[sl].tolist() for s in sections},
                         "F_bc_in_W": {s: r[f"F_bc_{s}"][sl].tolist() for s in sections},
                         "bc_cumulative_J": {k: cb[k][sl].tolist() for k in ("P", "C", "F", "dE", "D")}}
    return out


# ---- caches (keyed on the source of the kernels; analysis_cache on LaCie) -----------------------

def _key(funcs=(), files=()) -> str:
    """sha256 (12 hex) of the source of the kernel functions and of whole library files: a change
    in the code that makes a cache misses it (thesis_contrast.cache_file convention)."""
    import hashlib
    import inspect
    h = hashlib.sha256()
    for f in funcs:
        h.update(inspect.getsource(f).encode())
    for p in files:
        h.update(Path(p).read_bytes())
    return h.hexdigest()[:12]


_LIB = Path(__file__).resolve().parent
KEYS = {
    "kf": lambda: _key(files=(_LIB / "budget.py", _LIB / "load.py")),
    "faces": lambda: _key((face_series,), (_LIB / "budget.py", _LIB / "load.py")),
    "surface": lambda: _key((wind_points, surface_series, wind_work, band_pass, _keep, _stations_of,
                             tau_from_diagnostics),
                            (_LIB / "load.py", _LIB / "validate.py", _LIB / "budget.py")),
    "energy": lambda: _key((level_energy, band_power, high_pass, _keep, _rfft_weights, _band_of_bins),
                           (_LIB / "budget.py", _LIB / "load.py")),
}


def cache_file(cache: Path, kind: str, run: str) -> Path:
    """Keyed cache file; the key is the one frozen when the module was imported (FROZEN_KEYS)."""
    key = FROZEN_KEYS[kind] if FROZEN_KEYS else KEYS[kind]()
    return Path(cache) / f"t7_{kind}_{run}_{key}.npz"


def _save(f: Path, d: dict) -> None:
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_name(f.stem + ".tmp.npz")
    np.savez(tmp, **d)
    tmp.replace(f)


def _load(f: Path) -> dict:
    return dict(np.load(f, allow_pickle=False))


# ---- geometry of a run ----------------------------------------------------------------------------

def _box(name: str, ds, runs_root: Path):
    """(mask, groups) of the run: the inner fjord box closed at Akselsundet and Mariasundet
    (budget.fjord_box; V4: the 200 m box mapped to 400 m, validate.coarse_box, Task 6), or the
    thesis bay box closed at the mouth (R1-R3). None for R3 (no bay)."""
    from exp4_arctic_bayiw.analysis import load
    if name == "V4":
        f1 = load.run(Path(runs_root) / "F1")
        m1, _ = budget.fjord_box(f1)
        return va.coarse_box(m1, ds["hFacC"].values[0] > 0, ds["XC"].values[0], ds["YC"].values[:, 0],
                             va._sections())
    if name in FJORD_RUNS:
        return budget.fjord_box(ds)
    try:
        return bay_box(ds)
    except ValueError:
        return None


def _stations_of(name: str, ds, mask, th: dict) -> dict:
    """{station: (j, i)} at T points. Fjord runs: the V2a shelf station and the V2b south-coast
    stations every 5 km (validate.shelf_station / coastal_stations with the Task 6 rules, fixed
    before the campaign was analysed); bay runs: the thesis stations (v1_thesis.STATIONS,
    picos_seiche_inercial.ESTACIONES) that are wet."""
    wet = ds["hFacC"].values[0] > 0
    if name in FJORD_RUNS:
        p = va._case_params(ds)
        rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, ds["hFacC"].values > 0)
        k = va.pycnocline_level(rho_b, ds["RC"].values, float(ds.attrs["rho0"]), float(ds.attrs["g"]))
        dx = float(np.median(np.diff(ds["XC"].values[0])))
        rules = va.station_rules(th, dx)
        deep = (ds["hFacC"].values[min(k + 1, ds.sizes["k"] - 1)] > 0) & \
            (ds["Depth"].values >= rules["min_depth_m"])
        out = {"shelf": va.shelf_station(wet, mask, p["ob_sides"], int(p["n_sponge"]))}
        st = va.coastal_stations(wet, deep, mask, ds["XC"].values[0], ds["YC"].values[:, 0],
                                 th["fjord_geometry"]["axis_angle_deg"], 5000.0,
                                 max_coast_cells=rules["max_coast_cells"])
        for q, c in enumerate(st["cells"]):
            out[f"coast{q:02d}"] = c
        return out
    from exp4_arctic_bayiw.analysis.v1_thesis import STATIONS
    return {n: c for n, c in STATIONS.items() if wet[c]}


def _station_points(st: dict) -> dict:
    """U and V points around each T-point station (west/east, south/north faces)."""
    ju, iu, jv, iv = [], [], [], []
    for j, i in st.values():
        ju += [j, j]; iu += [i, i + 1]
        jv += [j, j + 1]; iv += [i, i]
    return {"u_ji": (np.array(ju), np.array(iu)), "v_ji": (np.array(jv), np.array(iv))}


# ---- compute one run --------------------------------------------------------------------------

def compute(name: str, runs_root: Path, cache: Path, th: dict, recompute: bool = False,
            progress: bool = True, run_dir: Path | None = None, label: str | None = None,
            store_root: Path | None = None, keep_store: bool = False) -> dict:
    """Compute (or read from the keyed cache) the kf_budget series, the face series, the surface
    series with the wind work by band, and (LEVEL_RUNS) the box energy by band of one run.
    The box fields are read from LaCie once into a local store (build_store; ``store_root``,
    default <system temp>/exp4_t7_store, on the local disk), removed afterwards unless
    ``keep_store``.
    Returns {kind: cache file}. ``run_dir``/``label``: another directory classified as ``name``
    (smoke tests on short runs), cached under ``label``."""
    from exp4_arctic_bayiw.analysis import load
    ds = load.run(Path(run_dir) if run_dir else Path(runs_root) / name)
    tag = label or name
    t = ds["time"].values.astype(float)
    if t[0] != 0.0:
        raise ValueError(f"{name}: t[0] = {t[0]} s; the record must start at 0 (0-T_fin)")
    box = _box(name, ds, runs_root)
    files = {}
    rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, ds["hFacC"].values > 0)
    ice = name in ICE_RUNS
    tau = tau_from_diagnostics(ds) if ice else (lambda tt: load.forcing_stress(ds, tt))
    tau_src = ("oceTAUX/oceTAUY diagnostics (ice run: the stress the ocean receives)" if ice
               else "forcing files (load.forcing_stress)")
    store = None
    if box is not None:
        mask, groups = box
        need = [k for k in ("kf", "faces") + (("energy",) if name in LEVEL_RUNS else ())
                if recompute or not cache_file(cache, k, tag).exists()]
        if need:
            print(f"[partition] {tag}: one LaCie pass into the local store ({need})", file=sys.stderr, flush=True)
            import tempfile
            root = Path(store_root) if store_root else Path(tempfile.gettempdir()) / "exp4_t7_store"
            store = build_store(ds, mask, root / tag, progress=progress)
            dl, ml, gl, crop = store_dataset(ds, store, mask, groups)
            taul = crop_tau(tau, crop)
        f = cache_file(cache, "kf", tag)
        if "kf" in need:
            print(f"[partition] {tag}: kf_budget", file=sys.stderr, flush=True)
            r = budget.kf_budget(dl, ml, tau=taul, rho_b=rho_b, face_groups=gl)
            d = {k: v for k, v in r.items() if isinstance(v, np.ndarray)}
            d.update(mask=mask, sections=np.array(sorted(groups)), tau_source=np.array(tau_src))
            _save(f, d)
        files["kf"] = f
        f = cache_file(cache, "faces", tag)
        if "faces" in need:
            print(f"[partition] {tag}: face_series", file=sys.stderr, flush=True)
            _save(f, face_series(dl, ml, gl, rho_b=rho_b, progress=progress))
        files["faces"] = f
    f = cache_file(cache, "surface", tag)
    if recompute or not f.exists():
        print(f"[partition] {name}: surface series and wind work", file=sys.stderr, flush=True)
        st = _stations_of(name, ds, box[0] if box is not None else None, th)
        sp = _station_points(st)
        d = {"station_names": np.array(list(st)), "station_ji": np.array(list(st.values()))}
        bands = run_bands(name, th)["bands_hz"]
        if box is not None:
            wp = wind_points(ds, box[0])
            nu, nv = wp["u_ji"][0].size, wp["v_ji"][0].size
            pts = {"u_ji": tuple(np.concatenate((a, b)) for a, b in zip(wp["u_ji"], sp["u_ji"])),
                   "v_ji": tuple(np.concatenate((a, b)) for a, b in zip(wp["v_ji"], sp["v_ji"]))}
        else:
            nu = nv = 0
            pts = sp
        u, v = surface_series(ds, pts, progress=progress)
        d["station_u"] = 0.5 * (u[:, nu::2] + u[:, nu + 1::2])            # T point, k = 1
        d["station_v"] = 0.5 * (v[:, nv::2] + v[:, nv + 1::2])
        d["t"] = t
        if box is not None:
            uw, vw = u[:, :nu], v[:, :nv]
            d["wind_all_J"] = wind_work(t, tau, wp, uw, vw)
            d["wind_band_names"] = np.array(list(bands))
            d["wind_band_J"] = np.array([wind_work(t, tau, wp, uw, vw, band=bands[b]) for b in bands])
            d["tau_source"] = np.array(tau_src)
        _save(f, d)
    files["surface"] = f
    if name in LEVEL_RUNS and box is not None:
        f = cache_file(cache, "energy", tag)
        if "energy" in need:
            print(f"[partition] {tag}: level_energy", file=sys.stderr, flush=True)
            e = level_energy(dl, ml, rho_b, run_bands(name, th)["bands_hz"], progress=progress)
            J0, _, I0, _ = crop
            e["cells"] = (e["cells"][0] + J0, e["cells"][1] + I0)          # back to model indices
            import json
            _save(f, {"jj": e["cells"][0], "ii": e["cells"][1], "area_m2": e["area_m2"],
                      "band_names": np.array(e["band_names"]), "ke_bc_band": e["ke_bc_band_J_m-2"],
                      "ape_exact_all": e["ape_exact_all_J_m-2"], "ape_exact_hp": e["ape_exact_hp_J_m-2"],
                      "ape_linear_band": e["ape_linear_band_J_m-2"], "low_cut_hz": e["low_cut_hz"],
                      "totals_json": np.array(json.dumps(e["totals"])),
                      "dist_land_m": (va._dist_land(ds["hFacC"].values[0] > 0)[e["cells"]] - 0.5)
                      * float(np.median(np.diff(ds["XC"].values[0])))})
        files["energy"] = f
    if store is not None and not keep_store:
        remove_store(store)
    return files


# ---- assembly: partition.json and budget.json -----------------------------------------------------

POST_HOC = "POST HOC descriptive (task-7 rulings, R24): not a verdict"


def _cph(b: dict) -> dict:
    return {k: [lo * 3600.0, hi * 3600.0] for k, (lo, hi) in b.items()}


def _al_cover(nyq_hz: float) -> dict:
    lo, hi = AL_BAND_CPH[0] / 3600.0, AL_BAND_CPH[1] / 3600.0
    return {"below": (0.0, lo), "al": (lo, hi), "above": (hi, nyq_hz)}


def _slice_faces(fs: dict, n: int) -> dict:
    return {k: (v[:n] if k in ("up", "pp", "un", "e1", "fbt", "abt", "t") else v) for k, v in fs.items()}


def _straits_block(fs: dict, rb: dict, kf: dict, sections) -> dict:
    """Band split of F_bc, A_bc, F_bc + A_bc per strait and total [J, + out]; the post hoc
    Arneborg & Liljebladh band; and the check that the face series reproduce kf_budget."""
    bands = rb["bands_hz"]
    out = band_flux_series(fs, bands)
    dt = float(fs["t"][1] - fs["t"][0])
    al = band_flux_series(fs, _al_cover(0.5 / dt))
    ser = section_series(fs)
    worst = 0.0
    for s in sections:
        for k in ("F_bc", "A_bc", "F_bt", "A_bt"):
            ref = kf[f"{k}_{s}"][:fs["t"].size]
            sc = max(np.abs(ref).max(), 1e-300)
            worst = max(worst, float(np.abs(ser[s][k] - ref).max() / sc))
    return {"bands": out,
            "check_face_series_vs_kf_budget_max_rel": worst,
            "arneborg_liljebladh_band_post_hoc": {
                "label": POST_HOC, "band_cph": list(AL_BAND_CPH),
                "reference": "Arneborg & Liljebladh (2001, JPO 31:2549) Gullmar fjord: net inward "
                             "baroclinic energy flux at the mouth, band-passed 0.01-0.07 cph",
                **{g: {"FA_bc_in_J": -v["FA_bc_out_J"]["al"], "F_bc_in_J": -v["F_bc_out_J"]["al"]}
                   for g, v in al.items()}}}


def _station_block(sf: dict, rb: dict, n: int | None = None) -> dict:
    names = [str(x) for x in sf["station_names"]]
    t = sf["t"][:n] if n else sf["t"]
    dt = float(t[1] - t[0])
    out, coast = {}, []
    for q, nm in enumerate(names):
        u = sf["station_u"][:t.size, q]
        v = sf["station_v"][:t.size, q]
        s = station_band_energy(u, v, dt, rb["bands_hz"])
        out[nm] = {"ji": [int(x) for x in sf["station_ji"][q]], "band_m2_s-2": s["band_m2_s-2"],
                   "fraction": s["fraction"], "total_m2_s-2": s["total_m2_s-2"]}
        if nm.startswith("coast"):
            coast.append(s)
        out[nm]["_psd"] = s
    spec = {}
    for nm in ("shelf", "Seiche", "Seiche2", "Omega"):
        if nm in out:
            s = out[nm]["_psd"]
            spec[nm] = {"freq_cph": (s["freq_hz"] * 3600.0).tolist(), "psd_m2_s-2_cph-1": (s["psd"] / 3600.0).tolist()}
    if coast:
        P = np.mean([c["psd"] for c in coast], 0)
        fr, df = coast[0]["freq_hz"], coast[0]["df_hz"]
        e = {b: float(P[s_].sum() * df) for b, s_ in _band_of_bins(fr, rb["bands_hz"]).items()}
        tot = float(P.sum() * df)
        out["coast_south_mean"] = {"n_stations": len(coast), "band_m2_s-2": e, "total_m2_s-2": tot,
                                   "fraction": {b: x / tot for b, x in e.items()}}
        spec["coast_south_mean"] = {"freq_cph": (fr * 3600.0).tolist(), "psd_m2_s-2_cph-1": (P / 3600.0).tolist()}
    for nm in names:
        out[nm].pop("_psd")
    return {"stations": out, "spectra": spec,
            "method": "surface (k = 1) u, v at the T point; one Hann window over the record "
                      "(spectra.psd); fjord: V2a shelf station and V2b south-coast stations (Task 6 "
                      "rules); bay: thesis stations (v1_thesis.STATIONS) that are wet"}


def _energy_block(ef: dict, name: str, th: dict) -> dict:
    import json
    tot = json.loads(str(ef["totals_json"]))
    names = [str(x) for x in ef["band_names"]]
    e = ef["ke_bc_band"] * ef["area_m2"][:, None]
    season = "vmf_summer" if name in SUMMER else "vmf_winter"
    tf = {}
    for tag, case in (("R1_Hmax", season), ("R1_Hmean", season + "_Hmean")):
        r1 = th["cases"][case]["rossby_radius_m"][0]
        tf[tag] = {"R1_m": r1, "fraction": trapped_fraction(e, ef["dist_land_m"], r1, names)}
    return {"totals": tot, "trapped_fraction_ke_bc": tf,
            "method": "time-mean KE_bc (Parseval per band) and exact APE of the inner-fjord box, "
                      "fields high-passed at 2 df for E_in_hp (H4b(i)); trapped fraction: share of "
                      "the depth-integrated KE_bc of each band in cells closer than 2 R1 to land "
                      "(distance = (EDT - 1/2) dx); linear APE is a comparison only"}


def assemble(cache: Path, th: dict, runs=FJORD_RUNS + BAY_RUNS, out_dir: Path | None = None):
    """partition.json and budget.json from the caches of ``runs`` (missing ones are listed)."""
    import json
    out_dir = Path(out_dir or Path(__file__).resolve().parents[1] / "figures")
    part = {"meta": {"sign": "*_out_J: + out of the box (inner fjord, or bay for R1-R3); *_in_J: + into it",
                     "flux": "Kang & Fringer (2012) on the box faces (budget.kf_budget); F_bc pressure "
                             "work, A_bc advection; FA_bc = F_bc + A_bc",
                     "band_split": "exact cospectral split of the time mean (validate.band_flux); band "
                                   "sums = -dt sum over the records",
                     "window": "0-T_fin, every record from t = 0; R1k_0_120h: R1k records 0-120 h"},
            "runs": {}, "missing": []}
    bud = {"meta": {"sign": "F, A + INTO the box; F_bc_pressure_out + out (Enmienda 2)",
                    "closure": "budget.closure: D = P -/+ C - dE + F (bt/bc), F = pressure work + advection",
                    "window": part["meta"]["window"]}, "runs": {}, "missing": []}
    names = list(runs) + (["R1k_0_120h"] if "R1k" in runs else [])
    for name in names:
        src = "R1k" if name == "R1k_0_120h" else name
        f = {k: cache_file(cache, k, src) for k in KEYS}
        if not f["surface"].exists():
            part["missing"].append(name)
            bud["missing"].append(name)
            continue
        rb = run_bands(name, th)
        sf = _load(f["surface"])
        n = None
        if name == "R1k_0_120h":
            n = int(np.flatnonzero(np.isclose(sf["t"], 120 * 3600.0))[0]) + 1
        pr = {"bands_cph": _cph(rb["bands_hz"]), "resolvable": rb["resolvable"],
              "n_bins": rb["n_bins"], "record": {"n": int(n or sf["t"].size), "dt_s": float(sf["t"][1] - sf["t"][0]),
                                                 "t_last_h": float(sf["t"][(n or sf["t"].size) - 1] / 3600.0)}}
        pr.update(_station_block(sf, rb, n))
        if f["kf"].exists() and f["faces"].exists():
            kf = _load(f["kf"])
            sections = [str(x) for x in kf["sections"]]
            fs = _load(f["faces"])
            if n:
                fs = _slice_faces(fs, n)
                kf = window(kf, 120 * 3600.0)
            blk = _straits_block(fs, rb, kf, sections)
            pr["straits"] = blk["bands"]
            pr["check_face_series_vs_kf_budget_max_rel"] = blk["check_face_series_vs_kf_budget_max_rel"]
            pr["arneborg_liljebladh_band_post_hoc"] = blk["arneborg_liljebladh_band_post_hoc"]
            pr["tau_source"] = str(kf["tau_source"])
            if "wind_all_J" in sf and not n:
                pr["wind_work_J"] = {"all": float(sf["wind_all_J"]),
                                     "bands": {str(b): float(v) for b, v in zip(sf["wind_band_names"], sf["wind_band_J"])},
                                     "tau_source": str(sf["tau_source"]),
                                     "check_all_vs_kf_P_J": float(budget.cumulative(kf["t"], kf["P_bt"] + kf["P_bc"])[-1])}
            bs = budget_summary(kf, sections, every=3 if name in FJORD_RUNS else 4)
            bs["box_cells"] = int(np.asarray(kf["mask"]).sum())
            bs["tau_source"] = str(kf["tau_source"])
            bud["runs"][name] = bs
        else:
            pr["straits"] = None
            pr["note"] = "no box: the bay is land in R3 (no flux, no budget)" if name == "R3" else "box not computed"
        if f["energy"].exists() and not n:
            pr["energy"] = _energy_block(_load(f["energy"]), name, th)
        part["runs"][name] = pr
    (out_dir / "partition.json").write_text(json.dumps(part, indent=1, default=va.jsonable) + "\n")
    (out_dir / "budget.json").write_text(json.dumps(bud, indent=1, default=va.jsonable) + "\n")
    return part, bud


# keys frozen at import: a later edit of this file cannot rename a running process's caches
FROZEN_KEYS: dict = {}
FROZEN_KEYS.update({k: f() for k, f in KEYS.items()})


def main(argv=None) -> int:
    import argparse
    from pfcommon.paths import EXP4_RUNS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute", help="fill the caches of the runs (reads LaCie)")
    c.add_argument("runs", nargs="+")
    c.add_argument("--recompute", action="store_true")
    c.add_argument("--store-root", type=Path, default=None, help="local disk for the box store")
    a_ = sub.add_parser("assemble", help="partition.json and budget.json from the caches")
    a_.add_argument("--out-dir", type=Path, default=None)
    for p in (c, a_):
        p.add_argument("--runs-root", type=Path, default=EXP4_RUNS)
        p.add_argument("--cache", type=Path, default=EXP4_RUNS / "analysis_cache")
    a = ap.parse_args(argv)
    th = va.load_theory()
    if a.cmd == "compute":
        for r in a.runs:
            print(r, compute(r, a.runs_root, a.cache, th, recompute=a.recompute,
                             store_root=a.store_root), flush=True)
    else:
        part, bud = assemble(a.cache, th, out_dir=a.out_dir)
        print("missing:", part["missing"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
