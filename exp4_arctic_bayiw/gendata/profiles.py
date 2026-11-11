"""Representative T/S profiles of Van Mijenfjorden from the UNIS Hydrographic Database (Exp. 4, Task 2).

Pipeline: ``load_unis_hd`` (one dict per CTD cast) -> ``select`` (polygon, months, years) -> ``qc``
(range checks, density-inversion fraction) -> ``representative`` (median per model level, extension
below the last datum, minimal static-stability adjustment under the model EOS).

Conventions: ``z`` is in metres, positive up (z <= 0 in the water), as in theory/modes.py; profiles
are ordered from the surface down. ``p`` is sea pressure in dbar. ``T`` is the in-situ temperature
(deg C, ITS-90) and ``S`` the practical salinity, as stored in UNIS HD (TEMP, PSAL).

Sources:
- UNIS HD (Skogseth, Ellingsen et al., https://doi.org/10.21334/unis-hydrography, CC BY 4.0); the
  data summary warns that "calibration has not been quality checked in all the data, so use with
  caution (salinity values in particular)" (global attribute ``summary`` of every file).
- Depth from pressure and the TEOS-10 density used by ``qc``: gsw (IOC, SCOR and IAPSO 2010,
  TEOS-10 manual) through ``gsw.z_from_p``, ``gsw.SA_from_SP``, ``gsw.CT_from_t``, ``gsw.sigma0``.
- Model EOS (``EOS_EXP4``): MITgcm linear EOS, rho = rhoNil (sBeta (S - sRef) - tAlpha (T - tRef))
  + (rhoNil - rhoConst) (model/src/find_rho.F:104-107), with the thesis coefficients tAlpha = 2e-4
  K^-1 and sBeta = 7.4e-4 (thesis_copy/BayIW_Rectan/run_expand/data.stage1:15-16).
"""
from __future__ import annotations

import datetime as dt
import warnings
from pathlib import Path

import gsw
import numpy as np
import shapely
from scipy.ndimage import uniform_filter1d
from shapely.geometry import Polygon

EOS_EXP4 = {"tAlpha": 2.0e-4, "sBeta": 7.4e-4}   # K^-1, (g/kg)^-1 ; data.stage1:15-16
QC_DEFAULTS = {"s_range": (0.0, 35.5), "t_range": (-2.1, 10.0), "smooth_m": 3.0,
               "max_inversion_fraction": 0.10, "inversion_tol_kg_m3": 4.0e-3, "min_levels": 5}


# ----------------------------------------------------------------------------- reading
def _read_one(fn: Path) -> dict | None:
    import netCDF4

    with netCDF4.Dataset(fn) as d:
        d.set_auto_mask(False)
        p = np.asarray(d["PRES"][:], float).ravel()
        T = np.asarray(d["TEMP"][:], float).ravel()
        has_s = "PSAL" in d.variables          # 68 of the 1547 casts are temperature-only
        S = np.asarray(d["PSAL"][:], float).ravel() if has_s else np.full_like(T, np.nan)
        lat = float(np.asarray(d["LATITUDE"][:]).ravel()[0])
        lon = float(np.asarray(d["LONGITUDE"][:]).ravel()[0])
        tv = d["TIME"]
        t = netCDF4.num2date(np.asarray(tv[:]).ravel()[0], tv.units,
                             getattr(tv, "calendar", "standard"),
                             only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        attrs = {k: str(getattr(d, k)) for k in ("instrument_type", "QC_indicator", "cruise_ID")
                 if k in d.ncattrs()}
        for name, v in (("T", T), ("S", S)):
            if name == "S" and not has_s:
                continue
            fv = getattr(d[{"T": "TEMP", "S": "PSAL"}[name]], "_FillValue", None)
            if fv is not None and np.isfinite(fv):
                v[v == fv] = np.nan
    ok = np.isfinite(p) & np.isfinite(T) & (np.isfinite(S) | (not has_s))
    if not ok.any():
        return None
    p, T, S = p[ok], T[ok], S[ok]
    o = np.argsort(p)
    p, T, S = p[o], T[o], S[o]
    z = gsw.z_from_p(p, lat)                       # m, negative below the sea surface
    return {"lat": lat, "lon": lon, "time": dt.datetime(t.year, t.month, t.day, t.hour, t.minute),
            "z": z, "p": p, "T": T, "S": S, "has_salinity": has_s, "file": fn.name, **attrs}


def load_unis_hd(path: Path) -> list[dict]:
    """Read UNIS HD single-profile NetCDF files (a directory or one file).

    Each profile: ``lat``, ``lon`` (deg), ``time`` (datetime, UTC), ``z`` (m, positive up, from
    ``gsw.z_from_p``), ``p`` (dbar), ``T`` (deg C), ``S`` (PSU), ``file`` and the global attributes
    instrument_type, QC_indicator, cruise_ID. Levels with a missing value are dropped; casts
    without PSAL keep S = NaN (``has_salinity`` False) and are rejected by ``qc``.
    """
    path = Path(path)
    files = sorted(path.glob("*.nc")) if path.is_dir() else [path]
    out = []
    for fn in files:
        p = _read_one(fn)
        if p is not None:
            out.append(p)
    return out


# ----------------------------------------------------------------------------- selection
def select(profs, polygon_lonlat, months: set[int], years: tuple[int, int] | None) -> list[dict]:
    """Profiles inside the lon/lat polygon (straight edges in lon/lat), in ``months`` and, if
    given, in the inclusive year range ``years``."""
    poly = Polygon([tuple(v) for v in polygon_lonlat])
    months = set(int(m) for m in months)
    out = []
    for p in profs:
        if p["time"].month not in months:
            continue
        if years is not None and not (years[0] <= p["time"].year <= years[1]):
            continue
        if shapely.contains_xy(poly, p["lon"], p["lat"]):
            out.append(p)
    return out


# ----------------------------------------------------------------------------- quality control
def sigma0(prof) -> np.ndarray:
    """TEOS-10 potential density anomaly sigma0 (kg m^-3) of a profile (gsw)."""
    p = prof.get("p")
    if p is None:
        p = gsw.p_from_z(prof["z"], prof["lat"])
    SA = gsw.SA_from_SP(prof["S"], p, prof["lon"], prof["lat"])
    CT = gsw.CT_from_t(SA, prof["T"], p)
    return gsw.sigma0(SA, CT)


def inversion_fraction(prof, smooth_m: float = 3.0, tol: float = 4.0e-3) -> float:
    """Fraction of the 1 m steps of the profile where sigma0 decreases with depth by more than
    ``tol`` (kg m^-3), after interpolation to 1 m and a running mean over ``smooth_m`` metres."""
    d = -np.asarray(prof["z"], float)             # depth, m, increasing
    sig = sigma0(prof)
    o = np.argsort(d)
    d, sig = d[o], sig[o]
    grid = np.arange(np.ceil(d[0]), np.floor(d[-1]) + 1e-9, 1.0)
    if grid.size < 3:
        return 0.0
    s1 = np.interp(grid, d, sig)
    w = max(1, int(round(smooth_m)))
    s1 = uniform_filter1d(s1, size=w, mode="nearest")
    ds = np.diff(s1)
    return float(np.mean(ds < -tol))


def qc(prof, *, s_range=(0.0, 35.5), t_range=(-2.1, 10.0), smooth_m: float = 3.0,
       max_inversion_fraction: float = 0.10, inversion_tol_kg_m3: float = 4.0e-3,
       min_levels: int = 5) -> dict | None:
    """Quality control of one profile. Returns a copy (NaN levels dropped, key ``qc`` added) or None.

    Rejects the profile if any S is outside ``s_range`` or any T outside ``t_range`` (the whole
    cast: an out-of-range value flags a bad sensor or calibration, the UNIS HD caveat), if fewer
    than ``min_levels`` valid levels remain, or if more than ``max_inversion_fraction`` of its 1 m
    steps are density inversions (TEOS-10 sigma0, after a ``smooth_m`` running mean; steps with
    a decrease larger than ``inversion_tol_kg_m3``; 0.004 kg m^-3 is the sigma0 change of a 0.005
    salinity error, the order of the CTD salinity accuracy, so sensor noise in well-mixed winter
    columns is not counted as an inversion).
    """
    T, S, z = (np.asarray(prof[k], float) for k in ("T", "S", "z"))
    ok = np.isfinite(T) & np.isfinite(S) & np.isfinite(z)
    q = dict(prof)
    for k in ("T", "S", "z", "p"):
        if k in prof and prof[k] is not None:
            q[k] = np.asarray(prof[k], float)[ok]
    if ok.sum() < min_levels:
        return None
    if np.any((q["S"] < s_range[0]) | (q["S"] > s_range[1])):
        return None
    if np.any((q["T"] < t_range[0]) | (q["T"] > t_range[1])):
        return None
    frac = inversion_fraction(q, smooth_m, inversion_tol_kg_m3)
    if frac > max_inversion_fraction:
        return None
    q["qc"] = {"inversion_fraction": frac, "n_levels": int(ok.sum())}
    return q


# ----------------------------------------------------------------------------- representative
def rho_prime(T, S, eos=EOS_EXP4):
    """Linear-EOS density anomaly divided by rhoNil, up to a constant: sBeta S - tAlpha T
    (model/src/find_rho.F:104-107). Dimensionless."""
    return eos["sBeta"] * np.asarray(S, float) - eos["tAlpha"] * np.asarray(T, float)


def _pava_increasing(y, w):
    """Weighted least-squares non-decreasing fit (pool-adjacent-violators algorithm)."""
    vals, wts, cnt = [], [], []
    for yi, wi in zip(y, w):
        vals.append(float(yi))
        wts.append(float(wi))
        cnt.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            wv = wts[-2] + wts[-1]
            v = (vals[-2] * wts[-2] + vals[-1] * wts[-1]) / wv
            n = cnt[-2] + cnt[-1]
            vals[-2:], wts[-2:], cnt[-2:] = [v], [wv], [n]
    return np.repeat(vals, cnt)


def _thickness(zc):
    f = np.empty(zc.size + 1)
    f[1:-1] = 0.5 * (zc[1:] + zc[:-1])
    f[0] = 0.0
    f[-1] = zc[-1] - (f[-2] - zc[-1])
    return f[:-1] - f[1:]


def stabilize(zc, T, S, eos=EOS_EXP4) -> tuple[np.ndarray, int]:
    """Minimal static-stability adjustment under the linear EOS.

    The density anomaly r = sBeta S - tAlpha T must not decrease with depth (N^2 >= 0 at every
    interface). r is replaced by its thickness-weighted least-squares non-decreasing fit (pool
    adjacent violators), and the change is put on S only, S' = S + (r' - r)/sBeta, so T keeps its
    median values. Returns the new S and the number of levels whose S changed.
    """
    zc = np.asarray(zc, float)
    if np.any(np.diff(zc) >= 0):
        raise ValueError("zc must be ordered from the surface down (strictly decreasing)")
    r = rho_prime(T, S, eos)
    r2 = _pava_increasing(r, _thickness(zc))
    S2 = np.asarray(S, float) + (r2 - r) / eos["sBeta"]
    changed = np.abs(S2 - S) > 1e-12
    return S2, int(changed.sum())


def _on_grid(prof, zc, surface_extend_m: float):
    z, T, S = prof["z"], prof["T"], prof["S"]
    o = np.argsort(z)                              # ascending for np.interp
    z, T, S = z[o], T[o], S[o]
    Ti = np.interp(zc, z, T, left=np.nan, right=np.nan)
    Si = np.interp(zc, z, S, left=np.nan, right=np.nan)
    top = z[-1]                                    # shallowest datum (largest z)
    if -top <= surface_extend_m:
        up = zc > top
        Ti[up], Si[up] = T[-1], S[-1]
    return Ti, Si


def representative(profs, z_grid, *, eos=EOS_EXP4, min_count: int = 5,
                   surface_extend_m: float = 5.0, info: dict | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Median T and S per model level, made statically stable.

    ``z_grid``: cell-centre heights (m, positive up, surface first). Each profile is interpolated
    linearly onto the levels within its data range (and repeated up to the surface if its
    shallowest datum is within ``surface_extend_m``). The median at a level uses the profiles that
    reach it, if there are at least ``min_count``; the levels below the last such level take the
    last value ("below the last datum the last value is repeated") and the levels above the first
    one take the first value. Finally ``stabilize`` makes N^2 >= 0 with the minimal change of S.
    ``info`` (optional dict) receives counts and the number of adjusted levels.
    """
    zc = np.asarray(z_grid, float)
    if not profs:
        raise ValueError("no profiles")
    if np.any(np.diff(zc) >= 0):
        raise ValueError("z_grid must be ordered from the surface down")
    TT = np.full((len(profs), zc.size), np.nan)
    SS = np.full_like(TT, np.nan)
    for i, p in enumerate(profs):
        TT[i], SS[i] = _on_grid(p, zc, surface_extend_m)
    count = np.sum(np.isfinite(TT) & np.isfinite(SS), axis=0)
    good = count >= min_count
    if not good.any():
        raise ValueError(f"no level has {min_count} profiles")
    with warnings.catch_warnings():               # all-NaN levels below the deepest cast
        warnings.simplefilter("ignore", RuntimeWarning)
        T = np.where(good, np.nanmedian(TT, axis=0), np.nan)
        S = np.where(good, np.nanmedian(SS, axis=0), np.nan)
    kg = np.flatnonzero(good)
    k0, k1 = kg[0], kg[-1]
    # interior gaps (should not happen with sorted casts) are filled linearly
    gap = ~good
    gap[:k0] = gap[k1 + 1:] = False
    if gap.any():
        T[gap] = np.interp(zc[gap][::-1], zc[good][::-1], T[good][::-1])[::-1]
        S[gap] = np.interp(zc[gap][::-1], zc[good][::-1], S[good][::-1])[::-1]
    T[:k0], S[:k0] = T[k0], S[k0]
    T[k1 + 1:], S[k1 + 1:] = T[k1], S[k1]
    S_med = S.copy()
    S, n_adj = stabilize(zc, T, S, eos)
    if info is not None:
        info.update({
            "n_profiles": len(profs),
            "count_per_level": count.tolist(),
            "last_data_z_m": float(zc[k1]),
            "first_data_z_m": float(zc[k0]),
            "n_levels_extended_below": int(zc.size - 1 - k1),
            "n_levels_extended_above": int(k0),
            "n_levels_interior_gap": int(gap.sum()),
            "n_levels_adjusted": n_adj,
            "max_abs_dS_adjust": float(np.max(np.abs(S - S_med))),
        })
    return T, S
