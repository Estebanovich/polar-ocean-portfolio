"""Validation V2-V4 of the Exp. 4 fjord runs against the frozen PD0 (notes/pd0_preregistro.md 4).

``THRESHOLDS`` is the transcription of the PD0 thresholds (sections 3-5); tests/test_validate.py
parses the PD0 text and fails if any number here differs. Nothing here redefines a threshold.

Criteria (PD0 section 4, verbatim in substance):

* **V2a** (F1): the peak of the horizontal kinetic energy spectrum inside the near-inertial band
  [0.85 f, 1.20 f) at the **shelf station** lies within 1 df of f at 78 N; df = 1/T_fin, the
  Rayleigh resolution of the analysed record.
* **V2b** (F1): the internal Kelvin wave phase speed, from the arrival time of the pycnocline
  displacement front at >= 3 successive coastal stations of the inner fjord over the along-coast
  distance, lies in [0.8 c1(H mean), 1.2 c1(H max)] = [0.28, 0.62] m s^-1 (theory.json).
* **V3** (tide only): (1) presence: M2 pycnocline displacement trapped at the coast with an
  offshore decay scale of 1-8 km; (2) cyclonic sense: coast on the right, eastward along the south
  coast of the inner fjord; (3) phase speed within a factor 2 of 0.42 m s^-1 (0.21-0.84). The
  reference is Stoylen & Fer (2014, NPG 21:87-100, sections 4.2-4.3, p. 93-94), read online for
  the PD0; it is NOT in the local libraries (portfolio/literatura, thesis Literatura) [verificar].
* **V4** (F1 at 200 m against V4 at 400 m): V4a the net outward baroclinic flux over 0-T_fin
  summed over both sounds differs by <= 25 %; V4b the flux fraction of each resolved band
  differs by <= 0.10 in absolute value. Per-sound fluxes are descriptive.

Method choices the PD0 does not fix (fixed here, Task 6 part A, 2026-11-12, BEFORE any campaign
output was analysed; they are labelled "method, not PD0" in validation.json):

* pycnocline level k* = the level of maximum background N^2 (budget.background_profile and
  n2_background of the t = 0 density); displacement zeta = -rho'/(d rho_b/dz) at k*
  (linear isopycnal displacement, Gill 1982 sec. 6.4 [verificar]) [m, + up];
* shelf station: the wet cell outside the inner fjord at >= 2 n_sponge cells from every open
  boundary that maximises min(distance to land, distance to the open boundaries) (V2a: "W >> R1");
* coastal stations: inner-fjord cells on the south side of the local fjord centreline, holding
  k* + 1, with depth >= 1/2 H mean (31 m) and at <= R1(H mean) = 2.47 km (12 cells) from land,
  nearest to land at targets every ``spacing`` along the fjord axis (V2b: 5 km; V3: 2 km);
  along-coast distance = shortest 8-connected path through the south coastal strip (cells
  <= 12 from land). The V3 amplitude map uses the same depth rule;
* V2a uses the surface level (k = 0) velocity at the T point, as the thesis shelf test
  (analysis/v1_thesis.py);
* V2b arrival time: first time |zeta| reaches half of its record maximum at the station (linear
  interpolation between records); c = least-squares slope of s against arrival time;
* V3: M2 harmonic fit (amplitude, phase lag) of zeta over the records after the one-cycle tide
  ramp; decay scale = -1/slope of the log of the median M2 amplitude per 200 m distance-to-land
  bin (cells of the south half of the inner fjord, the coast of the stations; bins with >= 5
  cells, up to 4 km); c = omega / (d phase / d s) along the south coast;
* V4: flux of Kang & Fringer (2012), the H4 definition after PD0 Enmiendas 1-2 (pressure work
  F_bc = sum u_n' p' dz on the sound faces, + out of the inner fjord); the band partition is the
  exact cospectral (Parseval) split of the time mean of u_n' p' (no window, so the bands add up to
  the mean, the low band holding the product of the means).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp4_arctic_bayiw.analysis import budget, ice_check, spectra  # noqa: E402

PKG = Path(__file__).resolve().parents[1]
THEORY = PKG / "figures" / "theory.json"

# ---- frozen thresholds (notes/pd0_preregistro.md; tests/test_validate.py checks the text) -----
THRESHOLDS = {
    "V1": {"exact_max_abs_deta_m": 1e-10,             # PD0 4, V1-exacto
           "spectral_psd_max_rel_diff": 0.10},        # PD0 4, V1-espectral (ii)
    "V2a": {"near_inertial_band_f": (0.85, 1.2),      # PD0 3, band [0.85 f, 1.20 f)
            "max_peak_offset_df": 1.0},               # PD0 4, "a <= 1 df de f"
    "V2b": {"c1_factor_mean_depth": 0.8, "c1_factor_max_depth": 1.2,
            "range_m_s": (0.28, 0.62),                # PD0 4, fitted EOS
            "min_stations": 3.0},
    "V3": {"decay_scale_km": (1.0, 8.0), "sense": "cyclonic: eastward along the south coast",
           "c_factor": 2.0, "c_ref_m_s": 0.42, "c_range_m_s": (0.21, 0.84)},
    "V4a": {"max_rel_diff": 0.25},
    "V4b": {"max_abs_diff": 0.10},
    "factor_min": 1.5,                                # PD0 5, regla general
    "H4a": {"phi_ni_support": 0.5, "phi_k_alternative": 0.5, "low_band_max_frac": 0.25},
    "H4b": {"energy_ratio_max": 0.5},
    "H4c": {"ratio_min": 1.5},
    "H4d": {"ratio_min": 1.5},
}
RESOLVED = ("kelvin_subinertial", "near_inertial", "superinertial", "seiche")
METHOD_NOTE = "method, not PD0: fixed in Task 6 part A before any campaign output was analysed"


def jsonable(o):
    """json.dumps default for numpy values."""
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(type(o))


def load_theory(path: Path | None = None) -> dict:
    return json.loads(Path(path or THEORY).read_text())


def v2b_range(th: dict) -> tuple[float, float]:
    """[0.8 c1(H mean), 1.2 c1(H max)] [m s^-1], summer, fitted EOS (PD0 4, V2b; theory.json)."""
    hmax, hmean = summer_cases(th)
    c_mean, c_max = hmean["c_n"][0], hmax["c_n"][0]
    T = THRESHOLDS["V2b"]
    return T["c1_factor_mean_depth"] * c_mean, T["c1_factor_max_depth"] * c_max


def bands_hz(th: dict, record_h: int = 280) -> dict:
    """PD0 section 3 bands of the fjord at 78 N for the record length (theory.json, cph -> Hz)."""
    b = th["bands"][f"fjord_78N_{record_h}h"]["bands"]
    return {x["name"]: (x["lo_cph"] / 3600.0, x["hi_cph"] / 3600.0) for x in b}


# ---- pure evaluators ---------------------------------------------------------------------------

def v2a_inertial_peak(u, v, t, f: float) -> dict:
    """V2a. u, v [m s^-1] at the shelf station, t [s] equally spaced, f [rad s^-1].
    HKE spectrum = P_u + P_v (spectra.psd, one Hann window, [m^2 s^-2 Hz^-1]); peak inside
    [0.85, 1.20) f/(2 pi); pass if |f_peak - f/(2 pi)| <= 1 df (PD0 4)."""
    t = np.asarray(t, float)
    dt = float(t[1] - t[0])
    fr, Pu, df = spectra.psd(u, dt)
    _, Pv, _ = spectra.psd(v, dt)
    P = Pu + Pv
    f_hz = f / (2 * np.pi)                                    # Hz, cycles s^-1
    lo, hi = THRESHOLDS["V2a"]["near_inertial_band_f"]
    sel = (fr >= lo * f_hz) & (fr < hi * f_hz)
    k = np.flatnonzero(sel)[np.argmax(P[sel])]
    off = abs(fr[k] - f_hz) / df
    return {"f_hz": f_hz, "f_peak_hz": float(fr[k]), "df_hz": df, "offset_df": float(off),
            "pass": bool(off <= THRESHOLDS["V2a"]["max_peak_offset_df"]),
            "spectrum": {"freq_hz": fr, "hke": P}, "band_f": [lo, hi]}


def arrival_time(z, t, frac: float = 0.5) -> float | None:
    """First time |z| reaches frac * max|z| (linear interpolation between records) [s]."""
    a = np.abs(np.asarray(z, float))
    if not np.isfinite(a).all() or a.max() <= 0:
        return None
    thr = frac * a.max()
    k = int(np.argmax(a >= thr))
    if k == 0:
        return float(t[0])
    return float(t[k - 1] + (thr - a[k - 1]) / (a[k] - a[k - 1]) * (t[k] - t[k - 1]))


def v2b_kelvin_speed(zeta, t, s, crange) -> dict:
    """V2b. zeta (nt, ns) pycnocline displacement [m] at successive coastal stations, s (ns,)
    along-coast distance [m], crange the frozen range [m s^-1]. c = slope of s on arrival time
    (least squares) [m s^-1]; pass if >= 3 stations have a front and c is inside crange."""
    z = np.asarray(zeta, float)
    s = np.asarray(s, float)
    ta = [arrival_time(z[:, k], t) for k in range(z.shape[1])]
    ok = [k for k, x in enumerate(ta) if x is not None]
    n = len(ok)
    res = {"s_m": s, "arrival_s": [np.nan if x is None else x for x in ta], "n_stations": n,
           "range_m_s": list(crange), "method": METHOD_NOTE}
    if n < THRESHOLDS["V2b"]["min_stations"]:
        return {**res, "c_m_s": None, "pass": False, "reason": "fewer than 3 stations with a front"}
    tt, ss = np.array([ta[k] for k in ok]), s[ok]
    if np.ptp(tt) == 0:
        return {**res, "c_m_s": np.inf, "pass": False, "reason": "simultaneous arrival"}
    c = float(np.polyfit(tt, ss, 1)[0])                      # m s^-1
    return {**res, "c_m_s": c, "pass": bool(crange[0] <= c <= crange[1])}


def harmonic(x, t, omega: float):
    """Least-squares fit x = a0 + A cos(omega t - phi): returns (A, phi) [units of x, rad]."""
    t = np.asarray(t, float)
    M = np.column_stack([np.ones_like(t), np.cos(omega * t), np.sin(omega * t)])
    c = np.linalg.lstsq(M, np.asarray(x, float), rcond=None)[0]
    return float(np.hypot(c[1], c[2])), float(np.arctan2(c[2], c[1]))


def decay_scale(amp, d) -> float:
    """Offshore e-folding scale [m] of amplitude amp(d): -1/slope of ln amp on d (least squares);
    inf if amp does not decrease."""
    amp, d = np.asarray(amp, float), np.asarray(d, float)
    ok = amp > 0
    slope = np.polyfit(d[ok], np.log(amp[ok]), 1)[0]
    return float(-1.0 / slope) if slope < 0 else float("inf")


def v3_kelvin(coast, s, offshore, d, t, omega: float, ramp_s: float = 0.0) -> dict:
    """V3. coast (nt, ns): zeta at the south-coast stations ordered eastward, s (ns,) [m];
    offshore: (nt, nd) series or (nd,) M2 amplitudes against distance to land d (nd,) [m];
    t [s]; omega the M2 frequency [rad s^-1]; records with t < ramp_s are dropped.
    Phase lag phi(s) of zeta = A cos(omega t - phi): eastward propagation <=> d phi/ds > 0, and
    c = omega / (d phi/ds) [m s^-1] (+ eastward). PD0 4: presence (decay 1-8 km), sense
    (cyclonic = eastward along the south coast, coast on the right), |c| in 0.21-0.84 m s^-1."""
    t = np.asarray(t, float)
    m = t >= ramp_s
    T = THRESHOLDS["V3"]
    co = np.asarray(coast, float)[m]
    fits = [harmonic(co[:, k], t[m], omega) for k in range(co.shape[1])]
    amp_c = np.array([a for a, _ in fits])
    ph = np.unwrap(np.array([p for _, p in fits]))
    kk = float(np.polyfit(np.asarray(s, float), ph, 1)[0])   # rad m^-1
    c = omega / kk if kk != 0 else np.inf
    off = np.asarray(offshore, float)
    amp_d = np.array([harmonic(off[m][:, k], t[m], omega)[0] for k in range(off.shape[1])]) \
        if off.ndim == 2 else off
    L = decay_scale(amp_d, d)
    lo, hi = T["decay_scale_km"]
    pres = {"decay_scale_km": L / 1e3, "d_m": d, "amplitude": amp_d,
            "pass": bool(lo <= L / 1e3 <= hi)}
    sense = {"dphase_ds_rad_per_m": kk, "direction": "east" if kk > 0 else "west",
             "pass": bool(kk > 0)}
    clo, chi = T["c_range_m_s"]
    speed = {"c_m_s": float(c), "s_m": s, "phase_rad": ph, "amplitude": amp_c,
             "pass": bool(clo <= abs(c) <= chi)}
    return {"presence": pres, "sense": sense, "speed": speed,
            "pass": bool(pres["pass"] and sense["pass"] and speed["pass"]),
            "reference": "Stoylen & Fer 2014 sec. 4.2-4.3 (PD0 4) [verificar: not in the local library]",
            "method": METHOD_NOTE}


def band_flux(up, pp, dz, w, dt: float, bands: dict) -> dict:
    """Exact band split of the time-mean flux sum_faces w sum_z dz <u' p'> [W] (+ into the box).
    up, pp (nt, nz, nf); dz (nz, nf) [m]; w (nf,) signed face width [m]; bands {name: (lo, hi)}
    [Hz], half-open, the last one closed at Nyquist. Parseval: <u p> = sum_k Re(U_k P_k*) / nt^2
    over all FFT bins k (no window, mean included: the product of the means falls in bin 0)."""
    up, pp = np.asarray(up, float), np.asarray(pp, float)
    nt = up.shape[0]
    X, Y = np.fft.fft(up, axis=0), np.fft.fft(pp, axis=0)
    co = (X * np.conj(Y)).real / nt ** 2                     # (nt, nz, nf), sums to <u p>
    per_bin = (co * (np.asarray(dz, float) * np.asarray(w, float)[None])[None]).sum(axis=(1, 2))
    fk = np.abs(np.fft.fftfreq(nt, d=dt))
    items = sorted(bands.items(), key=lambda kv: kv[1][0])
    out = {}
    for q, (name, (lo, hi)) in enumerate(items):
        last = q == len(items) - 1
        sel = (fk >= lo) & ((fk <= hi + 1e-15) if last else (fk < hi))
        out[name] = float(per_bin[sel].sum())
    if abs(sum(out.values()) - per_bin.sum()) > 1e-9 * max(abs(per_bin).sum(), 1e-300):
        raise ValueError("the bands do not cover every frequency bin")
    return out


def v4_resolution(f1: dict, v4: dict) -> dict:
    """V4 from the flux summaries of F1 (200 m) and V4 (400 m): {F_bc_out_J, band_out_J}.
    V4a |F_V4 - F_F1| / |F_F1| <= 0.25; V4b |phi_b(V4) - phi_b(F1)| <= 0.10 for every resolved
    band, phi_b = band / sum of the four resolved bands (as H4a)."""
    rel = abs(v4["F_bc_out_J"] - f1["F_bc_out_J"]) / abs(f1["F_bc_out_J"])

    def frac(r):
        tot = sum(r["band_out_J"][b] for b in RESOLVED if b in r["band_out_J"])
        return {b: r["band_out_J"][b] / tot for b in RESOLVED if b in r["band_out_J"]}

    p1, p4 = frac(f1), frac(v4)
    diff = {b: abs(p4[b] - p1[b]) for b in p1}
    a = {"F1_out_J": f1["F_bc_out_J"], "V4_out_J": v4["F_bc_out_J"], "rel_diff": rel,
         "pass": bool(rel <= THRESHOLDS["V4a"]["max_rel_diff"])}
    b = {"phi_F1": p1, "phi_V4": p4, "abs_diff": diff,
         "pass": bool(max(diff.values()) <= THRESHOLDS["V4b"]["max_abs_diff"])}
    desc = {k: {"F1": f1.get("per_section_out_J", {}).get(k), "V4": v4.get("per_section_out_J", {}).get(k)}
            for k in set(f1.get("per_section_out_J", {})) | set(v4.get("per_section_out_J", {}))}
    return {"V4a": a, "V4b": b, "per_section_descriptive": desc, "pass": bool(a["pass"] and b["pass"])}


# ---- geometry: stations ------------------------------------------------------------------------

def pycnocline_level(rho_b, z, rho0: float, g: float = 9.81) -> int:
    """Level of maximum background N^2 = -(g/rho0) d rho_b/dz (budget.n2_background) [index]."""
    return int(np.argmax(budget.n2_background(rho_b, z, rho0, g)))


def displacement(rho, rho_b, z, k: int):
    """zeta = -(rho - rho_b(k)) / (d rho_b/dz at k) [m, + up] (linear isopycnal displacement)."""
    drdz = np.gradient(np.asarray(rho_b, float), np.asarray(z, float))[k]   # kg m^-4 (< 0)
    return -(np.asarray(rho, float) - rho_b[k]) / drdz


def _dist_land(wet):
    from scipy import ndimage
    return ndimage.distance_transform_edt(np.asarray(wet, bool))


def shelf_station(wet, inner, ob_sides, n_sponge: int, depth=None, min_depth: float = 0.0):
    """(j, i) of the wet cell outside ``inner`` at >= 2 n_sponge cells from every open boundary,
    deeper than min_depth, that maximises min(distance to land, distance to the open boundaries)
    [cells]: as far as possible from both coasts and boundary sponges."""
    wet = np.asarray(wet, bool)
    dob = ice_check.ob_distance_field(wet, ob_sides)
    cand = wet & ~np.asarray(inner, bool) & (dob >= 2 * n_sponge)
    if depth is not None:
        cand &= np.asarray(depth) >= min_depth
    dl = np.where(cand, np.minimum(_dist_land(wet), dob), -1.0)
    j, i = np.unravel_index(int(np.argmax(dl)), dl.shape)
    if dl[j, i] < 0:
        raise ValueError("no shelf-station candidate")
    return int(j), int(i)


def coastal_stations(wet, deep, inner, xc, yc, axis_deg: float, spacing_m: float,
                     side: str = "south", max_coast_cells: int = 5) -> dict:
    """Stations along one coast of the inner fjord, ordered along the axis (west to east for
    axis_deg in (-90, 90)). deep: cells deep enough for the analysis level. Returns {cells: [(j,
    i)], s_m: along-coast distance from the first station [m], a_m: along-axis coordinate}."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    wet, deep, inner = (np.asarray(x, bool) for x in (wet, deep, inner))
    X, Y = np.meshgrid(np.asarray(xc, float), np.asarray(yc, float))
    th = np.deg2rad(axis_deg)
    A = X * np.cos(th) + Y * np.sin(th)                       # along axis [m]
    C = -X * np.sin(th) + Y * np.cos(th)                      # across axis [m], + to the left
    dx = float(np.median(np.diff(xc)))
    fj = wet & inner
    bins = np.floor(A / dx).astype(int)
    mid = np.full(A.shape, np.nan)
    for b in np.unique(bins[fj]):
        sel = fj & (bins == b)
        mid[bins == b] = C[sel].mean()
    sidem = (C < mid) if side == "south" else (C > mid)
    dl = _dist_land(wet)
    strip = fj & sidem & (dl <= max_coast_cells)
    cand = strip & deep
    if not cand.any():
        raise ValueError("no coastal-station candidate")
    # targets span the coastal strip (not only the candidates), so a target without a candidate
    # is a gap in the "successive" stations of the PD0 (V2b) and is reported
    a0, a1 = A[strip].min(), A[strip].max()
    cells, gaps = [], []
    for q, at in enumerate(np.arange(a0 + 0.5 * spacing_m, a1, spacing_m)):
        sel = cand & (np.abs(A - at) <= 0.5 * spacing_m)
        if not sel.any():
            gaps.append({"target_index": q, "a_m": float(at)})
            continue
        jj, ii = np.nonzero(sel)
        key = np.lexsort((np.abs(A[jj, ii] - at), dl[jj, ii]))[0]
        cells.append((int(jj[key]), int(ii[key])))
    # along-coast distance: shortest 8-connected path through the strip
    idx = -np.ones(wet.shape, int)
    jj, ii = np.nonzero(strip)
    idx[jj, ii] = np.arange(jj.size)
    rows, cols, wts = [], [], []
    ny, nx = wet.shape
    for dj, di in ((0, 1), (1, 0), (1, 1), (1, -1)):
        j2, i2 = jj + dj, ii + di
        ok = (j2 >= 0) & (j2 < ny) & (i2 >= 0) & (i2 < nx)
        j2, i2, j1, i1 = j2[ok], i2[ok], jj[ok], ii[ok]
        ok2 = strip[j2, i2]
        rows += list(idx[j1[ok2], i1[ok2]])
        cols += list(idx[j2[ok2], i2[ok2]])
        wts += [dx * np.hypot(dj, di)] * int(ok2.sum())
    G = coo_matrix((wts, (rows, cols)), shape=(jj.size, jj.size)).tocsr()
    dist = dijkstra(G, directed=False, indices=idx[cells[0]])
    s = np.array([dist[idx[c]] for c in cells])
    if not np.isfinite(s).all():
        raise ValueError("coastal stations not connected through the coastal strip")
    return {"cells": cells, "s_m": s, "a_m": np.array([A[c] for c in cells]),
            "side": side, "max_coast_cells": max_coast_cells, "method": METHOD_NOTE,
            "gaps": [g for g in gaps if A[cand].min() < g["a_m"] < A[cand].max()],
            "successive": not any(A[cand].min() < g["a_m"] < A[cand].max() for g in gaps),
            "side_mask": fj & sidem}


def _public(st: dict) -> dict:
    """Station record for validation.json (without the side mask array)."""
    return {k: v for k, v in st.items() if k != "side_mask"}


# ---- flux on the sound faces (Kang & Fringer 2012, as budget.kf_budget) -------------------------

def face_profiles(ds, mask, face_groups: dict, rho_b=None, time=slice(None)) -> dict:
    """u_n' and p' on the faces of ``face_groups`` (the open faces of the box ``mask``), as
    budget.kf_budget computes the baroclinic pressure work (KF12 eqs. 15-16; Nash 2005 eq. 3 for
    p~): returns up, pp (nt, nz, nf) [m s^-1, Pa], dz (nz, nf) [m], w (nf,) signed face width
    [m] (+ when +x/+y points into the box), group (nf,), t (nt,) [s]. Only a small rectangle around
    each group is read per record."""
    rho0, g = float(ds.attrs["rho0"]), float(ds.attrs.get("g", budget.G))
    drF = ds["drF"].values
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
    up, pp = np.zeros((nt, nz, nf)), np.zeros((nt, nz, nf))
    dz = np.zeros((nz, nf))
    w = np.zeros(nf)
    for q, (kind, j, i, s, _) in enumerate(sel):
        dz[:, q] = drF * (hS[:, j, i] if kind == "y" else hW[:, j, i])
        w[q] = s * (dxG[j, i] if kind == "y" else dyG[j, i])
    # one rectangle per group (+1 halo)
    groups = {}
    for q, (kind, j, i, s, n) in enumerate(sel):
        groups.setdefault(n, []).append(q)
    for n, qs in groups.items():
        if np.any(dz[:, qs][int(np.flatnonzero((dz[:, qs] > 0).any(1)).max()) + 1:] != 0):
            raise AssertionError("face thickness below the last wet level")
        js = [sel[q][1] for q in qs]
        is_ = [sel[q][2] for q in qs]
        j0, j1, i0, i1 = max(min(js) - 1, 0), min(max(js) + 2, ny), max(min(is_) - 1, 0), min(max(is_) + 2, nx)
        # only the levels that are wet somewhere in the rectangle are read (I/O: the strait
        # rectangles are small, but each level is a separate read); the dry levels below carry
        # no face thickness and do not enter p~ above them, so the result is unchanged
        nk = int(np.flatnonzero((hC[:, j0:j1, i0:i1] > 0).any(axis=(1, 2))).max()) + 1
        dzC = drF[:nk, None, None] * hC[:nk, j0:j1, i0:i1]
        wetC = dzC > 0
        for r, n_rec in enumerate(tidx):
            win = dict(time=int(n_rec), k=slice(0, nk), j=slice(j0, j1), i=slice(i0, i1))
            U = ds["U"].isel(**win).values
            V = ds["V"].isel(**win).values
            R = ds["RHO"].isel(**win).values
            E = ds["ETA"].isel(time=int(n_rec), j=slice(j0, j1), i=slice(i0, i1)).values
            _progress(f"face_profiles {n}", r, nt)
            rp = np.where(wetC, R - rho_b[:nk, None, None], 0.0)
            pt = budget.pressure_tilde(rp, E, dzC, rho0, g)            # Pa, Nash 2005 eq. 3
            for q in qs:
                kind, j, i, _, _ = sel[q]
                a, b = j - j0, i - i0
                if kind == "y":
                    un = V[:, a, b] * (hS[:nk, j, i] > 0)
                    p = 0.5 * (pt[:, a - 1, b] + pt[:, a, b])
                else:
                    un = U[:, a, b] * (hW[:nk, j, i] > 0)
                    p = 0.5 * (pt[:, a, b - 1] + pt[:, a, b])
                dzq = dz[:nk, q]
                H = dzq.sum()
                Hs = H if H > 0 else 1.0
                ub, pb = (un * dzq).sum() / Hs, (p * dzq).sum() / Hs
                up[r, :nk, q], pp[r, :nk, q] = un - ub, p - pb           # KF12: u', p'
    return {"up": up, "pp": pp, "dz": dz, "w": w, "group": np.array([x[4] for x in sel]),
            "t": ds["time"].values[time].astype(float)}


def flux_summary(fp: dict, bands: dict | None = None) -> dict:
    """Outward KF baroclinic pressure work integrated over the record [J] (budget.cumulative,
    trapezoid, as budget.f_bc_pressure_out), per section, and (with ``bands``) per band [J] =
    -(band share of the time mean) x record length n dt.

    The record must start at t = 0 (the initial state at rest of a load.run dataset): the PD0
    integrates over 0-T_fin, and a record starting later would silently drop the first interval
    (ValueError). The band total, -dt sum F, differs from the trapezoid integral by the end terms
    dt (F_0 + F_n)/2 only; F_0 = 0 for a run that starts at rest."""
    Fq = (fp["up"] * fp["pp"] * fp["dz"][None]).sum(1) * fp["w"][None]     # (nt, nf) W, + in
    t = fp["t"]
    if t[0] != 0.0:
        raise ValueError(f"flux_summary: t[0] = {t[0]} s, the record must start at 0 (0-T_fin)")
    out = {"F_bc_out_J": float(-budget.cumulative(t, Fq.sum(1))[-1]),
           "per_section_out_J": {str(n): float(-budget.cumulative(t, Fq[:, fp["group"] == n].sum(1))[-1])
                                 for n in np.unique(fp["group"])}}
    if bands is not None:
        dt = float(t[1] - t[0])
        mean_in = band_flux(fp["up"], fp["pp"], fp["dz"], fp["w"], dt, bands)
        out["band_out_J"] = {b: -v * t.size * dt for b, v in mean_in.items()}
    return out


# ---- on load.run datasets (part B) ----------------------------------------------------------

def _progress(what: str, r: int, n: int, every: int = 120) -> None:
    if r % every == 0 or r == n - 1:
        print(f"[validate] {what}: record {r + 1}/{n}", file=sys.stderr, flush=True)


def level_values(ds, var: str, k: int, cells, rect=None):
    """(nt, ncells) of ds[var] at level k on the cells [(j, i)], and, with rect = (j0, j1, i0, i1),
    also the (nt, j1 - j0, i1 - i0) block. One contiguous read of the whole level per record
    (MDS layout: a level of one field is contiguous), instead of one read per cell and record."""
    jj = np.array([c[0] for c in cells], int)
    ii = np.array([c[1] for c in cells], int)
    nt = ds.sizes["time"]
    out = np.empty((nt, jj.size))
    blk = None if rect is None else np.empty((nt, rect[1] - rect[0], rect[3] - rect[2]))
    for n in range(nt):
        a = ds[var].isel(time=n, k=k).values
        out[n] = a[jj, ii]
        if blk is not None:
            blk[n] = a[rect[0]:rect[1], rect[2]:rect[3]]
        _progress(f"{var} level {k}", n, nt)
    return (out, blk) if rect is not None else out

def open_run(run_dir: Path):
    """load.run: require_ok on every stage, and the ice-OB check for ice runs."""
    from exp4_arctic_bayiw.analysis import load
    return load.run(Path(run_dir))


def _case_params(ds) -> dict:
    return json.loads((Path(ds.attrs["run_dir"]) / "input" / "case.json").read_text())["params"]


def _case_actual(ds) -> dict:
    return json.loads((Path(ds.attrs["run_dir"]) / "input" / "case.json").read_text())["actual"]


def station_rules(th: dict, dx: float) -> dict:
    """Coastal-station rules (method, not PD0): depth >= 1/2 H mean of the inner fjord and at most
    R1(H mean) from land, so that every station is inside the Kelvin trapping scale and in water
    where the pycnocline displacement is that of the mode-1 wave (theory.json, summer, H mean)."""
    return {"min_depth_m": 0.5 * th["fjord_geometry"]["mean_depth_m"],
            "max_coast_cells": int(summer_cases(th)[1]["rossby_radius_m"][0] // dx)}


def _geometry(ds, th: dict):
    wet = ds["hFacC"].values[0] > 0
    mask, groups = budget.fjord_box(ds)
    rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, ds["hFacC"].values > 0)
    k = pycnocline_level(rho_b, ds["RC"].values, float(ds.attrs["rho0"]), float(ds.attrs["g"]))
    dx = float(np.median(np.diff(ds["XC"].values[0])))
    rules = station_rules(th, dx)
    deep = (ds["hFacC"].values[min(k + 1, ds.sizes["k"] - 1)] > 0) & \
        (ds["Depth"].values >= rules["min_depth_m"])
    return wet, mask, groups, rho_b, k, deep, rules


def v2_from_run(ds, th: dict) -> dict:
    """V2a and V2b on F1 (0-T_fin, every record)."""
    p = _case_params(ds)
    wet, inner, _, rho_b, k, deep, rules = _geometry(ds, th)
    f = summer_cases(th)[0]["f_s-1"]
    t = ds["time"].values.astype(float)
    js, is_ = shelf_station(wet, inner, p["ob_sides"], int(p["n_sponge"]))
    uu = level_values(ds, "U", 0, [(js, is_), (js, is_ + 1)])
    vv = level_values(ds, "V", 0, [(js, is_), (js + 1, is_)])
    u, v = 0.5 * (uu[:, 0] + uu[:, 1]), 0.5 * (vv[:, 0] + vv[:, 1])          # T point, k = 0
    v2a = v2a_inertial_peak(u, v, t, f)
    v2a["station_ji"] = [js, is_]
    axis = th["fjord_geometry"]["axis_angle_deg"]
    st = coastal_stations(wet, deep, inner, ds["XC"].values[0], ds["YC"].values[:, 0], axis, 5000.0,
                          max_coast_cells=rules["max_coast_cells"])
    z = displacement(level_values(ds, "RHO", k, st["cells"]), rho_b, ds["RC"].values, k)
    v2b = v2b_kelvin_speed(z, t, st["s_m"], v2b_range(th))
    v2b.update(stations=_public(st), pycnocline_level=k, station_gaps=st["gaps"],
               stations_successive=st["successive"])
    return {"V2a": v2a, "V2b": v2b}


def v3_from_run(ds) -> dict:
    """V3 on the tide-only run, after the tide ramp of case.json (actual.tide.ramp_s)."""
    p = _case_params(ds)
    ramp = float(_case_actual(ds)["tide"]["ramp_s"])
    th = load_theory()
    wet, inner, _, rho_b, k, deep, rules = _geometry(ds, th)
    omega = 2 * np.pi / float(p["tide"]["period_s"])          # M2 [rad s^-1] (case.json)
    t = ds["time"].values.astype(float)
    st = coastal_stations(wet, deep, inner, ds["XC"].values[0], ds["YC"].values[:, 0],
                          th["fjord_geometry"]["axis_angle_deg"], 2000.0,
                          max_coast_cells=rules["max_coast_cells"])
    RC = ds["RC"].values
    # amplitude map at k* over the south side of the inner fjord (the coast of the stations),
    # binned by distance to land
    jj, ii = np.nonzero(inner & deep & st["side_mask"])
    j0, j1, i0, i1 = jj.min(), jj.max() + 1, ii.min(), ii.max() + 1
    cs, lev = level_values(ds, "RHO", k, st["cells"], rect=(j0, j1, i0, i1))     # one pass
    coast = displacement(cs, rho_b, RC, k)
    m = t >= ramp
    dx = float(np.median(np.diff(ds["XC"].values[0])))
    dl = _dist_land(wet)[jj, ii] * dx                                             # m
    amp = np.array([harmonic(displacement(lev[:, j - j0, i - i0], rho_b, RC, k)[m], t[m], omega)[0]
                    for j, i in zip(jj, ii)])
    edges = np.arange(0.0, 4000.0 + dx, dx) + 0.5 * dx
    cen, med = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        s_ = (dl >= a) & (dl < b)
        if s_.sum() >= 5:
            cen.append(0.5 * (a + b))
            med.append(float(np.median(amp[s_])))
    r = v3_kelvin(coast, st["s_m"], np.array(med), np.array(cen), t, omega, ramp_s=ramp)
    r.update(stations=_public(st), pycnocline_level=k, ramp_s=ramp)
    return r


def _seg_dist(px, py, p0, p1) -> float:
    """Distance [m] from (px, py) to the segment p0-p1."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    u = np.clip(((px - p0[0]) * d[0] + (py - p0[1]) * d[1]) / (d @ d), 0.0, 1.0)
    return float(np.hypot(px - p0[0] - u * d[0], py - p0[1] - u * d[1]))


def coarse_box(fine_mask, wet_coarse, xc, yc, sections: dict, max_dist_m: float = 2000.0):
    """The fine (200 m) inner-fjord box mapped onto a coarser grid (V4, 400 m): a coarse cell is in
    the box if it is wet and at least half of its children are in the fine box. Each open face
    (budget.box_faces with the coarse wet mask) is assigned to the nearest section line; a face
    farther than ``max_dist_m`` from every section is an error (the box would leak elsewhere).

    Why (Task 6 part B, method, not PD0; decided after budget.fjord_box failed on the V4 grid and
    BEFORE any V4 flux was computed): at 400 m the Akselsundet line (from bathy_report.json) ends
    west of the coarsened tip of Akseloya, so the line-cut box of fjord_box leaks along the west
    flank of the island for any extension of 1-3 cells. Mapping the 200 m box gives the same box
    (532.2 km2 against 531.6 km2) with 10 open faces, all at the two sounds, so F1 and V4 compare
    the flux through the same sections."""
    wc = np.asarray(wet_coarse, bool)
    nyc, nxc = wc.shape
    mf = np.asarray(fine_mask, bool)
    fy, fx = mf.shape[0] // nyc, mf.shape[1] // nxc
    if fy * nyc != mf.shape[0] or fx * nxc != mf.shape[1]:
        raise ValueError("the fine grid is not an integer refinement of the coarse grid")
    frac = mf.reshape(nyc, fy, nxc, fx).mean(axis=(1, 3))
    mask = (frac >= 0.5) & wc
    faces = budget.box_faces(mask, wet=wc)
    dx, dy = float(np.median(np.diff(xc))), float(np.median(np.diff(yc)))
    groups = {n: {"x": np.zeros((nyc, nxc), bool), "y": np.zeros((nyc, nxc), bool)} for n in sections}
    for kind in ("x", "y"):
        for j, i, _ in faces[kind]:
            px = xc[i] - 0.5 * dx if kind == "x" else xc[i]
            py = yc[j] if kind == "x" else yc[j] - 0.5 * dy
            d = {n: _seg_dist(px, py, *v) for n, v in sections.items()}
            n = min(d, key=d.get)
            if d[n] > max_dist_m:
                raise RuntimeError(f"open {kind} face ({j}, {i}) is {d[n]:.0f} m far from every section")
            groups[n][kind][j, i] = True
    return mask, groups


def _sections():
    rep = json.loads((PKG / "figures" / "bathy_report.json").read_text())
    return {n: (v["p0"], v["p1"]) for n, v in rep["sounds"].items()}


def flux_from_run(ds, th: dict, box=None) -> dict:
    mask, groups = box if box is not None else budget.fjord_box(ds)
    fp = face_profiles(ds, mask, groups)
    out = flux_summary(fp, bands_hz(th))
    out["box_cells"] = int(np.asarray(mask).sum())
    out["open_faces"] = {n: int(g["x"].sum() + g["y"].sum()) for n, g in groups.items()}
    return out


def v4_from_runs(ds_f1, ds_v4) -> dict:
    """V4: F1 with budget.fjord_box; V4 with that box mapped onto the 400 m grid (coarse_box)."""
    th = load_theory()
    m1, g1 = budget.fjord_box(ds_f1)
    m4, g4 = coarse_box(m1, ds_v4["hFacC"].values[0] > 0, ds_v4["XC"].values[0],
                        ds_v4["YC"].values[:, 0], _sections())
    r = v4_resolution(flux_from_run(ds_f1, th, (m1, g1)), flux_from_run(ds_v4, th, (m4, g4)))
    r["box_method"] = ("F1: budget.fjord_box; V4: coarse_box (the F1 box mapped to 400 m; "
                       "method, not PD0, fixed in Task 6 part B before any V4 flux was computed)")
    return r


# ---- PD2 blocks (Esteban, PD2 2026-11-13): rebuilt from the post hoc JSONs, never by hand -------

def summer_cases(th: dict) -> tuple[dict, dict]:
    """(H max, H mean) theory cases of the fjord summer profile, found by content rather than by
    key: profile 'vmf_summer' with the fitted EOS; H mean is the case whose mode depth equals the
    Merian (mean) depth, H max the deeper one (theory_numbers.py, PD0 2.2)."""
    cs = [c for c in th["cases"].values()
          if c.get("profile") == "vmf_summer" and c.get("eos") == "fitted_vmf"]
    mean = [c for c in cs if np.isclose(c["H_modes_m"], c["H_merian_m"])]
    deep = [c for c in cs if c["H_modes_m"] > c["H_merian_m"] * (1 + 1e-9)]
    if len(mean) != 1 or len(deep) != 1:
        raise ValueError(f"theory.json: expected one summer H-mean and one summer H-max case, found "
                         f"{len(mean)} and {len(deep)}")
    return deep[0], mean[0]


def _poincare_width_note(th: dict) -> str:
    """omega/f = (1 + (c1 pi / (W f))^2)^(1/2) (transverse Poincare mode, Gill 1982 sec. 10.4
    [verificar]) for the 60 km box width, and the width that 1.09-1.14 f would need."""
    hmax, hmean = summer_cases(th)
    f = hmax["f_s-1"]
    cs = (hmean["c_n"][0], hmax["c_n"][0])
    r60 = [float(np.sqrt(1 + (np.pi * c / (60e3 * f)) ** 2)) for c in cs]
    need = {q: [float(np.pi * c / f / np.sqrt(q ** 2 - 1) / 1e3) for c in cs] for q in (1.09, 1.14)}
    return (f"The 60 km width alone does not explain it: omega/f = (1 + (c1 pi/(W f))^2)^(1/2) = "
            f"{r60[0]:.3f}-{r60[1]:.3f} for W = 60 km and c1 = {cs[0]:.2f}-{cs[1]:.2f} m/s; 1.09 f would need "
            f"W = {min(need[1.09]):.0f}-{max(need[1.09]):.0f} km and 1.14 f W = {min(need[1.14]):.0f}-"
            f"{max(need[1.14]):.0f} km (Gill 1982 sec. 10.4 [verificar]).")


def _wind_pulse() -> dict:
    """centre_h and sigma_h of the fjord wind pulse (cases.yaml configs.fjord_vmf.wind, PD0 7)."""
    import yaml
    w = yaml.safe_load((PKG / "cases.yaml").read_text())["configs"]["fjord_vmf"]["wind"]
    return {"center_h": float(w["center_h"]), "sigma_h": float(w["sigma_h"])}


def _v2b_text(v2b: dict, pulse: dict) -> str:
    """V2b class (a) text with its numbers (fitted c, stations that cross half-maximum before the
    pulse is over, i.e. before centre + 3 sigma, and the range of their arrival times)."""
    arr = np.asarray(v2b["arrival_s"], float) / 3600.0
    end = pulse["center_h"] + 3 * pulse["sigma_h"]
    early = arr[arr <= end]
    return ("Our error in the definition of the pre-registered metric: with a spatially uniform wind "
            "every coastal station responds locally during the pulse (half-maximum reached at "
            f"{early.min():.1f}-{early.max():.1f} h at {early.size} of {arr.size} stations, before the "
            f"pulse centre + 3 sigma = {end:.1f} h; all arrivals {arr.min():.1f}-{arr.max():.1f} h), so "
            "the first-arrival time measures the forced response, not a propagating Kelvin front, and "
            f"the fitted c = {v2b['c_m_s']:.2f} m/s is not a phase speed.")


def _v4_text(v4a: dict, d4: dict) -> str:
    """V4 class (c) text with its numbers, from validation.json V4a and pd2_v4_diagnostics.json."""
    f1, v4 = d4["F1"], d4["V4"]
    low = sum(x["band_out_J"]["low"] for x in f1.values())
    res = sum(sum(v for k, v in x["band_out_J"].items() if k != "low") for x in f1.values())
    net = sum(x["net_out_J"] for x in f1.values())
    per = "; ".join(
        f"{s}: {f1[s]['gross_J'] / 1e9:.2f} -> {v4[s]['gross_J'] / 1e9:.2f} GJ "
        f"(x{v4[s]['gross_J'] / f1[s]['gross_J']:.2f}), {f1[s]['wet_levels_max']} -> "
        f"{v4[s]['wet_levels_max']} wet levels, {f1[s]['open_faces']} -> {v4[s]['open_faces']} faces"
        for s in f1)
    return ("Numerical choice (resolution): the net outward baroclinic pressure work changes sign "
            f"between 200 m ({v4a['F1_out_J'] / 1e9:+.2f} GJ) and 400 m ({v4a['V4_out_J'] / 1e9:+.2f} GJ). "
            f"At 200 m the net is {100 * low / net:.1f} % in the unresolved low band (time-mean exchange); "
            f"the four resolved bands add up to {res / 1e9:+.3f} GJ, so the V4b fractions are "
            f"ill-conditioned. Gross strait flux, 200 -> 400 m: {per}. 400 m does not resolve the strait "
            "exchange. PD0 4: if V4 fails, the limitation is declared.")


def pd2_merge(val: dict, figures: Path | None = None, th: dict | None = None,
              pulse: dict | None = None) -> dict:
    """Rebuild the PD2 blocks of validation.json (classification, pd2_decision, column_test,
    descriptive_free_wave, note_on_lag_correlation) from figures/pd2_column_test.json,
    pd2_v2_diagnostics.json, pd2_v2b_free_wave.json, pd2_v4_diagnostics.json and
    pd2_v4_descriptive.json, with fixed classification text. Verdicts ("pass") are never touched.
    The documented path after any `evaluate` (which does not hold these blocks) is
    `validate.py pd2-merge`."""
    figures = Path(figures or PKG / "figures")
    th = th or load_theory()
    rd = lambda n: json.loads((figures / n).read_text())  # noqa: E731
    col, d2, fw = rd("pd2_column_test.json"), rd("pd2_v2_diagnostics.json"), rd("pd2_v2b_free_wave.json")
    d4 = rd("pd2_v4_diagnostics.json"); rd("pd2_v4_descriptive.json")        # must exist
    pulse = pulse or _wind_pulse()
    out = json.loads(json.dumps(val))
    pm, rot = d2["V2a_peak_map"], d2["V2a_station_rotary"]
    pc = pm["peak_over_f_percentiles"]
    out["V2a"]["classification"] = {
        "class": "(c)",
        "text": ("(a) in the inertial frequency excluded by the PD2 column test (doubly periodic f-plane "
                 f"slab, F1 configuration, peak at {col['peak_over_f']:.3f} f, {col['offset_df']:.2f} df), "
                 "which tests the Coriolis term and time stepping only: with useOBCS=.FALSE. it does not "
                 "test the open-boundary code, the sponge or lateral viscosity. Class (c) by the rule "
                 "fixed before the test. The cause of the domain-wide blue shift of F1 (median "
                 f"{pc['50']:.2f} f, 5-95 % {pc['5']:.2f}-{pc['95']:.2f} f, clockwise, shelf peak "
                 f"{rot['cw_peak_over_f']:.3f} f) is not identified. Candidates, none tested: bounded "
                 "geometry and lateral modes; the open boundaries and sponge; background vorticity or "
                 "mean flow. " + _poincare_width_note(th)),
        "evidence": "figures/pd2_v2_diagnostics.json, figures/pd2_column_test.json (post hoc)",
        "post_hoc_evidence": True,
        "column_test": {**{k: col[k] for k in ("peak_over_f", "offset_df", "within_1df") if k in col},
                        "rule": "fixed before running (PD2, Esteban 2026-11-13): <= 1 df -> (c); > 1 df -> (a)",
                        "source": "figures/pd2_column_test.json; scripts/column_test.py (post hoc)"}}
    out["V2b"]["classification"] = {
        "class": "(a)",
        "text": _v2b_text(val["V2b"], pulse),
        "pd2_decision": ("Esteban, PD2 2026-11-13: verdict stays FAIL, class (a): the pre-registered "
                         "front-arrival metric is not defined under spatially uniform forcing. Gill (2016) "
                         "sec. 10.10, pp. 399-401, eq. 10.10.21 (thesis library Gill_2016 .md l.15449-15486): "
                         "with a uniform wind of finite duration the coastal response is 1 - cos(pi t), "
                         "independent of x, until the signal from the end of the coast arrives. Yao (1986), "
                         "doi:10.1080/07055900.1986.9649249 (Yao_1986 .md l.553-558): the forced response has "
                         "about half the free-wave alongshore phase gradient, so a free wave cannot be told "
                         "from the forced response by phase alone."),
        "note_on_lag_correlation": ("The earlier post hoc lag correlation (0.64 m/s) spans the forced "
                                    "period and is therefore expected to be biased high (Yao 1986)."),
        "evidence": "figures/pd2_v2_diagnostics.json (post hoc)", "post_hoc_evidence": True}
    out["V2b"]["descriptive_free_wave"] = {
        **{k: fw[k] for k in ("extremum", "arrival_h", "n_stations", "t_min_h", "band_hz", "c_m_s",
                              "c_stderr_m_s", "reference") if k in fw},
        "label": ("POST HOC, DESCRIPTIVE, not a validation (PD2 2026-11-13); method of Antenucci et al. "
                  "2000, Hodges et al. 2000 (l.798-800), Beletsky et al. 1997 (l.424-428); pycnocline "
                  "displacement at k* as a mode-1 proxy, not a modal projection"),
        "rules": ("Kelvin/subinertial band of PD0 section 3; record cut to t >= 6 h (pulse centre 2.5 h "
                  "+ 3 sigma) before filtering; first crest/trough with magnitude >= 1/2 of the station "
                  "maximum; least-squares slope with standard error"),
        "source": "figures/pd2_v2b_free_wave.json; analysis/pd2_v2_diagnostics.py --free-wave"}
    out["V3"]["classification"] = "not needed (pass)"
    out["V4"]["classification"] = {
        "class": "(c)",
        "text": _v4_text(val["V4"]["V4a"], d4),
        "pd2_decision": ("limitation declared (PD0 4 rule; Esteban, PD2 2026-11-13); post hoc descriptive "
                         "per strait in figures/pd2_v4_descriptive.json"),
        "evidence": "figures/pd2_v4_diagnostics.json, figures/pd2_v4_descriptive.json (post hoc)",
        "post_hoc_evidence": True}
    return out


# ---- CLI -----------------------------------------------------------------------------------------

def evaluate(f1=None, v3=None, v4=None, previous: dict | None = None, out: Path | None = None) -> dict:
    """Evaluate the criteria whose runs are given. ``previous`` (--resume): criteria already in it
    are kept and not recomputed. With ``out``, the JSON is rewritten after each criterion, so a
    failure later does not lose the hours of reading before it."""
    th = load_theory()
    res = dict(previous or {})
    res.update(thresholds=THRESHOLDS, pd0="notes/pd0_preregistro.md sections 3-4, Enmiendas 1-2")
    res.setdefault("runs", {})
    res["runs"].update({k: str(v) for k, v in (("F1", f1), ("V3", v3), ("V4", v4)) if v})

    def save():
        if out is not None:
            Path(out).parent.mkdir(parents=True, exist_ok=True)
            Path(out).write_text(json.dumps(res, indent=1, default=jsonable))

    ds_f1 = open_run(f1) if f1 and not ({"V2a", "V2b", "V4"} <= set(res)) else None
    ds_v3 = open_run(v3) if v3 and "V3" not in res else None
    ds_v4 = open_run(v4) if v4 and "V4" not in res else None
    if ds_f1 is not None and not ({"V2a", "V2b"} <= set(res)):
        res.update(v2_from_run(ds_f1, th))
        save()
    if ds_v3 is not None:
        res["V3"] = v3_from_run(ds_v3)
        save()
    if ds_f1 is not None and ds_v4 is not None:
        res["V4"] = v4_from_runs(ds_f1, ds_v4)
        save()
    done = [k for k in ("V2a", "V2b", "V3", "V4") if k in res]
    for k in done:
        res[k].setdefault("classification", "pending")     # (a)-(d), filled when a criterion fails
    failed = [k for k in done if not res[k]["pass"]]
    missing = [k for k in ("V2a", "V2b", "V3", "V4") if k not in res]
    # undecided (None) while a criterion is not evaluated, unless one already failed
    res["summary"] = {"evaluated": done, "failed": failed,
                      "go_to_PD2": True if failed else (None if missing else False),
                      "not_evaluated": missing}
    save()
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("evaluate", help="V2 (F1), V3 (tide-only run), V4 (F1 and V4)")
    e.add_argument("--f1", type=Path)
    e.add_argument("--v3", type=Path)
    e.add_argument("--v4", type=Path)
    e.add_argument("--out", type=Path, default=PKG / "figures" / "validation.json")
    e.add_argument("--resume", action="store_true",
                   help="keep the criteria already in --out and evaluate only the missing ones")
    c = sub.add_parser("ice-check", help="mandatory ice-OB distance check (pd1_seaice.md 8)")
    c.add_argument("runs", nargs="+", type=Path)
    c.add_argument("--stamps", type=Path, default=None)
    m = sub.add_parser("pd2-merge", help="rebuild the PD2 blocks of validation.json (after evaluate)")
    m.add_argument("--validation", type=Path, default=PKG / "figures" / "validation.json")
    m.add_argument("--figures", type=Path, default=PKG / "figures")
    a = ap.parse_args(argv)
    if a.cmd == "pd2-merge":
        res = pd2_merge(json.loads(a.validation.read_text()), a.figures)
        a.validation.write_text(json.dumps(res, indent=1, default=jsonable))
        return 0
    if a.cmd == "ice-check":
        bad = 0
        for r in a.runs:
            res = ice_check.check_and_stamp(r, a.stamps)
            print(f"{r}: min distance {res['min_cells']} cells (required {res['required_cells']}), "
                  f"{res['records']} records: {'PASS' if res['passes'] else 'FAIL'}")
            bad += not res["passes"]
        return 1 if bad else 0
    prev = json.loads(a.out.read_text()) if a.resume and a.out.exists() else None
    res = evaluate(a.f1, a.v3, a.v4, previous=prev, out=a.out)
    print(json.dumps(res["summary"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
