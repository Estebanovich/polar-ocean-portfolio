"""POST HOC diagnostics of the V2 failures of F1 (Task 6 part B, written after V2a and V2b failed).

Nothing here is a criterion and nothing changes a verdict of validate.py: these numbers back the
(a)-(d) classification of the failures in notes/pd2_validacion.md.

1. V2a: map of the near-inertial peak frequency of the surface HKE (k = 0, T points, every
   ``step`` cells) over the whole wet domain, with the same spectral estimate as V2a
   (spectra.psd, one Hann window over 0-T_fin) and the rotary split (spectra.rotary) at the shelf
   station: is the 1.14 f peak local to the station or domain wide, and is it clockwise?
   The bounded-basin reference is the transverse Poincare frequency
   omega_n = (f^2 + c_n^2 (pi/W)^2)^(1/2) [rad s^-1] (Gill 1982 sec. 10.3, eq. 10.3.5, with c_n for
   (gH)^(1/2); thesis library Gill_2016 .md l. 14593-14600).
2. V2b: pycnocline displacement at the V2b stations: time of the forced response (wind pulse,
   centre 2.5 h, sigma 1 h; PD0 7) against the half-maximum arrival time, and the lag of the
   maximum cross-correlation of the Kelvin-band (PD0 3 band) displacement between successive
   stations, giving a descriptive along-coast speed [m s^-1].

Usage: pd2_v2_diagnostics.py --f1 RUN --validation JSON [--out figures/pd2_v2_diagnostics.json]
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

from exp4_arctic_bayiw.analysis import spectra, validate as va  # noqa: E402

PKG = Path(__file__).resolve().parents[1]


def peak_map(U, V, t, f, band=(0.85, 1.2)):
    """Near-inertial peak frequency / f of P_u + P_v at every column of U, V (nt, n) [-]."""
    dt = float(t[1] - t[0])
    fr, Pu, df = spectra.psd(U[:, 0], dt)
    fh = f / (2 * np.pi)
    sel = (fr >= band[0] * fh) & (fr < band[1] * fh)
    Uw = (U - U.mean(0)) * np.hanning(U.shape[0])[:, None]
    Vw = (V - V.mean(0)) * np.hanning(V.shape[0])[:, None]
    P = np.abs(np.fft.rfft(Uw, axis=0)) ** 2 + np.abs(np.fft.rfft(Vw, axis=0)) ** 2
    k = np.flatnonzero(sel)[np.argmax(P[sel], axis=0)]
    return fr[k] / fh, df / fh


def band_pass(x, dt, lo, hi):
    """FFT band-pass of x (nt, ...) along time between lo and hi [Hz] (no window)."""
    X = np.fft.rfft(x - x.mean(0), axis=0)
    fr = np.fft.rfftfreq(x.shape[0], dt)
    X[(fr < lo) | (fr >= hi)] = 0
    return np.fft.irfft(X, n=x.shape[0], axis=0)


def lag_speed(z, t, s, max_lag_s=48 * 3600.0):
    """Descriptive speed from the lag of maximum cross-correlation between station 0 and each
    other station (least-squares slope of s on lag) [m s^-1]."""
    dt = float(t[1] - t[0])
    lags = []
    for q in range(z.shape[1]):
        a, b = z[:, 0], z[:, q]
        n = int(max_lag_s // dt)
        cc = [np.dot(a[:len(a) - L], b[L:]) if L >= 0 else np.dot(a[-L:], b[:len(b) + L])
              for L in range(-n, n + 1)]
        lags.append((int(np.argmax(cc)) - n) * dt)
    lags = np.array(lags)
    c = float(np.polyfit(lags, s, 1)[0]) if np.ptp(lags) > 0 else float("inf")
    return c, lags


def _first_extremum(x, t, t_min, kind, frac: float = 0.5):
    """Time [s] of the first local maximum with x >= frac max|x| ('crest') or minimum with
    x <= -frac max|x| ('trough') at t >= t_min; None if there is none. The amplitude floor keeps
    filter ripples from counting as the crest (fixed on the synthetic test, before F1)."""
    thr = frac * np.abs(x).max()
    for n in range(1, x.size - 1):
        if t[n] < t_min:
            continue
        if kind == "crest" and x[n] >= thr and x[n - 1] < x[n] >= x[n + 1]:
            return float(t[n])
        if kind == "trough" and x[n] <= -thr and x[n - 1] > x[n] <= x[n + 1]:
            return float(t[n])
    return None


def free_wave_speed(z, t, s, band, t_min: float = 6 * 3600.0, cut_before_filter: bool = True,
                    frac: float = 0.5) -> dict:
    """PD2 V2b closure (Esteban, 2026-11-13): POST HOC and DESCRIPTIVE free-wave speed, the
    literature method (Antenucci et al. 2000; Hodges et al. 2000, filtered isotherm displacement;
    Beletsky et al. 1997, tracing the displacement maximum and minimum along the coast).

    The signal is the pycnocline displacement at the level k* of maximum background N^2, used as a
    proxy for the mode-1 displacement: it is NOT a modal projection.

    Rules fixed BEFORE computing (2026-11-13): z (nt, ns) = pycnocline displacement [m] at the
    V2b south-coast stations ordered eastward, band-passed in the Kelvin/subinertial band of the
    PD0 section 3 (band = (lo, hi) [Hz]) after cutting the record to t >= t_min = 6 h (pulse centre 2.5 h + 3 sigma,
    sigma = 1 h); the extremum kind is the first extremum (crest or trough) after t_min at the
    first (westmost) station, and its first occurrence after t_min is tracked at every station
    (an extremum counts if its magnitude is >= 1/2 of the station's maximum |filtered z|);
    c = least-squares slope of s on arrival time [m s^-1] with its standard error (np.polyfit
    covariance). Not a validation: V2b stays FAIL, class (a)."""
    dt = float(t[1] - t[0])
    # the record is cut at t_min BEFORE filtering: band-passing the whole record spreads the
    # forced pulse past t_min (filter ringing of a spatially uniform signal; found on the
    # synthetic test before any F1 value was computed)
    # (cut_before_filter=False and frac=0 exist only to test why each rule is needed)
    t = np.asarray(t, float)
    if cut_before_filter:
        keep = t >= t_min
        t = t[keep]
        zk = band_pass(np.asarray(z, float)[keep], dt, *band)
    else:
        zk = band_pass(np.asarray(z, float), dt, *band)
    tc = _first_extremum(zk[:, 0], t, t_min, "crest", frac)
    tt = _first_extremum(zk[:, 0], t, t_min, "trough", frac)
    kind = "crest" if (tt is None or (tc is not None and tc <= tt)) else "trough"
    arr = [_first_extremum(zk[:, q], t, t_min, kind, frac) for q in range(zk.shape[1])]
    ok = [q for q, a in enumerate(arr) if a is not None]
    res = {"extremum": kind, "arrival_h": [None if a is None else a / 3600 for a in arr],
           "n_stations": len(ok), "t_min_h": t_min / 3600, "band_hz": list(band),
           "post_hoc": "descriptive, not a validation (PD2 2026-11-13)"}
    if len(ok) < 4:
        return {**res, "c_m_s": None, "c_stderr_m_s": None}
    p, cov = np.polyfit(np.array([arr[q] for q in ok]), np.asarray(s, float)[ok], 1, cov=True)
    return {**res, "c_m_s": float(p[0]), "c_stderr_m_s": float(np.sqrt(cov[0, 0]))}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--f1", type=Path, required=True)
    ap.add_argument("--validation", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=PKG / "figures" / "pd2_v2_diagnostics.json")
    ap.add_argument("--step", type=int, default=4)
    ap.add_argument("--free-wave", action="store_true",
                    help="only the PD2 V2b free-wave speed -> figures/pd2_v2b_free_wave.json")
    a = ap.parse_args(argv)
    if a.free_wave:
        return free_wave_main(a)
    val = json.loads(a.validation.read_text())
    th = va.load_theory()
    f = th["cases"]["vmf_summer"]["f_s-1"]
    ds = va.open_run(a.f1)
    t = ds["time"].values.astype(float)
    wet = ds["hFacC"].values[0] > 0
    st = val["V2b"]["stations"]["cells"]
    k = int(val["V2b"]["pycnocline_level"])
    sl = (slice(None, None, a.step), slice(None, None, a.step))
    nt = t.size
    U = np.empty((nt,) + wet[sl].shape); V = np.empty_like(U)
    jj = np.array([c[0] for c in st]); ii = np.array([c[1] for c in st])
    R = np.empty((nt, len(st)))
    for n in range(nt):
        u = ds["U"].isel(time=n, k=0).values
        v = ds["V"].isel(time=n, k=0).values
        uc = 0.5 * (u + np.roll(u, -1, axis=1))                  # T point (last column unused)
        vc = 0.5 * (v + np.roll(v, -1, axis=0))
        U[n], V[n] = uc[sl], vc[sl]
        R[n] = ds["RHO"].isel(time=n, k=k).values[jj, ii]
        va._progress("pd2 diagnostics", n, nt)
    w = wet[sl]
    ratio, df_f = peak_map(U[:, w], V[:, w], t, f)
    js, is_ = val["V2a"]["station_ji"]
    fr, cw, ccw, df = spectra.rotary(U[:, js // a.step, is_ // a.step], V[:, js // a.step, is_ // a.step],
                                     float(t[1] - t[0]))
    fh = f / (2 * np.pi)
    sel = (fr >= 0.85 * fh) & (fr < 1.2 * fh)
    from exp4_arctic_bayiw.analysis import budget
    rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, ds["hFacC"].values > 0)
    z = va.displacement(R, rho_b, ds["RC"].values, k)
    s = np.asarray(val["V2b"]["s_m"], float)
    bands = va.bands_hz(th)
    zk = band_pass(z, float(t[1] - t[0]), *bands["kelvin_subinertial"])
    c_lag, lags = lag_speed(zk, t, s)
    t_absmax = t[np.argmax(np.abs(z), axis=0)]
    pulse = th.get("wind", {})
    out = {
        "post_hoc": "diagnostics written after V2a and V2b failed; not criteria",
        "V2a_peak_map": {"step_cells": a.step, "n_columns": int(w.sum()),
                         "peak_over_f_percentiles": {str(p): float(np.percentile(ratio, p))
                                                     for p in (5, 25, 50, 75, 95)},
                         "fraction_within_1df_of_f": float(np.mean(np.abs(ratio - 1) <= df_f)),
                         "df_over_f": float(df_f)},
        "V2a_station_rotary": {"cw_near_inertial": float(cw[sel].sum() * df),
                               "ccw_near_inertial": float(ccw[sel].sum() * df),
                               "cw_peak_over_f": float(fr[sel][np.argmax(cw[sel])] / fh)},
        "V2b": {"arrival_h": [x / 3600 for x in val["V2b"]["arrival_s"]],
                "time_of_max_abs_zeta_h": (t_absmax / 3600).tolist(),
                "zeta_max_abs_m": np.abs(z).max(0).tolist(),
                "kelvin_band_lag_h": (lags / 3600).tolist(),
                "kelvin_band_lag_speed_m_s": c_lag, "wind": pulse},
    }
    a.out.write_text(json.dumps(out, indent=1, default=va.jsonable))
    print(json.dumps({k: v for k, v in out.items() if k != "V2b"}, indent=1, default=va.jsonable))
    print(json.dumps(out["V2b"], default=va.jsonable))
    return 0


def free_wave_main(a) -> int:
    from exp4_arctic_bayiw.analysis import budget
    val = json.loads(a.validation.read_text())
    th = va.load_theory()
    ds = va.open_run(a.f1)
    t = ds["time"].values.astype(float)
    k = int(val["V2b"]["pycnocline_level"])
    R = va.level_values(ds, "RHO", k, val["V2b"]["stations"]["cells"])
    rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, ds["hFacC"].values > 0)
    z = va.displacement(R, rho_b, ds["RC"].values, k)
    r = free_wave_speed(z, t, np.asarray(val["V2b"]["s_m"], float), va.bands_hz(th)["kelvin_subinertial"])
    r["reference"] = {"V3_c_m_s": val["V3"]["speed"]["c_m_s"],
                      "theory_c1_m_s": {"H_mean": th["cases"]["vmf_summer_Hmean"]["c_n"][0],
                                        "H_max": th["cases"]["vmf_summer"]["c_n"][0]}}
    out = a.out.parent / "pd2_v2b_free_wave.json"
    out.write_text(json.dumps(r, indent=1, default=va.jsonable))
    print(json.dumps(r, indent=1, default=va.jsonable))
    return 0


if __name__ == "__main__":
    sys.exit(main())
