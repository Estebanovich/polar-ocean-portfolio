"""Spectra for Exp. 4: one-sided and rotary periodograms with a single Hann window, band energy.

Conventions
-----------
* ``dt`` in seconds gives frequencies in Hz (cycles s^-1); any other unit of ``dt`` gives
  cycles per that unit. ``df = 1 / (n dt)`` is the Rayleigh resolution of the record.
* The mean is removed and one Hann window covers the whole record (no segment averaging, as the
  thesis periodograms, v1_thesis.psd_hann). The normalisation is the window-corrected
  periodogram, P = |FFT(x w)|^2 / (fs sum w^2) (Bendat & Piersol 2010 sec. 11.5 [verificar]), so
  that sum(P) df = sum((x w)^2) / sum(w^2): the windowed variance (Parseval). Units: x^2 / Hz.
* Rotary spectra (Gonella 1972 [verificar]) use w = u + i v: positive frequencies are
  counter-clockwise, negative ones clockwise. In the northern hemisphere inertial motion
  (u = cos f t, v = -sin f t) is clockwise.
"""
from __future__ import annotations

import numpy as np


def _window(n: int) -> np.ndarray:
    return np.hanning(n)


def psd(x, dt: float):
    """One-sided periodogram of the real series x. Returns (freq >= 0, P [x^2/Hz], df)."""
    x = np.asarray(x, float)
    x = x - x.mean()
    n = x.size
    w = _window(n)
    fs = 1.0 / dt
    P = np.abs(np.fft.rfft(x * w)) ** 2 / (fs * np.sum(w ** 2))
    P[1:-1 if n % 2 == 0 else None] *= 2.0             # fold negative frequencies (not DC, Nyquist)
    return np.fft.rfftfreq(n, d=dt), P, 1.0 / (n * dt)


def rotary(u, v, dt: float):
    """Rotary periodogram of (u, v). Returns (freq >= 0, P_cw, P_ccw, df); both P in
    (u^2 + v^2)/Hz with sum(P_cw + P_ccw) df = windowed variance of u plus that of v. The DC
    (and, for even n, the Nyquist) bin is split evenly between the two rotations."""
    u = np.asarray(u, float) - np.mean(u)
    v = np.asarray(v, float) - np.mean(v)
    n = u.size
    w = _window(n)
    fs = 1.0 / dt
    Pfull = np.abs(np.fft.fft((u + 1j * v) * w)) ** 2 / (fs * np.sum(w ** 2))
    f = np.fft.fftfreq(n, d=dt)
    nf = n // 2 + 1
    freq = np.abs(f[:nf])
    ccw = Pfull[:nf].copy()
    cw = Pfull[(-np.arange(nf)) % n].copy()             # P(-f)
    cw[0] = ccw[0] = 0.5 * Pfull[0]
    if n % 2 == 0:                                      # Nyquist bin -fs/2, shared
        cw[-1] = ccw[-1] = 0.5 * Pfull[n // 2]
    return freq, cw, ccw, 1.0 / (n * dt)


def check_bands(bands: dict, df: float) -> None:
    """ValueError if two bands overlap or a band is narrower than 2 df (PD0 section 3)."""
    items = sorted(bands.items(), key=lambda kv: kv[1][0])
    for name, (lo, hi) in items:
        if not hi > lo:
            raise ValueError(f"band {name}: upper edge {hi} <= lower edge {lo}")
        if hi - lo < 2.0 * df:
            raise ValueError(f"band {name} ({lo:.3g}-{hi:.3g}) is narrower than 2 df = {2 * df:.3g}")
    for (n1, (_, h1)), (n2, (l2, _)) in zip(items[:-1], items[1:]):
        if l2 < h1:
            raise ValueError(f"bands {n1} and {n2} overlap ({l2:.3g} < {h1:.3g})")


def band_energy(freq, P, bands: dict, df: float) -> dict:
    """{band: sum(P df) over lo <= freq < hi} [x^2]; bands = {name: (lo, hi)} in the units of freq.
    Half-open bins so adjacent bands that share an edge do not count a bin twice."""
    check_bands(bands, df)
    freq, P = np.asarray(freq, float), np.asarray(P, float)
    return {k: float(P[(freq >= lo) & (freq < hi)].sum() * df) for k, (lo, hi) in bands.items()}
