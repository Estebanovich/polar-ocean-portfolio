#!/opt/anaconda3/envs/pyocean/bin/python
r"""Four DESCRIPTIVE GIF animations of Exp. 4 (Task 9): plan views of the mode-1 isopycnal
displacement a_1 (raw) and of its inertially filtered part a_1' for two run pairs.

Usage::

    animate.py fjord-ice [--runs F3 F2]   -> figures/anim_fjord_ice_{raw,filt}.gif
    animate.py bay-fjord [--runs R1 F1]   -> figures/anim_bay_vs_fjord_{raw,filt}.gif
        common options: [--out-dir DIR] [--t-end-h 96] [--step-s S] [--stride N] [--dpi D]
                        [--max-mb 5] [--frame-png PREFIX] [--frame-h H ...] [--save-frames NPZ]
                        [--fields raw filt] [--check-scales] [--runs-root EXP4_RUNS]

One pass over the run records writes both GIFs of a pair. No number of these animations enters
the README claims; no hypothesis or criterion uses them.

Field (choice and justification). a_1(x, y, t) = the mode-1 amplitude of the isopycnal
displacement zeta = -rho'/(d rho_b/dz) [m, + up], from each 3-D density snapshot:

    a_1 = sum_k N_b^2 zeta W_1 dz / sum_k N_b^2 W_1^2 dz = (g/rho0) sum_k rho' W_1 dz / sum_k N_b^2 W_1^2 dz

(N_b^2-weighted projection, the Sturm-Liouville weight of the W modes, Gill 1982 sec. 6.11
[verificar]; N_b^2 zeta = (g/rho0) rho'). rho' = rho - rho_b with rho_b the FIXED reference
profile of the energy budget (PD0 Enmienda 1: budget.background_profile of the t = 0 density);
W_1 = theory.modes.vertical_modes of rho_b over the wet thickness of each column (dz = drF hFacC,
the same dz in both sums), scaled to max |W_1| = 1, so a_1 is the displacement where the mode-1
displacement peaks. Why not the single-level displacement at the level k* of maximum N_b^2
(validate.displacement, used by V2b/V3): it was tried first and is unusable in winter. The winter
profile is a staircase with N_b^2 ~ 0 between steps; its k* is the top cell (z = -0.75 m), whose
0.03 kg m^-3 surface step the pulse mixes away within 3 h, so F2 and F3 showed a uniform +1.6 to
+1.8 m "displacement" (median over the wet cells) that is mixing, not waves; deeper levels either
have d rho_b/dz ~ 0 or drift by metres over 96 h. The projection needs no division by
d rho_b/dz and weights the surface (W_1 = 0 at the lid) little. Caveats: W_1 is that of a
flat-bottom column of the local depth (topographic coupling ignored); rho' changes that are not
vertical advection (KPP mixing, the profile drift, ice-ocean fluxes in F2) still project on mode
1, as same-sign features that persist from frame to frame (coastal fringes, a shelf strip, the
open-boundary sponges; under the F2 ice a slow monotonic negative a_1). The raw GIFs say so in
their footer.

Filtered field (Esteban, fix round 1). a_1' = a_1 minus its centred running mean over exactly one
inertial period T_f = 2 pi/|f| of each run (inertial_filter; 12.23 h fjord, 25.04 h thesis bay).
The boxcar removes a constant and a linear trend exactly and the slowly varying same-sign part
(wind mixing, drift) mostly, and keeps oscillations: those near f are kept whole, slower ones
(sub-inertial Kelvin and seiche signals) are attenuated. Frames within T_f/2 of the record ends
have no full window and are dropped, not padded: the record read extends to t_end + max T_f/2 so
the filtered GIFs reach t_end, and they start at max T_f/2 over the pair.

Animations.

* ``fjord-ice`` (Anim A): open water against landfast ice. The brief named "F1 vs F1c"; in
  cases.yaml F1c is the thesis closure (constant mixing), not an ice case, and the landfast-ice
  case is F2 (winter profile). The comparison that changes only the ice is F2 against F3 (same
  winter profile, grid and wind; the pair of H4b and Figure 4), so the default is F3 | F2. Shared
  colour scale (same profile, grid and pulse). F2 shows the outline of the initial ice cover
  (AREA.bin >= 0.5, the actual model input), with the strait sections drawn above it.
* ``bay-fjord`` (Anim B): the thesis bay R1 against Van Mijenfjorden F1, each on its own axes.
  Separate fixed scales, stated in the title: different profiles, depths, f (6.97e-5 against
  1.43e-4 s^-1) and wind footprints (a Gaussian patch of radius ~300 km, 60 cells, against the
  fjord-wide pulse), so the amplitudes are not physically comparable (decided from the
  configurations before any frame was drawn). The bay view (bay_crop) is the bay (bay_mask;
  60.6 km x 121 km in figures/theory.json) widened by one bay width on each side and extended one
  bay length onto the shelf; the thesis grid stretches to ~2800 km. The R1 limit is set from the
  bay pixels, so the subject does not saturate. Footer: the bridge is not resolution-converged
  and the fixed reference drifts.

Conventions (as exp1_fjord_plumes/analysis/animate_c.py and animate_b.py): land light grey; water
columns shallower than MIN_DEPTH_M = 10 m (or < 3 wet levels) a lighter grey (no a_1); FIXED
symmetric colour limits for all frames of a GIF, L = the ceiling on the NICE ladder (1, 1.5, 2,
3, 4, 5, 6, 8 x 10^n) of the 99.5th percentile of |field| over the limit region (every wet pixel;
the bay for R1) of every frame, pooled over the panels when shared; ``--check-scales`` reports
the saturated fractions (shown pixels, limit region) and the share of the saturated pixels in
columns shallower than 30 m. Strait sections (figures/bathy_report.json) on the fjord panels; a
clock in hours after the pulse onset and in inertial periods of each run; a strip with the peak
wind stress of the forcing files (load.forcing_stress; under the landfast ice of F2 the ocean
does not receive it), sampled every TAU_STEP_S = 600 s on its own times (the 3600 s frame step
flattened the 1 h ramp to 0.085 of the 0.1 N m^-2 peak), with a marker interpolated to the frame.

Time. Pulse onset t_on = centre - 2.5 sigma of the half-Gaussian ramp (gendata/forcing.pulse and
make_case.thesis_wind_patch: Gaussian rise, zero after the centre), clipped at 0 (fjord and
thesis: centre 2.5 h, sigma 1 h, so t_on = 0 = the model start, where tau is 4 % of its peak).
Raw frames from t_on to t_on + 96 h by default (7.8 fjord inertial periods; 3.8 of the thesis f).
Anim A reads every 1200 s snapshot (before the GIF stride); Anim B subsamples both runs to the
common step 3600 s (bay 900 s, fjord 1200 s snapshots). GIF: Pillow, the stride doubles until
the file fits --max-mb (animate_b.write_within_budget); exit code 1 if one is still over budget.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

from exp1_fjord_plumes.analysis.animate_b import write_within_budget  # noqa: E402
from pfcommon import style  # noqa: E402

PKG = REPO / "exp4_arctic_bayiw"
FIG_DIR = PKG / "figures"
T_END_H = 96.0            # frames from the onset to onset + 96 h
PCTL = 99.5               # colour limit: NICE ceiling of this percentile of |a_1|
ONSET_SIGMAS = 2.5
MIN_LEVELS = 3            # columns with fewer wet levels have no mode-1 amplitude
MIN_DEPTH_M = 10.0        # nor columns shallower than this (only the pulse-mixed surface layer)
SHALLOW_REPORT_M = 30.0   # check_scales: share of the saturated pixels in columns shallower than this
TAU_STEP_S = 600.0        # the wind strip samples tau on its own times, not on the frame times
FPS = 10
FOOTER_WIDTH = 190        # characters per footer line
LAND_GREY, SHALLOW_GREY = "0.80", "0.92"
A1_DEF = ("Descriptive only. $a_1$ = N$^2$-weighted projection of the isopycnal displacement on the local mode-1 "
          "structure (max = 1), fixed t = 0 reference profile (energy-budget $\\rho_b$).")
RAW_CAVEAT = ("Same-sign features that persist from frame to frame (coastal fringes, shelf strip, "
              "sponge) may be wind mixing or reference drift, not internal waves")
NO_ICE_WIND = "; there is no wind under the ice in F2."
BRIDGE = ("The bay-fjord bridge is not resolution-converged, and the fixed t = 0 reference profile "
          "drifts during the runs.")


FILTER_NOTE = ("The 1-T$_f$ running-mean filter attenuates sub-inertial signal (Kelvin waves, slow "
               "seiches) and has ripple above f; $a_1'$ shows mainly near- and super-inertial "
               "oscillations.")
SHALLOW_NOTE_FACTOR = 2.0  # footer note when shallow columns hold >= this x their area share


def saturation_note(scales: dict, factor: float = SHALLOW_NOTE_FACTOR) -> str:
    """Footer sentence, computed from check_scales at render time: the runs whose saturated pixels
    lie in columns shallower than SHALLOW_REPORT_M at least ``factor`` times as often as those
    columns occur among the shown pixels (saturation concentrated there); "" if none."""
    key = f"share_of_saturated_in_columns_shallower_than_{SHALLOW_REPORT_M:g}m"
    area = f"fraction_of_shown_pixels_in_columns_shallower_than_{SHALLOW_REPORT_M:g}m"
    hits = [(n, v[key], v[area]) for n, v in scales.items()
            if isinstance(v, dict) and v.get(key) is not None and v.get(area)
            and v[key] >= factor * v[area]]
    return " ".join(f"{n}: saturation concentrates in columns shallower than {SHALLOW_REPORT_M:g} m "
                    f"({100 * sh:.0f} % of its saturated pixels, {100 * ar:.0f} % of its area)."
                    for n, sh, ar in hits)


def filt_definition(tf_h: dict) -> str:
    """Footer line that defines the filtered field a_1' ({run: T_f [h]})."""
    per = ", ".join(f"{n} {v:.2f} h" for n, v in tf_h.items())
    return ("$a_1'$ = $a_1$ minus its centred running mean over exactly one inertial period "
            f"T$_f$ = 2$\\pi$/|f| of each run ({per}); frames within T$_f$/2 of the record ends are "
            "dropped (no padding).")


# ---- pure pieces --------------------------------------------------------------------------------

def mode1_shape(n2, dz):
    """Mode-1 vertical-displacement structure W_1 (rigid lid, flat bottom; theory.modes.
    vertical_modes, Gill 1982 sec. 6.11 [verificar]) of one column, scaled to max |W_1| = 1
    (positive lobe). n2 (kb,) N_b^2 of the wet levels [s^-2], dz (kb,) their wet thicknesses
    drF hFacC [m]: the column spans 0 .. -H, H = sum dz, with the level nodes at the centres of
    the wet thicknesses (a partial bottom cell moves its node and the bottom); W = 0 at both ends.
    None if kb < MIN_LEVELS or H < MIN_DEPTH_M: the winter surface step (top ~5-7 m) is mixed
    away by the pulse within 3 h, and columns of only that layer gave the largest |a_1| (up to
    19 m, one-cell coastal fringes) with no wave content."""
    from exp4_arctic_bayiw.theory import modes
    dz = np.asarray(dz, float)
    H = float(dz.sum())
    if dz.size < MIN_LEVELS or H < MIN_DEPTH_M:
        return None
    n2 = np.maximum(np.asarray(n2, float), modes.N2_FLOOR)
    zc = -(np.cumsum(dz) - dz / 2.0)
    _, phi = modes.vertical_modes(np.r_[0.0, zc, -H], np.r_[n2[0], n2, n2[-1]], H, n_modes=1)
    w = phi[0, 1:-1]
    return w / np.max(np.abs(w))


def mode1_weights(rho_b, RC, drF, hfac, rho0: float, g: float = 9.81) -> tuple:
    """(W (nz, ny, nx), den (ny, nx)) of mode1_amplitude: W_1 of every column (columns with the
    same wet-thickness profile share one eigen-solve) and den = sum_k N_b^2 W_1^2 drF hFacC
    (NaN where mode1_shape is None). Depends on the grid and rho_b only (computed once per run)."""
    from exp4_arctic_bayiw.analysis import budget
    from exp4_arctic_bayiw.theory import modes
    hfac = np.asarray(hfac, float)
    nz, ny, nx = hfac.shape
    n2 = np.maximum(budget.n2_background(rho_b, RC, rho0, g), modes.N2_FLOOR)
    dz = np.asarray(drF, float)[:, None, None] * hfac
    cols = np.round(dz.reshape(nz, -1).T, 6)
    uniq, inv = np.unique(cols, axis=0, return_inverse=True)
    W = np.zeros((nz, ny * nx))
    den = np.full(ny * nx, np.nan)
    for u, dzc in enumerate(uniq):
        kb = int(np.sum(dzc > 0))
        if kb == 0 or not np.all(dzc[:kb] > 0):
            continue
        w = mode1_shape(n2[:kb], dzc[:kb])
        if w is None:
            continue
        m = np.ravel(inv) == u
        W[:kb, m] = w[:, None]
        den[m] = np.sum(n2[:kb] * w * w * dzc[:kb])
    return W.reshape(nz, ny, nx), den.reshape(ny, nx)


def mode1_amplitude(rho, rho_b, RC, drF, hfac, rho0: float, g: float = 9.81,
                    weights: tuple | None = None) -> np.ndarray:
    """Mode-1 displacement amplitude a(x, y) [m, + up] of the 3-D density rho (nz, ny, nx):
    the N_b^2-weighted projection of the isopycnal displacement zeta = -rho'/(d rho_b/dz) on the
    local W_1 (max |W_1| = 1, so a is the displacement where the mode-1 displacement peaks),

        a = sum_k N_b^2 zeta W_1 dz / sum_k N_b^2 W_1^2 dz
          = (g / rho0) sum_k rho' W_1 dz / sum_k N_b^2 W_1^2 dz,

    using N_b^2 zeta = (g / rho0) rho' (N_b^2 = -(g/rho0) d rho_b/dz), so no division by
    d rho_b/dz is needed where the profile is unstratified (N_b^2 = 0 layers of the winter
    profile). rho' = rho - rho_b with the fixed reference rho_b (Enmienda 1); dz = drF hFacC in the
    numerator, in the denominator and in the depth W_1 is solved on (mode1_weights). NaN where
    mode1_shape is None (too few levels or too shallow). ``weights`` = mode1_weights(...) to
    reuse across frames."""
    rho = np.asarray(rho, float)
    hfac = np.asarray(hfac, float)
    nz = rho.shape[0]
    W, den = weights if weights is not None else mode1_weights(rho_b, RC, drF, hfac, rho0, g)
    rp = np.where(hfac > 0, rho - np.asarray(rho_b, float)[:nz, None, None], 0.0)
    dz = np.asarray(drF, float)[:, None, None] * hfac
    return g / rho0 * np.sum(rp * W * dz, axis=0) / den


def inertial_filter(x, t_s, T: float) -> tuple:
    """(t_kept, x') with x' = x minus its centred running mean over exactly one period T [s]:

        x'(t) = x(t) - (1/T) int_{t - T/2}^{t + T/2} x dt',

    the integral exact for the piecewise-linear interpolant of the uniformly sampled x (nt, ...),
    so T need not be a multiple of the step. Only the frames whose whole window lies inside the
    record are kept (t - T/2 >= t_0 and t + T/2 <= t_end; no padding). A centred boxcar removes a
    constant and a linear trend exactly and a sinusoid of period T to O((dt/T)^2)."""
    t = np.asarray(t_s, float)
    x = np.asarray(x)
    dt = float(t[1] - t[0])
    if not np.allclose(np.diff(t), dt):
        raise ValueError("inertial_filter needs uniformly spaced times")
    eps = 1e-6 * dt
    keep = (t - T / 2 >= t[0] - eps) & (t + T / 2 <= t[-1] + eps)
    if not keep.any():
        raise ValueError(f"record of {t[-1] - t[0]:.0f} s is shorter than the window {T:.0f} s")
    xd = x.astype(np.float64)
    C = np.concatenate([np.zeros((1,) + xd.shape[1:]),
                        np.cumsum(0.5 * (xd[1:] + xd[:-1]) * dt, axis=0)])

    def cint(tq):                                   # int_{t0}^{tq} of the linear interpolant
        k = int(np.clip(np.floor((tq - t[0]) / dt + 1e-9), 0, len(t) - 2))
        sx = tq - t[k]
        return C[k] + xd[k] * sx + (xd[k + 1] - xd[k]) * sx * sx / (2 * dt)

    idx = np.flatnonzero(keep)
    out = np.stack([xd[n] - (cint(t[n] + T / 2) - cint(t[n] - T / 2)) / T for n in idx])
    return t[idx], out.astype(x.dtype if np.issubdtype(x.dtype, np.floating) else np.float64)


def pulse_onset_s(center_h: float, sigma_h: float) -> float:
    """Onset of the half-Gaussian wind ramp = centre - 2.5 sigma, clipped at the model start [s]."""
    return max(0.0, (center_h - ONSET_SIGMAS * sigma_h) * 3600.0)


def clock(t_s: float, onset_s: float, f: float) -> tuple[float, float]:
    """(hours after the onset, inertial periods 2 pi/|f| after the onset)."""
    dt = float(t_s) - float(onset_s)
    return dt / 3600.0, dt * abs(f) / (2.0 * np.pi)


def clock_text(t_s: float, onset_s: float, f: float) -> str:
    h, n = clock(t_s, onset_s, f)
    return f"t = {h:4.1f} h = {n:4.2f} T_f"


def common_times(records: list, step_s: float, t0_s: float, t1_s: float) -> np.ndarray:
    """Frame times t0, t0 + step, ... <= t1 [s]; ValueError if one is missing in any record."""
    t = t0_s + step_s * np.arange(int(math.floor((t1_s - t0_s) / step_s + 1e-9)) + 1)
    for r in records:
        r = np.asarray(r, float)
        miss = [x for x in t if not np.any(np.abs(r - x) < 1e-6)]
        if miss:
            raise ValueError(f"{len(miss)} frame times are not records (first {miss[0]} s)")
    return t


NICE = (1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)


def nice_ceil(x: float) -> float:
    """Smallest m x 10^n >= x with m in NICE (a 1-2-5 ladder left the F1 scale 1.7 x too wide)."""
    e = math.floor(math.log10(x))
    for m in NICE:
        if m * 10.0 ** e >= x * (1 - 1e-12):
            return m * 10.0 ** e
    return 10.0 ** (e + 1)


def colour_limits(frames: dict, shared: bool, pctl: float = PCTL) -> dict:
    """{run: (-L, L)} fixed for all frames; L = nice_ceil(pctl-th percentile of |a_1|), pooled
    over the runs when ``shared``."""
    def lim(a):
        a = np.abs(np.asarray(a, float))
        return nice_ceil(float(np.nanpercentile(a[np.isfinite(a)], pctl)))
    if shared:
        L = lim(np.concatenate([np.ravel(v) for v in frames.values()]))
        return {r: (-L, L) for r in frames}
    return {r: (-lim(v), lim(v)) for r, v in frames.items()}


def saturated_fraction(eta, lim: float) -> float:
    """Fraction of the finite pixels with |eta| >= lim."""
    a = np.asarray(eta, float)
    a = a[np.isfinite(a)]
    return float(np.mean(np.abs(a) >= lim)) if a.size else 0.0


def edges(centres, widths) -> np.ndarray:
    """Cell edges from the centres and widths of a (possibly stretched) grid."""
    c, w = np.asarray(centres, float), np.asarray(widths, float)
    return np.r_[c[0] - w[0] / 2.0, c + w / 2.0]


# ---- figure -------------------------------------------------------------------------------------

@dataclass
class Panel:
    name: str
    title: str
    x_km: np.ndarray            # (nx + 1) edges
    y_km: np.ndarray            # (ny + 1) edges
    land: np.ndarray            # (ny, nx) bool
    shallow: np.ndarray         # (ny, nx) bool, water columns without a mode-1 amplitude
    eta: np.ndarray             # (nt, ny, nx) mode-1 amplitude [m], NaN on land / shallow
    lim: tuple                  # (-L, L)
    f: float                    # Coriolis parameter [s^-1]
    lines: list = field(default_factory=list)       # [(x_km array, y_km array, label)]
    ice: np.ndarray | None = None                   # (ny, nx) ice concentration for the outline
    xlabel: str = "x (km)"
    ylabel: str = "y (km)"
    tau: np.ndarray | None = None                   # (nt,) peak wind stress [N m^-2]
    mesh: object = None
    clock: object = None


def build_figure(panels: list, t_s, onset_s: float, suptitle: str, footer: str,
                 tau_t_s=None, cbar_label: str = "$a_1$") -> tuple:
    """Figure with one plan-view panel per run (fixed colour limits), a clock per panel and, if any
    panel has ``tau``, a wind-stress strip (tau sampled at ``tau_t_s``, default the frame times;
    the moving marker is interpolated to the frame time). Returns (fig, draw(n))."""
    style.apply()
    import textwrap
    footer = "\n".join(textwrap.fill(par, FOOTER_WIDTH) for par in footer.split("\n")) if footer else ""
    t_s = np.asarray(t_s, float)
    n = len(panels)
    shared = len({p.lim for p in panels}) == 1
    has_tau = any(p.tau is not None for p in panels)
    widths = []
    for p in panels:
        widths.append((p.x_km[-1] - p.x_km[0]) / (p.y_km[-1] - p.y_km[0]))
    fig = plt.figure(figsize=(5.2 * n + 0.8, 5.6 if has_tau else 4.6), facecolor="white")
    hr = [1.0, 0.28] if has_tau else [1.0]
    gs = fig.add_gridspec(len(hr), n, height_ratios=hr, width_ratios=widths, hspace=0.42,
                          wspace=0.28, left=0.06, right=0.90 if shared else 0.95, top=0.84,
                          bottom=(0.16 if has_tau else 0.12) + 0.03 * footer.count("\n"))
    cmap = style.CMAPS["diverging"].copy()
    cmap.set_bad((0, 0, 0, 0))
    bg = ListedColormap([(1, 1, 1, 0), SHALLOW_GREY, LAND_GREY])
    for q, p in enumerate(panels):
        ax = fig.add_subplot(gs[0, q])
        code = np.where(p.land, 2, np.where(p.shallow, 1, 0))
        ax.pcolormesh(p.x_km, p.y_km, code, cmap=bg, vmin=-0.5, vmax=2.5, shading="flat")
        p.mesh = ax.pcolormesh(p.x_km, p.y_km, np.ma.masked_invalid(p.eta[0]), cmap=cmap,
                               vmin=p.lim[0], vmax=p.lim[1], shading="flat")
        if p.ice is not None:
            xc = 0.5 * (p.x_km[1:] + p.x_km[:-1])
            yc = 0.5 * (p.y_km[1:] + p.y_km[:-1])
            ax.contour(xc, yc, p.ice, levels=[0.5], colors="#7a3fa0", linewidths=1.2, zorder=3)
            ax.plot([], [], color="#7a3fa0", lw=1.2, label="landfast ice (initial)")
            ax.legend(loc="lower right", fontsize=7, frameon=True)
        for xs, ys, lab in p.lines:          # sections above the ice outline
            ax.plot(xs, ys, color="k", lw=1.8, solid_capstyle="butt", zorder=5)
            ax.annotate(lab, (xs[-1], ys[-1]), xytext=(4, 2), textcoords="offset points",
                        fontsize=7, zorder=5)
        ax.set_aspect("equal")
        ax.set(xlim=(p.x_km[0], p.x_km[-1]), ylim=(p.y_km[0], p.y_km[-1]),
               xlabel=p.xlabel, ylabel=p.ylabel)
        ax.set_title(p.title, fontsize=9.5)
        p.clock = ax.text(0.02, 0.97, "", transform=ax.transAxes, ha="left", va="top",
                          family="monospace", fontsize=8,
                          bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6", alpha=0.85))
        if not shared:
            cb = fig.colorbar(p.mesh, ax=ax, orientation="vertical", fraction=0.04, pad=0.02,
                              extend="both")
            cb.set_label(f"{cbar_label} (m, + up)")
    if shared:
        cax = fig.add_axes([0.915, 0.30 if has_tau else 0.16, 0.012, 0.50 if has_tau else 0.64])
        cb = fig.colorbar(panels[0].mesh, cax=cax, extend="both")
        cb.set_label(f"{cbar_label} (m, + up)")
    marks = []
    tt = np.asarray(tau_t_s if tau_t_s is not None else t_s, float)
    if has_tau:
        axt = fig.add_subplot(gs[1, :])
        th = (tt - onset_s) / 3600.0
        for q, p in enumerate(panels):
            if p.tau is None:
                continue
            c = style.OKABE_ITO[q % len(style.OKABE_ITO)]
            axt.plot(th, p.tau, color=c, lw=1.2, label=p.name)
            mk, = axt.plot([], [], "o", color=c, ms=5)
            marks.append((mk, p))
        axt.set(xlim=(th[0], th[-1]), xlabel="hours after the pulse onset",
                ylabel=r"max $|\tau|$ (N m$^{-2}$)")
        axt.set_ylim(0, max(float(np.nanmax(p.tau)) for p in panels if p.tau is not None) * 1.15)
        if len(marks) > 1:
            axt.legend(loc="upper right", fontsize=7)
    fig.suptitle(suptitle, x=0.06, ha="left", fontsize=10.5)
    if footer:
        fig.text(0.06, 0.012, footer, ha="left", va="bottom", fontsize=7.5, color="0.25")

    def draw(k: int):
        for p in panels:
            p.mesh.set_array(np.ma.masked_invalid(p.eta[k]).ravel())
            p.clock.set_text(clock_text(t_s[k], onset_s, p.f))
        for mk, p in marks:
            mk.set_data([(t_s[k] - onset_s) / 3600.0], [float(np.interp(t_s[k], tt, p.tau))])

    draw(0)
    return fig, draw


def write_gif(fig, draw: Callable, n_frames: int, out: Path, max_mb: float = 5.0,
              dpi: int = 80, stride: int = 1, fps: int = FPS) -> dict:
    """Pillow GIF of frames 0, stride, 2 stride, ... (and the last); the stride doubles until the
    file fits ``max_mb`` (animate_b.write_within_budget)."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    nsel = [0]

    def write(st: int) -> float:
        sel = list(range(0, n_frames, st))
        if sel[-1] != n_frames - 1:
            sel.append(n_frames - 1)
        w = PillowWriter(fps=fps)
        with w.saving(fig, str(out), dpi=dpi):
            for k in sel:
                draw(k)
                w.grab_frame()
        nsel[0] = len(sel)
        return out.stat().st_size / 1e6

    st, mb, over = write_within_budget(write, n_frames, stride, max_mb)
    return {"gif": str(out), "mb": mb, "over_budget": over, "stride": st, "frames": nsel[0],
            "dpi": dpi, "fps": fps}


# ---- model runs ---------------------------------------------------------------------------------

def _case(run_dir: Path) -> dict:
    return json.loads((Path(run_dir) / "case.json").read_text())


def _f0(run_dir: Path) -> float:
    from exp4_arctic_bayiw.analysis import load
    act = _case(run_dir).get("actual", {})
    if act.get("f0_s-1") is not None:
        return float(act["f0_s-1"])
    return float(load._num(load._nml(Path(run_dir)).get("PARM01", {}).get("f0")))


def _onset(run_dir: Path) -> float:
    w = _case(run_dir)["params"]["wind"]
    c = w.get("center_h", w.get("center_t_h"))
    s = w.get("sigma_h", w.get("sigma_t_h"))
    return pulse_onset_s(float(c), float(s))


def is_bay_run(case: dict) -> bool:
    """The one rule for bay vs fjord runs: fjord_vmf cases carry the UTM box ``bbox_utm``."""
    return "bbox_utm" not in case["params"]


def _tau_peak(ds, t_s) -> np.ndarray:
    """max over the grid of |tau| (forcing files, components averaged to the cell centres) [N m^-2]."""
    from exp4_arctic_bayiw.analysis import load
    out = []
    for t in t_s:
        tx, ty = load.forcing_stress(ds, float(t))
        txc = 0.5 * (tx[:, :-1] + tx[:, 1:])          # (ny, nx - 1), U points -> centres
        tyc = 0.5 * (ty[:-1] + ty[1:])                # (ny - 1, nx), V points -> centres
        out.append(float(np.max(np.hypot(txc[:-1], tyc[:, :-1]))))
    return np.array(out)


def read_run(run_dir: Path) -> dict:
    """Grid, fixed reference rho_b, mode-1 weights, f and pulse onset of one run (load.run;
    LaCie read-only)."""
    from exp4_arctic_bayiw.analysis import budget, load
    ds = load.run(run_dir)
    hfac = ds["hFacC"].values
    wet3 = hfac > 0
    rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, wet3)   # Enmienda 1
    rho0, g = float(ds.attrs["rho0"]), float(ds.attrs["g"])
    W, den = mode1_weights(rho_b, ds["RC"].values, ds["drF"].values, hfac, rho0, g)
    land = ~wet3.any(axis=0)
    H = np.sum(ds["drF"].values[:, None, None] * hfac, axis=0)               # wet depth [m]
    return {"ds": ds, "rho_b": rho_b, "hfac": hfac, "weights": (W, den), "land": land,
            "shallow": ~land & np.isnan(den), "H": H, "f": _f0(run_dir),
            "onset_s": _onset(run_dir), "case": _case(run_dir), "run_dir": Path(run_dir)}


def eta_frames(r: dict, t_s) -> np.ndarray:
    """(nt, ny, nx) mode-1 displacement amplitude at the frame times [m], float32 (one 3-D RHO
    snapshot per frame, load.run)."""
    ds = r["ds"]
    tm = ds["time"].values
    rho0, g = float(ds.attrs["rho0"]), float(ds.attrs["g"])
    RC, drF = ds["RC"].values, ds["drF"].values
    out = []
    for t in t_s:
        n = int(np.argmin(np.abs(tm - t)))
        if abs(tm[n] - t) > 1e-6:
            raise ValueError(f"{r['run_dir']}: no record at t = {t} s")
        rho = ds["RHO"].isel(time=n).values
        out.append(mode1_amplitude(rho, r["rho_b"], RC, drF, r["hfac"], rho0, g,
                                   weights=r["weights"]).astype(np.float32))
    return np.stack(out)


def _sections_km() -> list:
    rep = json.loads((FIG_DIR / "bathy_report.json").read_text())
    return [(np.array([v["p0"][0], v["p1"][0]]) / 1e3, np.array([v["p0"][1], v["p1"][1]]) / 1e3, n)
            for n, v in rep["sounds"].items()]


def _ice(r: dict):
    ice = r["case"].get("actual", {}).get("ice")
    if not ice or ice.get("mode") in (None, "none"):
        return None
    ny, nx = r["land"].shape
    return np.fromfile(r["run_dir"] / "AREA.bin", ">f8").reshape(ny, nx)


def bay_mask(land) -> np.ndarray:
    """Cells of the thesis bay: the water north of the coast row, the first row (from the south)
    that is mostly land (the bay opens southward onto the shelf; thesis BayIW_Rectan geometry)."""
    land = np.asarray(land, bool)
    rows = np.flatnonzero(land.mean(axis=1) > 0.5)
    if rows.size == 0:
        raise ValueError("no coast row: not a bay grid")
    m = np.zeros_like(land)
    m[rows[0]:] = ~land[rows[0]:]
    return m


def bay_crop(land) -> tuple:
    """(j0, j1, i0, i1) of the bay view: the bay box (bay_mask) widened by one bay width on each
    side and extended one bay length south of the mouth, up to the north edge. The stretched
    outer grid (thousands of km) and most of the shelf are left out, so the bay fills the panel."""
    js, is_ = np.nonzero(bay_mask(land))
    ny, nx = np.asarray(land).shape
    jb0, ib0, ib1 = js.min(), is_.min(), is_.max() + 1
    width, length = ib1 - ib0, js.max() + 1 - jb0
    return (max(0, int(jb0 - length)), ny, max(0, int(ib0 - width)), min(nx, int(ib1 + width)))


def _crop_mask(shape, crop) -> np.ndarray:
    j0, j1, i0, i1 = crop
    m = np.zeros(shape, bool)
    m[j0:j1, i0:i1] = True
    return m


def make_panel(r: dict, eta: np.ndarray, lim, title: str, tau=None) -> Panel:
    ds = r["ds"]
    x = edges(ds["XC"].values[0], ds["dxF"].values[0]) / 1e3
    y = edges(ds["YC"].values[:, 0], ds["dyF"].values[:, 0]) / 1e3
    if not is_bay_run(r["case"]):
        return Panel(name=r["case"]["params"]["case"], title=title, x_km=x, y_km=y,
                     land=r["land"], shallow=r["shallow"], eta=eta, lim=lim, f=r["f"],
                     lines=_sections_km(), ice=_ice(r), xlabel="easting UTM 33N (km)",
                     ylabel="northing (km)", tau=tau)
    j0, j1, i0, i1 = bay_crop(r["land"])
    return Panel(name=r["case"]["params"].get("case", "R1"), title=title,
                 x_km=x[i0:i1 + 1] - x[i0], y_km=y[j0:j1 + 1] - y[j0], land=r["land"][j0:j1, i0:i1],
                 shallow=r["shallow"][j0:j1, i0:i1], eta=eta[:, j0:j1, i0:i1], lim=lim, f=r["f"],
                 xlabel="x (km, from the crop edge)", ylabel="y (km)", tau=tau)


def _limit_region(r: dict) -> np.ndarray:
    """Pixels that set the fixed colour limit: the bay itself for the thesis bay (the subject of
    the panel must not saturate), every wet pixel for the fjord."""
    return bay_mask(r["land"]) if is_bay_run(r["case"]) else ~r["land"]


def _view_mask(r: dict) -> np.ndarray:
    """Pixels shown in the panel (the bay crop for the thesis bay)."""
    if is_bay_run(r["case"]):
        return _crop_mask(r["land"].shape, bay_crop(r["land"]))
    return np.ones(r["land"].shape, bool)


def tf_s(f: float) -> float:
    """Inertial period 2 pi / |f| [s]."""
    return 2.0 * np.pi / abs(f)


def prepare(kind: str, runs: list, root: Path, t_end_h: float = T_END_H, step_s=None) -> dict:
    """Read the runs and build the raw and filtered frames and their fixed limits.

    Raw frames: onset .. onset + t_end_h. The record read extends to onset + t_end_h + max T_f/2
    (when the run has it), so the filtered frames reach t_end_h; filtered frames start at the
    largest T_f/2 of the runs (both panels need a full window). Limits per field and GIF."""
    rr = {n: read_run(root / n) for n in runs}
    onsets = {n: r["onset_s"] for n, r in rr.items()}
    if len(set(onsets.values())) != 1:
        raise ValueError(f"pulse onsets differ: {onsets}")
    t0 = next(iter(onsets.values()))
    recs = [r["ds"]["time"].values for r in rr.values()]
    if step_s is None:
        steps = [float(np.median(np.diff(t))) for t in recs]
        step_s = float(np.lcm.reduce([int(round(s)) for s in steps]))
    half = max(tf_s(r["f"]) for r in rr.values()) / 2
    n_want = math.ceil((t_end_h * 3600.0 + half) / step_s - 1e-9)
    n_have = math.floor((min(float(t[-1]) for t in recs) - t0) / step_s + 1e-9)
    t_all = common_times(recs, step_s, t0, t0 + min(n_want, n_have) * step_s)
    raw_keep = t_all <= t0 + t_end_h * 3600.0 + 1e-6
    raw, filt, tk = {}, {}, None
    for n, r in rr.items():
        e = eta_frames(r, t_all)
        e = np.where(_view_mask(r), e, np.nan).astype(np.float32)
        raw[n] = e[raw_keep]
        tkn, ef = inertial_filter(e, t_all, tf_s(r["f"]))
        filt[n] = (tkn, ef)
    # common filtered times: full window in every run, and <= t_end
    tset = None
    for tkn, _ in filt.values():
        tset = set(np.round(tkn, 3)) if tset is None else tset & set(np.round(tkn, 3))
    tk = np.array(sorted(x for x in tset if x <= t0 + t_end_h * 3600.0 + 1e-3))
    filt = {n: ef[np.isin(np.round(tkn, 3), np.round(tk, 3))] for n, (tkn, ef) in filt.items()}
    shared = kind == "fjord-ice"
    out = {"runs": rr, "onset_s": t0, "shared": shared, "step_s": step_s,
           "t_read_end_h": (t_all[-1] - t0) / 3600.0}
    for fld, t, frames in (("raw", t_all[raw_keep], raw), ("filt", tk, filt)):
        lim_src = {n: np.where(_limit_region(rr[n]), e, np.nan) for n, e in frames.items()}
        out[fld] = {"t": t, "eta": frames, "lims": colour_limits(lim_src, shared=shared)}
    return out


def check_scales(prep: dict, fld: str) -> dict:
    """Fixed limits and their saturation for one field ('raw' or 'filt'): over the shown pixels,
    over the limit region (the bay for R1), and the share of the saturated pixels that lie in
    columns shallower than SHALLOW_REPORT_M."""
    P = prep[fld]
    out = {}
    for n, e in P["eta"].items():
        r = prep["runs"][n]
        L = P["lims"][n][1]
        reg = _limit_region(r)
        shal = r["H"] < SHALLOW_REPORT_M
        sat = [saturated_fraction(e[k], L) for k in range(len(e))]
        sat_reg = [saturated_fraction(np.where(reg, e[k], np.nan), L) for k in range(len(e))]
        hits = np.abs(e) >= L
        nh = int(hits.sum())
        a = np.abs(e[np.isfinite(e)])
        out[n] = {"lim_m": L, "p99_5_abs_m_limit_region": float(np.nanpercentile(
                      np.abs(np.where(reg, e, np.nan)), PCTL)),
                  "max_abs_m": float(a.max()),
                  "saturated_mean": float(np.mean(sat)), "saturated_max": max(sat),
                  "saturated_limit_region_mean": float(np.mean(sat_reg)),
                  "saturated_limit_region_max": max(sat_reg),
                  "limit_region": "bay" if is_bay_run(r["case"]) else "all wet cells",
                  f"share_of_saturated_in_columns_shallower_than_{SHALLOW_REPORT_M:g}m":
                      (float(hits[:, shal].sum() / nh) if nh else None),
                  f"fraction_of_shown_pixels_in_columns_shallower_than_{SHALLOW_REPORT_M:g}m":
                      float(np.mean(shal[np.isfinite(e[0])])) if np.isfinite(e[0]).any() else None,
                  "frames": len(e)}
        if fld == "filt":
            out[n]["T_f_h"] = tf_s(r["f"]) / 3600.0
    out["shared"] = prep["shared"]
    out["t_first_h"] = float((P["t"][0] - prep["onset_s"]) / 3600)
    out["t_last_h"] = float((P["t"][-1] - prep["onset_s"]) / 3600)
    out["step_s"] = prep["step_s"]
    return out


def figure_for(kind: str, prep: dict, fld: str):
    rr = prep["runs"]
    P = prep[fld]
    t = P["t"]
    t_tau = np.arange(prep["onset_s"], t[-1] + 1e-6, TAU_STEP_S)    # the pulse is always shown
    tau = {n: _tau_peak(r["ds"], t_tau) for n, r in rr.items()}
    sym = "$a_1$" if fld == "raw" else "$a_1'$"
    tf_h = {n: tf_s(r["f"]) / 3600.0 for n, r in rr.items()}
    panels = []
    if kind == "fjord-ice":
        for n, r in rr.items():
            lab = "landfast ice" if _ice(r) is not None else "open water"
            panels.append(make_panel(r, P["eta"][n], P["lims"][n], f"{n}: {lab}, winter profile",
                                     tau=tau[n]))
        names = list(rr)
        if np.allclose(tau[names[0]], tau[names[-1]]):     # same forcing: one curve in the strip
            for p in panels[1:]:
                p.tau = None
            panels[0].name = "wind forcing (both runs)"
        L = P["lims"][names[0]][1]
        what = ("mode-1 isopycnal displacement $a_1$" if fld == "raw" else
                "inertially filtered mode-1 displacement $a_1'$ (running 1-T$_f$ mean removed)")
        sup = (f"DESCRIPTIVE: {what}\nopen water vs landfast ice, same winter profile and wind "
               f"pulse; shared fixed scale $\\pm${L:g} m")
        foot = (A1_DEF + "\n" + (RAW_CAVEAT + NO_ICE_WIND if fld == "raw" else
                                  " ".join(x for x in (filt_definition(tf_h), FILTER_NOTE,
                                                       saturation_note(check_scales(prep, fld))) if x)))
    else:
        for n, r in rr.items():
            bay = is_bay_run(r["case"])
            lab = "thesis rectangular bay" if bay else "Van Mijenfjorden, summer profile"
            panels.append(make_panel(r, P["eta"][n], P["lims"][n],
                                     f"{n}: {lab}; {sym} $\\pm${P['lims'][n][1]:g} m"
                                     + (" (bay)" if bay else ""), tau=tau[n]))
        what = ("mode-1 isopycnal displacement $a_1$" if fld == "raw" else
                "inertially filtered $a_1'$ (running 1-T$_f$ mean removed)")
        sup = (f"DESCRIPTIVE: thesis bay vs Van Mijenfjorden after a wind pulse; {what}\nSEPARATE "
               "fixed colour scales (different profiles, depths, f and wind footprints: amplitudes "
               "not comparable between panels)")
        foot = (A1_DEF + " " + BRIDGE + "\n" +
                (RAW_CAVEAT + "." if fld == "raw" else
                 " ".join(x for x in (filt_definition(tf_h), FILTER_NOTE,
                                      saturation_note(check_scales(prep, fld))) if x)))
    return build_figure(panels, t, prep["onset_s"], sup, foot, tau_t_s=t_tau, cbar_label=sym)


DEFAULTS = {"fjord-ice": (["F3", "F2"], "anim_fjord_ice"),
            "bay-fjord": (["R1", "F1"], "anim_bay_vs_fjord")}


def main(argv: list | None = None) -> int:
    from pfcommon.paths import EXP4_RUNS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("kind", choices=sorted(DEFAULTS))
    ap.add_argument("--runs", nargs=2, default=None)
    ap.add_argument("--runs-root", type=Path, default=EXP4_RUNS)
    ap.add_argument("--out-dir", type=Path, default=FIG_DIR,
                    help="writes <name>_raw.gif and <name>_filt.gif here")
    ap.add_argument("--t-end-h", type=float, default=T_END_H)
    ap.add_argument("--step-s", type=float, default=None,
                    help="frame step [s] (default: lcm of the record spacings)")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--dpi", type=int, default=80)
    ap.add_argument("--max-mb", type=float, default=5.0)
    ap.add_argument("--frame-png", type=Path, default=None,
                    help="prefix of still PNGs: <prefix>_{raw,filt}_<h>h.png")
    ap.add_argument("--frame-h", type=float, nargs="+", default=[24.0],
                    help="time(s) of the --frame-png still(s) [h]")
    ap.add_argument("--save-frames", type=Path, default=None,
                    help="also write the raw a_1 frames, times and land masks to this .npz")
    ap.add_argument("--fields", nargs="+", choices=("raw", "filt"), default=["raw", "filt"],
                    help="which GIFs to write (default both)")
    ap.add_argument("--check-scales", action="store_true",
                    help="only report the fixed scales and their saturation (all frames); no GIF")
    a = ap.parse_args(argv)
    runs, name = DEFAULTS[a.kind]
    runs = a.runs or runs
    prep = prepare(a.kind, runs, a.runs_root, a.t_end_h, a.step_s)
    res = {"t_read_end_h": prep["t_read_end_h"]}
    if a.save_frames is not None:          # for checks outside the GIF (not a figure output)
        np.savez_compressed(a.save_frames, t_s=prep["raw"]["t"], onset_s=prep["onset_s"],
                            **{f"a1_{n}": e for n, e in prep["raw"]["eta"].items()},
                            **{f"land_{n}": r["land"] for n, r in prep["runs"].items()})
    over = False
    for fld in a.fields:
        cs = check_scales(prep, fld)
        if a.check_scales:
            res[fld] = {"scales": cs}
            continue
        fig, draw = figure_for(a.kind, prep, fld)
        P = prep[fld]
        r = write_gif(fig, draw, len(P["t"]), a.out_dir / f"{name}_{fld}.gif", a.max_mb, a.dpi,
                      a.stride)
        over |= r["over_budget"]
        if a.frame_png is not None:
            r["frame_png"] = []
            for h in a.frame_h:
                k = int(np.argmin(np.abs((P["t"] - prep["onset_s"]) / 3600 - h)))
                draw(k)
                png = a.frame_png.with_name(f"{a.frame_png.name}_{fld}_{h:g}h.png")
                fig.savefig(png, dpi=a.dpi * 2)
                r["frame_png"].append({"png": str(png),
                                       "t_h": float((P["t"][k] - prep["onset_s"]) / 3600)})
        plt.close(fig)
        r["scales"] = cs
        res[fld] = r
    print(json.dumps(res, indent=2))
    return 1 if over else 0


if __name__ == "__main__":
    sys.exit(main())
