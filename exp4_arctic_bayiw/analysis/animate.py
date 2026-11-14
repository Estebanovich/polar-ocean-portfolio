#!/opt/anaconda3/envs/pyocean/bin/python
r"""Two DESCRIPTIVE GIF animations of Exp. 4 (Task 9): plan views of the mode-1 isopycnal displacement.

Usage::

    animate.py fjord-ice [--runs F3 F2] [--out exp4_arctic_bayiw/figures/anim_fjord_ice.gif]
    animate.py bay-fjord [--runs R1 F1] [--out exp4_arctic_bayiw/figures/anim_bay_vs_fjord.gif]
        common options: [--t-end-h 96] [--step-s S] [--stride N] [--dpi D] [--max-mb 5]
                        [--frame-png PATH] [--frame-h H] [--check-scales] [--runs-root EXP4_RUNS]

No number of these animations enters the README claims; no hypothesis or criterion uses them.

Field (choice and justification). a_1(x, y, t) = the mode-1 amplitude of the isopycnal
displacement zeta = -rho'/(d rho_b/dz) [m, + up], from each 3-D density snapshot:

    a_1 = sum_k N_b^2 zeta W_1 dz / sum_k N_b^2 W_1^2 dz = (g/rho0) sum_k rho' W_1 dz / sum_k N_b^2 W_1^2 dz

(N_b^2-weighted projection, the Sturm-Liouville weight of the W modes, Gill 1982 sec. 6.11
[verificar]; N_b^2 zeta = (g/rho0) rho'). rho' = rho - rho_b with rho_b the FIXED reference
profile of the energy budget (PD0 Enmienda 1: budget.background_profile of the t = 0 density);
W_1 = theory.modes.vertical_modes of rho_b over the local number of wet levels, scaled to
max |W_1| = 1, so a_1 is the displacement where the mode-1 displacement peaks. Why not the
single-level displacement at the level k* of maximum N_b^2 (validate.displacement, used by
V2b/V3): it was tried first and is unusable in winter. The winter profile is a staircase with
N_b^2 ~ 0 between steps; its k* is the top cell (z = -0.75 m), whose 0.03 kg m^-3 surface step the
pulse mixes away within 3 h, so F2 and F3 showed a uniform +1.6 to +1.8 m "displacement" (median
over the wet cells) that is mixing, not waves; deeper levels either have d rho_b/dz ~ 0 (eta
undefined) or drift by metres over 96 h. The projection needs no division by d rho_b/dz, weights
the surface (W_1 = 0 at the lid) little, and gave wet-cell medians of -0.1 to -0.2 m in F2/F3 over
3-96 h against a spatial spread of 0.7-1.6 m. Caveats: W_1 is that of a flat-bottom column of the
local depth (topographic coupling ignored); rho' changes that are not vertical advection (KPP
mixing, the profile drift, ice-ocean fluxes in F2) still project on mode 1: under the F2 ice the
inner fjord shows a slow, monotonic negative a_1 with no oscillation (drift, not waves), and the
open-boundary corners (NW, S sponge) carry noise. a_1 is descriptive. Checked on the frames
(Hovmoeller along the inner-fjord coasts, --save-frames): in F1 the signal travels east along the
south coast and west along the north coast, coast on the right (cyclonic, as an internal Kelvin
wave in the northern hemisphere); in F3 the near-inertial oscillation dominates the coasts.

Animations.

* ``fjord-ice`` (Anim A): open water against landfast ice. The brief names "F1 vs F1c"; in
  cases.yaml F1c is the thesis closure (constant mixing), not an ice case, and the landfast-ice
  case is F2 (winter profile). The comparison that changes only the ice is F2 against F3 (same
  winter profile, grid and wind; the pair of H4b and Figure 4), so the default is F3 | F2. Shared
  colour scale (same profile, grid and pulse, so a_1 is directly comparable). F2 shows the outline of the
  initial ice cover (AREA.bin >= 0.5, the actual model input, not the nominal polygon).
* ``bay-fjord`` (Anim B): the thesis bay R1 against Van Mijenfjorden F1, each on its own axes.
  Separate fixed scales, stated in the title: the two runs have different profiles, depths
  (R1 bay 163 m, open part 600 m; F1 mostly < 120 m), f (6.97e-5 against 1.43e-4 s^-1) and pulse
  footprints (a 300 km thesis patch against the fjord-wide pulse), so the amplitudes are not
  physically comparable (decided from the configurations before any frame was drawn). The bay
  view is cropped to the wind patch (bay_crop): the thesis grid stretches to ~2800 km.
  Descriptive footer: the bridge is not resolution-converged and the fixed reference drifts.

Conventions (as exp1_fjord_plumes/analysis/animate_c.py and animate_b.py): land light grey; water
columns shallower than MIN_DEPTH_M = 10 m (or < 3 wet levels) a lighter grey (no a_1); FIXED symmetric colour
limits for all frames, L = the ceiling on the NICE ladder (1, 1.5, 2, 3, 4, 5, 6, 8 x 10^n) of the 99.5th percentile of |a_1| over every pixel of every
frame (pooled over the panels when shared), so at most ~0.5 % of the pixels saturate
(``--check-scales`` reports the fraction); the strait sections (figures/bathy_report.json) on the
fjord panels; a clock in hours after the pulse onset and in inertial periods 2 pi / |f| of each
run; a strip with the peak wind stress of the forcing files (load.forcing_stress, the stress the
model reads; under the landfast ice of F2 the ocean does not receive it) with a moving marker.

Time. Pulse onset t_on = centre - 2.5 sigma of the half-Gaussian ramp (gendata/forcing.pulse and
make_case.thesis_wind_patch: Gaussian rise, zero after the centre), clipped at 0 (fjord and thesis:
centre 2.5 h, sigma 1 h, so t_on = 0 = the model start, where tau is 4 % of its peak). Frames from t_on to
t_on + 96 h by default (4 days, 7.8 fjord inertial periods of 12.23 h; 3.8 of the thesis f). Anim
A uses every 1200 s snapshot (before the GIF stride); Anim B subsamples both runs to the common step 3600 s (the bay has
900 s snapshots, the fjord 1200 s). GIF: Pillow, the stride doubles until the file fits --max-mb
(animate_b.write_within_budget); exit code 1 if it is still over budget.
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
MIN_DEPTH_M = 10.0        # nor columns shallower than this (only the pulse-mixed surface layer)        # pulse onset = centre - 2.5 sigma
FPS = 10
LAND_GREY, SHALLOW_GREY = "0.80", "0.92"
FOOTER_B = ("Descriptive only: the bay-fjord bridge is not resolution-converged, and the fixed t = 0 "
            "reference profile drifts during the runs.")
FOOTER_A = ("Descriptive only. $a_1$ = N$^2$-weighted projection of the isopycnal displacement on the "
            "local mode-1 structure (max = 1), fixed t = 0 reference profile (energy-budget $\\rho_b$).")


# ---- pure pieces --------------------------------------------------------------------------------

def mode1_shape(rho_b, RC, drF, kb: int, rho0: float, g: float = 9.81):
    """Mode-1 vertical-displacement structure W_1 (rigid lid, flat bottom; theory.modes.
    vertical_modes, Gill 1982 sec. 6.11 [verificar]) of the background profile in a column of kb
    wet levels, on the level centres RC[:kb], scaled to max |W_1| = 1 (positive lobe), and
    N_b^2 [s^-2] on the same levels. The column spans 0 .. -sum(drF[:kb]); W = 0 at both ends.
    None if kb < MIN_LEVELS or the column is shallower than MIN_DEPTH_M: the winter surface step
    (top ~5-7 m) is mixed away by the pulse within 3 h, and columns of only that layer gave the
    largest |a_1| (up to 19 m, one-cell coastal fringes) with no wave content."""
    from exp4_arctic_bayiw.analysis import budget
    from exp4_arctic_bayiw.theory import modes
    H = float(np.sum(drF[:kb]))
    if kb < MIN_LEVELS or H < MIN_DEPTH_M:
        return None
    n2 = budget.n2_background(rho_b, RC, rho0, g)[:kb]
    z = np.r_[0.0, np.asarray(RC[:kb], float), -H]
    n2e = np.r_[n2[0], n2, n2[-1]]
    _, phi = modes.vertical_modes(z, n2e, H, n_modes=1)
    w = phi[0, 1:-1]
    return w / np.max(np.abs(w)), np.maximum(n2, modes.N2_FLOOR)


def mode1_amplitude(rho, rho_b, RC, drF, hfac, rho0: float, g: float = 9.81) -> np.ndarray:
    """Mode-1 displacement amplitude a(x, y) [m, + up] of the 3-D density rho (nz, ny, nx):
    the N_b^2-weighted projection of the isopycnal displacement zeta = -rho'/(d rho_b/dz) on the
    local W_1 (max |W_1| = 1, so a is the displacement where the mode-1 displacement peaks),

        a = sum_k N_b^2 zeta W_1 dz / sum_k N_b^2 W_1^2 dz
          = (g / rho0) sum_k rho' W_1 dz / sum_k N_b^2 W_1^2 dz,

    using N_b^2 zeta = (g / rho0) rho' (N_b^2 = -(g/rho0) d rho_b/dz), so no division by
    d rho_b/dz is needed where the profile is unstratified (N_b^2 = 0 layers of the winter
    profile). rho' = rho - rho_b with the fixed reference rho_b (Enmienda 1); dz = drF hFacC;
    W_1 of the column depth (number of wet levels kb, theory.modes). NaN where mode1_shape is None (too few levels or too shallow)."""
    rho = np.asarray(rho, float)
    hfac = np.asarray(hfac, float)
    nz = rho.shape[0]
    kb = np.sum(hfac > 0, axis=0)
    W = np.zeros_like(rho)
    den = np.full(kb.shape, np.nan)
    for n in np.unique(kb):
        sh = mode1_shape(rho_b, RC, drF, int(n), rho0, g)
        if sh is None:
            continue
        w, n2 = sh
        m = kb == n
        W[:n, m] = w[:, None]
        den[m] = np.sum(n2 * w * w * np.asarray(drF[:n], float))
    rp = np.where(hfac > 0, rho - np.asarray(rho_b, float)[:nz, None, None], 0.0)
    dz = np.asarray(drF, float)[:, None, None] * hfac
    return g / rho0 * np.sum(rp * W * dz, axis=0) / den


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
                 tau_t_s=None) -> tuple:
    """Figure with one plan-view panel per run (fixed colour limits), a clock per panel and, if any
    panel has ``tau``, a wind-stress strip. Returns (fig, draw(n))."""
    style.apply()
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
                          bottom=0.16 if has_tau else 0.12)
    cmap = style.CMAPS["diverging"].copy()
    cmap.set_bad((0, 0, 0, 0))
    bg = ListedColormap([(1, 1, 1, 0), SHALLOW_GREY, LAND_GREY])
    for q, p in enumerate(panels):
        ax = fig.add_subplot(gs[0, q])
        code = np.where(p.land, 2, np.where(p.shallow, 1, 0))
        ax.pcolormesh(p.x_km, p.y_km, code, cmap=bg, vmin=-0.5, vmax=2.5, shading="flat")
        p.mesh = ax.pcolormesh(p.x_km, p.y_km, np.ma.masked_invalid(p.eta[0]), cmap=cmap,
                               vmin=p.lim[0], vmax=p.lim[1], shading="flat")
        for xs, ys, lab in p.lines:
            ax.plot(xs, ys, color="k", lw=1.6, solid_capstyle="butt")
            ax.annotate(lab, (xs[-1], ys[-1]), xytext=(4, 2), textcoords="offset points",
                        fontsize=7)
        if p.ice is not None:
            xc = 0.5 * (p.x_km[1:] + p.x_km[:-1])
            yc = 0.5 * (p.y_km[1:] + p.y_km[:-1])
            ax.contour(xc, yc, p.ice, levels=[0.5], colors="#7a3fa0", linewidths=1.2)
            ax.plot([], [], color="#7a3fa0", lw=1.2, label="landfast ice (initial)")
            ax.legend(loc="lower right", fontsize=7, frameon=True)
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
            cb.set_label(r"$a_1$ (m, + up)")
    if shared:
        cax = fig.add_axes([0.915, 0.30 if has_tau else 0.16, 0.012, 0.50 if has_tau else 0.64])
        cb = fig.colorbar(panels[0].mesh, cax=cax, extend="both")
        cb.set_label(r"$a_1$ (m, + up)")
    marks = []
    if has_tau:
        axt = fig.add_subplot(gs[1, :])
        tt = np.asarray(tau_t_s if tau_t_s is not None else t_s, float)
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
            mk.set_data([(t_s[k] - onset_s) / 3600.0], [p.tau[k]])

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
    """Grid, fixed reference rho_b, f and pulse onset of one run (load.run; LaCie read-only)."""
    from exp4_arctic_bayiw.analysis import budget, load
    ds = load.run(run_dir)
    hfac = ds["hFacC"].values
    wet3 = hfac > 0
    rho_b = budget.background_profile(ds["RHO"].isel(time=0).values, wet3)   # Enmienda 1
    kb = wet3.sum(axis=0)
    H = np.r_[0.0, np.cumsum(ds["drF"].values)][kb]           # depth of the wet levels [m]
    return {"ds": ds, "rho_b": rho_b, "hfac": hfac, "land": kb == 0,
            "shallow": (kb > 0) & ((kb < MIN_LEVELS) | (H < MIN_DEPTH_M)), "f": _f0(run_dir),
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
        out.append(mode1_amplitude(rho, r["rho_b"], RC, drF, r["hfac"], rho0, g).astype(np.float32))
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


def make_panel(r: dict, eta: np.ndarray, lim, title: str, tau=None, fjord: bool = True) -> Panel:
    ds = r["ds"]
    x = edges(ds["XC"].values[0], ds["dxF"].values[0]) / 1e3
    y = edges(ds["YC"].values[:, 0], ds["dyF"].values[:, 0]) / 1e3
    if fjord:
        return Panel(name=r["case"]["params"]["case"], title=title, x_km=x, y_km=y,
                     land=r["land"], shallow=r["shallow"], eta=eta, lim=lim, f=r["f"],
                     lines=_sections_km(), ice=_ice(r), xlabel="easting UTM 33N (km)",
                     ylabel="northing (km)", tau=tau)
    j0, j1, i0, i1 = bay_crop(r["case"]["params"]["wind"], *r["land"].shape)
    cx, cy = x[i0], y[j0]
    return Panel(name=r["case"]["params"].get("case", "R1"), title=title, x_km=x[i0:i1 + 1] - cx,
                 y_km=y[j0:j1 + 1] - cy, land=r["land"][j0:j1, i0:i1],
                 shallow=r["shallow"][j0:j1, i0:i1], eta=eta[:, j0:j1, i0:i1], lim=lim, f=r["f"],
                 xlabel="x (km, from the crop edge)", ylabel="y (km)", tau=tau)


def bay_crop(wind: dict, ny: int, nx: int, margin: float = 1.0) -> tuple:
    """(j0, j1, i0, i1) of the thesis-bay view: the wind patch (centre cx, cy, radius in cells,
    thesis wind_patch.py arguments in case.json) widened by ``margin``, up to the bay head (the
    north edge). The stretched outer grid (thousands of km) is left out."""
    r = margin * wind["radius"]
    return (max(0, int(wind["cy"] - r)), ny, max(0, int(wind["cx"] - r)),
            min(nx, int(wind["cx"] + r) + 1))


def prepare(kind: str, runs: list, root: Path, t_end_h: float = T_END_H, step_s=None) -> dict:
    """Read the runs and build the frame arrays, fixed limits and titles of one animation."""
    rr = {n: read_run(root / n) for n in runs}
    onsets = {n: r["onset_s"] for n, r in rr.items()}
    if len(set(onsets.values())) != 1:
        raise ValueError(f"pulse onsets differ: {onsets}")
    t0 = next(iter(onsets.values()))
    recs = [r["ds"]["time"].values for r in rr.values()]
    if step_s is None:
        steps = [float(np.median(np.diff(t))) for t in recs]
        step_s = float(np.lcm.reduce([int(round(s)) for s in steps]))
    t = common_times(recs, step_s, t0, t0 + t_end_h * 3600.0)
    eta = {}
    for n, r in rr.items():
        e = eta_frames(r, t)
        if "bbox_utm" not in r["case"]["params"]:          # thesis bay: the cropped view only
            j0, j1, i0, i1 = bay_crop(r["case"]["params"]["wind"], *r["land"].shape)
            e = np.where(_crop_mask(e.shape[1:], (j0, j1, i0, i1)), e, np.nan).astype(np.float32)
        eta[n] = e
    shared = kind == "fjord-ice"
    lims = colour_limits(eta, shared=shared)
    return {"runs": rr, "t": t, "onset_s": t0, "eta": eta, "lims": lims, "shared": shared,
            "step_s": step_s}


def _crop_mask(shape, crop) -> np.ndarray:
    j0, j1, i0, i1 = crop
    m = np.zeros(shape, bool)
    m[j0:j1, i0:i1] = True
    return m


def check_scales(prep: dict) -> dict:
    out = {}
    for n, e in prep["eta"].items():
        L = prep["lims"][n][1]
        sat = [saturated_fraction(e[k], L) for k in range(len(e))]
        a = np.abs(e[np.isfinite(e)])
        out[n] = {"lim_m": L, "p99_5_abs_eta_m": float(np.percentile(a, PCTL)), "max_abs_eta_m": float(a.max()),
                  "saturated_max": max(sat), "saturated_mean": float(np.mean(sat)),
                  "frames": len(e)}
    out["shared"] = prep["shared"]
    out["t_first_h"] = float((prep["t"][0] - prep["onset_s"]) / 3600)
    out["t_last_h"] = float((prep["t"][-1] - prep["onset_s"]) / 3600)
    out["step_s"] = prep["step_s"]
    return out


def figure_for(kind: str, prep: dict):
    rr, t = prep["runs"], prep["t"]
    tau = {n: _tau_peak(r["ds"], t) for n, r in rr.items()}
    panels = []
    if kind == "fjord-ice":
        for n, r in rr.items():
            ice = _ice(r) is not None
            lab = "landfast ice" if ice else "open water"
            panels.append(make_panel(r, prep["eta"][n], prep["lims"][n],
                                     f"{n}: {lab}, winter profile", tau=tau[n]))
        # same wind forcing in both runs: one curve in the strip
        names = list(rr)
        if np.allclose(tau[names[0]], tau[names[-1]]):
            for p in panels[1:]:
                p.tau = None
            panels[0].name = "wind forcing (both runs)"
        L = prep["lims"][names[0]][1]
        sup = (f"Mode-1 isopycnal displacement $a_1$: open water vs landfast ice, same winter "
               f"profile and wind pulse; shared fixed scale $\\pm${L:g} m")
        foot = FOOTER_A
    else:
        for n, r in rr.items():
            fj = "ice" in r["case"]["params"] and "bbox_utm" in r["case"]["params"]
            lab = ("Van Mijenfjorden, summer profile" if fj else "thesis rectangular bay")
            panels.append(make_panel(r, prep["eta"][n], prep["lims"][n],
                                     f"{n}: {lab}; scale $\\pm${prep['lims'][n][1]:g} m",
                                     tau=tau[n], fjord=fj))
        sup = ("DESCRIPTIVE: thesis bay vs Van Mijenfjorden after a wind pulse; mode-1 isopycnal "
               "displacement $a_1$, SEPARATE fixed colour scales\n(different profiles, depths, "
               "f and pulse footprints: amplitudes not comparable between panels)")
        foot = FOOTER_B
    return build_figure(panels, t, prep["onset_s"], sup, foot)


DEFAULTS = {"fjord-ice": (["F3", "F2"], "anim_fjord_ice.gif"),
            "bay-fjord": (["R1", "F1"], "anim_bay_vs_fjord.gif")}


def main(argv: list | None = None) -> int:
    from pfcommon.paths import EXP4_RUNS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("kind", choices=sorted(DEFAULTS))
    ap.add_argument("--runs", nargs=2, default=None)
    ap.add_argument("--runs-root", type=Path, default=EXP4_RUNS)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--t-end-h", type=float, default=T_END_H)
    ap.add_argument("--step-s", type=float, default=None,
                    help="frame step [s] (default: lcm of the record spacings)")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--dpi", type=int, default=80)
    ap.add_argument("--max-mb", type=float, default=5.0)
    ap.add_argument("--frame-png", type=Path, default=None)
    ap.add_argument("--frame-h", type=float, nargs="+", default=[24.0],
                    help="time(s) of the --frame-png still(s) [h]; several add a _<h>h suffix")
    ap.add_argument("--save-frames", type=Path, default=None,
                    help="also write the a_1 frames, times and land masks to this .npz")
    ap.add_argument("--check-scales", action="store_true",
                    help="only report the fixed scales and their saturation (all frames); no GIF")
    a = ap.parse_args(argv)
    runs, name = DEFAULTS[a.kind]
    runs = a.runs or runs
    prep = prepare(a.kind, runs, a.runs_root, a.t_end_h, a.step_s)
    cs = check_scales(prep)
    if a.save_frames is not None:          # for checks outside the GIF (not a figure output)
        np.savez_compressed(a.save_frames, t_s=prep["t"], onset_s=prep["onset_s"],
                            **{f"a1_{n}": e for n, e in prep["eta"].items()},
                            **{f"land_{n}": r["land"] for n, r in prep["runs"].items()})
    if a.check_scales:
        print(json.dumps(cs, indent=2))
        return 0
    fig, draw = figure_for(a.kind, prep)
    res = write_gif(fig, draw, len(prep["t"]), a.out or FIG_DIR / name, a.max_mb, a.dpi, a.stride)
    if a.frame_png is not None:
        res["frame_png"] = []
        for h in a.frame_h:
            k = int(np.argmin(np.abs((prep["t"] - prep["onset_s"]) / 3600 - h)))
            draw(k)
            png = a.frame_png if len(a.frame_h) == 1 else \
                a.frame_png.with_name(f"{a.frame_png.stem}_{h:g}h{a.frame_png.suffix}")
            fig.savefig(png, dpi=a.dpi * 2)
            res["frame_png"].append({"png": str(png),
                                     "t_h": float((prep["t"][k] - prep["onset_s"]) / 3600)})
    plt.close(fig)
    res["scales"] = cs
    print(json.dumps(res, indent=2))
    return 1 if res["over_budget"] else 0


if __name__ == "__main__":
    sys.exit(main())
