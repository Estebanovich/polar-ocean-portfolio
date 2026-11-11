"""Figure 0 of Experiment 3: the protocol, in the style and method colours of Figures 1-4.

(a) The training data: the ISOMIP+ ice shelf on its centre line (the real Ocean1 input geometry, Cornford and
    Asay-Davis 2016), the ocean coloured by one ramp profile of the Exp. 2 design, and the melt map that MITgcm returns.
(b) The 29 runs of Exp. 2 in the (Tbot, z_c) plane as designed, with the held-out sets of pd0_preregistro.md and the
    hypothesis each one serves.
(c) The eight methods between physics and data, and the factorisation of the hybrid M7.
(d) What H3c asks of the ensemble: a spread that grows with the error.

No result of Exp. 3 is drawn. The melt maps, the power law and the ensemble members are an illustration computed here
from the input geometry and one input profile (melt = TF_+^2 x a geometric pattern); they are not output of MITgcm or
of any method.

Usage: fig0_protocol.py   (writes figures/fig0_protocol.png and .pdf; needs data_external/Ocean1_input_geom_v1.01.nc)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
REPO = FIG.parent.parent
sys.path.insert(0, str(REPO))

import cmocean  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4  # noqa: E402
from matplotlib.colors import ListedColormap, PowerNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

from exp2_isomip_cavity.gendata.profiles import ramp  # noqa: E402
from exp3_melt_emulator.analysis.figures import METHOD_COLOURS  # noqa: E402
from exp3_melt_emulator.data import freezing_temperature  # noqa: E402
from pfcommon import style  # noqa: E402

GEOM = REPO / "exp2_isomip_cavity" / "data_external" / "Ocean1_input_geom_v1.01.nc"
OUT = FIG / "fig0_protocol"
FW, FH = 8.0, 6.25                         # figure size (in); every position below is in inches
TOP_SHIFT, TOP_ROW = 0.25, 2.7             # the top row, (a) and (b), drawn above y = 2.7 in, is raised by 0.25 in

# The design of the runs (exp3 dataset; the case names encode Tbot, z_c and s): inputs, not results.
FACTORIAL = [(t, z) for z in (250, 400, 550) for t in (-0.5, 0.0, 0.5, 1.0, 1.5)]
GEOMETRY = [(0.75, 527, 0.8), (1.43, 295, 0.8), (-0.05, 328, 0.8), (-0.42, 402, 0.8), (0.34, 444, 0.8),
            (1.42, 481, 1.2), (0.02, 291, 1.2), (-0.39, 503, 1.2), (0.97, 386, 1.2), (0.66, 311, 1.2)]
FINAL = [(2.0, 250), (2.0, 400), (2.0, 550)]
OCEAN0 = (1.0, 360)                       # ISOMIP+ Ocean0 linear profile (z_c 360 m, d 720 m): profile test
EXAMPLE = dict(Tbot=1.0, z_c=400.0)       # the profile drawn in (a), one factorial run

TRAIN_C, GEOM_C, TEST_C, GREY = "#0072B2", "#56B4E9", "#C0392B", "0.35"
ICE_C, ICE_EDGE, BED_C = "#eef4f8", "#7d9cb5", "#9c8670"
X0, X1 = 420.0, 680.0                     # km shown in the maps


def load_geometry():
    with netCDF4.Dataset(GEOM) as d:
        x = np.asarray(d["x"][:], dtype=float) / 1e3
        y = np.asarray(d["y"][:], dtype=float) / 1e3
        g = {k: np.asarray(d[k][:], dtype=float) for k in
             ("upperSurface", "lowerSurface", "bedrockTopography", "floatingMask", "groundedMask")}
    return x, y, g


def illustration(x, y, g):
    """Toy melt = TF_+^2 x pattern on the floating cells, for the drawing only (not a model run)."""
    T, S = ramp(g["lowerSurface"], **EXAMPLE)
    tf = np.clip(T - freezing_temperature(S, g["lowerSurface"]), 0, None)
    fl = g["floatingMask"] > 0.5
    gl = np.array([x[np.argmax(fl[j])] if fl[j].any() else np.nan for j in range(len(y))])
    dist = x[None, :] - gl[:, None]
    pattern = (0.35 + np.exp(-dist / 45.0)) * (1.0 + 0.45 * (y[:, None] - 40.0) / 40.0)
    pattern = np.where(fl, pattern / np.nanmax(np.where(fl, pattern, np.nan)), np.nan)
    melt = np.where(fl, tf ** 2 * pattern, np.nan)
    return np.where(fl, tf, np.nan), pattern, melt / np.nanmax(melt)


def plan(ax, x, y, field, g, cmap, gamma=1.0):
    sel = (x >= X0) & (x <= X1)
    ax.set_facecolor("#cfe0ec")                                   # open ocean
    grounded = np.where(g["groundedMask"][:, sel] > 0.5, 1.0, np.nan)
    ax.pcolormesh(x[sel], y, grounded, cmap=ListedColormap([ICE_C]), shading="nearest", rasterized=True)
    ax.pcolormesh(x[sel], y, field[:, sel], cmap=cmap, norm=PowerNorm(gamma, vmin=0, vmax=1), shading="nearest",
                  rasterized=True)
    ax.set_xlim(X0, X1)
    ax.set_ylim(0, 80)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.5)
        s.set_color("0.5")


def up(yin):
    return yin + TOP_SHIFT if yin > TOP_ROW else yin


def letter(fig, xin, yin, k, title):
    yin = up(yin)
    fig.text(xin / FW, yin / FH, f"({k})", fontsize=10, fontweight="bold", va="baseline")
    fig.text((xin + 0.3) / FW, yin / FH, title, fontsize=9.5, va="baseline")


def axes(fig, l, b, w, h):
    return fig.add_axes([l / FW, up(b) / FH, w / FW, h / FH])


def ftext(fig, xin, yin, s, **kw):
    return fig.text(xin / FW, up(yin) / FH, s, **kw)


def arrow(fig, p, q, text=None, dy=0.07, color="0.25"):
    fig.patches.append(FancyArrowPatch((p[0] / FW, up(p[1]) / FH), (q[0] / FW, up(q[1]) / FH),
                                       transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=9, lw=1.0, color=color))
    if text:
        ftext(fig, (p[0] + q[0]) / 2, (p[1] + q[1]) / 2 + dy, text, ha="center", va="bottom", fontsize=7.2,
              color=color)


# ----------------------------------------------------------------------------------------------- (a) training data
def panel_a(fig, x, y, g, melt):
    letter(fig, 0.05, 5.78, "a", "The training data: the MITgcm runs of Exp. 2")
    ax = axes(fig, 0.42, 3.55, 2.75, 1.95)
    j = int(np.argmin(np.abs(y - 40.0)))
    sel = (x >= 380) & (x <= 700)
    xs, top, base, bed = x[sel], g["upperSurface"][j, sel], g["lowerSurface"][j, sel], g["bedrockTopography"][j, sel]
    fl = g["floatingMask"][j, sel] > 0.5
    zz = np.linspace(-760, 0, 200)
    T, _ = ramp(zz, **EXAMPLE)
    Z = np.repeat(zz[:, None], xs.size, axis=1)
    ceiling = np.where(fl, base, np.where(xs > 600, 0.0, -9e9))
    ocean = (Z < ceiling[None, :]) & (Z > bed[None, :])
    pm = ax.pcolormesh(xs, zz, np.where(ocean, np.repeat(T[:, None], xs.size, axis=1), np.nan),
                       cmap=cmocean.cm.thermal, vmin=-2, vmax=2, shading="nearest", rasterized=True)
    has_ice = fl | (xs < 600)
    ax.fill_between(xs, np.where(fl, base, bed), np.minimum(top, 60), where=has_ice, color=ICE_C, ec=ICE_EDGE,
                    lw=0.6)
    ax.fill_between(xs, -800, bed, color=BED_C, lw=0)
    for s in (0.8, 1.2):
        ok = fl & (base * s > bed + 5)
        ax.plot(xs[ok], base[ok] * s, color="w", lw=0.8, ls=(0, (3, 2)))
    ax.plot([560, 700], [-EXAMPLE["z_c"]] * 2, color="w", lw=0.7, ls=":")
    ax.text(697, -EXAMPLE["z_c"] + 14, "$z_c$", color="w", ha="right", fontsize=8.5)
    ax.text(697, -690, "$T_{bot}$", color="k", ha="right", fontsize=8.5)
    ax.text(697, -70, "cold", color="w", ha="right", fontsize=7.5)
    ax.text(575, -95, "ice shelf", ha="center", fontsize=7.5, color="0.2")
    ax.text(410, -260, "grounded\nice", ha="center", fontsize=7, color="0.2")
    k = int(np.argmin(np.abs(xs - 560)))
    ax.annotate("draft × $s$\n($s$ = 0.8, 1, 1.2)", xy=(xs[k], base[k] * 1.2), xytext=(470, -600), fontsize=7,
                color="w", arrowprops=dict(arrowstyle="-", lw=0.5, color="w"))
    ax.set_xlim(380, 700)
    ax.set_ylim(-760, 60)
    ax.set_xlabel("distance along the centre line (km)", fontsize=7.5, labelpad=1)
    ax.set_ylabel("depth (m)", fontsize=7.5, labelpad=1)
    ax.set_yticks([-600, -300, 0], ["600", "300", "0"])
    ax.tick_params(labelsize=7, length=2)
    cax = axes(fig, 3.24, 3.62, 0.06, 0.55)
    cb = fig.colorbar(pm, cax=cax, ticks=[-2, 0, 2])
    cb.ax.tick_params(labelsize=6.5, length=2, pad=1)
    cb.outline.set_linewidth(0.4)
    ftext(fig, 3.2, 4.24, "T (°C)", fontsize=6.8)
    ftext(fig, 0.42, 3.15, "inputs: $T_{bot}$, $z_c$ and $s$, three numbers per run", fontsize=7.5, color=GREY)
    # MITgcm -> melt map
    mx = axes(fig, 3.82, 4.55, 1.03, 0.317)
    plan(mx, x, y, melt, g, cmocean.cm.amp, gamma=0.5)
    ftext(fig, 3.82 + 0.515, 4.97, "basal melt rate", ha="center", fontsize=7.5)
    ftext(fig, 3.82, 4.42, "grounding line", fontsize=6.4, color=GREY, va="top")
    ftext(fig, 4.85, 4.42, "front", fontsize=6.4, color=GREY, ha="right", va="top")
    arrow(fig, (3.27, 4.71), (3.77, 4.71), "MITgcm", dy=0.04)
    ftext(fig, 3.52, 4.64, "hours\nper run", ha="center", va="top", fontsize=6.6, color=GREY)
    ftext(fig, 3.82 + 0.515, 4.12, "an emulator learns this\nmapping from the runs and\npredicts a map in < 1 s",
          ha="center", va="top", fontsize=7, color=GREY)


# ----------------------------------------------------------------------------------------------- (b) tests
def panel_b(fig):
    letter(fig, 4.95, 5.78, "b", "Where each method is tested")
    ax = axes(fig, 5.45, 3.72, 2.45, 1.78)
    ax.add_patch(Rectangle((-0.68, 228), 2.28, 344, fc="#e8f0f7", ec="none", zorder=0))
    ax.add_patch(Rectangle((1.36, 228), 0.24, 344, fc="#f6d5cf", ec="none", zorder=0))
    for zc in (325, 475):
        ax.plot([-0.68, 1.6], [zc, zc], color="0.45", lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.scatter(*zip(*FACTORIAL), s=16, color=TRAIN_C, lw=0, zorder=3)
    ax.scatter([t for t, _, _ in GEOMETRY], [z for _, z, _ in GEOMETRY], s=18, marker="^", color=GEOM_C,
               lw=0, zorder=3)
    ax.scatter(*zip(*FINAL), s=70, marker="*", color=TEST_C, lw=0, zorder=3)
    ax.scatter(*OCEAN0, s=24, marker="D", color="k", lw=0, zorder=3)
    ax.scatter([0.5], [400], s=110, facecolor="none", edgecolor="k", lw=0.8, ls=(0, (2, 1.5)), zorder=4)
    kw = dict(fontsize=6.8, zorder=5)
    ax.text(0.57, 388, "loro", va="bottom", **kw)
    ax.text(-0.62, 471, "zc_out: hide a third of $z_c$", color="0.3", va="bottom", **kw)
    ax.text(1.48, 222, "warm_out", color=TEST_C, ha="center", va="bottom", **kw)
    ax.annotate("Ocean0: no\nthermocline\n(profile test)", xy=(1.03, 357), xytext=(2.15, 300), **kw,
                ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="k", shrinkB=3))
    ax.text(2.18, 470, "final test,\n$T_{bot}$ = 2.0 °C,\npredicted\nonce  [H3a]", color=TEST_C, va="center",
            **kw)

    ax.set_xlim(-0.75, 3.2)
    ax.set_ylim(575, 200)
    ax.set_xticks([-0.5, 0, 0.5, 1, 1.5, 2])
    ax.set_yticks([250, 400, 550])
    ax.tick_params(labelsize=7, length=2)
    ax.set_xlabel("$T_{bot}$, deep-water temperature (°C)", fontsize=7.5, labelpad=1)
    ax.set_ylabel("$z_c$, thermocline depth (m)", fontsize=7.5, labelpad=1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ftext(fig, 5.45, 3.15, "shaded: the training range (25 runs). loro hides\neach run in turn and predicts it, "
          "inside the range\n[H3b]; the other tests extrapolate.", fontsize=6.8, color=GREY, va="top")
    hs = [Line2D([], [], ls="", marker="o", ms=4, color=TRAIN_C, label="15 runs with $s$ = 1"),
          Line2D([], [], ls="", marker="^", ms=4.5, color=GEOM_C, label="10 with $s$ = 0.8 or 1.2")]
    fig.legend(handles=hs, loc="lower left", bbox_to_anchor=(5.35 / FW, up(3.2) / FH), ncol=2, frameon=False,
               fontsize=7, handletextpad=0.2, columnspacing=1.2)


# ----------------------------------------------------------------------------------------------- (c) methods
def panel_c(fig, x, y, g, pattern, melt):
    letter(fig, 0.05, 2.62, "c", "Eight methods, from physics to data")
    ax = axes(fig, 0.42, 1.72, 4.15, 0.68)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.annotate("", xy=(1.0, 0.2), xytext=(0.0, 0.2), arrowprops=dict(arrowstyle="<|-|>", lw=0.8, color="0.4"))
    ax.text(0.0, 0.08, "physics: a formula, one fitted constant", fontsize=6.8, color=GREY, va="top")
    ax.text(1.0, 0.08, "data: learned from the runs", fontsize=6.8, color=GREY, va="top", ha="right")
    items = [(0.06, "M1", "linear"), (0.16, "M2", "quadratic"), (0.26, "M3", "nonlocal"),
             (0.44, "M7", "hybrid"), (0.62, "M4", "regression\nper cell"), (0.76, "M5", "Gaussian\nprocess"),
             (0.92, "M6", "deep\nensemble")]
    for xp, m, lab in items:
        ax.plot([xp, xp], [0.15, 0.25], color="0.4", lw=0.8)
        ax.text(xp, 0.32, lab, ha="center", va="bottom", fontsize=6.8, color="0.25")
        ax.text(xp, 0.32 + (0.3 if "\n" in lab else 0.16), m, ha="center", va="bottom", fontsize=8.5,
                fontweight="bold", color=METHOD_COLOURS[m])
    ax.text(0.0, 0.98, "baseline M0: the mean of the training maps", fontsize=6.8, color=METHOD_COLOURS["M0"],
            va="top")
    # the hybrid
    ftext(fig, 0.42, 1.38, "The hybrid M7:", fontsize=8, fontweight="bold")
    ftext(fig, 1.42, 1.38, r"melt = K TF$_+^{\,n}$ × g$_\theta$(geometry)", fontsize=8.5)
    a1 = axes(fig, 0.6, 0.42, 0.9, 0.72)
    t = np.linspace(0, 1, 100)
    tr = t <= 0.62
    a1.axvspan(0, 0.62, color="#e8f0f7", lw=0)
    a1.plot(t[tr], t[tr] ** 2.5, color="k", lw=1.4)
    a1.plot(t[~tr], t[~tr] ** 2.5, color="k", lw=1.4, ls=(0, (3, 1.5)))
    a1.text(0.03, 0.97, "training", fontsize=6.4, color=TRAIN_C, va="top")
    a1.set_xlim(0, 1)
    a1.set_ylim(0, 1.02)
    a1.set_xticks([])
    a1.set_yticks([])
    a1.set_xlabel("thermal forcing TF", fontsize=6.8, labelpad=1)
    a1.set_ylabel("melt", fontsize=6.8, labelpad=1)
    for s in ("top", "right"):
        a1.spines[s].set_visible(False)
    ftext(fig, 1.05, 0.26, "physics, how much:\nK and n fitted; continues\ninto warmer water",
          ha="center", va="top", fontsize=6.5, color=GREY)
    ftext(fig, 1.64, 0.78, "×", fontsize=13, ha="center", va="center")
    a2 = axes(fig, 1.78, 0.62, 1.03, 0.317)
    plan(a2, x, y, pattern, g, cmocean.cm.dense)
    ftext(fig, 1.78 + 0.515, 0.5, "network, where:\nsees only the geometry,\nnever a temperature",
          ha="center", va="top", fontsize=6.5, color=GREY)
    ftext(fig, 2.95, 0.78, "=", fontsize=13, ha="center", va="center")
    a3 = axes(fig, 3.08, 0.62, 1.03, 0.317)
    plan(a3, x, y, melt, g, cmocean.cm.amp, gamma=0.5)
    ftext(fig, 3.08 + 0.515, 0.5, "melt map", ha="center", va="top", fontsize=6.5, color=GREY)
    ftext(fig, 4.25, 0.78, "M6 instead takes\ntemperature as an\ninput, so it must\nextrapolate it",
          fontsize=6.5, color=METHOD_COLOURS["M6"], va="center")


# ----------------------------------------------------------------------------------------------- (d) H3c
def panel_d(fig):
    letter(fig, 4.95, 2.62, "d", "What H3c asks of the ensemble")
    s = np.linspace(0, 1, 60)
    truth = 0.25 + 0.75 * np.exp(-s / 0.3)
    rng = np.random.default_rng(3)
    cases = [("honest", "spread as large\nas the error", 1.15, 0.16),
             ("confidently wrong", "members agree,\nall far off", 1.55, 0.025)]
    for k, (title, note, bias, spread) in enumerate(cases):
        ax = axes(fig, 5.45 + k * 1.3, 0.72, 1.1, 1.5)
        mem = np.array([truth * bias * (1 + spread * rng.normal() + 0.5 * spread * rng.normal() * s)
                        for _ in range(5)])
        ax.fill_between(s, mem.min(0), mem.max(0), color=METHOD_COLOURS["M6"], alpha=0.18, lw=0)
        for m in mem:
            ax.plot(s, m, color=METHOD_COLOURS["M6"], lw=0.7)
        ax.plot(s, truth, color="k", lw=1.5)
        ax.set_title(title, fontsize=7.5, pad=3)
        ax.set_ylim(0, 1.9)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("along the ice shelf", fontsize=6.6, labelpad=1)
        if k == 0:
            ax.set_ylabel("melt", fontsize=6.8, labelpad=1)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(0.97, 0.97, note, transform=ax.transAxes, ha="right", va="top", fontsize=6.4, color=GREY)
    ftext(fig, 5.45, 0.42, "black: MITgcm;  orange: the 5 members of M6", fontsize=6.6, color=GREY)
    ftext(fig, 5.45, 0.24, "H3c fails if any extrapolation case has\nspread < 5 % and integrated error > 20 %",
          fontsize=6.6, color="k", va="top")


def main() -> None:
    style.apply()
    plt.rcParams.update({"font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "axes.linewidth": 0.6})
    x, y, g = load_geometry()
    _, pattern, melt = illustration(x, y, g)
    fig = plt.figure(figsize=(FW, FH))
    panel_a(fig, x, y, g, melt)
    panel_b(fig)
    panel_c(fig, x, y, g, pattern, melt)
    panel_d(fig)
    out = [OUT.with_suffix(".png"), OUT.with_suffix(".pdf")]
    fig.savefig(out[0], dpi=250)
    fig.savefig(out[1])
    for p in out:
        print("wrote", p)


if __name__ == "__main__":
    main()
