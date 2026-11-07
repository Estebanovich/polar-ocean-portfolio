"""Data panels for Figure 0 (the schematic): cut-outs of our own model runs, the initial profiles and
the forcing time series. Each panel is a small vector PDF that ``fig0_schematic.py`` places with
VectorCraft (``file.place``), so the composite stays vector.

Panels (written to figures/src/panels/):
  p_plume2d.pdf  config (a): T and w next to the ice, a_line2d_twolayer_q0.03 at t = 48 h
  p_face_w.pdf   config (b): window-mean (2-3 h) w on the ice face, b_3d_twolayer_x0_80 (two plumes)
  p_profiles.pdf initial T and S: two-layer (interface 200 m) and Sermilik-like (centre 200 m), H = 500 m
  p_forcing.pdf  config (c) forcing: plume Q0, wind stress pulse, shelf pycnocline depth, days 0-12

Usage: fig0_panels.py [--runs-root DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from exp1_fjord_plumes.analysis.plume_diag import _snapshots  # noqa: E402
from exp1_fjord_plumes.gendata import profiles  # noqa: E402
from pfcommon import style  # noqa: E402
from pfcommon.mitgcm_io import read_bin, read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path(__file__).resolve().parent / "panels"
RC = {"font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
      "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42}


def _save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name.replace(".pdf", ".png"), dpi=300, transparent=False, facecolor="white")
    plt.close(fig)


def plume2d(root: Path):
    run = root / "a_line2d_twolayer_q0.03"
    require_ok(run)
    p = json.loads((run / "case.json").read_text())["params"]
    T = read_mds(run / "T")[0][0]
    W = read_mds(run / "W")[0][0]
    nz, nx = T.shape
    dx = read_bin(run / "delX.bin", (nx,))
    xf = np.concatenate([[0.0], np.cumsum(dx)])
    zf = -np.arange(nz + 1) * p["dz"]
    n = int(np.searchsorted(xf, 80.0))
    xc, zc = 0.5 * (xf[1:] + xf[:-1]), 0.5 * (zf[1:] + zf[:-1])
    Wc = 0.5 * (W + np.vstack([W[1:], np.zeros((1, nx))]))       # w at cell centres [m/s]
    fig, ax = plt.subplots(figsize=(2.25, 2.3))
    pc = ax.pcolormesh(xf[: n + 1], zf, Wc[:, :n], cmap=style.CMAPS["diverging"], vmin=-0.6, vmax=0.6,
                       shading="flat", rasterized=True)
    cs = ax.contour(xc[:n], zc, T[:, :n], levels=[1.25], colors="k", linewidths=0.8)
    ax.text(62, -158, "1.25 °C", fontsize=6.5, ha="center", va="bottom")
    ax.set(xlim=(0, 80), ylim=(-500, 0), xlabel="distance from ice (m)", ylabel="depth (m)")
    cb = fig.colorbar(pc, ax=ax, shrink=0.8, pad=0.03, extend="both")
    cb.set_label("w (m/s)", labelpad=1)
    fig.tight_layout(pad=0.2)
    _save(fig, "p_plume2d.pdf")


def face_w(root: Path):
    run = root / "b_3d_twolayer_x0_80"
    require_ok(run)
    meta = json.loads((run / "case.json").read_text())
    p = meta["params"]
    snaps = _snapshots(run, "WVEL", float(p["deltaT"]), 2.0, 3.0)
    Wm = np.mean([snaps[i][..., 0] for i in sorted(snaps)], axis=0)       # (nz, ny), i = 0
    nz, ny = Wm.shape
    yf, zf = np.arange(ny + 1) * p["dy"], -np.arange(nz + 1) * p["dz"]
    fig, ax = plt.subplots(figsize=(2.5, 2.3))
    pc = ax.pcolormesh(yf, zf, Wm, cmap=style.CMAPS["diverging"], vmin=-0.6, vmax=0.6, shading="flat",
                       rasterized=True)
    for (_, j, _) in meta["source_cells"]:
        ax.plot((j + 0.5) * p["dy"], -p["H"] + 4, "k^", ms=3.5, clip_on=False)
    yc = 600.0
    ax.set(xlim=(yc - 250, yc + 250), ylim=(-p["H"], 0), xlabel="along the ice face (m)", ylabel="depth (m)")
    ax.annotate("", xy=((meta["source_cells"][1][1] + 0.5) * p["dy"], -p["H"] + 25),
                xytext=((meta["source_cells"][0][1] + 0.5) * p["dy"], -p["H"] + 25),
                arrowprops=dict(arrowstyle="<->", lw=0.7))
    ax.text(yc, -p["H"] + 45, f"x0 = {meta['x0_actual']:.0f} m", ha="center", va="bottom", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    cb = fig.colorbar(pc, ax=ax, shrink=0.8, pad=0.03, extend="both")
    cb.set_label("w (m/s)", labelpad=1)
    fig.tight_layout(pad=0.2)
    _save(fig, "p_face_w.pdf")


def profile_panel():
    z = -np.linspace(0, 500, 501)
    T2, S2 = profiles.two_layer(z, 200.0, profiles.T_PW, profiles.S_PW, profiles.T_AW, profiles.S_AW, 10.0)
    Tc, Sc = profiles.continuous_sermilik(z, z_c=200.0)
    fig, ax = plt.subplots(figsize=(1.9, 2.3))
    ax.plot(T2, z, color=style.OKABE_ITO[1], lw=1.2, label="two-layer")
    ax.plot(Tc, z, color=style.OKABE_ITO[1], lw=1.2, ls="--", label="Sermilik-like")
    ax.set(xlabel="T (°C)", ylabel="depth (m)", ylim=(-500, 0), xlim=(-2, 4.5))
    ax.xaxis.label.set_color(style.OKABE_ITO[1])
    ax2 = ax.twiny()
    ax2.plot(S2, z, color=style.OKABE_ITO[0], lw=1.2)
    ax2.plot(Sc, z, color=style.OKABE_ITO[0], lw=1.2, ls="--")
    ax2.set(xlim=(32.6, 34.9))
    ax2.set_xlabel("S (psu)", color=style.OKABE_ITO[0], labelpad=2)
    ax.text(0.9, -235, "AW", fontsize=7, ha="center")
    ax.text(0.9, -150, "PW", fontsize=7, ha="center")
    ax.legend(frameon=False, fontsize=6.5, loc="lower left", handlelength=1.6)
    fig.tight_layout(pad=0.2)
    _save(fig, "p_profiles.pdf")


def forcing_panel():
    t = np.linspace(0, 12, 1201)
    fig, axs = plt.subplots(3, 1, figsize=(2.2, 2.3), sharex=True)
    axs[0].plot(t, np.full_like(t, 100.0), color=style.OKABE_ITO[5], lw=1.2)
    axs[0].set(ylim=(0, 130), ylabel="Q0\n(m³/s)")
    axs[1].plot(t, np.where(t < 0.5, 0.2, 0.0), color="k", lw=1.2, drawstyle="steps-post")
    axs[1].set(ylim=(-0.02, 0.25), ylabel="τ\n(N/m²)")
    axs[2].plot(t, 200 - 30 * np.sin(2 * np.pi * t / 6.0), color=style.OKABE_ITO[2], lw=1.2)
    axs[2].set(ylim=(240, 160), ylabel="z_c\n(m)", xlabel="time (days)", xlim=(0, 12))
    for a in axs:
        a.yaxis.set_label_coords(-0.25, 0.5)
    fig.tight_layout(pad=0.2, h_pad=0.3)
    _save(fig, "p_forcing.pdf")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-root", type=Path, default=Path("/Volumes/LaCie/portfolio_runs/exp1"))
    a = ap.parse_args()
    style.apply()
    plt.rcParams.update(RC)
    plume2d(a.runs_root)
    face_w(a.runs_root)
    profile_panel()
    forcing_panel()
    print("wrote", ", ".join(sorted(f.name for f in OUT.glob("p_*.png"))))


if __name__ == "__main__":
    main()
