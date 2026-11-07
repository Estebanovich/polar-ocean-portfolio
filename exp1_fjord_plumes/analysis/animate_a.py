"""GIF animation of one configuration (a) case (2-D line plume against the ice face; ruling R29).

Frames are the run's THETA/WVEL snapshots (production runs: hourly over 48 h; the dedicated
animation rerun: every 5 min over 6 h) with the hour-mean ice-front melt (icefrontDiag) ending at
each snapshot. Panels:
  (a) T in the x-z section next to the ice (first X_SHOW metres), with w contours (0.1, 0.3 m/s);
  (b) w in the cell touching the ice against depth, with the detachment depth (plume_diag.detach_depth,
      5 % of the maximum);
  (c) ice melt rate against depth (m/day; ICFfwFlx / rho_ice, plume_diag.melt_rate_m_per_s).
Colour and axis limits are FIXED (thesis make_animation.py convention, R29) and checked with
--check-scales, which reports the fraction of saturated pixels / clipped values; scales are not
changed silently. GIFs are kept <= --max-mb (stride doubled until they fit, animate_b.write_within_budget).

Usage: animate_a.py RUN_DIR [--out GIF] [--stride N] [--dpi 80] [--max-mb 8] [--frame-png PNG] [--check-scales]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from exp1_fjord_plumes.analysis.animate_b import saturated_fraction, write_within_budget  # noqa: E402
from exp1_fjord_plumes.analysis.plume_diag import detach_depth, melt_rate_m_per_s  # noqa: E402
from pfcommon import style  # noqa: E402
from pfcommon.mitgcm_io import read_bin, read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "exp1_fjord_plumes" / "figures"
X_SHOW = 200.0          # m from the ice face shown in (a)
T_LIM = (-1.6, 4.1)     # degC, as animate_c
W_LIM = 0.8             # m/s, x-limit of (b)
MELT_LIM = 1.5          # m/day, x-limit of (c) (max 1.0 in the q0 = 0.03 runs, --check-scales)
DAY = 86400.0


def _iters(run_dir: Path, name: str) -> list[int]:
    return sorted(int(p.name.split(".")[1]) for p in run_dir.glob(f"{name}.[0-9]*.meta"))


def frames(run_dir: Path) -> list[int]:
    """Iterations with THETA and WVEL snapshots (t > 0)."""
    w = set(_iters(run_dir, "WVEL"))
    return [i for i in _iters(run_dir, "THETA") if i in w and i > 0]


def _grid(run_dir: Path):
    p = json.loads((run_dir / "case.json").read_text())["params"]
    nx, nz = p["nx"], p["nz"]
    dx = read_bin(run_dir / "delX.bin", (nx,))
    xf = np.concatenate([[0.0], np.cumsum(dx)])
    zf = -np.arange(nz + 1) * p["dz"]
    return p, xf, zf


def _load(run_dir: Path, it: int):
    T = read_mds(run_dir / "THETA", it)[0][0]               # (nz, nx)
    W = read_mds(run_dir / "WVEL", it)[0][0]                # w at the top face of each cell
    melt = None
    if (run_dir / f"icefrontDiag.{it:010d}.meta").exists():
        fw = read_mds(run_dir / "icefrontDiag", it)[0][0]    # record ICFfwFlx, (nz, nx)
        melt = melt_rate_m_per_s(fw[:, 0]) * DAY             # face column, m/day
    return T, W, melt


def check_scales(run_dir: Path) -> dict:
    p, xf, zf = _grid(run_dir)
    n = int(np.searchsorted(xf, X_SHOW))
    sat_T, wmax, mmax = [], 0.0, 0.0
    for it in frames(run_dir):
        T, W, melt = _load(run_dir, it)
        sat_T.append(saturated_fraction(T[:, :n] - 0.5 * (T_LIM[0] + T_LIM[1]), 0.5 * (T_LIM[1] - T_LIM[0]),
                                        np.ones_like(T[:, :n], bool)))
        wmax = max(wmax, float(W[:, 0].max()))
        if melt is not None:
            mmax = max(mmax, float(np.nanmax(melt)))
    return {"a_T_saturated": {"max": max(sat_T), "mean": float(np.mean(sat_T)), "frames": len(sat_T)},
            "b_w_face_max_m_s": wmax, "b_over_lim": wmax / W_LIM,
            "c_melt_max_m_day": mmax, "c_over_lim": mmax / MELT_LIM,
            "limits": {"T_LIM": list(T_LIM), "W_LIM": W_LIM, "MELT_LIM": MELT_LIM}}


def make_gif(run_dir: Path, out: Path, stride: int = 1, dpi: int = 80, max_mb: float = 8.0,
             frame_png: Path | None = None) -> dict:
    run_dir = Path(run_dir)
    require_ok(run_dir)
    p, xf, zf = _grid(run_dir)
    its = frames(run_dir)
    if not its:
        raise FileNotFoundError(f"no THETA/WVEL snapshots in {run_dir}")
    dt, nz, nx = p["deltaT"], p["nz"], p["nx"]
    n = int(np.searchsorted(xf, X_SHOW))
    xc, zc = 0.5 * (xf[1:] + xf[:-1]), 0.5 * (zf[1:] + zf[:-1])
    strat = {"twolayer": "two-layer", "sermilik": "Sermilik-like"}.get(p.get("strat"), p.get("strat"))

    style.apply()
    cmap = style.CMAPS["temperature"].copy()
    fig = plt.figure(figsize=(10.0, 5.4), facecolor="white")
    gs = fig.add_gridspec(2, 4, width_ratios=[2.4, 0.06, 1, 1], height_ratios=[1, 0.09], wspace=0.35,
                          hspace=0.45, left=0.07, right=0.97, top=0.86, bottom=0.17)
    axT, cax, axW, axM = (fig.add_subplot(gs[0, k]) for k in range(4))
    axC = fig.add_subplot(gs[1, :])
    T0, W0, _ = _load(run_dir, its[0])
    mT = axT.pcolormesh(xf[: n + 1], zf, T0[:, :n], cmap=cmap, vmin=T_LIM[0], vmax=T_LIM[1], shading="flat")
    fig.colorbar(mT, cax=cax, extend="both", label="T (°C)")
    axT.set(xlim=(0, X_SHOW), ylim=(zf[-1], 0), xlabel="distance from the ice face (m)", ylabel="depth (m)",
            title="(a) temperature, with w contours 0.1 and 0.3 m/s")
    axT.axvline(0, color="#6b8fa3", lw=3)
    cont = [None]
    lw_, = axW.plot([], [], color=style.OKABE_ITO[0], lw=1.6)
    zdl = axW.axhline(np.nan, color="k", ls="--", lw=0.9)
    axW.set(xlim=(0, W_LIM), ylim=(zf[-1], 0), xlabel="w at the ice (m/s)", title="(b) plume velocity")
    axW.text(0.97, 0.03, "dashed: detachment", transform=axW.transAxes, ha="right", fontsize=8)
    lm_, = axM.plot([], [], color=style.OKABE_ITO[1], lw=1.6)
    axM.set(xlim=(0, MELT_LIM), ylim=(zf[-1], 0), xlabel="melt rate (m/day)", title="(c) ice melt (1-h mean)")
    for a in (axW, axM):
        a.set_yticklabels([])
    t_end = p["endTime"]
    axC.set(xlim=(0, t_end / 3600.0), ylim=(0, 1), yticks=[], xlabel="model time (h)")
    for s_ in ("top", "right", "left"):
        axC.spines[s_].set_visible(False)
    bar = axC.barh(0.5, 0.0, height=0.6, color=style.OKABE_ITO[0])[0]
    fig.suptitle(f"{p['case']}: line plume, q0 = {p['q0']} m²/s, {strat} ambient", x=0.07, ha="left", fontsize=11)
    clock = fig.text(0.97, 0.955, "", ha="right", va="top", family="monospace", fontsize=11)
    fig.text(0.07, 0.01, "Boundaries: the ice face (x = 0) is closed and melts; the bottom is closed and the top is a free "
             "surface.\nThe far end (x = 20 km, off the panel) is open: Orlanski radiation with a velocity sponge. "
             "2-D: no variation along the ice face.",
             ha="left", va="bottom", fontsize=7.5, color="0.25")

    def draw(it):
        T, W, melt = _load(run_dir, it)
        mT.set_array(T[:, :n].ravel())
        if cont[0] is not None:
            cont[0].remove()
        Wc = 0.5 * (W + np.vstack([W[1:], np.zeros((1, nx))]))
        cont[0] = axT.contour(xc[:n], zc, Wc[:, :n], levels=[0.1, 0.3], colors="w", linewidths=[0.6, 1.2])
        lw_.set_data(W[:, 0], zf[:-1])
        zd = detach_depth(zf[:-1], W[:, 0])
        zdl.set_ydata([zd, zd])
        if melt is not None:
            lm_.set_data(melt, zc)
        ts = it * dt
        bar.set_width(ts / 3600.0)
        hh, rem = divmod(int(round(ts)), 3600)
        clock.set_text(f"t = {hh:d}:{rem // 60:02d} h")

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    nsel = [0]

    def write(st: int) -> float:
        sel = its[::st]
        if sel[-1] != its[-1]:
            sel = sel + [its[-1]]
        writer = PillowWriter(fps=6)
        with writer.saving(fig, str(out), dpi=dpi):
            for it in sel:
                draw(it)
                writer.grab_frame()
        nsel[0] = len(sel)
        return out.stat().st_size / 1e6

    stride, mb, over = write_within_budget(write, len(its), stride, max_mb)
    if frame_png is not None:
        draw(its[-1])
        fig.savefig(frame_png, dpi=dpi * 2)
    plt.close(fig)
    res = {"gif": str(out), "mb": mb, "over_budget": over, "frames": nsel[0], "stride": stride,
           "t_first_h": its[0] * dt / 3600, "t_last_h": its[-1] * dt / 3600}
    print(json.dumps(res))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--dpi", type=int, default=80)
    ap.add_argument("--max-mb", type=float, default=8.0)
    ap.add_argument("--frame-png", type=Path, default=None)
    ap.add_argument("--check-scales", action="store_true")
    a = ap.parse_args(argv)
    if a.check_scales:
        print(json.dumps(check_scales(a.run_dir), indent=2))
        return 0
    case = json.loads((a.run_dir / "case.json").read_text())["params"]["case"]
    make_gif(a.run_dir, a.out or DEFAULT_OUT / f"fig1_anim_{case}.gif", a.stride, a.dpi, a.max_mb, a.frame_png)
    return 0


if __name__ == "__main__":
    sys.exit(main())
