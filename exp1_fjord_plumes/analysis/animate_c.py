"""GIF animation of one config (c) case (ruling R29).

Usage::

    animate_c.py <run_dir> [--out exp1_fjord_plumes/figures/anim/<case>.gif]
                 [--stride N] [--dpi D] [--max-mb 5] [--frame-png PATH]
                 [--check-scales]

Panels (after verification/BayIW_Rectan/make_animation.py and animate_b.py:
FIXED colour/axis scales shared by ALL four cases, land light grey, clock):

  (a) plan view of T at the plume outflow level k = 18 (z = -185 m, the cell
      of z_neutral = -181 m), with (u, v) arrows (fixed reference speed);
  (b) along-fjord section of T at the plume row j = ny//2, with the Q_ex
      section (i = 100) marked and the OBCS sponge (i >= 188) shaded;
  (c) hourly Q_ex(t) at i = 100 (``c_runs.exchange_series``, R31) with a
      moving marker at the frame time.

Frames: the 3-h means ``mean3h`` (UVEL, THETA, full 3-D) and ``animV``
(VVEL at level 19, 1-based; data.diagnostics list 4, added for this script).
3-h means rather than snapshots: smooth enough for a 30-day GIF (240 frames)
and already needed for Fig. 3, so only animV is extra (~12 MB per case).
Two-stage cases are read from both stage directories (``c_runs.case_parts``,
which also refuses runs not OK). Frame time = centre of the 3-h window.

GIF size: Pillow; the stride doubles until the GIF fits ``--max-mb``
(``animate_b.write_within_budget``); exit code 1 if it is still over budget. ``--check-scales`` writes no GIF; it
prints the fraction of wet pixels outside the fixed T scale in (a) and (b)
(max and mean over frames), max Q_ex / QEX_LIM and the max speed in (a).
Run it on the first finished case and widen the limits for all cases alike.
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
from exp1_fjord_plumes.analysis import c_runs  # noqa: E402
from exp1_fjord_plumes.analysis.animate_b import write_within_budget  # noqa: E402
from exp1_fjord_plumes.analysis.fig3_exchange_partition import I_SECTION, T_REF  # noqa: E402
from pfcommon import style  # noqa: E402
from pfcommon.mitgcm_io import read_mds  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402

# Fixed scales for ALL (c) cases (set a priori; check with --check-scales)
T_LIM = (-1.6, 4.1)      # degC: PW -1.5 .. AW 4.0 (Sciascia 2013 layer values)
QEX_LIM = 1.2e5          # m3/s: slab estimate of the shelf-forced layer flux ~8e4
U_REF = 0.2              # m/s, arrow reference speed
K_PLAN = 18              # 0-based outflow level (z = -185 m)
DEFAULT_OUT = REPO / "exp1_fjord_plumes" / "figures" / "anim"


def frame_index(run_dir: Path) -> list[tuple[int, Path]]:
    return c_runs.mds_iters(c_runs.case_parts(run_dir), "mean3h")


def _load(d: Path, it: int, k_plan: int):
    a = read_mds(d / "mean3h", it)[0]                  # (2, nz, ny, nx)
    u, T = a[0], a[1]
    v = None
    if (d / f"animV.{it:010d}.meta").exists():
        v = read_mds(d / "animV", it)[0][0]
        v = v[0] if v.ndim == 3 else v                 # (ny, nx)
    return u, T, v


def _out_of_scale(T, mask) -> float:
    t = np.asarray(T)[np.asarray(mask, bool)]
    return float(np.mean((t < T_LIM[0]) | (t > T_LIM[1]))) if t.size else 0.0


def check_scales(run_dir: Path, i_section: int = I_SECTION, k_plan: int = K_PLAN) -> dict:
    run_dir = Path(run_dir)
    p = c_runs.params(run_dir)["params"]
    dep = c_runs.depth(run_dir)
    wet2 = dep > 0
    zb = (np.arange(p["nz"]) + 1) * p["dz"]
    j0 = p["ny"] // 2
    wet_sec = zb[:, None] <= dep[None, j0, :] + 1e-9
    a, b, spd = [], [], []
    for it, d in frame_index(run_dir):
        u, T, v = _load(d, it, k_plan)
        a.append(_out_of_scale(T[k_plan], wet2))
        b.append(_out_of_scale(T[:, j0, :], wet_sec))
        uu = u[k_plan][wet2]
        vv = v[wet2] if v is not None else 0.0
        spd.append(float(np.max(np.hypot(uu, vv))) if uu.size else 0.0)
    q = c_runs.exchange_series(run_dir, i_section, T_REF)["Q_ex"]   # hourly
    return {"a_T_plan": {"max": max(a), "mean": float(np.mean(a)), "frames": len(a)},
            "b_T_section": {"max": max(b), "mean": float(np.mean(b)), "frames": len(b)},
            "c_Q_ex": {"max_m3_s": float(q.max()), "max_over_lim": float(q.max() / QEX_LIM)},
            "speed_plan_max_m_s": max(spd),
            "limits": {"T_LIM": list(T_LIM), "QEX_LIM": QEX_LIM, "U_REF": U_REF}}


def make_gif(run_dir: Path, out: Path, stride: int = 1, dpi: int = 80, max_mb: float = 5.0,
             frame_png: Path | None = None, i_section: int = I_SECTION,
             k_plan: int = K_PLAN) -> dict:
    run_dir = Path(run_dir)
    frames = frame_index(run_dir)
    if not frames:
        raise FileNotFoundError(f"no mean3h output in {run_dir}")
    p = c_runs.params(run_dir)["params"]
    nx, ny, nz, dx, dy, dz, dt = (p["nx"], p["ny"], p["nz"], p["dx"], p["dy"], p["dz"],
                                  p["deltaT"])
    dep = c_runs.depth(run_dir)
    land = dep <= 0
    j0 = ny // 2
    zb = (np.arange(nz) + 1) * dz
    land_sec = ~(zb[:, None] <= dep[None, j0, :] + 1e-9)
    xf, yf, zf = np.arange(nx + 1) * dx / 1e3, np.arange(ny + 1) * dy / 1e3, -np.arange(nz + 1) * dz
    xc, yc = (np.arange(nx) + 0.5) * dx / 1e3, (np.arange(ny) + 0.5) * dy / 1e3
    qs = c_runs.exchange_series(run_dir, i_section, T_REF)
    t_end = p["endTime"]
    sx, sy = max(1, nx // 25), max(1, ny // 8)

    style.apply()
    cmap = style.CMAPS["temperature"].copy()
    cmap.set_bad("0.82")
    fig = plt.figure(figsize=(10.0, 7.0), facecolor="white")
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 0.03], height_ratios=[0.8, 1.2, 0.8],
                          hspace=0.55, wspace=0.05, left=0.08, right=0.92, top=0.9, bottom=0.13)
    axA, axB, axC = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[2, 0])
    cax = fig.add_subplot(gs[0:2, 1])
    zA = np.ma.masked_array(np.zeros((ny, nx)), land)
    mA = axA.pcolormesh(xf, yf, zA, cmap=cmap, vmin=T_LIM[0], vmax=T_LIM[1], shading="flat")
    Q = axA.quiver(xc[::sx], yc[::sy], np.zeros((len(yc[::sy]), len(xc[::sx]))),
                   np.zeros((len(yc[::sy]), len(xc[::sx]))), scale=U_REF * 12, width=0.002)
    axA.quiverkey(Q, 0.9, 1.08, U_REF, f"{U_REF} m/s", labelpos="E")
    axA.set(xlabel="x from ice face (km)", ylabel="y (km)",
            title=f"(a) T at z = {-(k_plan + 0.5) * dz:.0f} m (plume outflow level)")
    zB = np.ma.masked_array(np.zeros((nz, nx)), land_sec)
    mB = axB.pcolormesh(xf, zf, zB, cmap=cmap, vmin=T_LIM[0], vmax=T_LIM[1], shading="flat")
    xs = i_section * dx / 1e3
    for ax in (axA, axB):
        ax.axvline(xs, color="k", lw=1.0, ls="--")
        ax.axvspan((nx - 2 - c_runs.SPONGE) * dx / 1e3, xf[-1], color="0.5", alpha=0.25, lw=0)
    axB.set(xlabel="x from ice face (km)", ylabel="z (m)",
            title=f"(b) T along the fjord at y = {yc[j0]:.1f} km (plume row); "
                  "dashed: Q_ex section, grey: OBCS sponge")
    fig.colorbar(mA, cax=cax, extend="both", label="T (°C)")
    axC.plot(qs["t_hourly_s"] / 86400, qs["Q_ex"] / 1e3, color=style.OKABE_ITO[0], lw=1.0)
    mk, = axC.plot([], [], "o", color=style.OKABE_ITO[1], ms=6)
    axC.set(xlim=(0, t_end / 86400), ylim=(0, QEX_LIM / 1e3), xlabel="time (days)",
            ylabel=r"$Q_{ex}$ (10$^3$ m$^3$/s)", title=f"(c) exchange flux at x = {xs:.0f} km")
    fig.suptitle(f"{p['case']}", x=0.06, ha="left", fontsize=11)
    # Boundary types stated explicitly (as in animate_b): only the fjord mouth is open.
    fig.text(0.08, 0.008,
             "Boundaries: ice face (x = 0) and fjord side walls are closed (no normal flow), so waves "
             "reflect there by design.\nThe only open boundary is the mouth (x ≈ %.0f km): Stevens (1990) "
             "radiation with interior velocity, plus the grey sponge." % (nx * dx / 1e3),
             ha="left", va="bottom", fontsize=7.5, color="0.25")
    clock = fig.text(0.97, 0.965, "", ha="right", va="top", family="monospace", fontsize=11)

    def draw(fr):
        it, d = fr
        u, T, v = _load(d, it, k_plan)
        mA.set_array(np.ma.masked_array(T[k_plan], land).ravel())
        mB.set_array(np.ma.masked_array(T[:, j0, :], land_sec).ravel())
        uc = u[k_plan][::sy, ::sx]
        vc = (v if v is not None else np.zeros((ny, nx)))[::sy, ::sx]
        msk = land[::sy, ::sx]
        Q.set_UVC(np.where(msk, 0, uc), np.where(msk, 0, vc))
        ts = it * dt - c_runs.MEAN_S / 2
        mk.set_data([ts / 86400], [np.interp(ts, qs["t_hourly_s"], qs["Q_ex"]) / 1e3])
        clock.set_text(f"t = {ts / 86400:5.2f} d (3-h mean)")

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    nsel = [0]

    def write(st: int) -> float:
        sel = frames[::st]
        if sel[-1][0] != frames[-1][0]:
            sel = sel + [frames[-1]]
        w = PillowWriter(fps=8)
        with w.saving(fig, str(out), dpi=dpi):
            for f in sel:
                draw(f)
                w.grab_frame()
        nsel[0] = len(sel)
        return out.stat().st_size / 1e6

    stride, mb, over = write_within_budget(write, len(frames), stride, max_mb)
    if frame_png is not None:
        draw(frames[-1])
        fig.savefig(frame_png, dpi=dpi * 2)
    plt.close(fig)
    res = {"gif": str(out), "mb": mb, "over_budget": over, "frames": nsel[0], "stride": stride,
           "dpi": dpi, "t_first_s": frames[0][0] * dt - c_runs.MEAN_S / 2,
           "t_last_s": frames[-1][0] * dt - c_runs.MEAN_S / 2}
    print(json.dumps(res))
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--dpi", type=int, default=80)
    ap.add_argument("--max-mb", type=float, default=5.0)
    ap.add_argument("--frame-png", type=Path, default=None)
    ap.add_argument("--check-scales", action="store_true",
                    help="only report saturation of the fixed scales (all frames); no GIF")
    a = ap.parse_args(argv)
    if a.check_scales:
        print(json.dumps(check_scales(a.run_dir), indent=2))
        return 0
    case = c_runs.params(a.run_dir)["params"]["case"]
    res = make_gif(a.run_dir, a.out or DEFAULT_OUT / f"{case}.gif", a.stride, a.dpi,
                   a.max_mb, a.frame_png)
    return 1 if res.get("over_budget") else 0


if __name__ == "__main__":
    sys.exit(main())
