"""GIF animation of one config (b) case (ruling R29).

Usage::

    animate_b.py <run_dir> [--out exp1_fjord_plumes/figures/anim/<case>.gif]
                 [--stride N] [--dpi D] [--max-mb 8] [--frame-png PATH] [--mp4]

Panels (after the thesis script verification/BayIW_Rectan/make_animation.py:
fixed colour scales shared by ALL cases so they can be compared frame by
frame, cmocean maps, land in light grey, a clock with the model time):

  (a) w on the ice face, the y-z section of the face column i = 0: the two
      plumes rising from the sources (triangles) and their merging;
  (b) T anomaly T - T_init(z) on a horizontal slice in the outflow layer
      (z = -75 m);
  (c) the same at mid-depth (z = -155 m).
  The anomaly (not T itself) is shown because at mid-depth the plume is 1-2
  cells wide and only ~0.2 degC colder than the ambient, invisible on a
  -1.6..4.1 degC scale (checked on the PD1 short run); it also makes the two
  stratifications comparable. Each panel has its own FIXED scale;
  (d) clock: model time and progress bar over endTime.

Frames: the R29 animation diagnostics ``animW`` (3-D WVEL) and ``animT``
(THETA at levels 8 and 16, 1-based) every 120 s. If a run has none (runs made
before R29, e.g. the 30-min short run), the 3-D snapshots ``WVEL``/``THETA``
and the state dumps ``W``/``T`` are used instead, whatever times exist.
The frame time is iteration * deltaT (diagnostics snapshots are written one
step late, pkg/diagnostics/diagnostics_write.F:62-72, so the label can be one
step early; irrelevant at 120 s).

Size: the GIF is written with Pillow; if it exceeds ``--max-mb`` the stride is
doubled and the GIF rewritten (the stride used is printed); a WARNING is
printed if it still exceeds the budget at the largest stride.

``--check-scales`` writes no GIF: it prints, per panel, the fraction of wet
pixels at or beyond the fixed colour limit (max and mean over all frames). The
limits were set from the 30-min short run; run this on the first 3 h case and
widen W_LIM / DT_LIM (for all cases alike) if saturation is large. Refuses runs that
did not finish OK (pfcommon.runctl.require_ok).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from pfcommon import style  # noqa: E402
from pfcommon.mitgcm_io import read_bin, read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FFMpegWriter, PillowWriter  # noqa: E402

# Fixed scales for ALL (b) cases (PD1 short run, x0 = 40 m at 30 min: w up to
# 1.28 m/s; T - T_init up to +5.4 degC at -75 m and -0.19 degC at -155 m)
W_LIM = 1.2              # m/s, symmetric
DT_LIM = (4.0, 0.2)      # degC, symmetric, panels (b) and (c)
K_T = (7, 15)            # 0-based levels of the T slices (animT levels 8, 16)
X_SHOW = 500.0           # m from the ice face shown in the horizontal slices
FFMPEG = "/opt/homebrew/bin/ffmpeg"
DEFAULT_OUT = REPO / "exp1_fjord_plumes" / "figures" / "anim"


def _iters(run_dir: Path, name: str) -> list[int]:
    out = []
    for f in run_dir.glob(f"{name}.[0-9]*.meta"):
        m = re.fullmatch(rf"{re.escape(name)}\.(\d+)\.meta", f.name)
        if m:
            out.append(int(m.group(1)))
    return sorted(out)


def frame_index(run_dir: Path) -> list[tuple[int, dict]]:
    """[(iteration, {"W": (file, it), "T": (file, it, levels)})] in time order."""
    run_dir = Path(run_dir)
    w_it, t_it = _iters(run_dir, "animW"), _iters(run_dir, "animT")
    if w_it and t_it:
        common = sorted(set(w_it) & set(t_it))
        return [(i, {"W": ("animW", i), "T": ("animT", i, (0, 1))}) for i in common]
    frames: dict[int, dict] = {}
    for wname, tname in (("WVEL", "THETA"), ("W", "T")):
        for i in sorted(set(_iters(run_dir, wname)) & set(_iters(run_dir, tname))):
            frames.setdefault(i, {"W": (wname, i), "T": (tname, i, K_T)})
    return sorted(frames.items())


def _load(run_dir: Path, spec: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    w = read_mds(run_dir / spec["W"][0], spec["W"][1])[0][0]  # (nz, ny, nx)
    tname, it, lev = spec["T"]
    t = read_mds(run_dir / tname, it)[0][0]                    # (nlev or nz, ny, nx)
    return w[:, :, 0], t[lev[0]], t[lev[1]]


def write_within_budget(write, n_frames: int, stride: int, max_mb: float):
    """Call ``write(stride) -> MB`` doubling the stride until the GIF fits ``max_mb``.

    Returns (stride, mb, over_budget). If even one frame in ``n_frames`` is too
    big, a WARNING goes to stderr and over_budget is True.
    """
    while True:
        mb = write(stride)
        if mb <= max_mb:
            return stride, mb, False
        if stride >= n_frames:
            print(f"WARNING: GIF is {mb:.1f} MB > {max_mb} MB even at stride {stride}; "
                  "lower --dpi", file=sys.stderr)
            return stride, mb, True
        stride *= 2


def saturated_fraction(values: np.ndarray, lim: float, mask: np.ndarray) -> float:
    """Fraction of the pixels where ``mask`` is True with |value| >= lim (colour saturated)."""
    v = np.asarray(values)[np.asarray(mask, bool)]
    return float(np.mean(np.abs(v) >= lim)) if v.size else 0.0


def check_scales(run_dir: Path) -> dict:
    """Saturated-pixel fraction per panel over ALL frames (max and mean over frames)."""
    run_dir = Path(run_dir)
    require_ok(run_dir)
    p = json.loads((run_dir / "case.json").read_text())["params"]
    nx, ny, nz, dx = p["nx"], p["ny"], p["nz"], p["dx"]
    wet = read_bin(run_dir / "topog.bin", (ny, nx)) < 0
    nxs = int(round(X_SHOW / dx))
    T0 = read_bin(run_dir / "T.init.bin", (nz, ny, nx))
    face = np.broadcast_to(wet[:, 0], (nz, ny))
    fr: dict[str, list[float]] = {"a_w_face": [], "b_dT_upper": [], "c_dT_mid": []}
    for _, spec in frame_index(run_dir):
        w, t1, t2 = _load(run_dir, spec)
        fr["a_w_face"].append(saturated_fraction(w, W_LIM, face))
        for key, t, k, lim in (("b_dT_upper", t1, K_T[0], DT_LIM[0]),
                               ("c_dT_mid", t2, K_T[1], DT_LIM[1])):
            fr[key].append(saturated_fraction((t - T0[k])[:, :nxs], lim, wet[:, :nxs]))
    res = {k: {"max": max(v), "mean": float(np.mean(v)), "frames": len(v)} for k, v in fr.items()}
    res["limits"] = {"W_LIM": W_LIM, "DT_LIM": list(DT_LIM)}
    return res


def make_gif(run_dir: Path, out: Path, stride: int = 1, dpi: int = 80,
             max_mb: float = 8.0, frame_png: Path | None = None, mp4: bool = False) -> dict:
    run_dir = Path(run_dir)
    require_ok(run_dir)
    meta = json.loads((run_dir / "case.json").read_text())
    p = meta["params"]
    nx, ny, nz, dx, dy, dz, dt = (p["nx"], p["ny"], p["nz"], p["dx"], p["dy"], p["dz"],
                                  p["deltaT"])
    frames = frame_index(run_dir)
    if not frames:
        raise FileNotFoundError(f"no W/T snapshots in {run_dir}")
    depth = -read_bin(run_dir / "topog.bin", (ny, nx))
    land = depth <= 0
    yf, xf, zf = np.arange(ny + 1) * dy, np.arange(nx + 1) * dx, -np.arange(nz + 1) * dz
    nxs = int(round(X_SHOW / dx))
    face_land = np.broadcast_to(land[:, 0], (nz, ny))
    src = sorted({c[1] for c in meta["source_cells"]})
    x0 = meta.get("x0_actual")
    label = "single source (Q0)" if x0 is None else f"two sources, x0 = {x0:.0f} m"
    if x0 == 0:
        label = "x0 = 0: one source of 2 Q0"

    style.apply()
    cmapW = style.CMAPS["diverging"].copy()
    cmapT = style.CMAPS["diverging"].copy()
    T0 = read_bin(run_dir / "T.init.bin", (nz, ny, nx))
    for c in (cmapW, cmapT):
        c.set_bad("0.82")  # land light grey (thesis convention)
    fig = plt.figure(figsize=(10.0, 5.6), facecolor="white")
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1, 1], height_ratios=[1, 0.13],
                          hspace=0.42, wspace=0.55, left=0.07, right=0.93, top=0.86, bottom=0.08)
    axW, axT1, axT2 = (fig.add_subplot(gs[0, i]) for i in range(3))
    axC = fig.add_subplot(gs[1, :])
    z0 = np.ma.masked_array(np.zeros((nz, ny)), face_land)
    mW = axW.pcolormesh(yf, zf, z0, cmap=cmapW, vmin=-W_LIM, vmax=W_LIM, shading="flat")
    fig.colorbar(mW, ax=axW, extend="both", label="w (m/s)", shrink=0.9)
    for j in src:
        axW.plot((j + 0.5) * dy, zf[-1], "k^", ms=6, clip_on=False)
    axW.set(xlabel="y along the ice face (m)", ylabel="z (m)", title="(a) w on the ice face")
    hz = np.ma.masked_array(np.zeros((ny, nxs)), land[:, :nxs])
    mT = []
    for ax, k, lim, tag in zip((axT1, axT2), K_T, DT_LIM, "bc"):
        m = ax.pcolormesh(xf[: nxs + 1], yf, hz, cmap=cmapT, vmin=-lim, vmax=lim, shading="flat")
        ax.set(xlabel="distance from ice face (m)", ylabel="y (m)",
               title=f"({tag}) T − T_init at z = {-(k + 0.5) * dz:.0f} m")
        for j in src:
            ax.plot(0, (j + 0.5) * dy, "k>", ms=5, clip_on=False)
        fig.colorbar(m, ax=ax, extend="both", label="ΔT (°C)", shrink=0.9)
        mT.append(m)
    t_end = p["endTime"]
    axC.set_xlim(0, t_end / 60.0)
    axC.set_ylim(0, 1)
    axC.set_yticks([])
    axC.set_xlabel("model time (min)")
    for s in ("top", "right", "left"):
        axC.spines[s].set_visible(False)
    bar = axC.barh(0.5, 0.0, height=0.6, color=style.OKABE_ITO[0])[0]
    fig.suptitle(f"{p['case']}: {label}, {p['strat']} stratification", x=0.06, ha="left",
                 fontsize=11)
    clock = fig.text(0.97, 0.955, "", ha="right", va="top", family="monospace", fontsize=11)

    def draw(spec_it):
        it, spec = spec_it
        w, t1, t2 = _load(run_dir, spec)
        mW.set_array(np.ma.masked_array(w, face_land).ravel())
        for m, t, k in zip(mT, (t1, t2), K_T):
            m.set_array(np.ma.masked_array((t - T0[k])[:, :nxs], land[:, :nxs]).ravel())
        ts = it * dt
        bar.set_width(ts / 60.0)
        hh, rem = divmod(int(round(ts)), 3600)
        clock.set_text(f"t = {hh:d}:{rem // 60:02d}:{rem % 60:02d}  ({ts / 60:5.1f} min)")

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    nsel = [0]

    def write(st: int) -> float:
        sel = frames[::st]
        if sel[-1][0] != frames[-1][0]:
            sel = sel + [frames[-1]]  # always end on the last state
        writer = PillowWriter(fps=8)
        with writer.saving(fig, str(out), dpi=dpi):
            for f in sel:
                draw(f)
                writer.grab_frame()
        nsel[0] = len(sel)
        return out.stat().st_size / 1e6

    stride, mb, over = write_within_budget(write, len(frames), stride, max_mb)
    if frame_png is not None:
        draw(frames[-1])
        fig.savefig(frame_png, dpi=dpi * 2)
    if mp4:
        plt.rcParams["animation.ffmpeg_path"] = FFMPEG  # before the writer is built
        mw = FFMpegWriter(fps=12, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "23"])
        with mw.saving(fig, str(out.with_suffix(".mp4")), dpi=120):
            for f in frames:
                draw(f)
                mw.grab_frame()
    plt.close(fig)
    res = {"gif": str(out), "mb": mb, "over_budget": over, "frames": nsel[0], "stride": stride,
           "dpi": dpi,
           "t_first_s": frames[0][0] * dt, "t_last_s": frames[-1][0] * dt}
    print(json.dumps(res))
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--dpi", type=int, default=80)
    ap.add_argument("--max-mb", type=float, default=8.0)
    ap.add_argument("--frame-png", type=Path, default=None)
    ap.add_argument("--mp4", action="store_true")
    ap.add_argument("--check-scales", action="store_true",
                    help="only report the saturated-pixel fraction per panel (all frames); no GIF")
    a = ap.parse_args(argv)
    if a.check_scales:
        print(json.dumps(check_scales(a.run_dir), indent=2))
        return 0
    case = json.loads((a.run_dir / "case.json").read_text())["params"]["case"]
    out = a.out or DEFAULT_OUT / f"{case}.gif"
    make_gif(a.run_dir, out, a.stride, a.dpi, a.max_mb, a.frame_png, a.mp4)
    return 0


if __name__ == "__main__":
    sys.exit(main())
