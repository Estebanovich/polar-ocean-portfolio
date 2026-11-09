#!/opt/anaconda3/envs/pyocean/bin/python
r"""GIF of the meltwater tracer of one Exp. 2 run (default: the 4 km Ocean0 reference).

Usage::

    animate_tracer.py <run_dir> [--out exp2_isomip_cavity/figures/fig3_anim_<run>.gif]
                      [--stride N] [--dpi D] [--max-mb 5] [--frame-png PATH]

Frames: the monthly-mean records of TRAC01 (C = meltwater volume fraction, pd1 sect. 4),
read through ``metrics.load_run`` (require_ok, melt_law guard, staged runs concatenated).
Panels, with FIXED colour scales over the whole GIF (max over all frames, reported):

  (a) plan view of the vertically integrated tracer sum_k C hFacC drF [m of meltwater], with
      the ice-front outline (shelf mask) in white;
  (b) x-z section of C along the centre row j = ny // 2, dry cells (ice, bed) grey.

Size: Pillow; the stride doubles until the GIF fits ``--max-mb`` (Exp. 1
``animate_b.write_within_budget``, as ``animate_c.py``); exit code 1 if still over budget.
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

from exp1_fjord_plumes.analysis.animate_b import write_within_budget  # noqa: E402
from exp2_isomip_cavity.analysis import metrics as mt  # noqa: E402
from pfcommon import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402

FIG_DIR = REPO / "exp2_isomip_cavity" / "figures"


def _frame(rd: mt.RunData, k: int):
    C = rd.field("dynDiag", "TRAC01", k)
    col = (C * rd.hfac * rd.drF[:, None, None]).sum(axis=0)
    return C, col


def make_gif(run_dir: Path, out: Path, stride: int = 1, dpi: int = 80, max_mb: float = 5.0,
             frame_png: Path | None = None) -> dict:
    rd = mt.load_run(run_dir)
    n = rd.n_records
    j0 = rd.ny // 2
    wet = rd.wet
    dry_sec = rd.hfac[:, j0, :] == 0
    vmax_col = vmax_c = 0.0
    for k in range(n):                                   # pass 1: the fixed scales
        C, col = _frame(rd, k)
        vmax_col = max(vmax_col, float(col[wet].max()))
        vmax_c = max(vmax_c, float(C[:, j0, :][~dry_sec].max()))
    vmax_col = vmax_col if vmax_col > 0 else 1.0
    vmax_c = vmax_c if vmax_c > 0 else 1.0
    xg = np.append(rd.xg[0], rd.xg[0, -1] + rd.dxG[0, -1]) / 1e3
    yg = np.append(rd.yg[:, 0], rd.yg[-1, 0] + rd.dyG[-1, 0]) / 1e3

    style.apply()
    cmap = plt.get_cmap("cmo.dense").copy()
    cmap.set_bad("0.8")
    fig = plt.figure(figsize=(7.0, 5.6), facecolor="white")
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 0.03], hspace=0.45, wspace=0.05,
                          left=0.1, right=0.9, top=0.9, bottom=0.1)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    z0 = np.ma.masked_array(np.zeros((rd.ny, rd.nx)), ~wet)
    mA = axA.pcolormesh(xg, yg, z0, cmap=cmap, vmin=0, vmax=vmax_col, shading="flat")
    axA.contour(rd.xc[0] / 1e3, rd.yc[:, 0] / 1e3, rd.shelf.astype(float), [0.5], colors="w",
                linewidths=0.7)
    axA.set(xlabel="x (km)", ylabel="y (km)",
            title="(a) vertically integrated meltwater tracer (white: ice-shelf edge)")
    fig.colorbar(mA, cax=fig.add_subplot(gs[0, 1]), label=r"$\int C\,dz$ (m)")
    s0 = np.ma.masked_array(np.zeros((rd.nz, rd.nx)), dry_sec)
    mB = axB.pcolormesh(xg, rd.rf, s0, cmap=cmap, vmin=0, vmax=vmax_c, shading="flat")
    axB.set(xlabel="x (km)", ylabel="z (m)",
            title=f"(b) C along y = {rd.yc[j0, 0] / 1e3:.0f} km (grey: ice / bed)")
    fig.colorbar(mB, cax=fig.add_subplot(gs[1, 1]), label="C (meltwater fraction)")
    fig.suptitle(f"{rd.params.get('case')}: meltwater tracer, monthly means", x=0.1,
                 ha="left", fontsize=10)
    clock = fig.text(0.97, 0.965, "", ha="right", va="top", family="monospace", fontsize=10)

    def draw(k):
        C, col = _frame(rd, k)
        mA.set_array(np.ma.masked_array(col, ~wet).ravel())
        mB.set_array(np.ma.masked_array(C[:, j0, :], dry_sec).ravel())
        clock.set_text(f"month {k + 1:3d} / {n}")

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    frames = list(range(n))
    nsel = [0]

    def write(st: int) -> float:
        sel = frames[::st]
        if sel[-1] != frames[-1]:
            sel = sel + [frames[-1]]
        w = PillowWriter(fps=4)
        with w.saving(fig, str(out), dpi=dpi):
            for k in sel:
                draw(k)
                w.grab_frame()
        nsel[0] = len(sel)
        return out.stat().st_size / 1e6

    stride, mb, over = write_within_budget(write, n, stride, max_mb)
    if frame_png is not None:
        draw(n - 1)
        fig.savefig(frame_png, dpi=dpi * 2)
    plt.close(fig)
    res = {"gif": str(out), "mb": mb, "over_budget": over, "frames": nsel[0],
           "stride": stride, "dpi": dpi, "n_records": n, "vmax_column": vmax_col,
           "vmax_C_section": vmax_c, "row_j": j0}
    print(json.dumps(res))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--dpi", type=int, default=80)
    ap.add_argument("--max-mb", type=float, default=5.0)
    ap.add_argument("--frame-png", type=Path, default=None)
    a = ap.parse_args(argv)
    out = a.out or FIG_DIR / f"fig3_anim_{a.run_dir.name}.gif"
    res = make_gif(a.run_dir, out, a.stride, a.dpi, a.max_mb, a.frame_png)
    return 1 if res["over_budget"] else 0


if __name__ == "__main__":
    sys.exit(main())
