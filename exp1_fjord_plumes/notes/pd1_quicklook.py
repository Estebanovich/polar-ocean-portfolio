"""PD1 quick look (Task 4): T at the end of the 2 h short run with w contours.

Usage: pd1_quicklook.py <run_dir> <out_stem>
Reads the final state dumps T/W (dumpInitAndLast; the hourly diagnostics
snapshots stop one step before the end) and delX/delZ from <run_dir>,
refuses runs that did not finish OK (pfcommon.runctl.require_ok).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pfcommon import style  # noqa: E402
from pfcommon.mitgcm_io import read_bin, read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

XMAX = 1500.0  # m shown from the ice face


def main(run_dir: Path, stem: Path) -> None:
    require_ok(run_dir)
    T, mt = read_mds(run_dir / "T")
    W, _ = read_mds(run_dir / "W")
    T, W = T[0], W[0]  # (nz, nx)
    nz, nx = T.shape
    dx = read_bin(run_dir / "delX.bin", (nx,))
    dz = read_bin(run_dir / "delZ.bin", (nz,))
    xf = np.concatenate([[0.0], np.cumsum(dx)])
    zf = -np.concatenate([[0.0], np.cumsum(dz)])
    xc, zc = 0.5 * (xf[1:] + xf[:-1]), 0.5 * (zf[1:] + zf[:-1])
    # WVEL sits on the top face of each cell; interpolate to centres for contouring
    Wc = 0.5 * (W + np.vstack([W[1:], np.zeros((1, nx))]))
    n = int(np.searchsorted(xf, XMAX))
    dt = json.loads((run_dir / "case.json").read_text())["params"]["deltaT"]
    t_h = mt["timeStepNumber"] * dt / 3600.0

    style.apply()
    fig, ax = plt.subplots()
    pc = ax.pcolormesh(xf[: n + 1], zf, T[:, :n], cmap=style.CMAPS["temperature"], shading="flat")
    fig.colorbar(pc, ax=ax, label="T (°C)")
    lv = [-0.2, -0.1, -0.05, -0.02, 0.02, 0.05, 0.1, 0.2]
    cs = ax.contour(xc[:n], zc, Wc[:, :n], levels=lv, colors="k", linewidths=0.6)
    ax.clabel(cs, fmt="%.2f", fontsize=6)
    ax.set_xlabel("distance from ice face (m)")
    ax.set_ylabel("z (m)")
    ax.set_title(f"{run_dir.name}: T and w (m/s, contours) at t = {t_h:.1f} h")
    style.save(fig, stem)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
