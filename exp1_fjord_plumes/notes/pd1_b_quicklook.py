"""PD1 quick look for config (b) (Task 6 Step 2): w sections at the end of the short run.

Usage: pd1_b_quicklook.py <run_dir> <out_stem>

Reads the final state dump W (dumpInitAndLast; diagnostics snapshots are
written one step late, so none exists at endTime) and case.json from
<run_dir>; refuses runs that did not finish OK (pfcommon.runctl.require_ok).
Panels: (1) y-z section of w in the face column i = 0 (through both sources),
(2) x-z section of w through the first source row, (3) x-y slice of w at
mid-depth of the lower layer. Prints a grid-scale noise indicator

    G = rms(w[n+1] - 2 w[n] + w[n-1]) / (4 rms(w))      [dimensionless]

along x, y and z (G = 1 for a pure 2-dx checkerboard, G ~ (pi k dx)^2/4 << 1
for a resolved wave of wavenumber k; defined here, not from a reference),
for the whole wet field and for the ambient region away from the plume.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pfcommon import style  # noqa: E402
from pfcommon.mitgcm_io import read_mds  # noqa: E402
from pfcommon.runctl import require_ok  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def grid_noise(w: np.ndarray, axis: int, mask: np.ndarray) -> float:
    """G indicator (module docstring) along ``axis`` over cells where ``mask``."""
    d2 = np.zeros_like(w)
    sl = [slice(None)] * w.ndim
    c, p, m = list(sl), list(sl), list(sl)
    c[axis], p[axis], m[axis] = slice(1, -1), slice(2, None), slice(None, -2)
    d2[tuple(c)] = w[tuple(p)] - 2 * w[tuple(c)] + w[tuple(m)]
    ok = mask.copy()
    edge = [slice(None)] * w.ndim
    for e in (0, -1):
        edge[axis] = e
        ok[tuple(edge)] = False
    rw = np.sqrt(np.mean(w[mask] ** 2))
    return float(np.sqrt(np.mean(d2[ok] ** 2)) / (4 * rw)) if rw > 0 else float("nan")


def main(run_dir: Path, stem: Path) -> dict:
    require_ok(run_dir)
    meta = json.loads((run_dir / "case.json").read_text())
    p = meta["params"]
    nz, ny, nx, dx, dy, dz = p["nz"], p["ny"], p["nx"], p["dx"], p["dy"], p["dz"]
    W, mt = read_mds(run_dir / "W")
    W = W[0]  # (nz, ny, nx), w on the top face of each cell
    t_min = mt["timeStepNumber"] * p["deltaT"] / 60.0
    wet = np.ones_like(W, dtype=bool)
    wet[:, [0, -1], :] = False
    wet[:, :, -1] = False
    src_j = sorted({c[1] for c in meta["source_cells"]})
    xf, yf, zf = np.arange(nx + 1) * dx, np.arange(ny + 1) * dy, -np.arange(nz + 1) * dz

    # ambient region: farther than 200 m from the face (plume and outflow are near it)
    amb = wet.copy()
    amb[:, :, : int(200 / dx)] = False
    res = {"t_min": t_min, "w_max": float(W[wet].max()), "w_min": float(W[wet].min()),
           "w_rms_ambient": float(np.sqrt(np.mean(W[amb] ** 2))),
           "w_absmax_ambient": float(np.abs(W[amb]).max())}
    for name, ax_ in (("z", 0), ("y", 1), ("x", 2)):
        res[f"G_{name}_all"] = grid_noise(W, ax_, wet)
        res[f"G_{name}_ambient"] = grid_noise(W, ax_, amb)
    # face column only (the plume itself), along y and z
    face = np.zeros_like(wet)
    face[:, 1:-1, 0] = True
    res["G_y_face"] = grid_noise(W, 1, face)
    res["G_z_face"] = grid_noise(W, 0, face)

    style.apply()
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    vm = max(abs(res["w_max"]), abs(res["w_min"]))
    kw = dict(cmap=style.CMAPS["diverging"], vmin=-vm, vmax=vm, shading="flat")
    j0 = src_j[0]
    pc = axs[0].pcolormesh(yf, zf, W[:, :, 0], **kw)
    for j in src_j:
        axs[0].plot((j + 0.5) * dy, zf[-1] + 0.5 * dz, "k^", ms=6)
    axs[0].set(xlabel="y (m)", ylabel="z (m)", title="w at the face column (i = 0)")
    axs[1].pcolormesh(xf[:41], zf, W[:, j0, :40], **kw)
    axs[1].set(xlabel="distance from ice face (m)", ylabel="z (m)",
               title=f"w through the source row j = {j0}")
    k = int(round(0.75 * nz))  # z ~ -225 m, mid lower layer
    axs[2].pcolormesh(xf[:41], yf, W[k, :, :40], **kw)
    axs[2].set(xlabel="distance from ice face (m)", ylabel="y (m)",
               title=f"w at z = {zf[k]:.0f} m (top face of level {k})")
    fig.colorbar(pc, ax=axs, label="w (m/s)", shrink=0.9)
    fig.suptitle(f"{run_dir.name}: w at t = {t_min:.0f} min")
    style.save(fig, stem)
    print(json.dumps(res, indent=2))
    return res


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
