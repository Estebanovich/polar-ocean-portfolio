"""DESCRIPTIVE GIF of the eddy spin-up of the 2.5-km model (F1, then F1b): plan view of the 5-day
mean meridional velocity v near the surface, and the monthly eddy kinetic energy of the same runs.
No number of the README comes from it and no criterion or hypothesis uses it.

Field: VVEL of the 5-day diagnostics jflux5d (record 1; data.diagnostics of chan2p5) at level k
(default 2, the cell centred at RC[k]), from the start of F1 (the F0b pickup plus member-0 theta
noise <= 1e-3 degC) to the end of F1b. EKE: crosscheck.json (EKE_spinup_m2s2, analysis/crosscheck.eke).

Usage: python -m exp5_adjoint_slope.analysis.animate [--k 2] [--vmax 0.15] [--dpi 90] [--max-mb 6]
       (writes figures/anim_fine_spinup.gif)"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pfcommon import style

style.apply()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402

FIG = Path(__file__).resolve().parents[1] / "figures"
DAY = 86400.0
KM = 1e-3


def frame_times(t_end_s, t0_s: float) -> np.ndarray:
    """Centre of each 5-day mean, in days since t0 (records are stamped at their end)."""
    return (np.asarray(t_end_s, float) - 2.5 * DAY - t0_s) / DAY


def eke_months(series_f1, series_f1b, f1_days: float = 365.0, month_d: float = 30.0):
    """Mid-month times (days since the start of F1) and EKE of F1 then F1b (F1b starts at f1_days)."""
    t1 = (np.arange(len(series_f1)) + 0.5) * month_d
    t2 = f1_days + (np.arange(len(series_f1b)) + 0.5) * month_d
    return np.concatenate([t1, t2]), np.asarray(list(series_f1) + list(series_f1b), float)


def _records(run: Path):
    its = sorted(int(f.name.split(".")[1]) for f in run.glob("jflux5d.*.data") if not f.name.startswith("._"))
    return its


def load(runs_root: Path, k: int):
    from pfcommon.mitgcm_io import read_bin, read_mds
    f1 = runs_root / "F1"
    case = json.loads((f1 / "input" / "case.json").read_text())
    nr, ny, nx = case["params"]["nz"], case["params"]["ny"], case["params"]["nx"]
    data = (f1 / "input" / "data").read_text()
    from pfcommon.campaign import nml_float, nml_int
    dt = nml_float(data, "deltaT")
    t0 = nml_int(data, "nIter0") * dt
    v, t = [], []
    for name in ("F1", "F1b"):
        run = runs_root / name
        for it in _records(run):
            v.append(np.asarray(read_mds(run / "jflux5d", it)[0])[1, k])
            t.append(it * dt)
    grid = {"xc": read_bin(f1 / "XC.data", (ny, nx))[0] * KM, "yc": read_bin(f1 / "YC.data", (ny, nx))[:, 0] * KM,
            "depth": read_bin(f1 / "Depth.data", (ny, nx)), "z": float(read_bin(f1 / "RC.data", (nr,))[k]),
            "y_j": float(read_bin(f1 / "YG.data", (ny, nx))[case["j_qb"], 0]) * KM}
    return np.array(v), frame_times(t, t0), grid


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=2)
    ap.add_argument("--vmax", type=float, default=0.15)
    ap.add_argument("--dpi", type=int, default=90)
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--max-mb", type=float, default=6.0)
    a = ap.parse_args()
    from pfcommon.paths import EXP5_RUNS
    v, t, g = load(Path(EXP5_RUNS), a.k)
    spin = json.loads((FIG / "crosscheck.json").read_text())["EKE_spinup_m2s2"]
    tm, e = eke_months(spin["F1"], spin["F1b"])
    try:
        import cmocean
        cmap = cmocean.cm.balance
    except ImportError:
        cmap = "RdBu_r"

    fig = plt.figure(figsize=(9.0, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.28, left=0.06, right=0.98, bottom=0.12, top=0.92)
    ax, ae = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    wet = np.where(g["depth"] > 0, 1.0, np.nan)
    im = ax.pcolormesh(g["xc"], g["yc"], v[0] * wet, cmap=cmap, vmin=-a.vmax, vmax=a.vmax, shading="auto")
    cs = ax.contour(g["xc"], g["yc"], g["depth"], levels=[600, 1500, 2500], colors="0.25", linewidths=0.6)
    ax.axhline(g["y_j"], color="k", lw=0.9, ls="--")
    ax.text(g["xc"][-1], g["y_j"] + 6, "J line", ha="right", va="bottom", fontsize=7)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_aspect("equal")
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.16, extend="both")
    cb.set_label(f"v at {abs(g['z']):.0f} m (m s$^{{-1}}$)")
    title = ax.set_title("", loc="left", fontsize=9)

    ae.semilogy(tm, e, "o-", color="0.35", ms=3, lw=1)
    ae.axvline(365, color="0.6", lw=0.6, ls=":")
    ae.text(368, e.max(), "extension", fontsize=7, va="top", color="0.4")
    ae.set_xlabel("days since the θ noise was added")
    ae.set_ylabel("eddy kinetic energy (m$^2$ s$^{-2}$)")
    ae.set_title("monthly EKE, volume mean", loc="left", fontsize=9)
    ae.set_xlim(0, tm[-1] + 15)
    mark = ae.axvline(t[0], color="C3", lw=1.2)
    fig.text(0.01, 0.01, "Descriptive. Grey lines: 600, 1500 and 2500 m isobaths (shelf to the south). The colour scale saturates in the burst of months 2-3.",
             fontsize=6.5, color="0.35")

    def draw(i):
        im.set_array((v[i] * wet).ravel())
        title.set_text(f"2.5-km model, day {t[i]:.0f}")
        mark.set_xdata([t[i], t[i]])

    out = FIG / "anim_fine_spinup.gif"
    w = PillowWriter(fps=a.fps)
    with w.saving(fig, out, a.dpi):
        for i in range(len(t)):
            draw(i)
            w.grab_frame()
    mb = out.stat().st_size / 2 ** 20
    print(json.dumps({"gif": str(out), "frames": len(t), "mb": round(mb, 2), "over_budget": mb > a.max_mb}))


if __name__ == "__main__":
    main()
