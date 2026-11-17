"""PD1 control figure (plan Task 1): initial T and S sections of chan10 with the J line,
the wind profile and the slope-aware GM factor. Writes figures/pd1_section.{png,pdf}."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pfcommon import style

from ..gendata import fields
from ..gendata.make_case import grid_of, load_case

FIG = Path(__file__).resolve().parents[1] / "figures"


def make(out_stem: Path = FIG / "pd1_section") -> list[Path]:
    style.apply()
    p = load_case("G0")
    g = grid_of(p)
    T, S = fields.initial_ts(g)
    d = fields.depth(g.yc)
    jq = fields.j_qb(g, p["y_qb_m"])
    yk = g.yc / 1e3
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.0), sharex=True,
                           gridspec_kw={"height_ratios": [2.2, 1]}, layout="constrained")
    for a, F, cm, lab in ((ax[0, 0], T, style.CMAPS["temperature"], "θ (°C)"),
                          (ax[0, 1], S, style.CMAPS["salinity"], "S")):
        sec = np.where(g.zc[:, None] <= d[None, :], F[:, :, 0], np.nan)
        pc = a.pcolormesh(yk, g.zc, sec, cmap=cm, shading="nearest")
        a.fill_between(yk, d, 3100, color="0.6")
        a.axvline(g.yg[jq] / 1e3, color="k", lw=0.8, ls="--")
        a.set_ylim(3050, 0)
        a.set_ylabel("depth (m)")
        fig.colorbar(pc, ax=a, label=lab)
    ax[0, 0].text(g.yg[jq] / 1e3 + 4, 2800, f"J line, y = {g.yg[jq] / 1e3:.0f} km", fontsize=7)
    ax[1, 0].plot(yk, fields.wind_tau_x(g.yc))
    ax[1, 0].axhline(0, color="0.5", lw=0.5)
    ax[1, 0].set_ylabel("τx (N m⁻²)")
    ax[1, 1].plot(yk, fields.gm_factor(g.yc))
    ax[1, 1].set_ylabel("GM factor")
    for a in ax[1]:
        a.set_xlabel("y (km)")
    out = style.save(fig, out_stem)
    plt.close(fig)
    return out


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__).parse_args()
    print(make())
