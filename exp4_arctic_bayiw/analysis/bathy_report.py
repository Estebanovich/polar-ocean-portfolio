"""Control report and figure of the model bathymetry of Van Mijenfjorden (Exp. 4, Task 1b).

Runs gendata.bathy.build on the fjord_vmf config of cases.yaml, writes figures/bathy_report.json and
figures/bathy_report.{png,pdf}, and (optionally) the MITgcm bathyFile. The binary is not committed
(> 1 MB); Task 3/4 make_case regenerates it.

Usage: bathy_report.py [--data DIR] [--bathy-out FILE] [--no-figure]
Data credits: IBCAO v5.2 (doi:10.17043/ibcao-5.2); coastline (c) Norwegian Polar Institute, S100
Kartdata, CC BY 4.0 (doi:10.21334/npolar.2014.645336c7).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import matplotlib.colors as mcolors
import numpy as np
import yaml
from pyproj import Transformer
from shapely.geometry import box

from exp4_arctic_bayiw.gendata import bathy  # noqa: E402
from pfcommon import style  # noqa: E402

PKG = Path(__file__).resolve().parents[1]
CASES = PKG / "cases.yaml"
FIG = PKG / "figures" / "bathy_report"


def section_profile(grid, src, sec, ds=10.0):
    """Distance from the island end (m), model depth and IBCAO depth (nearest 100 m point) along a section."""
    P, step = bathy._samples(sec["p0"], sec["p1"], ds)
    xc, yc, dx = grid["xc"], grid["yc"], grid["dx_m"]
    i = np.clip(((P[:, 0] - (xc[0] - dx / 2)) // dx).astype(int), 0, xc.size - 1)
    j = np.clip(((P[:, 1] - (yc[0] - dx / 2)) // dx).astype(int), 0, yc.size - 1)
    dg = np.where(grid["wet"][j, i], grid["depth"][j, i], np.nan)
    px, py = Transformer.from_crs(bathy.CRS_MODEL, bathy.CRS_IBCAO, always_xy=True).transform(P[:, 0], P[:, 1])
    ii = np.abs(src["x"][None] - px[:, None]).argmin(axis=1)
    jj = np.abs(src["y"][None] - py[:, None]).argmin(axis=1)
    z = src["z"][jj, ii]
    s = (np.arange(len(P)) + 0.5) * step
    return s, dg, np.where(z < 0, -z, np.nan)


def _coast(ax, land, bb, **kw):
    g = land.intersection(box(bb[0], bb[2], bb[1], bb[3]))
    for p in getattr(g, "geoms", [g]):
        if p.geom_type != "Polygon":
            continue
        for ring in [p.exterior, *p.interiors]:
            x, y = ring.xy
            ax.plot(np.asarray(x) / 1e3, np.asarray(y) / 1e3, **kw)


def figure(src, grid, rep, stem: Path):
    import matplotlib.pyplot as plt

    style.apply()
    xc, yc, dx = grid["xc"] / 1e3, grid["yc"] / 1e3, grid["dx_m"] / 1e3
    xe = np.r_[xc - dx / 2, xc[-1] + dx / 2]
    ye = np.r_[yc - dx / 2, yc[-1] + dx / 2]
    d = np.where(grid["wet"], grid["depth"], np.nan)
    land, secs, cfg = grid["land"], grid["sections"], grid["cfg"]
    fig = plt.figure(figsize=(7.0, 8.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 1.0, 0.8], hspace=0.38, wspace=0.28)

    # (a) model depth with the NPI coastline, sponge band and sound sections
    ax = fig.add_subplot(gs[0, :])
    pc = ax.pcolormesh(xe, ye, d, cmap="cmo.deep", vmin=0, vmax=240, rasterized=True)
    bb = cfg["bbox_utm"]
    _coast(ax, land, bb, color="k", lw=0.35)
    n = cfg["n_sponge"] * dx
    x0, x1, y0, y1 = (v / 1e3 for v in bb)
    ax.plot([x0 + n, x0 + n, x1], [y1 - n, y0 + n, y0 + n], color="#D55E00", lw=0.8, ls="--")
    ax.plot([x0 + n, x1], [y1 - n, y1 - n], color="#D55E00", lw=0.8, ls="--", label="sponge (10 cells)")
    for name, s in secs.items():
        ax.plot([s["p0"][0] / 1e3, s["p1"][0] / 1e3], [s["p0"][1] / 1e3, s["p1"][1] / 1e3], color="#E69F00", lw=1.5)
    ax.add_patch(plt.Rectangle((488.5, 8620.5), 9.0, 12.5, fill=False, ec="#CC79A7", lw=0.9))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xlabel("UTM 33N easting (km)")
    ax.set_ylabel("northing (km)")
    ax.set_title("(a) Model depth, 200 m (IBCAO v5.2 area mean, NPI S100 coastline in black)", loc="left")
    ax.legend(loc="lower right", frameon=False)
    fig.colorbar(pc, ax=ax, label="depth (m)", shrink=0.8, pad=0.01)

    # (b) zoom on Akseloya and its two sounds, grid cells visible
    ax = fig.add_subplot(gs[1, 0])
    zb = (488500, 497500, 8620500, 8633000)
    ax.pcolormesh(xe, ye, d, cmap="cmo.deep", vmin=0, vmax=120, edgecolors="none", rasterized=True)
    _coast(ax, land, zb, color="k", lw=0.6)
    for name, s in secs.items():
        ax.plot([s["p0"][0] / 1e3, s["p1"][0] / 1e3], [s["p0"][1] / 1e3, s["p1"][1] / 1e3], color="#E69F00", lw=1.5)
        ax.annotate(name, (s["p1"][0] / 1e3, s["p1"][1] / 1e3), xytext=(5, -2) if name == "Mariasundet" else (5, 6),
                    textcoords="offset points", fontsize=7,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
    ax.set_xlim(zb[0] / 1e3, zb[1] / 1e3)
    ax.set_ylim(zb[2] / 1e3, zb[3] / 1e3)
    ax.set_aspect("equal")
    ax.set_title("(b) Akseløya sounds (0-120 m)", loc="left")
    ax.set_xlabel("easting (km)")

    # (c) data type of the model sea cells (GEBCO TID of the most frequent IBCAO point)
    ax = fig.add_subplot(gs[1, 1])
    t = grid["tid_mode"]
    cls = np.full(t.shape, np.nan)
    cls[np.isin(t, list(bathy.TID_CLASSES["measured"]))] = 0
    cls[np.isin(t, list(bathy.TID_CLASSES["indirect"]))] = 1
    cls[np.isin(t, list(bathy.TID_CLASSES["unknown"]) + [0, 100])] = 2
    cls[~grid["wet"]] = np.nan
    cmap = mcolors.ListedColormap(["#0072B2", "#E69F00", "#999999"])
    ax.pcolormesh(xe, ye, cls, cmap=cmap, vmin=-0.5, vmax=2.5, rasterized=True)
    _coast(ax, land, bb, color="k", lw=0.25)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    f = rep["tid_fraction"]["cells_inner"]
    ax.set_title(f"(c) IBCAO TID: inner fjord {100 * f['measured']:.0f} % measured", loc="left")
    for k, (lab, c) in enumerate([("measured (TID 10-17)", "#0072B2"), ("interpolated (40-47)", "#E69F00"),
                                 ("land code / unknown", "#999999")]):
        ax.plot([], [], "s", color=c, label=lab)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False, fontsize=6,
              handlelength=1, columnspacing=0.8)
    ax.set_xlabel("easting (km)", labelpad=1)

    # (d) hypsometry of the inner fjord
    ax = fig.add_subplot(gs[2, 0])
    h = rep["hypsometry_inner"]
    e = np.asarray(h["bin_edges_m"])
    m = e[:-1] < 130
    ax.barh(e[:-1][m] + 5, np.asarray(h["area_ibcao_m2"])[m] / 1e6, height=4.2, color="#56B4E9", label="IBCAO 100 m")
    ax.barh(e[:-1][m] + 5, np.asarray(h["area_grid_m2"])[m] / 1e6, height=2.0, color="#D55E00", label="model 200 m")
    ax.invert_yaxis()
    ax.set_xlabel("area (km²) per 10 m bin")
    ax.set_ylabel("depth (m)")
    v = rep["inner_fjord"]
    ax.set_title(f"(d) Inner fjord hypsometry; volume ratio {v['volume_ratio_grid_to_ibcao']:.3f}", loc="left")
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.92))

    # (e) sections across the two sounds
    ax = fig.add_subplot(gs[2, 1])
    for (name, s), c in zip(secs.items(), ("#0072B2", "#009E73")):
        dist, dg, di = section_profile(grid, src, s)
        ax.plot(dist, di, color=c, lw=0.8, ls=":", label=f"{name} IBCAO")
        ax.plot(dist, dg, color=c, lw=1.4, drawstyle="steps-mid", label=f"{name} model")
    ax.invert_yaxis()
    ax.set_xlabel("distance from Akseløya along the section (m)")
    ax.set_ylabel("depth (m)")
    ax.set_title("(e) Sections across the sounds", loc="left")
    ax.legend(frameon=False, fontsize=6.5)
    fig.text(0.01, 0.005, "Data: IBCAO v5.2; coastline © Norwegian Polar Institute (S100 Kartdata, CC BY 4.0).",
             fontsize=6)
    return style.save(fig, stem)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, default=PKG / "data_external")
    ap.add_argument("--bathy-out", type=Path, default=PKG / "data_external" / "derived" / "bathy_fjord_vmf.bin")
    ap.add_argument("--no-figure", action="store_true")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(CASES.read_text())["configs"]["fjord_vmf"]
    src, grid = bathy.build(cfg, a.data)
    rep = bathy.report(src, grid)
    FIG.parent.mkdir(parents=True, exist_ok=True)
    FIG.with_suffix(".json").write_text(json.dumps(rep, indent=1) + "\n")
    if a.bathy_out:
        a.bathy_out.parent.mkdir(parents=True, exist_ok=True)
        bathy.write_bathy(a.bathy_out, grid["depth"], grid["wet"])
        print(f"bathyFile -> {a.bathy_out} ({grid['depth'].shape[0]} x {grid['depth'].shape[1]})")
    if not a.no_figure:
        for p in figure(src, grid, rep, FIG):
            print(p)
    for c in rep["controls"]:
        print(f"{c['name']:34s} {c['value']:12.4g} ref {c['reference']:12.4g} tol {c['tolerance']:10.3g} "
              f"{'PASS' if c['pass'] else 'FAIL'} {c['classification'] or ''}")


if __name__ == "__main__":
    main()
