"""Figure 0d of Experiment 2: the ISOMIP+ geometry in perspective, annotated with VectorCraft.

Pipeline (each step can be run on its own):
  1. fig0_geometry_export.py      (pyocean Python) geometry file -> WORK/geom.npz
  2. fig0_geometry_blender.py     (blender -b --python ... -- --work WORK) -> WORK/geometry_render.png and
                                  WORK/geometry_points.json (image positions of the labelled 3-D points)
  3. fig0_geometry.py             (this script) places the render, draws the labels at those positions, the
                                  dimensions, process 7 and an inset of the ramp profile family, and exports
                                  figures/fig0_geometry.png with the VectorCraft CLI.

Usage: fig0_geometry.py [--work DIR] [--cli PATH_TO_vectorcraft-cli]   (or env VECTORCRAFT_CLI)
       fig0_geometry.py --all [--blender PATH] ...   runs steps 1-3
No model output enters this figure: the render is of the input geometry; the arrows are a schematic of published
circulation features (Yung et al., 2026) and the inset shows the forcing design, not results.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import fig0_schematic as S
from fig0_schematic import INK, WARM, COOL, PURP, GREY, ICE_LINE, arrow, c, curve_arrow, fill, line, marker, rect, text

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
W, H = 1100, 610
RX, RY, RW = 40.0, 40.0, 1020.0          # the render on the page


def build(work: Path) -> list[dict]:
    meta = json.loads((work / "geometry_points.json").read_text())
    rw, rh = meta["res"]
    RH = RW * rh / rw
    def P(name, dx=0.0, dy=0.0):
        fx_, fy_ = meta["points"][name]
        return RX + fx_ * RW + dx, RY + fy_ * RH + dy

    S.cmds.clear()
    c("file.new", width=W, height=H, units="Points", backgroundContents="white", name="fig0_geometry")
    c("file.place", path=str(work / "geometry_render.png"), rect=[RX, RY, RW, RH], link=False)
    text(20, 26, "(d)", 13, style="Bold")
    text(46, 26, "The ISOMIP+ geometry in perspective: the real ice draft and bed of the input file", 11)

    def callout(xy, s, dx, dy, color=INK, size=9.5, justify=None, box=True):
        x, y = P(xy) if isinstance(xy, str) else xy
        line(x, y, x + dx, y + dy, color=color, weight=0.8)
        c("shape.ellipse", x=x - 2.2, y=y - 2.2, width=4.4, height=4.4)
        fill(color=color)
        j = justify or ("left" if dx >= 0 else "right")
        tx = x + dx + (3 if j == "left" else -3)
        for k, s_ in enumerate(s.split("\n")):
            ty = y + dy + 3.5 + k * (size + 2.5) - (0 if dy >= 0 or "\n" not in s else (size + 2.5) * s.count("\n"))
            if box:
                w_ = 0.5 * size * len(s_)
                rect(tx - (w_ if j == "right" else 0) - 2, ty - size + 1, w_ + 4, size + 3, color="#ffffff", opacity=80)
            text(tx, ty, s_, size, color=color, justify=j)

    # the geometry
    callout("ice_grounded", "grounded ice", 40, -26)
    callout("wall_far_top", "side wall y = 80 km (closed, no-slip)", -22, -30)
    callout("ice_shelf_top", "ice shelf (floating; drawn translucent)", 30, -70)
    callout("front_far_top", "calving front\nx = 640 km", 22, -36)
    callout("gl_cut", "grounding line", -44, -46)
    callout("bed_grounded_cut", "bed", -30, -26, color="#5f4f40")
    callout("tracer_cut", "meltwater tracer source: the ice base", 8, -34, color=PURP)
    callout("restore_top", "restoring zone\nx = 790–800 km", -26, -46, color="#c4622f")
    x, y = P("cavity_cut")
    text(x + 30, y + 4, "cavity", 10.5, color="#ffffff", style="Bold")
    x, y = P("open_ocean")
    text(x + 30, y + 30, "open ocean", 10, color="#1f4e79", style="Bold")
    # dimensions: 480 km along the near bottom edge, 80 km across the end, 720 m deep
    (x0, y0), (x1, y1) = P("c_0_0_bot2", 0, 14), P("c_480_0_bot", 0, 14)
    arrow(x0, y0, x1, y1, color=INK, weight=1.1, both=True)
    text((x0 + x1) / 2 + 20, (y0 + y1) / 2 + 24, "480 km (x = 320–800 km)", 10, style="Bold", justify="center")
    (x0, y0), (x1, y1) = P("c_480_0_bot", 14, 8), P("c_480_80_bot", 14, 8)
    arrow(x0, y0, x1, y1, color=INK, weight=1.1, both=True)
    text(x0 + 8, (y0 + y1) / 2 + 8, "80 km", 10, style="Bold")
    (x0, y0), (x1, y1) = P("c_0_0_sl", -16, 0), P("c_0_0_bot2", -16, 0)
    arrow(x0, y0, x1, y1, color=INK, weight=1.1, both=True)
    text(min(x0, x1) - 6, (y0 + y1) / 2 + 4, "720 m", 10, style="Bold", justify="right")
    # process 7: rotation (f-plane, 75 S) and the cyclonic cell along the far wall
    xo, yo = P("outflow_end")
    marker(xo + 18, yo - 4, 7)
    xi, yi = P("inflow_mid")
    S.label(xi - 20, yi + 26, "inflow at depth (seen through the cut)", 8.5, color=WARM)
    callout("outflow_open", "outflow along the y = 80 km wall", 26, -40, color=PURP, size=8.5)
    # inset: the ramp profile family (forcing design)
    inset_profile(834, 50)
    # notes
    ny = RY + RH - 26
    marker(49, ny + 14, 7)
    text(66, ny + 18, "Rotation: the domain is an f-plane at 75° S (f < 0). In most ISOMIP+ Ocean1 models the cavity holds a "
         "cyclonic (clockwise from above) circulation,", 9.5)
    text(66, ny + 31, "intensified along the y = 80 km wall (Yung et al., 2026); the arrows sketch that published pattern, "
         "not our output.", 9.5)
    text(40, ny + 52, "Rendered with Blender from the ISOMIP+ input geometry (Cornford and Asay-Davis, 2016) on its 1 km "
         "grid; vertical exaggeration ×100. The geometry is symmetric about y = 40 km:", 9, color=GREY)
    text(40, ny + 64, "the near half (y = 0–40 km) is cut away along the centreline and drawn in outline. Arrows and the "
         "tracer layer are enlarged.", 9, color=GREY)
    return S.cmds


def inset_profile(x0, y0, w=226, h=176):
    """The ramp family: T = -1.9 C above z_c - d/2, T_bot below z_c + d/2, linear across a thermocline d thick."""
    rect(x0, y0, w, h, color="#ffffff", stroke=ICE_LINE, weight=1.0, opacity=94)
    text(x0 + 8, y0 + 15, "forcing design: ramp profiles", 9, style="Bold")
    ax, ay, aw, ah = x0 + 40, y0 + 34, w - 70, h - 62         # plot box: T to the right, depth down
    line(ax, ay, ax + aw, ay, color=INK, weight=0.8)
    line(ax, ay, ax, ay + ah, color=INK, weight=0.8)
    arrow(ax, ay + ah - 2, ax, ay + ah + 8, color=INK, weight=0.8)
    text(ax - 4, ay + 4, "0", 8, justify="right")
    text(ax - 4, ay + ah + 4, "720 m", 8, justify="right")
    text(ax + aw / 2, ay - 4, "temperature T", 8, justify="center")
    tc = ax + 8                                               # -1.9 C
    tb = ax + aw - 12                                         # T_bot
    for zc, op, wgt in ((0.30, 50, 1.0), (0.70, 50, 1.0), (0.50, 100, 2.0)):
        za, zb_ = ay + ah * (zc - 0.21), ay + ah * (zc + 0.21)
        S.stroke_path(f"M {tc} {ay} L {tc} {za} L {tb} {zb_} L {tb} {ay + ah}", WARM, wgt, opacity=op)
    za, zm, zb_ = ay + ah * 0.29, ay + ah * 0.5, ay + ah * 0.71
    line(ax, zm, ax + aw, zm, color=INK, weight=0.7, dash=[2, 2])
    text(ax + aw + 4, zm + 3.5, "z_{c}", 10, style="Bold")
    arrow(ax + aw - 2, za, ax + aw - 2, zb_, color=GREY, weight=0.8, both=True)
    text(ax + aw + 4, za + 8, "d", 9, color=GREY)
    text(tc + 7, ay + 14, "−1.9 °C", 8, color=COOL)
    text(tb - 2, ay + ah - 6, "T_{bot}", 9, color=WARM, justify="right")
    text(x0 + 8, y0 + h - 14, "T_{bot}, z_{c} varied; d = 300 m; S ramps", 8, color=GREY)
    text(x0 + 8, y0 + h - 4, "the same way (33.8 to 34.7)", 8, color=GREY)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", default=str(Path(tempfile.gettempdir()) / "fig0_geometry"))
    ap.add_argument("--cli", default=os.environ.get("VECTORCRAFT_CLI"))
    ap.add_argument("--all", action="store_true", help="run the export and the Blender render first")
    ap.add_argument("--blender", default=os.environ.get("BLENDER", "blender"))
    a = ap.parse_args()
    work = Path(a.work)
    if a.all:
        subprocess.run([sys.executable, str(HERE / "fig0_geometry_export.py"), "--work", str(work)], check=True)
        subprocess.run([a.blender, "-b", "--python", str(HERE / "fig0_geometry_blender.py"), "--", "--work", str(work)],
                       check=True, stdout=subprocess.DEVNULL)
    cmd_file = HERE / "fig0_geometry.commands.json"
    cmd_file.write_text(json.dumps(build(work), indent=1))
    if not a.cli:
        print(f"wrote {cmd_file}; pass --cli to render")
        return
    png = FIG / "fig0_geometry.png"
    r = subprocess.run([a.cli, "run", "--commands", str(cmd_file), "--ppi", "144", "--export", str(png)],
                       capture_output=True, text=True, env={**os.environ, "VECTORCRAFT_NO_PREFS": "1"})
    if r.returncode or any('"error"' in l for l in r.stdout.splitlines()):
        raise SystemExit(f"vectorcraft-cli failed (rc {r.returncode}):\n" + r.stdout[-2000:] + r.stderr[-2000:])
    # README copy: an adaptive 256-colour palette keeps the PNG under 1 MB (the render has no long smooth gradients)
    from PIL import Image
    im = Image.open(png).convert("RGB")
    im.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.FLOYDSTEINBERG).save(png, optimize=True)
    print("rendered", png.name, f"{png.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
