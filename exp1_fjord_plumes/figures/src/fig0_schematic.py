"""Figure 0: schematic of the phenomenon (along-fjord section), drawn with VectorCraft.

Writes the VectorCraft command list (``fig0_schematic.commands.json``) and, if the VectorCraft CLI is
available, renders ``figures/fig0_schematic.{png,svg,pdf}`` with its headless batch mode:

    vectorcraft-cli run --commands fig0_schematic.commands.json --export fig0_schematic.png ...

Tool: VectorCraft (https://github.com/storytold/vectorcraft), CLI release v0.22.0, macOS universal.
Usage: fig0_schematic.py [--cli PATH_TO_vectorcraft-cli]   (or env VECTORCRAFT_CLI)

Coordinates are points, y down. The drawing is a schematic, not to scale; it follows the set-up of
the experiment (ice at the left, as x = 0 in the model; fjord mouth and shelf at the right).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
W, H = 1000, 560
SEA, PYC, BED = 110, 260, 420          # sea surface, pycnocline centre, sea floor (y, points)
FACE, MOUTH = 170, 900                 # ice face and fjord mouth (x, points)

ICE, ICE_LINE = "#e8f3f9", "#6b8fa3"
PW, AW = "#a9cfe6", "#f2b48c"          # Polar Water (cold), Atlantic Water (warm)
PLUME = "#b9a7d8"                      # mixed plume water
INK, WARM, COOL, BEDC = "#222222", "#c0392b", "#2c6ea3", "#a39282"
FONT = "Helvetica"

cmds: list[dict] = []


def c(cmd: str, **params) -> None:
    cmds.append({"command": cmd, "params": params})


def fill(color=None, gradient=None, stroke=None, weight=1.0, opacity=None):
    c("paint.setFill", **({"gradient": gradient} if gradient else {"color": color} if color else {"none": True}))
    if stroke:
        c("paint.setStroke", color=stroke)
        c("stroke.set", weight=weight)
    else:
        c("paint.setStroke", none=True)
    if opacity is not None:
        c("object.setProps", opacity=opacity)


def rect(x, y, w, h, **kw):
    c("shape.rectangle", x=x, y=y, width=w, height=h)
    fill(**kw)


def path(d, **kw):
    c("path.create", d=d)
    fill(**kw)


def arrow(x1, y1, x2, y2, color=INK, weight=2.0, both=False, dash=None):
    c("shape.line", x1=x1, y1=y1, x2=x2, y2=y2)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    st = {"weight": weight, "endArrow": "Triangle", "cap": "round"}
    if both:
        st["startArrow"] = "Triangle"
    if dash:
        st["dash"] = dash
    c("stroke.set", **st)


def line(x1, y1, x2, y2, color=INK, weight=1.0, dash=None):
    c("shape.line", x1=x1, y1=y1, x2=x2, y2=y2)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, **({"dash": dash} if dash else {}))


def text(x, y, s, size=12, color=INK, style=None, justify=None, rotation=None):
    p = {"x": x, "y": y, "text": s, "size": size, "font": FONT, "color": color}
    if style:
        p["style"] = style
    if justify:
        p["justify"] = justify
    if rotation is not None:
        p["rotation"] = rotation
    c("text.create", **p)


def build() -> list[dict]:
    c("file.new", width=W, height=H, units="Points", backgroundContents="white", name="fig0_schematic")
    # sky, water masses, pycnocline band, sea floor
    rect(0, 0, W, SEA, color="#f6f9fb")
    rect(FACE, SEA, W - FACE, PYC - SEA, color=PW)
    rect(FACE, PYC, W - FACE, BED - PYC, color=AW)
    rect(FACE, PYC - 18, W - FACE, 36,
         gradient={"kind": "linear", "angle": 90, "stops": [{"offset": 0, "color": AW}, {"offset": 1, "color": PW}]})
    rect(0, BED, W, H - BED, color=BEDC)
    line(FACE, SEA, W, SEA, color=COOL, weight=1.5)
    # glacier (ice above and below the sea surface, grounded on the bed)
    path(f"M 0 52 L 60 40 L 120 46 L {FACE} 34 L {FACE} {BED} L 0 {BED} Z", color=ICE, stroke=ICE_LINE, weight=1.5)
    # surface meltwater draining to the bed and out at the grounding line
    line(105, 44, 105, BED - 8, color=COOL, weight=1.6, dash=[5, 4])
    arrow(105, BED - 8, FACE - 4, BED - 8, color=COOL, weight=1.6, dash=[5, 4])
    # buoyant plume (half cone against the face) and its outflow at neutral buoyancy
    path(f"M {FACE} {BED} L {FACE + 10} {BED} C {FACE + 25} 350 {FACE + 40} 300 {FACE + 58} 262 "
         f"C 300 266 450 262 560 246 C 450 230 300 226 {FACE} 236 Z", color=PLUME, opacity=90)
    arrow(FACE + 14, BED - 16, FACE + 24, 278, color="#ffffff", weight=3)
    arrow(560, 246, 700, 246, color="#6a4c9c", weight=3)
    # entrainment of ambient water into the plume
    for y in (395, 350, 305):
        arrow(FACE + 95, y, FACE + 50 - (BED - y) * 0.08, y - 4, color=WARM, weight=1.6)
    # warm inflow at depth
    arrow(840, 370, 470, 370, color=WARM, weight=3)
    # wind over the fjord
    for x0 in (300, 480, 660):
        arrow(x0, 82, x0 + 120, 82, color=INK, weight=2.2)
    # shelf pycnocline heaving at the mouth
    line(MOUTH, SEA, MOUTH, BED, color="#555555", weight=1.0, dash=[3, 3])
    arrow(945, PYC - 40, 945, PYC + 40, color=INK, weight=2.2, both=True)
    path(f"M {FACE + 260} {PYC} C 500 {PYC - 14} 600 {PYC + 14} 700 {PYC} C 780 {PYC - 14} 840 {PYC + 14} {MOUTH} {PYC}",
         stroke="#555555", weight=1.2)
    # labels
    text(30, 24, "Glacier", 15, style="Bold")
    text(520, 24, "Fjord", 15, style="Bold", justify="center")
    text(945, 24, "Shelf", 15, style="Bold", justify="center")
    text(46, 300, "surface meltwater drains to the bed", 10.5, color=COOL, rotation=-90)
    text(150, 395, "submarine melting", 11, color=WARM, rotation=-90)
    text(FACE + 6, BED + 18, "subglacial discharge leaves the bed at the grounding line", 11, color="#ffffff")
    text(FACE + 76, 228, "neutral buoyancy: the plume detaches and spreads", 11)
    text(FACE + 100, 330, "entrainment", 11, color=WARM)
    text(FACE + 62, 285, "buoyant plume rises along the ice", 10.5, color="#4b3577")
    text(575, 236, "outflow of plume water", 11, color="#6a4c9c")
    text(480, 390, "warm Atlantic Water flows toward the glacier", 11, color=WARM)
    text(560, 70, "along-fjord wind", 11, justify="center")
    text(890, 160, "Polar Water: cold, fresh", 11, justify="right")
    text(890, 410, "Atlantic Water: warm, salty", 11, justify="right")
    text(960, PYC + 58, "pycnocline", 10.5, justify="center")
    text(960, PYC + 72, "heaving (days)", 10.5, justify="center")
    # which part of the experiment covers what
    line(FACE, 484, 300, 484, color="#ffffff", weight=2)
    text(FACE, 500, "(a), (b): plume at the ice face, resolved (metres)", 11, color="#ffffff")
    line(FACE, 516, MOUTH, 516, color="#ffffff", weight=2)
    text(FACE, 532, "(c): whole fjord with plume, wind and shelf forcing; plume parameterised (40 km)", 11, color="#ffffff")
    text(W - 10, H - 8, "schematic, not to scale", 9, color="#ffffff", justify="right")
    return cmds


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cli", default=os.environ.get("VECTORCRAFT_CLI"))
    a = ap.parse_args()
    cmd_file = HERE / "fig0_schematic.commands.json"
    cmd_file.write_text(json.dumps(build(), indent=1))
    if not a.cli:
        print(f"wrote {cmd_file}; pass --cli to render")
        return
    outs = [FIG / f"fig0_schematic.{e}" for e in ("png", "svg", "pdf")]
    args = [a.cli, "run", "--commands", str(cmd_file), "--ppi", "200"]
    for o in outs:
        args += ["--export", str(o)]
    r = subprocess.run(args, capture_output=True, text=True, env={**os.environ, "VECTORCRAFT_NO_PREFS": "1"})
    bad = [l for l in r.stdout.splitlines() if '"error"' in l]
    if r.returncode or bad:
        raise SystemExit(f"vectorcraft-cli failed (rc {r.returncode}):\n" + "\n".join(bad) + r.stderr[-2000:])
    print("rendered", ", ".join(o.name for o in outs))


if __name__ == "__main__":
    main()
