"""Figure 0 of Experiment 3: a guide to the emulator protocol, drawn with VectorCraft.

Writes the VectorCraft command list (``fig0_protocol.commands.json``) and, if the VectorCraft CLI is available,
renders ``figures/fig0_protocol.png`` and ``.pdf`` with its headless batch mode, with the drawing primitives of
exp2_isomip_cavity/figures/src/fig0_schematic.py:

    vectorcraft-cli run --commands fig0_protocol.commands.json --export fig0_protocol.png ...

Tool: VectorCraft (https://github.com/storytold/vectorcraft), CLI release v0.22.0, macOS universal.
Usage: fig0_protocol.py [--cli PATH_TO_vectorcraft-cli]   (or env VECTORCRAFT_CLI)

Content (no result of Exp. 3 anywhere; melt maps, curves and ensemble members are sketches):
(1) what an emulator does: three numbers per run (Tbot, z_c, s) -> MITgcm melt map, or -> a fast emulator;
(2) where each method is tested: the 29 Exp. 2 runs in the (Tbot, z_c) plane as designed, with the held-out sets of
    pd0_preregistro.md; (3) the 8 methods on a physics-to-data axis and the idea of the hybrid M7; (4) the three
    pre-registered hypotheses H3a-H3c. The ice-shelf outline in the map sketches is the ISOMIP+ grounding line
    (exp2_isomip_cavity/figures/src/isomip_outline.json). Coordinates are points, y down.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
REPO = FIG.parent.parent
W, H = 1100, 945
FONT = "Helvetica"

INK, WARM, COOL, PURP, GREY = "#222222", "#b8322a", "#2c6ea3", "#5b3f8f", "#555555"
PAPER, PANEL_LINE, BOX_LINE = "#f7f4ee", "#9a8f80", "#b5ab9c"
TRAIN, GEOM, TEST, REF = "#2c6ea3", "#5aa0d8", "#b8322a", "#222222"
PHYS, STAT, NN, BASE = "#cfe3f2", "#d6ecd9", "#f6dcc8", "#e4e4e4"
ICE, ICE_DEEP, ICE_LINE, BED = "#e6f0f7", "#b9d2e6", "#7d9cb5", "#8a7560"
MELT = ((0, "#7f1d1d"), (0.35, "#d0573a"), (0.7, "#f2b07a"), (1, "#fbe7c6"))   # high -> low, a sketch

# The design of the runs (exp3 dataset; the case names encode Tbot, z_c and s): inputs, not results.
FACTORIAL = [(t, z) for z in (250, 400, 550) for t in (-0.5, 0.0, 0.5, 1.0, 1.5)]
GEOMETRY = [(0.75, 527, 0.8), (1.43, 295, 0.8), (-0.05, 328, 0.8), (-0.42, 402, 0.8), (0.34, 444, 0.8),
            (1.42, 481, 1.2), (0.02, 291, 1.2), (-0.39, 503, 1.2), (0.97, 386, 1.2), (0.66, 311, 1.2)]
FINAL = [(2.0, 250), (2.0, 400), (2.0, 550)]
OCEAN0 = (1.0, 360)          # ISOMIP+ Ocean0 linear profile (d = 720 m): profile test only

OUTLINE = json.loads((REPO / "exp2_isomip_cavity/figures/src/isomip_outline.json").read_text())

cmds: list[dict] = []


# ----------------------------------------------------------------------------------- primitives (as Exp. 2 Fig. 0)
def c(cmd: str, **params) -> None:
    cmds.append({"command": cmd, "params": params})


def fx(effect, **params):
    c("effect.apply", effect=effect, params=params)


def fill(color=None, gradient=None, stroke=None, weight=1.0, opacity=None, dash=None):
    c("paint.setFill", **({"gradient": gradient} if gradient else {"color": color} if color else {"none": True}))
    if stroke:
        c("paint.setStroke", color=stroke)
        c("stroke.set", weight=weight, **({"dash": dash} if dash else {}))
    else:
        c("paint.setStroke", none=True)
    if opacity is not None:
        c("object.setProps", opacity=opacity)


def grad(*stops, angle=90):
    """Linear gradient; with angle 90 offset 0 is at the bottom of the object, with angle 0 at its left."""
    return {"kind": "linear", "angle": angle, "stops": [{"offset": o, "color": col} for o, col in stops]}


def rect(x, y, w, h, radius=None, **kw):
    c("shape.rectangle", x=x, y=y, width=w, height=h, **({"radius": radius} if radius else {}))
    fill(**kw)


def poly(pts, close=True, **kw):
    c("path.create", d="M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in pts) + (" Z" if close else ""))
    fill(**kw)


def stroke_path(pts, color=INK, weight=1.0, dash=None, opacity=None, end_arrow=False):
    c("path.create", d="M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in pts))
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    st = dict(weight=weight, cap="round", join="round", **({"dash": dash} if dash else {}))
    if end_arrow:
        st["endArrow"] = "Triangle"
    c("stroke.set", **st)
    if opacity is not None:
        c("object.setProps", opacity=opacity)


def line(x1, y1, x2, y2, color=INK, weight=1.0, dash=None, opacity=None):
    stroke_path([(x1, y1), (x2, y2)], color, weight, dash, opacity)


def arrow(x1, y1, x2, y2, color=INK, weight=1.6, dash=None, both=False):
    c("shape.line", x1=x1, y1=y1, x2=x2, y2=y2)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    st = {"weight": weight, "endArrow": "Triangle", "cap": "round", **({"dash": dash} if dash else {})}
    if both:
        st["startArrow"] = "Triangle"
    c("stroke.set", **st)


def text(x, y, s, size=10, color=INK, style=None, justify=None, rotation=None):
    """Text with simple markup: _{..} is a subscript, ^{..} a superscript (set with text.setRangeStyle)."""
    plain, ranges, i = "", [], 0
    while i < len(s):
        if s[i] in "_^" and s[i + 1:i + 2] == "{":
            j = s.index("}", i)
            a = len(plain.encode())
            plain += s[i + 2:j]
            ranges.append((a, len(plain.encode()), s[i]))
            i = j + 1
        else:
            plain += s[i]
            i += 1
    p = {"x": x, "y": y, "text": plain, "size": size, "font": FONT, "color": color}
    if style:
        p["style"] = style
    if justify:
        p["justify"] = justify
    if rotation is not None:
        p["rotation"] = rotation
    c("text.create", **p)
    for a, b, k in ranges:
        c("text.setRangeStyle", start=a, end=b, size=round(size * 0.7, 2),
          baselineShift=round(-0.22 * size if k == "_" else 0.38 * size, 2))


def para(x, y, s, width, size=8.6, color=INK, style=None, lead=None):
    """Left-aligned paragraph wrapped at ``width`` characters; returns the y below its last line."""
    lead = lead or size * 1.32
    for k, ln in enumerate(textwrap.wrap(s, width)):
        text(x, y + k * lead, ln, size, color=color, style=style)
    return y + len(textwrap.wrap(s, width)) * lead


def marker(x, y, n, r=10):
    c("shape.ellipse", x=x - r, y=y - r, width=2 * r, height=2 * r)
    fill(color="#1d1d1d", stroke="#ffffff", weight=1.2)
    text(x, y + 4.4, str(n), 11, color="#ffffff", style="Bold", justify="center")


def dot(x, y, r, color, stroke="#ffffff", shape="circle"):
    if shape == "circle":
        c("shape.ellipse", x=x - r, y=y - r, width=2 * r, height=2 * r)
    elif shape == "triangle":
        c("path.create", d=f"M {x} {y - r * 1.2} L {x + r * 1.1} {y + r * 0.8} L {x - r * 1.1} {y + r * 0.8} Z")
    elif shape == "diamond":
        c("path.create", d=f"M {x} {y - r * 1.3} L {x + r * 1.1} {y} L {x} {y + r * 1.3} L {x - r * 1.1} {y} Z")
    elif shape == "star":
        c("shape.star", cx=x, cy=y, radius1=r * 1.5, radius2=r * 0.65, points=5)
    fill(color=color, stroke=stroke, weight=0.8)


def panel(x, y, w, h, n, title, sub=None):
    rect(x + 3, y + 4, w, h, radius=8, color="#000000", opacity=10)               # soft shadow
    rect(x, y, w, h, radius=8, color=PAPER, stroke=PANEL_LINE, weight=1.0)
    fx("texture.texturizer", texture="canvas", scaling=70, relief=2, lightDirection="topLeft")
    marker(x + 18, y + 20, n)
    text(x + 36, y + 25, title, 12.5, style="Bold")
    if sub:
        text(x + 36, y + 40, sub, 8.8, color=GREY)


def pill(x, y, tag, color=PURP, w=34):
    rect(x, y - 11, w, 15, radius=3, color=color)
    text(x + w / 2, y, tag, 8.6, color="#ffffff", style="Bold", justify="center")


# ----------------------------------------------------------------------------------- sketches
def shelf_map(x, y, sc, stops=MELT, question=False, label=None):
    """Plan view of the ISOMIP+ ice shelf (x = 440-660 km, y = 0-80 km) with a sketched melt pattern."""
    X = lambda km: x + (km - 440) * sc
    Y = lambda km: y + (80 - km) * sc
    gl = list(zip(OUTLINE["grounding_line_x_km"], OUTLINE["y_km"]))
    front = OUTLINE["front_x_km"]
    rect(X(440), Y(80), 220 * sc, 80 * sc, color="#a9c7de", stroke="#6f8fa8", weight=0.6)    # open ocean
    poly([(X(440), Y(0))] + [(X(g), Y(yk)) for g, yk in gl] + [(X(440), Y(80))],
         color="#f2f5f8", stroke="#9aaab8", weight=0.6)                                     # grounded ice
    poly([(X(g), Y(yk)) for g, yk in gl] + [(X(front), Y(80)), (X(front), Y(0))],
         gradient=grad(*stops, angle=0), stroke="#6b5a50", weight=0.6, opacity=55 if question else None)
    fx("texture.grain", intensity=8, contrast=40, grainType="soft")
    if question:
        text(X(565), Y(40) + 9, "?", 26, color=INK, style="Bold", justify="center")
    if label:
        text(x + 110 * sc, y + 80 * sc + 13, label, 8.6, justify="center")


def profile_sketch(x, y, w, h):
    """Ramp temperature profile: cold above, Tbot below, thermocline centred at z_c."""
    rect(x, y, w, h, color="#ffffff", stroke=BOX_LINE, weight=0.7)
    T = lambda t: x + 8 + (t + 2.2) / 4.6 * (w - 16)
    Z = lambda z: y + 6 + z / 720 * (h - 12)
    zc, d = 330, 220
    stroke_path([(T(-1.9), Z(0)), (T(-1.9), Z(zc - d / 2)), (T(1.3), Z(zc + d / 2)), (T(1.3), Z(720))],
                color=WARM, weight=2.0)
    stroke_path([(T(-1.9), Z(zc - d / 2)), (T(2.1), Z(zc + d / 2)), (T(2.1), Z(720))], color=WARM, weight=0.9,
                dash=[2, 2], opacity=70)
    line(x, Z(zc), x + w, Z(zc), color=GREY, weight=0.7, dash=[3, 2])
    text(x + w - 3, Z(zc) - 3, "z_{c}", 9, color=INK, style="Bold", justify="right")
    arrow(T(1.3), Z(640), T(2.1), Z(640), color=WARM, weight=1.0, both=True)
    text(T(1.7), Z(640) - 5, "T_{bot}", 8.6, color=WARM, style="Bold", justify="center")
    text(T(-1.9) + 3, Z(40) + 3, "cold", 7.5, color=COOL)
    text(x + w / 2, y - 5, "temperature profile", 8.2, color=GREY, justify="center")
    text(x - 4, y + h / 2, "depth", 7.5, color=GREY, justify="center", rotation=-90)


def draft_sketch(x, y, w, h):
    """Section of the ice shelf: the draft (ice base) scaled by s = 0.8, 1, 1.2."""
    rect(x, y, w, h, gradient=grad((0, "#1f5480"), (1, "#8fbbdc")), stroke=BOX_LINE, weight=0.7)
    sea = y + 18
    gl, fr = x + 14, x + w * 0.74
    base = lambda u, s: sea + s * (62 - 52 * u)        # draft below sea level, u = 0 at the grounding line
    pts = lambda s: [(gl + (fr - gl) * k / 20, base(k / 20, s)) for k in range(21)]
    poly([(x, y), (x + w * 0.74, y + 12), (fr, sea - 4)] + pts(1.0)[::-1] + [(x, base(0, 1.0) + 8)],
         gradient=grad((0, "#d3e3ef"), (1, "#fbfdff")), stroke="#4f6f88", weight=1.0)
    fx("texture.texturizer", texture="canvas", scaling=60, relief=2, lightDirection="topLeft")
    for s in (0.8, 1.2):
        stroke_path(pts(s), color="#ffffff", weight=1.1, dash=[2.5, 2])
    text(x + 22, y + 30, "ice shelf", 8, color=INK, style="Bold")
    text(x + w - 5, y + h - 30, "ocean", 8, color="#ffffff", style="Bold", justify="right")
    poly([(x, y + h), (x, base(0, 1.0) + 8), (gl + 10, base(0, 1.3) + 6), (x + w, y + h - 10), (x + w, y + h)],
         color=BED, stroke="#5e4d3d", weight=0.6)
    fx("texture.texturizer", texture="sandstone", scaling=80, relief=3, lightDirection="top")
    text(gl + 34, base(0.3, 1.2) + 12, "draft × s", 7.6, color="#ffffff", style="Bold")
    text(x + w / 2, y - 5, "ice-shelf section", 8.2, color=GREY, justify="center")


def net_icon(x, y, w, h, color=INK):
    cols = [(0, 3), (0.5, 4), (1, 2)]
    pts = [[(x + f * w, y + h * (k + 0.5) / n) for k in range(n)] for f, n in cols]
    for a, b in zip(pts, pts[1:]):
        for p in a:
            for q in b:
                line(*p, *q, color=color, weight=0.35, opacity=60)
    for col in pts:
        for px, py in col:
            dot(px, py, 2.2, color, stroke="#ffffff")


def power_icon(x, y, w, h, n=2.0, color=COOL):
    pts = [(x + w * k / 16, y + h - h * (k / 16) ** n) for k in range(17)]
    stroke_path(pts, color=color, weight=1.5)


# ----------------------------------------------------------------------------------- (1) what an emulator does
def emulator_panel(x0=30, y0=70, w=1040, h=262):
    panel(x0, y0, w, h, 1, "What an emulator does",
          "The 29 runs of Exp. 2 are the data: each pairs three numbers with one map of the basal melt rate.")
    # inputs
    profile_sketch(x0 + 30, y0 + 78, 92, 118)
    draft_sketch(x0 + 140, y0 + 78, 130, 118)
    text(x0 + 150, y0 + 218, "Inputs: three numbers per run", 9.6, style="Bold", justify="center")
    text(x0 + 150, y0 + 233, "T_{bot}, deep-water temperature;  z_{c}, thermocline depth;  s, draft scale", 8.2,
         color=GREY, justify="center")
    # two routes
    bx, bw, bh = x0 + 330, 180, 56
    ya, yb = y0 + 82, y0 + 158
    arrow(x0 + 280, y0 + 130, bx - 4, ya + bh / 2, color=PANEL_LINE, weight=2.2)
    arrow(x0 + 280, y0 + 150, bx - 4, yb + bh / 2, color=PANEL_LINE, weight=2.2)
    rect(bx, ya, bw, bh, radius=6, color="#ffffff", stroke=COOL, weight=1.4)
    text(bx + bw / 2, ya + 22, "MITgcm ocean model", 10.5, style="Bold", color=COOL, justify="center")
    text(bx + bw / 2, ya + 38, "solves the flow in the cavity", 8.4, color=GREY, justify="center")
    text(bx + bw / 2, ya + 50, "hours per run on 5 cores", 8.4, color=GREY, justify="center")
    rect(bx, yb, bw, bh, radius=6, color="#fff6ee", stroke=WARM, weight=1.4, dash=[4, 2.5])
    text(bx + bw / 2, yb + 22, "emulator (M0–M7)", 10.5, style="Bold", color=WARM, justify="center")
    text(bx + bw / 2, yb + 38, "fitted on the training runs only", 8.4, color=GREY, justify="center")
    text(bx + bw / 2, yb + 50, "a fraction of a second per map", 8.4, color=GREY, justify="center")
    mx, sc = bx + bw + 52, 0.78
    arrow(bx + bw + 4, ya + bh / 2, mx - 6, ya + bh / 2, color=COOL, weight=2.0)
    arrow(bx + bw + 4, yb + bh / 2, mx - 6, yb + bh / 2, color=WARM, weight=2.0, dash=[4, 2.5])
    shelf_map(mx, ya - 3, sc)
    text(mx + 110 * sc, ya - 9, "simulated melt map (the target)", 8.6, justify="center")
    shelf_map(mx, yb - 3, sc, question=True, label="predicted melt map")
    # melt key
    kx = mx + 30
    rect(kx, y0 + 244, 70, 7, gradient=grad(*[(1 - o, col) for o, col in MELT][::-1], angle=0))
    text(kx - 4, y0 + 251, "low", 7.2, color=GREY, justify="right")
    text(kx + 74, y0 + 251, "high melt (sketch)", 7.2, color=GREY)
    # compare
    cx = mx + 220 * sc + 22
    arrow(cx, ya + 30, cx, yb + 26, color=INK, weight=1.3, both=True)
    text(cx + 12, y0 + 92, "Compared on runs the", 9.4, style="Bold")
    text(cx + 12, y0 + 105, "emulator has not seen:", 9.4, style="Bold")
    for k, s in enumerate(("error of the whole map (RMSE)",
                           "error of the total melt of the shelf",
                           "error near the grounding line, where",
                           "    the melt is strongest",
                           "for M5–M7, also: is the stated",
                           "    uncertainty honest? (coverage, CRPS)")):
        text(cx + 14, y0 + 124 + 14 * k, ("• " if not s.startswith(" ") else "  ") + s.strip(), 8.6)
    text(cx + 12, y0 + 222, "Metrics of Burgard et al. (2022), Eq. 32–35.", 7.8, color=GREY)


# ----------------------------------------------------------------------------------- (2) where each method is tested
def tests_panel(x0=30, y0=347, w=520, h=345):
    panel(x0, y0, w, h, 2, "Where each method is tested",
          "Each test hides some runs, fits on the rest and predicts the hidden ones.")
    px, py, pw, ph = x0 + 58, y0 + 66, 230, 186
    tmin, tmax, zmin, zmax = -0.8, 2.35, 215, 590
    X = lambda t: px + (t - tmin) / (tmax - tmin) * pw
    Y = lambda z: py + (z - zmin) / (zmax - zmin) * ph
    rect(px, py, pw, ph, color="#ffffff", stroke=BOX_LINE, weight=0.8)
    rect(X(-0.66), Y(232), X(1.6) - X(-0.66), Y(570) - Y(232), radius=10, color="#e1ecf6", stroke="#9fbcd6",
         weight=0.8)
    text(X(-0.12), Y(508), "training range", 7.6, color=COOL, style="Bold")
    rect(X(1.36), Y(232), X(1.6) - X(1.36), Y(570) - Y(232), color="#f2c6bb", opacity=75)
    for zc in (325, 475):
        line(X(-0.66), Y(zc), X(1.6), Y(zc), color="#7d6f6a", weight=0.8, dash=[3, 3])
    for t in (-0.5, 0.0, 0.5, 1.0, 1.5, 2.0):
        line(X(t), py + ph, X(t), py + ph + 4, color=INK, weight=0.8)
        text(X(t), py + ph + 14, f"{t:g}", 8.2, justify="center")
    for z in (250, 400, 550):
        line(px - 4, Y(z), px, Y(z), color=INK, weight=0.8)
        text(px - 7, Y(z) + 3, f"{z}", 8.2, justify="right")
    text(px + pw / 2, py + ph + 28, "T_{bot}: deep-water temperature (°C)", 8.8, justify="center")
    text(px - 36, py + ph / 2, "z_{c}: thermocline depth (m)", 8.8, justify="center", rotation=-90)
    for t, z in FACTORIAL:
        dot(X(t), Y(z), 4.2, TRAIN)
    for t, z, s in GEOMETRY:
        dot(X(t), Y(z), 4.0, GEOM, shape="triangle")
    c("shape.ellipse", x=X(0.5) - 9, y=Y(400) - 9, width=18, height=18)                # loro: one run hidden
    fill(stroke=INK, weight=1.0, dash=[2, 1.6])
    arrow(X(1.62), Y(400), X(1.92), Y(400), color=TEST, weight=1.4)
    for t, z in FINAL:
        dot(X(t), Y(z), 4.4, TEST, shape="star")
    dot(X(OCEAN0[0]), Y(OCEAN0[1]), 4.2, REF, shape="diamond")
    # marker legend
    ly = py + ph + 46
    dot(x0 + 30, ly, 4.0, TRAIN)
    text(x0 + 40, ly + 3, "15 runs on a grid, s = 1", 8.2)
    dot(x0 + 30, ly + 15, 4.0, GEOM, shape="triangle")
    text(x0 + 40, ly + 18, "10 runs with the draft scaled, s = 0.8 or 1.2", 8.2)
    # the tests, in plain words
    tx, ty = x0 + 312, y0 + 70
    text(tx, ty, "Inside the training range", 9.2, style="Bold", color=COOL)
    rows_in = [("ring", "loro: hide one run, predict it; repeat for each of the 25 runs.")]
    rows_out = [("dash", "zc_out: hide the shallowest, middle or deepest third of z_{c}."),
                ("band", "warm_out: fit on T_{bot} ≤ 1.0 °C, predict the 5 warmest runs."),
                ("diamond", "profile_test: Ocean0, the ISOMIP+ standard run: a straight-line profile, "
                            "no thermocline."),
                ("star", "final: 3 runs at T_{bot} = 2.0 °C, warmer than every training run; predicted once, "
                         "with a frozen setup.")]

    def row(kind, s, y):
        ix = tx + 6
        if kind == "ring":
            c("shape.ellipse", x=ix - 6, y=y - 9, width=12, height=12)
            fill(stroke=INK, weight=1.0, dash=[2, 1.6])
        elif kind == "dash":
            line(ix - 7, y - 3, ix + 7, y - 3, color="#7d6f6a", weight=1.0, dash=[3, 2])
        elif kind == "band":
            rect(ix - 6, y - 9, 12, 12, color="#f2c6bb")
        elif kind == "diamond":
            dot(ix, y - 3, 4.0, REF, shape="diamond")
        elif kind == "star":
            dot(ix, y - 3, 4.4, TEST, shape="star")
        name, rest = s.split(": ", 1)
        text(tx + 18, y, name, 8.6, style="Bold", color=TEST if kind == "star" else INK)
        return para(tx + 18, y + 12, rest, 37, 8.2) + 6

    y = ty + 17
    for k, s in rows_in:
        y = row(k, s, y)
    y += 4
    text(tx, y, "Outside it: extrapolation", 9.2, style="Bold", color=TEST)
    y += 17
    for k, s in rows_out:
        y = row(k, s, y)


# ----------------------------------------------------------------------------------- (3) the methods
def methods_panel(x0=565, y0=347, w=505, h=345):
    panel(x0, y0, w, h, 3, "Eight methods, from pure physics to pure data")
    ax0, ax1, ay = x0 + 22, x0 + w - 22, y0 + 64
    arrow(ax0 + 70, ay, ax1 - 70, ay, color=PANEL_LINE, weight=1.6, both=True)
    text(ax0, ay + 3, "physics", 9, style="Bold", color=COOL)
    text(ax1, ay + 3, "data", 9, style="Bold", color="#b06a2a", justify="right")
    chips = [("M1–M3", PHYS, "melt ∝ TF or TF^{2}", "one fitted constant"),
             ("M7", NN, "hybrid:", "physics × network"),
             ("M4", STAT, "regression per", "cell on T_{bot}, z_{c}, s"),
             ("M5", STAT, "Gaussian process", "on map patterns"),
             ("M6", NN, "5 neural networks", "10 local inputs")]
    cw, gap = 87, (ax1 - ax0 - 5 * 87) / 4
    for k, (name, col, a, b) in enumerate(chips):
        x = ax0 + k * (cw + gap)
        rect(x, ay + 14, cw, 62, radius=5, color=col, stroke=WARM if name == "M7" else BOX_LINE,
             weight=1.6 if name == "M7" else 0.7)
        text(x + cw / 2, ay + 30, name, 10, style="Bold", color=WARM if name == "M7" else INK, justify="center")
        text(x + cw / 2, ay + 47, a, 7.8, justify="center")
        text(x + cw / 2, ay + 59, b, 7.8, justify="center")
    text(ax0, ay + 94, "and a baseline, M0: the mean of the training maps, the score any method must beat.", 8.4,
         color=GREY)
    # the hybrid idea
    hx, hy, hw, hh = x0 + 16, y0 + 168, w - 32, 162
    rect(hx, hy, hw, hh, radius=6, color="#ffffff", stroke=WARM, weight=1.0)
    text(hx + 12, hy + 18, "The idea behind the hybrid M7", 10, style="Bold", color=WARM)
    text(hx + hw / 2, hy + 42, "melt  =  K · TF_{+}^{n}  ×  g_{θ}(geometry)", 14, style="Bold", justify="center")
    # left: the power law extends past the training range
    lx, ly, lw, lh = hx + 16, hy + 58, 80, 46
    rect(lx, ly, lw * 0.62, lh, color="#e1ecf6")
    line(lx, ly + lh, lx + lw, ly + lh, color=INK, weight=0.8)
    line(lx, ly, lx, ly + lh, color=INK, weight=0.8)
    pts = [(lx + lw * k / 20, ly + lh - lh * (k / 20) ** 2.5) for k in range(21)]
    stroke_path(pts[:13], color=COOL, weight=1.8)
    stroke_path(pts[12:], color=COOL, weight=1.8, dash=[3, 2])
    text(lx + lw / 2, ly + lh + 10, "TF", 7.4, color=GREY, justify="center")
    text(lx + 3, ly + 9, "trained", 6.8, color=COOL)
    para(lx + lw + 10, ly + 8, "Physics sets how much: a power law in the local thermal forcing TF, with K and n "
         "fitted. It continues smoothly into warmer water.", 30, 8.0)
    # right: the network sees geometry only
    rx = hx + hw / 2 + 12
    shelf_map(rx, ly + 2, 0.38, stops=((0, "#5b3f8f"), (0.5, "#a995d1"), (1, "#ece6f7")))
    para(rx + 220 * 0.38 + 10, ly + 8, "A network sets where: it sees only the geometry (draft, distances, "
         "slope, s), never a temperature.", 25, 8.0)
    para(hx + 12, hy + hh - 26, "So the network is never asked about a temperature it has not seen. M6, by contrast, "
         "takes temperature as an input and must extrapolate it.", 92, 8.2, color=INK, style="Bold")


# ----------------------------------------------------------------------------------- (4) the hypotheses
def hypotheses_panel(x0=30, y0=707, w=1040, h=205):
    panel(x0, y0, w, h, 4, "Three questions, fixed before any hidden run was predicted",
          "Thresholds in notes/pd0_preregistro.md. Which method ranks first is reported, not hypothesised.")
    cw, cy, ch = 330, y0 + 56, 120
    cards = [("H3a", "Does the physics help the hybrid", "extrapolate to warmer water?",
              "On the 3 final-test runs, M7 has a smaller error in the total melt than M6, on average and in at "
              "least 2 of the 3 runs."),
             ("H3b", "Does the hybrid lose little", "inside the training range?",
              "In loro, the map error (RMSE) of M7 is at most 1.10 times that of M6."),
             ("H3c", "Does the network ensemble know", "when it is wrong?",
              "In the 33 extrapolation cases, never a spread below 5 % together with an error above 20 % "
              "(a replication of Ockenden et al., 2026).")]
    for k, (tag, q1, q2, crit) in enumerate(cards):
        x = x0 + 15 + k * (cw + 12)
        rect(x, cy, cw, ch, radius=6, color="#ffffff", stroke=BOX_LINE, weight=0.8)
        pill(x + 10, cy + 20, tag)
        text(x + 52, cy + 20, q1, 9.6, style="Bold")
        text(x + 52, cy + 33, q2, 9.6, style="Bold")
        para(x + 12, cy + 56, crit, 36 if tag == "H3c" else 62, 8.4, color=GREY)
        if tag == "H3a":
            sx, sy = x + 24, cy + 88
            rect(sx, sy, 120, 10, gradient=grad((0, "#7ab0d9"), (1, "#e7a07f"), angle=0))
            text(sx, sy + 22, "–0.5 °C", 7.4, color=GREY, justify="center")
            text(sx + 120, sy + 22, "1.5 °C", 7.4, color=GREY, justify="center")
            line(sx + 122, sy + 5, sx + 145, sy + 5, color=TEST, weight=1.0, dash=[2, 2])
            dot(sx + 150, sy + 5, 4.6, TEST, shape="star")
            text(sx + 150, sy + 22, "2.0 °C", 7.4, color=TEST, justify="center")
            text(sx + 162, sy + 9, "the final test, outside it", 8, color=TEST)
            text(sx + 60, sy - 4, "training runs", 7.2, color=GREY, justify="center")
        if tag == "H3b":
            text(x + 12, cy + 96, "The price of the physics constraint:", 8.2, color=INK)
            text(x + 12, cy + 109, "the hybrid may fit the training range a bit worse.", 8.2, color=INK)
        if tag == "H3c":
            for j, (honest, lab) in enumerate(((True, "honest"), (False, "confidently wrong"))):
                sx, sy, sw, sh = x + 190 + j * 70, cy + 52, 60, 40
                rect(sx, sy, sw, sh, color="#faf8f4", stroke=BOX_LINE, weight=0.5)
                truth = [(sx + sw * u / 10, sy + sh * (0.62 - 0.32 * (u / 10) ** 1.4)) for u in range(11)]
                for m in range(5):
                    off = (m - 2) * (5.5 if honest else 1.1) + (0 if honest else -11)
                    stroke_path([(px_, py_ + off * (0.4 + 0.6 * i / 10)) for i, (px_, py_) in enumerate(truth)],
                                color="#d08a4a", weight=0.8, opacity=85)
                stroke_path(truth, color=INK, weight=1.6)
                text(sx + sw / 2, sy + sh + 11, lab, 7.4, color=INK, justify="center")
                text(sx + sw - 4, sy + 11, "✓" if honest else "✗", 10, color="#2e7d32" if honest else TEST,
                     style="Bold", justify="right")
            text(x + 190, cy + 112, "black: model;  orange: 5 members", 6.8, color=GREY)
    text(x0 + 15, y0 + h - 12, "Not read: differences below the noise floor of the runs themselves "
         "(0.157 m a^{-1} for the map, 1.3 % for the total melt).", 8.4, color=GREY)


def build() -> list[dict]:
    c("file.new", width=W, height=H, units="Points", backgroundContents="white", name="fig0_protocol")
    text(30, 34, "Experiment 3: can a fast emulator replace the ocean model, and can it be trusted outside its "
                 "training data?", 15, style="Bold")
    text(30, 52, "A guide to the protocol: what is predicted from what (1), how it is tested (2), by which methods "
                 "(3), and the questions fixed in advance (4).", 9.6, color=GREY)
    emulator_panel()
    tests_panel()
    methods_panel()
    hypotheses_panel()
    text(30, H - 14, "A schematic: no result of Experiment 3 is drawn. The melt maps, curves and ensemble members "
                     "are sketches; the ice-shelf outline is the ISOMIP+ grounding line.", 8.2, color=GREY)
    return cmds


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cli", default=os.environ.get("VECTORCRAFT_CLI"))
    a = ap.parse_args()
    cmd_file = HERE / "fig0_protocol.commands.json"
    cmd_file.write_text(json.dumps(build(), indent=1))
    if not a.cli:
        print(f"wrote {cmd_file}; pass --cli to render")
        return
    outs = [FIG / f"fig0_protocol.{e}" for e in ("png", "pdf")]
    args = [a.cli, "run", "--commands", str(cmd_file), "--ppi", "120"]
    for o in outs:
        args += ["--export", str(o)]
    r = subprocess.run(args, capture_output=True, text=True, env={**os.environ, "VECTORCRAFT_NO_PREFS": "1"})
    bad = [l for l in r.stdout.splitlines() if '"error"' in l or '"warning"' in l]
    if r.returncode or any('"error"' in l for l in bad):
        raise SystemExit(f"vectorcraft-cli failed (rc {r.returncode}):\n" + "\n".join(bad) + r.stderr[-2000:])
    for l in bad:
        print("warning:", l[:300])
    print("rendered", ", ".join(o.name for o in outs))


if __name__ == "__main__":
    main()
