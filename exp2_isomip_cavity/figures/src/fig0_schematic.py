"""Figure 0 (a-c) of Experiment 2: the setting and a process schematic of an ice-shelf cavity, drawn with VectorCraft.

Writes the VectorCraft command list (``fig0_schematic.commands.json``) and, if the VectorCraft CLI is available,
renders ``figures/fig0_schematic.pdf`` and a JPEG copy for the README (``.jpg``) with its headless batch mode:

    vectorcraft-cli run --commands fig0_schematic.commands.json --export fig0_schematic.png ...

Tool: VectorCraft (https://github.com/storytold/vectorcraft), CLI release v0.22.0, macOS universal.
Usage: fig0_schematic.py [--cli PATH_TO_vectorcraft-cli]   (or env VECTORCRAFT_CLI)

(a) and (b) are satellite images (fetch_sat.py); the place labels are positioned by projecting their latitude and
longitude with pyproj, in the projection of each image. (c) is a schematic, not to scale, of a section along the
ice flow (grounded ice at the left, as x = 320 km in the model; the open ocean at the right). No model output is
drawn anywhere in this figure. Coordinates are points, y down.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
SAT = HERE / "sat"
W, H = 1100, 1140
FONT = "Helvetica"

ICE, ICE_LINE, ICE_DEEP = "#eef6fb", "#6b8fa3", "#bcd8ea"
INK, WARM, COOL, PURP, GREY = "#222222", "#b8322a", "#2c6ea3", "#5b3f8f", "#555555"
MELT, MELT_D = "#b9a7d8", "#6f58a8"
ROCK = "#8f7f6d"

cmds: list[dict] = []


# ----------------------------------------------------------------------------------------------- primitives
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
    """Linear gradient; with angle 90 offset 0 is at the bottom of the object."""
    return {"kind": "linear", "angle": angle, "stops": [{"offset": o, "color": col} for o, col in stops]}


def rect(x, y, w, h, **kw):
    c("shape.rectangle", x=x, y=y, width=w, height=h)
    fill(**kw)


def path(d, **kw):
    c("path.create", d=d)
    fill(**kw)


def poly(pts, close=True, **kw):
    d = "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in pts) + (" Z" if close else "")
    path(d, **kw)


def line(x1, y1, x2, y2, color=INK, weight=1.0, dash=None, opacity=None):
    c("shape.line", x1=x1, y1=y1, x2=x2, y2=y2)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, cap="round", **({"dash": dash} if dash else {}))
    if opacity is not None:
        c("object.setProps", opacity=opacity)


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


def curve_arrow(d, color, weight, dash=None, opacity=None, both=False):
    c("path.create", d=d)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    st = dict(weight=weight, endArrow="Triangle", cap="round", join="round")
    if both:
        st["startArrow"] = "Triangle"
    if dash:
        st["dash"] = dash
    c("stroke.set", **st)
    if opacity is not None:
        c("object.setProps", opacity=opacity)


def stroke_path(d, color, weight, dash=None, opacity=None):
    c("path.create", d=d)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, cap="round", join="round", **({"dash": dash} if dash else {}))
    if opacity is not None:
        c("object.setProps", opacity=opacity)


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


def label(x, y, s, size=9.5, color=INK, style=None, pad=3, width=None, opacity=80):
    """Text on a translucent white box (for labels over imagery)."""
    w = width or 0.52 * size * len(s)
    rect(x - pad, y - size - pad + 2, w + 2 * pad, size + 2 * pad, color="#ffffff", opacity=opacity)
    text(x, y, s, size, color=color, style=style)


def marker(x, y, n, r=9):
    c("shape.ellipse", x=x - r, y=y - r, width=2 * r, height=2 * r)
    fill(color="#1d1d1d", stroke="#ffffff", weight=1.2)
    text(x, y + 4.2, str(n), 10.5, color="#ffffff", style="Bold", justify="center")


def north_arrow(x, y, ux, uy, color="#ffffff", L=26):
    """North arrow at (x, y) pointing along the unit vector (ux, uy) (page coordinates)."""
    arrow(x - ux * L / 2, y - uy * L / 2, x + ux * L / 2, y + uy * L / 2, color=color, weight=2.0)
    text(x + ux * (L / 2 + 9), y + uy * (L / 2 + 9) + 3.5, "N", 9.5, color=color, style="Bold", justify="center")


# ----------------------------------------------------------------------------------------------- (a), (b)
PIG_CRS = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=112.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
BBOX_ASE = (-2050000, -860000, -1200000, -160000)      # as fetch_sat.py
BBOX_PIG = (806000, -1459000, 986000, -1319000)
BBOX_ANT = (-3400000, -3400000, 3400000, 3400000)
# Place names: the US-ACAN (USA) entries of the SCAR Composite Gazetteer of Antarctica (placenames.aq, checked
# 2026-11-10): Pine Island Glacier 75.1667 S, 100.0 W; Thwaites Glacier 75.5 S, 106.75 W. Both are nominal points
# on the grounded glacier, not positions of the ice front.
PLACES = {"Pine Island Glacier": (-75.167, -100.0), "Thwaites Glacier": (-75.5, -106.75)}


def satellite_row(y0=64.0, hgt=400.0):
    from pyproj import Transformer
    to3031 = Transformer.from_crs(4326, 3031, always_xy=True)
    toPIG = Transformer.from_crs(4326, PIG_CRS, always_xy=True)
    pig2ant = Transformer.from_crs(PIG_CRS, 3031, always_xy=True)

    # (a) MODIS, Amundsen Sea Embayment: 1215 x 1000 px placed at 486 x 400 pt
    ax0, aw = 40.0, hgt * 1215 / 1000
    def A(x, y):                                          # EPSG:3031 metres -> page
        return (ax0 + (x - BBOX_ASE[0]) / (BBOX_ASE[2] - BBOX_ASE[0]) * aw,
                y0 + (BBOX_ASE[3] - y) / (BBOX_ASE[3] - BBOX_ASE[1]) * hgt)
    def Al(lat, lon):
        return A(*to3031.transform(lon, lat))
    c("file.place", path=str(SAT / "ase_modis_2025-12-01.jpg"), rect=[ax0, y0, aw, hgt], link=False)
    # footprint of (b): the four corners of its box, through PIG_CRS -> EPSG:3031
    bx = [(BBOX_PIG[0], BBOX_PIG[3]), (BBOX_PIG[2], BBOX_PIG[3]), (BBOX_PIG[2], BBOX_PIG[1]), (BBOX_PIG[0], BBOX_PIG[1])]
    poly([A(*pig2ant.transform(*p)) for p in bx], stroke="#ffd400", weight=1.8)
    lx, ly = A(*pig2ant.transform(BBOX_PIG[2], BBOX_PIG[3]))
    text(lx + 4, ly + 2, "(b)", 9.5, color="#ffd400", style="Bold")
    for name, (lat, lon) in PLACES.items():
        px, py = Al(lat, lon)
        c("shape.ellipse", x=px - 2.6, y=py - 2.6, width=5.2, height=5.2)
        fill(color=INK, stroke="#ffffff", weight=0.8)
        if name.startswith("Pine"):
            label(px + 8, py - 6, name, 9)
        else:
            label(px + 8, py + 12, name, 9)
    text(*A(-2000000, -500000), "Amundsen Sea", 11, color="#3c5f80", style="Bold")
    text(*A(-2000000, -522000), "(sea ice)", 9, color="#3c5f80")
    p0, p1 = A(-1745000, -385000), A(-1790000, -505000)       # leader from the polynya to its label
    line(*p0, *p1, color="#5aa0d8", weight=1.1)
    label(p1[0] - 60, p1[1] + 12, "Pine Island Bay: open water (polynya)", 8.5, color=COOL)
    label(*A(-1450000, -700000), "West Antarctic Ice Sheet", 9.5)
    # the second open-water area: 73.0-74.2 S, 110-118 W, north of the Dotson Ice Shelf (SCAR gazetteer, US-ACAN:
    # 74.4 S, 112.37 W) and the eastern Getz Ice Shelf, where Jenkins et al. (2018, Nat. Geosci.) place the
    # Amundsen Sea polynya, "the most reliable area of summer open water" of the region
    q0, q1 = A(-1655000, -655000), A(-1565000, -655000)
    line(*q0, *q1, color="#5aa0d8", weight=1.1)
    label(q1[0] + 2, q1[1] + 4, "Amundsen Sea Polynya: open water", 8.5, color=COOL)
    # scale bar 100 km and north arrow (north = away from the pole, here toward the page left)
    sx, sy = A(-1420000, -830000)
    sk = aw / ((BBOX_ASE[2] - BBOX_ASE[0]) / 1e3)                 # pt per km
    line(sx, sy, sx + 100 * sk, sy, color=INK, weight=2.2)
    text(sx, sy - 5, "100 km", 9)
    n1, n2 = Al(-74.0, -100.0), Al(-73.9, -100.0)
    d = math.hypot(n2[0] - n1[0], n2[1] - n1[1])
    north_arrow(*A(-1260000, -810000), (n2[0] - n1[0]) / d, (n2[1] - n1[1]) / d, color=INK)
    # the ISOMIP+ domain at the same scale (480 x 80 km), from the geometry file: grounded ice, ice shelf, open ocean
    ol = json.loads((HERE / "isomip_outline.json").read_text())
    dx0, dy0 = ax0 + 14, y0 + hgt - 80 * sk - 36
    rect(dx0 - 6, dy0 - 6, 480 * sk + 12, 80 * sk + 30, color="#ffffff", opacity=80)
    rect(dx0, dy0, 480 * sk, 80 * sk, color="#dfe6ea")                                   # grounded ice
    gl = [(dx0 + (g - 320) * sk, dy0 + (80 - yk) * sk) for g, yk in zip(ol["grounding_line_x_km"], ol["y_km"])]
    fr = dx0 + (ol["front_x_km"] - 320) * sk
    gl = sorted(gl, key=lambda p: p[1])
    poly([(gl[0][0], dy0)] + gl + [(gl[-1][0], dy0 + 80 * sk), (fr, dy0 + 80 * sk), (fr, dy0)], color="#a9cfe6")
    rect(fr, dy0, (800 - ol["front_x_km"]) * sk, 80 * sk, color="#2f6f9f")              # open ocean
    rect(dx0, dy0, 480 * sk, 80 * sk, stroke=INK, weight=1.0)
    text(dx0, dy0 + 80 * sk + 13, "ISOMIP+ domain, same scale: 480 × 80 km", 8.5)
    text(dx0 + 4, dy0 + 80 * sk / 2 + 3, "grounded", 7.5, color=INK)
    text((gl[40][0] + fr) / 2 + 6, dy0 + 80 * sk / 2 + 3, "ice shelf", 7.5, color=INK, justify="center")
    text((fr + dx0 + 480 * sk) / 2, dy0 + 80 * sk / 2 + 3, "ocean", 7.5, color="#ffffff", justify="center")
    # locator: Antarctica (Blue Marble with bathymetry), the box of (a) in red
    lsz, lx0, ly0 = 122.0, ax0 + 8, y0 + 8
    c("file.place", path=str(SAT / "antarctica_bluemarble.jpg"), rect=[lx0, ly0, lsz, lsz], link=False)
    rect(lx0, ly0, lsz, lsz, stroke="#ffffff", weight=1.2)
    def L_(x, y):
        return (lx0 + (x - BBOX_ANT[0]) / (BBOX_ANT[2] - BBOX_ANT[0]) * lsz,
                ly0 + (BBOX_ANT[3] - y) / (BBOX_ANT[3] - BBOX_ANT[1]) * lsz)
    p1, p2 = L_(BBOX_ASE[0], BBOX_ASE[3]), L_(BBOX_ASE[2], BBOX_ASE[1])
    rect(p1[0], p1[1], p2[0] - p1[0], p2[1] - p1[1], stroke="#ff3b30", weight=1.4)
    text(lx0 + lsz / 2, ly0 + lsz - 6, "Antarctica", 8, color="#ffffff", justify="center")
    text(ax0, y0 - 12, "(a)", 13, style="Bold")
    text(ax0 + 26, y0 - 12, "Amundsen Sea Embayment, West Antarctica", 10.5)

    # (b) Landsat 9, Pine Island Glacier, 1000 x 778 px placed at 514 x 400 pt; ice flows to the right
    bx0, bw = ax0 + aw + 20, hgt * 1000 / 778
    sb = bw / 1000
    Q = lambda px, py: (bx0 + sb * px, y0 + sb * py)
    c("file.place", path=str(SAT / "pig_landsat9_2025-12-01.jpg"), rect=[bx0, y0, bw, hgt], link=False)
    rect(bx0, y0, bw, hgt, stroke="#ffd400", weight=1.6)
    # line of the section (c), along the central flow line
    line(*Q(40, 282), *Q(985, 282), color="#ffffff", weight=1.4, dash=[6, 4])
    label(*Q(700, 300), "line of the section in (c)", 8.5, color=INK)
    arrow(*Q(120, 252), *Q(400, 252), color=COOL, weight=2.8)
    label(*Q(70, 232), "Pine Island Glacier: ice flows toward the sea", 9.5, color=COOL)
    label(*Q(420, 200), "ice shelf (floating)", 9.5)
    front = [(598, 148), (615, 170), (620, 198), (607, 222), (614, 250), (608, 280), (622, 305), (643, 332)]
    stroke_path("M " + " L ".join("%.1f %.1f" % Q(x + 6, y) for x, y in front), WARM, 2.2, dash=[5, 3])
    label(*Q(560, 375), "calving front (ice front)", 9.5, color=WARM)
    text(*Q(745, 120), "Pine Island Bay", 11, color="#ffffff", style="Bold")
    text(*Q(745, 145), "open water", 9.5, color="#ffffff")
    # scale bar 20 km (180 m per pixel) and north arrow, from two points of the same meridian
    sx, sy = Q(40, 740)
    line(sx, sy, sx + 20 / 0.18 * sb, sy, color=INK, weight=2.2)
    text(sx, sy - 5, "20 km", 9)
    def Ql(lat, lon):
        x, y = toPIG.transform(lon, lat)
        return Q((x - BBOX_PIG[0]) / 180.0, (BBOX_PIG[3] - y) / 180.0)
    n1, n2 = Ql(-75.0, -101.0), Ql(-74.9, -101.0)
    d = math.hypot(n2[0] - n1[0], n2[1] - n1[1])
    north_arrow(*Q(900, 690), (n2[0] - n1[0]) / d, (n2[1] - n1[1]) / d, color="#ffffff")
    text(bx0, y0 - 12, "(b)", 13, style="Bold")
    text(bx0 + 26, y0 - 12, "Pine Island Glacier ice shelf, Landsat 9, 30 m", 10.5)
    text(40, y0 + hgt + 16, "Imagery: (a) MODIS Terra true colour, 1 Dec 2025, NASA Worldview/GIBS; inset: NASA Blue Marble "
         "with shaded relief and bathymetry. (b) Landsat 9 OLI-2, 1 Dec 2025, courtesy of the U.S. Geological Survey; "
         "rotated so that the ice flows to the right.", 8.0, color=GREY)


# ----------------------------------------------------------------------------------------------- (c)
O = 520                      # top of the section
SL = O + 72                  # sea level
PC = O + 222                 # thermocline centre (open ocean)
BOT = O + 470                # bottom of the drawing
GLX, FRX = 300, 700          # grounding line, ice front (x)
FRB = O + 138                # ice-front draft


def base_y(x):
    """Ice-shelf base: steep near the grounding line, flatter toward the front (as the ISOMIP+ draft)."""
    f = ((FRX - x) / (FRX - GLX)) ** 2.2
    return FRB + (O + 392 - FRB) * f


def bed_y(x):
    pts = [(40, O + 384), (200, O + 386), (GLX, O + 393), (380, O + 412), (520, O + 410), (700, O + 398),
           (840, O + 372), (915, O + 356), (945, O + 364), (985, O + 430), (1015, O + 462), (1060, O + 468)]
    for (x1, y1), (x2, y2) in zip(pts, pts[1:]):
        if x1 <= x <= x2:
            t = (x - x1) / (x2 - x1)
            t = t * t * (3 - 2 * t)
            return y1 + t * (y2 - y1)
    return pts[-1][1]


def top_y(x):
    return SL - 28 + (x - GLX) / (FRX - GLX) * 18


def boundary_layer_inset(ix, iy, w=244, h=214):
    """Inset: the ice-ocean boundary layer and the three-equation melt law (AD16 Eq. 24-27)."""
    rect(ix, iy, w, h, color="#ffffff", stroke=ICE_LINE, weight=1.2)
    text(ix + 8, iy + 15, "ice–ocean boundary layer (AD16 Eq. 24–27)", 8.6, style="Bold")
    yi, yb, yo = iy + 62, iy + 84, iy + h - 6          # interface, base of the boundary layer, bottom
    rect(ix + 1, iy + 22, w - 2, yi - iy - 22, gradient=grad((0, ICE_DEEP), (1, ICE)))
    fx("texture.texturizer", texture="canvas", scaling=70, relief=2, lightDirection="topLeft")
    text(ix + 80, iy + 37, "ice: melts at the rate m_{w}, with", 8.3, color="#3d6075")
    text(ix + 80, iy + 50, "ρ_{fw}L m_{w} = the heat flux (24)", 8.3, color="#3d6075")
    rect(ix + 1, yi, w - 2, yb - yi, color="#d9cdef")
    rect(ix + 1, yb, w - 2, yo - yb, gradient=grad((0, "#f0a37a"), (1, "#f8d3b8")))
    fx("texture.grain", intensity=10, contrast=50, grainType="soft")
    line(ix + 1, yi, ix + w - 1, yi, color="#3d6075", weight=1.4)
    text(ix + 80, yi + 15, "interface: T_{b} = T_{f}(S_{b}, p_{b})  (25)", 8.3, color=PURP)
    # velocity profile: zero at the ice, sheared across the boundary layer (friction velocity u*)
    vx = ix + 12
    line(vx, yi + 1, vx, yo - 4, color=INK, weight=0.8)
    for k, L in enumerate((5, 12, 19, 24, 27, 29, 30)):
        y = yi + 6 + k * 17
        arrow(vx, y, vx + L, y, color=INK, weight=1.0)
    # heat up into the ice (red) and fresh meltwater down into the ocean (blue)
    arrow(ix + 54, yo - 10, ix + 54, yi - 3, color=WARM, weight=2.4)
    arrow(ix + 66, yi + 2, ix + 66, yb + 40, color=COOL, weight=1.8)
    tx = ix + 80
    text(tx, yb + 15, "u_{*}^{2} = C_{d}(u_{w}^{2} + u_{tidal}^{2})  (27)", 8.3)
    text(tx, yb + 31, "heat: ρc u_{*}Γ_{T}(T_{w} − T_{b})  (24)", 8.3, color=WARM)
    text(tx, yb + 47, "salt: u_{*}Γ_{S}(S_{w} − S_{b})  (26)", 8.3, color=COOL)
    text(tx, yb + 66, "thermal driving T_{w} − T_{f}", 8.3)
    text(tx, yb + 82, "T_{w}, S_{w}, u_{w}: the ocean", 8.3)
    text(tx, yb + 94, "just below the layer", 8.3)
    text(vx + 2, yo - 6, "u", 8.3, style="Bold")


def process_section():
    # sky
    rect(40, O, 1020, SL - O, gradient=grad((0, "#e3edf4"), (1, "#ffffff")))
    # water: cold, fresh above; Circumpolar Deep Water below; the thermocline between (fine grain)
    h = BOT - SL
    tc = (BOT - PC) / h
    rect(GLX - 4, SL, 1060 - GLX + 4, h, gradient=grad(
        (0, "#d9764c"), (0.35, "#e88f63"), (tc - 0.075, "#f1ad86"), (tc, "#d5cdca"), (tc + 0.07, "#a7cde4"),
        (1, "#d2e8f4")))
    fx("texture.grain", intensity=12, contrast=50, grainType="soft")
    # thermocline: a band (dotted edges) around PC, and its centre, the lever z_c
    for yy in (PC - 18, PC + 18):
        line(430, yy, 1060, yy, color="#7d6f6a", weight=0.7, dash=[1.5, 3], opacity=70)
    stroke_path(f"M 430 {PC} C 600 {PC - 6} 760 {PC + 6} 900 {PC} S 1000 {PC - 4} 1060 {PC}", "#4d4d4d", 1.3)
    # sea-surface ripples and a few sea-ice floes beyond the ice front
    c("path.create", d=f"M {FRX} {SL} L 1060 {SL}")
    c("paint.setFill", none=True)
    c("paint.setStroke", color="#3b7fb3")
    c("stroke.set", weight=1.2)
    fx("distort.zigZag", size=1.4, relative=False, ridges=36, points="smooth")
    for (x0, w) in ((FRX + 150, 34), (FRX + 196, 18), (FRX + 236, 46), (FRX + 330, 24)):
        path(f"M {x0} {SL} L {x0 + 3} {SL - 3} L {x0 + w - 3} {SL - 3} L {x0 + w} {SL} L {x0 + w - 2} {SL + 4} "
             f"L {x0 + 2} {SL + 4} Z", color="#f7fbfd", stroke="#8fb3c9", weight=0.6)
        fx("distort.roughen", size=0.8, relative=False, detail=20, points="corner", seed=int(x0))
    # sea floor: the continental shelf (deepening toward the ice, as in the Amundsen Sea), shelf break at right
    xs = list(range(40, 1061, 10))
    poly([(x, bed_y(x)) for x in xs] + [(1060, BOT + 8), (40, BOT + 8)], color=ROCK)
    fx("texture.texturizer", texture="sandstone", scaling=120, relief=4, lightDirection="top")
    fx("distort.roughen", size=2.2, relative=False, detail=16, points="corner", seed=7)
    # ice: grounded ice sheet and the floating ice shelf in one body
    top = [(40, O + 8)] + [(x, O + 8 + (top_y(GLX) - O - 8) * (1 - (1 - (x - 40) / (GLX - 40)) ** 2.0))
                           for x in range(50, GLX + 1, 10)]
    top += [(x, top_y(x)) for x in range(GLX + 20, FRX + 1, 20)]
    basep = [(FRX, FRB)] + [(x, base_y(x)) for x in range(FRX - 10, GLX - 1, -10)]
    bedp = [(x, bed_y(x)) for x in range(GLX - 10, 39, -10)]
    poly(top + basep + bedp, gradient=grad((0, ICE), (1, ICE_DEEP)), stroke=ICE_LINE, weight=1.4)
    fx("texture.texturizer", texture="canvas", scaling=80, relief=2, lightDirection="topLeft")
    for k in range(7):                                      # flow bands in the ice shelf
        x = GLX + 40 + k * 50
        y_top, y_bot = top_y(x) + 6, base_y(x) - 4
        line(x, y_top + (y_bot - y_top) * 0.35, x + 18, y_top + (y_bot - y_top) * 0.38, color="#c3dbea", weight=1.0)
    # (4) meltwater plume: a thin buoyant layer along the ice base, thickening toward the front
    up = [(x, base_y(x)) for x in range(GLX + 6, FRX + 1, 6)]
    dn = [(x, base_y(x) + 4 + 13 * (x - GLX) / (FRX - GLX)) for x in range(FRX, GLX + 5, -6)]
    poly(up + dn, gradient=grad((0, "#7a63b8"), (1, "#d5c9ef"), angle=0), opacity=85)
    fx("distort.roughen", size=1.6, relative=False, detail=24, points="smooth", seed=11)
    for x1, x2 in ((GLX + 24, GLX + 64), (GLX + 110, GLX + 180), (GLX + 240, GLX + 330)):
        d1, d2 = 3 + 9 * (x1 - GLX) / 400, 3 + 9 * (x2 - GLX) / 400
        curve_arrow(f"M {x1} {base_y(x1) + d1} C {x1 + (x2 - x1) * 0.4} {base_y(x1 + (x2 - x1) * 0.4) + d1} "
                    f"{x1 + (x2 - x1) * 0.7} {base_y(x1 + (x2 - x1) * 0.7) + d2} {x2} {base_y(x2) + d2}", "#ffffff", 2.0)
    # (3) heat to the ice base where it is in the warm water: small red arrows up into the base
    for x in (GLX + 16, GLX + 46, GLX + 82):
        yb = base_y(x)
        curve_arrow(f"M {x + 14} {yb + 34} C {x + 6} {yb + 26} {x + 3} {yb + 18} {x + 4} {yb + 7}", WARM, 1.6)
    # (1) CDW floods onto the shelf at depth below the thermocline and (2) reaches the grounding line
    curve_arrow(f"M 1044 {BOT - 16} C 1032 {BOT - 70} 990 {O + 338} 930 {O + 336} C 860 {O + 340} 720 {O + 364} "
                f"560 {O + 380}", WARM, 7)
    curve_arrow(f"M 548 {O + 384} C 460 {O + 392} 380 {O + 394} {GLX + 14} {O + 380}", WARM, 3.2)
    # (5) melt-driven overturning: inflow at depth, outflow along the base and out at the front
    cx5, cy5, r5 = 560, O + 300, 34              # clockwise on the page: outflow above, inflow below
    pts = [(cx5 + r5 * 1.5 * math.cos(math.radians(a)), cy5 + r5 * math.sin(math.radians(a))) for a in range(200, 500, 10)]
    curve_arrow("M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in pts), "#4d4d4d", 1.5, dash=[4, 3])
    # (6) the meltwater leaves the cavity: (i) detaches near its neutral-buoyancy depth, or (ii) rises to the surface
    path(f"M {FRX} {FRB + 4} C {FRX + 60} {FRB + 10} {FRX + 110} {PC - 46} {FRX + 170} {PC - 44} "
         f"C {FRX + 240} {PC - 42} {FRX + 300} {PC - 46} {FRX + 330} {PC - 40} C {FRX + 300} {PC - 24} "
         f"{FRX + 230} {PC - 22} {FRX + 170} {PC - 24} C {FRX + 100} {PC - 26} {FRX + 50} {FRB + 30} {FRX} {FRB + 20} Z",
         gradient=grad((0, "#7a63b8"), (1, "#d5c9ef"), angle=0), opacity=55)
    fx("distort.roughen", size=2.4, relative=False, detail=20, points="smooth", seed=5)
    curve_arrow(f"M {FRX + 150} {PC - 34} L {FRX + 290} {PC - 34}", PURP, 2.6, dash=[6, 3])
    curve_arrow(f"M {FRX + 4} {FRB - 4} C {FRX + 24} {FRB - 20} {FRX + 30} {SL + 18} {FRX + 60} {SL + 10} "
                f"L {FRX + 128} {SL + 10}", PURP, 2.6)
    text(FRX + 136, SL + 14, "(ii)", 9, color=PURP, style="Bold")
    text(FRX + 298, PC - 30, "(i)", 9, color=PURP, style="Bold")
    # design parameters: z_c (the lever) and T_bot
    zx = 975
    arrow(zx, SL + 3, zx, PC - 2, color=INK, weight=1.4, both=True)
    text(zx + 8, SL + 52, "z_{c}", 13, style="Bold")
    text(zx + 8, SL + 66, "thermocline", 8.5)
    text(zx + 8, SL + 78, "depth:", 8.5)
    text(zx + 8, SL + 90, "the lever", 8.5)
    text(850, PC + 4 - 10, "thermocline", 9, color="#4d4d4d")
    text(780, O + 300, "Circumpolar Deep Water (CDW): warm, salty;", 10, color="#7a2a14")
    text(780, O + 314, "temperature T_{bot}, varied", 10, color="#7a2a14")
    text(780, SL + 40, "cold, fresh surface water", 10, color="#1f4e79")
    # ice labels
    text(48, O + 46, "Grounded ice sheet", 12, style="Bold")
    text(GLX + 150, SL - 40, "Ice shelf (floating)", 12, style="Bold")
    text(FRX + 180, O + 22, "Continental shelf sea", 12, style="Bold")
    line(GLX, bed_y(GLX) + 2, GLX, bed_y(GLX) + 14, color="#ffffff", weight=1.6)
    text(GLX, bed_y(GLX) + 26, "grounding line", 8.5, color="#ffffff", justify="center")
    text(FRX + 6, SL - 22, "ice front", 8.5, color=ICE_LINE)
    text(430, O + 352, "cavity", 10, color="#7a2a14", style="Bold")
    text(918, bed_y(918) + 22, "shelf break", 8.5, color="#ffffff", justify="center")
    # inset: the ice-ocean boundary layer, linked to a point on the deep ice base
    ix, iy = 50, O + 136
    lxp = GLX + 70
    line(ix + 244, iy + 70, lxp, base_y(lxp) + 2, color="#455a66", weight=0.8, dash=[2, 2])
    c("shape.ellipse", x=lxp - 4, y=base_y(lxp) - 2, width=8, height=8)
    fill(stroke="#455a66", weight=1.0)
    boundary_layer_inset(ix, iy)
    # numbered processes
    for (x, y, n) in ((900, O + 322, 1), (GLX + 40, O + 405 - 20, 2), (GLX + 108, base_y(GLX + 108) + 32, 3),
                      (GLX + 210, base_y(GLX + 210) + 28, 4), (560, O + 300, 5), (FRX + 100, PC - 64, 6),
                      ):
        marker(x, y, n)
    text(20, O - 8, "(c)", 13, style="Bold")
    text(46, O - 8, "Processes beneath an ice shelf (a section along the ice flow)", 11)
    text(40, BOT + 22, "(c) is a schematic: not to scale (the cavity, a few hundred metres thick, is several tens of "
         "kilometres long). Process 7, rotation, is drawn in Figure 0d.", 9, color=GREY)
    # legend
    lx, ly = 40, BOT + 48
    items = (
        "warm, salty Circumpolar Deep Water flows onto the continental shelf at depth, below the thermocline",
        "the CDW reaches the grounding line, where the ice base is deepest and the freezing point lowest",
        "basal melting: heat and salt cross a thin boundary layer at the ice base (inset)",
        "the meltwater is buoyant and rises as a plume along the sloping ice base",
        "the melt drives an overturning: inflow at depth, outflow along the ice base",
        "the meltwater leaves the cavity: (i) it detaches near its neutral-buoyancy depth, or (ii) it rises to the surface",
        "Earth's rotation steers the circulation into currents along the cavity walls (Figure 0d)",
    )
    for k, s_ in enumerate(items):
        col, row = k % 2, k // 2
        x, y = lx + col * 520, ly + row * 24
        marker(x + 9, y, k + 1)
        text(x + 26, y + 4, s_, 9.5)


def build() -> list[dict]:
    c("file.new", width=W, height=H, units="Points", backgroundContents="white", name="fig0_schematic")
    satellite_row()
    process_section()
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
    outs = [FIG / f"fig0_schematic.{e}" for e in ("png", "pdf")]
    args = [a.cli, "run", "--commands", str(cmd_file), "--ppi", "150"]
    for o in outs:
        args += ["--export", str(o)]
    r = subprocess.run(args, capture_output=True, text=True, env={**os.environ, "VECTORCRAFT_NO_PREFS": "1"})
    bad = [l for l in r.stdout.splitlines() if '"error"' in l or '"warning"' in l]
    if r.returncode or any('"error"' in l for l in bad):
        raise SystemExit(f"vectorcraft-cli failed (rc {r.returncode}):\n" + "\n".join(bad) + r.stderr[-2000:])
    for l in bad:
        print("warning:", l[:300])
    # PDF: the embedded satellite images and textures are re-encoded as JPEG (200 ppi) with Ghostscript when it is
    # installed: 0.65 MB instead of 1.6 MB. Only image streams change (a /printer preset dropped the translucent
    # shapes, so it is not used). An SVG export (2 MB, the images inline in base64) is not written.
    import shutil
    if shutil.which("gs"):
        tmp = outs[1].with_suffix(".gs.pdf")
        subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7",
                        "-dDownsampleColorImages=true", "-dColorImageResolution=200", "-dColorImageDownsampleType=/Bicubic",
                        "-dAutoFilterColorImages=false", "-dColorImageFilter=/DCTEncode", "-dJPEGQ=85",
                        f"-sOutputFile={tmp}", str(outs[1])], check=True)
        tmp.replace(outs[1])
    from PIL import Image
    jpg = FIG / "fig0_schematic.jpg"
    Image.open(outs[0]).convert("RGB").save(jpg, quality=86, subsampling=0, optimize=True)
    outs[0].unlink()
    print("rendered", ", ".join(o.name for o in outs[1:] + [jpg]))


if __name__ == "__main__":
    main()
