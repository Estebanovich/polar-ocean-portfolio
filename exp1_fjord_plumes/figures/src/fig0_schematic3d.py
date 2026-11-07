"""Figure 0 (3-D): block diagram of a tidewater glacier and its fjord, drawn with VectorCraft.

Companion to fig0_schematic.py (same processes 1-7, same colours and textures). The scene is built in
3-D world coordinates and projected orthographically, so every face, arrow and plume is geometrically
consistent; faces are drawn back to front (painter's order) and the water is translucent so the plumes
on the ice face can be seen through it. Schematic: not to scale (vertical exaggeration, enlarged plumes).

World axes: x along the fjord (0 at the ice face, L at the mouth), y across it (0 = cut-away near side,
Wd = far wall), z up (0 = sea level). Units are arbitrary drawing units.

Usage: fig0_schematic3d.py [--cli PATH_TO_vectorcraft-cli]   (or env VECTORCRAFT_CLI)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fig0_schematic as F  # noqa: E402  (drawing helpers; they append to F.cmds)

FIG = HERE.parent
CW, CH = 1180, 1000                      # canvas (points)

# scene (drawing units)
L, WD, H = 780.0, 300.0, 260.0           # fjord length, width, depth
XG = -300.0                              # inland end of the drawn glacier
HC = 46.0                                # ice cliff above the waterline
ZP = -112.0                              # pycnocline (Polar Water above, Atlantic Water below)
ZO = -122.0                              # neutral-buoyancy (outflow) depth
SLAB = 48.0                              # thickness of the rock slab under the bed
YP = (0.36 * WD, 0.64 * WD)              # across-fjord positions of the two discharge portals
ALPHA = 0.30                             # half-cone spreading rate (schematic)

# camera: direction from the scene to the viewer
_e = (0.68, -0.62, 0.34)
_n = math.sqrt(sum(a * a for a in _e))
E = tuple(a / _n for a in _e)
_r = (-E[1], E[0], 0.0)
_rn = math.hypot(_r[0], _r[1])
R = (_r[0] / _rn, _r[1] / _rn, 0.0)
UP = (E[1] * R[2] - E[2] * R[1], E[2] * R[0] - E[0] * R[2], E[0] * R[1] - E[1] * R[0])
S, U0, V0 = 1.0, 0.0, 0.0               # scale and offset: set by fit() so the scene fills the frame

ICE, ICE_DEEP, ICE_LINE = F.ICE, F.ICE_DEEP, F.ICE_LINE
MW, MW_D = "#4a9fd0", "#1f5f99"
ROCK = "#8f7f6d"


def fit(box=(30, 75, CW - 90, CH - 140)):
    """Choose S, U0, V0 so the scene's bounding points fill `box` (left, top, right, bottom)."""
    global S, U0, V0
    S, U0, V0 = 1.0, 0.0, 0.0
    pts = [P(x, y, z) for x in (XG, L) for y in (0, WD) for z in (-H - SLAB, zs(XG) if x == XG else HC)]
    pts += [P(XG, WD + 150, 150)]
    us, vs = [q[0] for q in pts], [q[1] for q in pts]
    S = min((box[2] - box[0]) / (max(us) - min(us)), (box[3] - box[1]) / (max(vs) - min(vs)))
    U0 = box[0] - S * min(us)
    V0 = box[1] - S * min(vs)


def P(x, y, z):
    """World point -> page point (orthographic projection, y down on the page)."""
    return (U0 + S * (x * R[0] + y * R[1] + z * R[2]), V0 - S * (x * UP[0] + y * UP[1] + z * UP[2]))


def d_of(pts3, close=True):
    pts = [P(*p) for p in pts3]
    s = "M " + " L ".join(f"{u:.2f} {v:.2f}" for u, v in pts)
    return s + (" Z" if close else "")


def poly(pts3, **kw):
    F.path(d_of(pts3), **kw)


def grad(angle, *stops):
    return {"kind": "linear", "angle": angle,
            "stops": [{"offset": o, "color": col} for o, col in stops]}


def zs(x):
    """Glacier surface elevation."""
    return HC + 0.24 * (-x) + 6 * math.sin(-x / 37.0) if x < 0 else HC


def hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3:
        return pts
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def ribbon(path3, width, normal_axis, color, opacity=None, head=2.2, stroke=None):
    """A flat 3-D arrow along a polyline; its width lies along `normal_axis` ('y' horizontal ribbon,
    'z' vertical ribbon). Head length = head * width."""
    pts = [tuple(p) for p in path3]
    left, right = [], []
    hw = width / 2
    off = (0, hw, 0) if normal_axis == "y" else (0, 0, hw)
    # shorten the shaft for the head
    hx, hy, hz = pts[-1]
    px, py, pz = pts[-2]
    seg = math.dist(pts[-1], pts[-2])
    t = max(0.0, 1 - head * width / seg)
    base = (px + t * (hx - px), py + t * (hy - py), pz + t * (hz - pz))
    shaft = pts[:-1] + [base]
    for p in shaft:
        left.append((p[0] + off[0], p[1] + off[1], p[2] + off[2]))
        right.append((p[0] - off[0], p[1] - off[1], p[2] - off[2]))
    hl = (base[0] + 2.1 * off[0], base[1] + 2.1 * off[1], base[2] + 2.1 * off[2])
    hr = (base[0] - 2.1 * off[0], base[1] - 2.1 * off[1], base[2] - 2.1 * off[2])
    poly(left + [hl, pts[-1], hr] + right[::-1], color=color, opacity=opacity,
         **({"stroke": stroke, "weight": 0.8} if stroke else {}))


def screen_label(pt, s, size=10, color=F.INK, style=None, box=True):
    u, v = pt
    (F.label if box else F.text)(u, v, s, size, color=color, style=style)


# ---------------------------------------------------------------------------------------------------
def background():
    F.rect(0, 0, CW, CH, color="#ffffff")
    F.rect(0, 40, CW, 560, gradient=grad(90, (0, "#ffffff"), (1, "#dce8f2")))


def back_wall():
    """Far fjord wall: mountains rising from the waterline at y = WD, with snow on the peaks."""
    rnd = random.Random(11)
    ridge = []
    for x in range(int(XG), int(L) + 1, 30):
        h = (100 + 55 * math.sin(x / 140.0) + rnd.uniform(-20, 30)) * min(1.0, 0.15 + (L - x) / 260.0)
        ridge.append((x, WD + 150, h))
    poly([(L, WD, 0)] + ridge[::-1] + [(XG, WD, 0)],
         gradient=grad(90, (0, "#8c8378"), (0.6, "#a79e93"), (1, "#c8c2b8")))
    F.fx("texture.texturizer", texture="sandstone", scaling=110, relief=4, lightDirection="topLeft")
    F.fx("distort.roughen", size=2.0, relative=False, detail=14, points="corner", seed=21)
    for (x, y, h) in ridge:
        if h > 135:
            poly([(x - 30, y - 20, h - 42), (x, y, h + 1), (x + 28, y - 18, h - 38), (x + 6, y - 30, h - 52)],
                 color="#fbfdff", opacity=92)
            F.fx("distort.roughen", size=2, relative=False, detail=12, points="corner", seed=int(x) % 97)


def bed_top():
    poly([(XG, 0, -H), (L, 0, -H), (L, WD, -H), (XG, WD, -H)], gradient=grad(0, (0, "#7d6e5e"), (1, "#a08f7c")))
    F.fx("texture.texturizer", texture="sandstone", scaling=90, relief=3, lightDirection="top")


def glacier():
    xs = [XG + k * (0 - XG) / 24 for k in range(25)]
    # ice face (x = 0), whole height; the part under water is seen through the water later
    poly([(0, 0, -H), (0, WD, -H), (0, WD, HC), (0, 0, HC)], gradient=grad(90, (0, "#9ec6de"), (1, "#eaf5fb")),
         stroke=ICE_LINE, weight=1.0)
    F.fx("texture.texturizer", texture="canvas", scaling=80, relief=2, lightDirection="topLeft")
    for k in range(1, 12):                                     # vertical seracs / fractures on the face
        y = k * WD / 12 + (5 if k % 2 else -4)
        F.path(d_of([(0, y, HC - 2), (0, y + 3, HC - 30), (0, y - 2, -10)], close=False), stroke="#c6dceb", weight=0.9)
    # top surface
    top = [(x, 0, zs(x)) for x in xs] + [(x, WD, zs(x) + 6) for x in xs[::-1]]
    poly(top, gradient=grad(0, (0, "#e3eef6"), (1, "#f8fcfe")), stroke=ICE_LINE, weight=0.8, opacity=93)
    F.fx("texture.texturizer", texture="canvas", scaling=70, relief=2, lightDirection="topLeft")
    # front cut face (y = 0): an ice section
    poly([(x, 0, zs(x)) for x in xs] + [(0, 0, -H), (XG, 0, -H)],
         gradient=grad(90, (0, ICE_DEEP), (1, ICE)), stroke=ICE_LINE, weight=1.2, opacity=92)
    F.fx("texture.texturizer", texture="canvas", scaling=80, relief=2, lightDirection="topLeft")
    for k in range(1, 9):                                      # annual layers bending toward the front
        z0 = -H + k * (H + HC) / 9
        F.path(d_of([(XG, 0, z0 + 0.10 * (-XG) * 0.24), (XG / 2, 0, z0 + 12), (-20, 0, z0 - 4)], close=False),
               stroke="#c3dbea", weight=0.9)
    # crevasses near the front, on the top surface
    rnd = random.Random(2)
    for k in range(26):
        x = -12 - rnd.uniform(0, 120)
        y = rnd.uniform(15, WD - 15)
        F.path(d_of([(x, y, zs(x) + 0.5), (x - 4, y + 22, zs(x) + 0.5)], close=False), stroke="#7e9fb4", weight=1.0)


LAKES = [(-250, 0.25 * WD), (-262, 0.72 * WD), (-160, 0.42 * WD)]


def _front(pts3):
    """Keep the part of a polyline that is not hidden behind the ice face (page u left of the face)."""
    u_face = P(0, 0, 0)[0]
    out = []
    for p in pts3:
        if P(*p)[0] > u_face - 1:
            break
        out.append(p)
    return out


def hydrology_inside():
    """Process 1, inside the ice (drawn before the glacier, seen through its translucent faces):
    moulins from the lakes to the bed, channels on the bed converging on the two portals, cavities."""
    for (lx, ly) in LAKES:
        top, bot = P(lx + 12, ly, zs(lx)), P(lx + 12, ly, -H)
        F.path(f"M {top[0] - 3} {top[1]} C {top[0] + 6} {top[1] + 80} {bot[0] - 6} {bot[1] - 80} {bot[0] - 3} {bot[1]} "
               f"L {bot[0] + 3} {bot[1]} C {bot[0] - 1} {bot[1] - 80} {top[0] + 12} {top[1] + 80} {top[0] + 3} {top[1]} Z",
               gradient=grad(90, (0, MW_D), (1, MW)), opacity=80)
    routes = [[(-238, 0.25 * WD), (-160, 0.30 * WD), (-70, 0.34 * WD), (0, YP[0])],
              [(-250, 0.72 * WD), (-150, 0.68 * WD), (-60, 0.65 * WD), (0, YP[1])],
              [(-148, 0.42 * WD), (-100, 0.33 * WD)]]
    for k, rt in enumerate(routes):
        dense = []
        for (xa, ya), (xb, yb) in zip(rt[:-1], rt[1:]):
            dense += [(xa + (xb - xa) * t / 12, ya + (yb - ya) * t / 12, -H + 2) for t in range(12)]
        vis = _front(dense + [(rt[-1][0], rt[-1][1], -H + 2)])
        if len(vis) < 2:
            continue
        F.path(d_of(vis, close=False), stroke=MW_D if k < 2 else MW, weight=5.0 if k < 2 else 2.8)
        F.c("object.setProps", opacity=85)
        F.path(d_of(vis, close=False), stroke="#bfe3f7", weight=1.4 if k < 2 else 0.9)
    for (lx, ly) in [(-290, 0.10 * WD), (-280, 0.90 * WD), (-200, 0.12 * WD)]:   # distributed cavities
        for j in range(3):
            c0 = P(lx + j * 22, ly + j * 6, -H + 1)
            F.c("shape.ellipse", x=c0[0] - 4, y=c0[1] - 2, width=8, height=4)
            F.fill(color=MW, opacity=85)


def hydrology_surface():
    """Process 1 at the surface: meandering streams into supraglacial lakes."""
    for (lx, ly) in LAKES:
        z = zs(lx) + 0.8
        poly([(lx + 26 * math.cos(t), ly + 18 * math.sin(t), z) for t in [k * math.pi / 14 for k in range(28)]],
             color=MW, stroke=MW_D, weight=0.7)
    rnd = random.Random(5)
    for (lx, ly) in LAKES:
        for k in range(3):
            x0, y0 = max(XG + 5, lx - 50 - rnd.uniform(0, 20)), min(WD - 5, max(5, ly + rnd.uniform(-60, 60)))
            pts = [(x0 + (lx - 20 - x0) * t, y0 + (ly - y0) * t + 7 * math.sin(9 * t + k), 0) for t in
                   [i / 10 for i in range(11)]]
            pts = [(x, y, zs(x) + 0.8) for x, y, _ in pts]
            F.curve_arrow(d_of(pts, close=False), MW, 1.0)


def cone_points(yp, z, r, n=20):
    return [(r * math.sin(t), yp + r * math.cos(t), z) for t in [k * math.pi / (n - 1) for k in range(n)]]


def plumes():
    """Processes 2-3: two half-cone plumes on the ice face, merging near the top (x0 < x0c), with jets,
    portals, edge billows and entrainment; meltwater streaks on the face between them."""
    ztop = ZO - 6
    for yp in YP:
        # portal on the ice face
        pts = [(0, yp + 9 * math.cos(t), -H + 13 * math.sin(t)) for t in [k * math.pi / 16 for k in range(17)]]
        poly(pts, color="#0f3f6b")
    for k, yp in enumerate(YP):
        samples = []
        for i in range(14):
            z = -H + (ztop + H) * i / 13
            samples += [P(*p) for p in cone_points(yp, z, 8 + ALPHA * (z + H))]
        hp = hull(samples)
        F.path("M " + " L ".join(f"{u:.2f} {v:.2f}" for u, v in hp) + " Z",
               gradient=grad(90, (0, "#6a52ad"), (1, "#c3b3ea")), stroke="#ffffff", weight=0.8, opacity=88)
        F.fx("distort.roughen", size=2.6, relative=False, detail=24, points="smooth", seed=31 + k)
    for yp in YP:                                              # cross-section rings: the cone shape
        for f in (0.35, 0.62, 0.88):
            z = -H + (ztop + H) * f
            F.path(d_of(cone_points(yp, z, 8 + ALPHA * (z + H)), close=False), stroke="#ffffff", weight=0.8)
            F.c("object.setProps", opacity=60)
    for yp in YP:                                              # fresh jets at the portals
        samples = []
        for i in range(6):
            z = -H + 46 * i / 5
            samples += [P(*p) for p in cone_points(yp, z, 6 + 0.18 * (z + H))]
        hp = hull(samples)
        F.path("M " + " L ".join(f"{u:.2f} {v:.2f}" for u, v in hp) + " Z",
               gradient=grad(90, (0, "#3fa6e0"), (1, "#b7a6e2")), stroke="#ffffff", weight=0.6, opacity=92)
    # billows on the outer (fjord-side) silhouette, growing with height
    for k, yp in enumerate(YP):
        for i, f in enumerate((0.32, 0.52, 0.72, 0.9)):
            z = -H + (ztop + H) * f
            r = 8 + ALPHA * (z + H)
            # rightmost projected point of the cross-section = outer silhouette
            u, v = max((P(*p) for p in cone_points(yp, z, r)), key=lambda q: q[0])
            rb = 0.20 * r + 3
            cx, cy = u - 0.1 * rb, v
            F.c("shape.ellipse", x=cx - rb * 1.05, y=cy - rb, width=2.1 * rb, height=2 * rb)
            F.fill(color="#b4a3df" if f > 0.5 else "#9c86cf", opacity=92)
            F.spiral(cx, cy, rb * 0.92, 1.8, math.pi, "#ffffff", 1.2)
            F.spiral(cx, cy, rb * 0.92, 1.8, 0.0, "#c4622f", 1.0)
            if k == 1 and i < 3:                               # entrainment into the billow underside
                F.curve_arrow(f"M {cx + 46} {cy + rb + 14} C {cx + 30} {cy + rb + 12} {cx + 14} {cy + rb + 6} "
                              f"{cx + 0.5 * rb} {cy + 0.9 * rb}", "#b8322a", 1.6)
    for yp in YP:                                              # rising arrows inside each cone
        a, b = P(4, yp, -H + 30), P(4, yp, ztop - 2)
        F.curve_arrow(f"M {a[0]} {a[1]} L {b[0]} {b[1]}", "#ffffff", 2.4)
    # meltwater boundary layer between and beside the plumes: thin rising streaks on the face
    for y in (0.12 * WD, 0.50 * WD, 0.88 * WD):
        F.path(d_of([(0.5, y, -H + 20), (0.5, y + 2, ZP + 30), (0.5, y - 1, -12)], close=False),
               stroke="#ffffff", weight=1.0)
        F.c("object.setProps", opacity=70)
    for (y, z) in ((0.10 * WD, -0.45 * H), (0.50 * WD, -0.30 * H), (0.86 * WD, -0.55 * H)):   # heat into the ice
        a, b = P(34, y, z), P(3, y, z)
        F.arrow(a[0], a[1], b[0], b[1], color=F.WARM, weight=1.4)


def interior_flows():
    """Processes 4-5 inside the water: outflow layer at neutral buoyancy, warm inflow at depth."""
    # outflow layer: a horizontal sheet spreading from the plumes down-fjord over the whole width
    poly([(0, 0, ZO + 12), (0.95 * L, 0, ZO + 12), (0.95 * L, WD, ZO + 12), (0, WD, ZO + 12)],
         gradient=grad(0, (0, "#7f66c0"), (0.35, "#b9a9e0"), (1, "#d9e0e8")), opacity=60)
    ribbon([(0.30 * WD, 0.5 * WD, ZO + 13), (0.45 * L, 0.5 * WD, ZO + 13)], 30, "y", "#5b3f8f", opacity=80)


def water():
    """Translucent water: top surface, the near cut face (y = 0) and the mouth (x = L) as sections."""
    sect = grad(90, (0, "#e98f62"), (0.40, "#f3b48e"), (0.54, "#cdd6dc"), (0.66, "#9ec9e2"), (1, "#c5e1f0"))
    xr = 0.34 * L                                             # near the ice: almost clear, to show the plumes
    nst = 30
    for k in range(nst):
        x0, x1 = xr * k / nst, xr * (k + 1) / nst
        poly([(x0, 0, 0), (x1, 0, 0), (x1, 0, -H), (x0, 0, -H)], gradient=sect, opacity=round(8 + 52 * (k + 0.5) / nst))
    poly([(xr, 0, 0), (L, 0, 0), (L, 0, -H), (xr, 0, -H)], gradient=sect, opacity=60)
    poly([(L, 0, 0), (L, WD, 0), (L, WD, -H), (L, 0, -H)], gradient=sect, opacity=90)
    F.fx("texture.grain", intensity=10, contrast=40, grainType="soft")
    for k in range(6):
        x0, x1 = xr * k / 6, xr * (k + 1) / 6
        poly([(x0, 0, 0), (x1 + 0.5, 0, 0), (x1 + 0.5, WD, 0), (x0, WD, 0)], color="#cfe6f3", opacity=int(10 + 32 * k / 6))
    poly([(xr, 0, 0), (L, 0, 0), (L, WD, 0), (xr, WD, 0)], gradient=grad(0, (0, "#cfe6f3"), (1, "#a9d1e8")), opacity=45)
    for k in range(1, 9):                                     # ripples on the surface
        y = k * WD / 9
        F.c("path.create", d=d_of([(30, y, 0), (L - 10, y, 0)], close=False))
        F.c("paint.setFill", none=True)
        F.c("paint.setStroke", color="#ffffff")
        F.c("stroke.set", weight=0.7)
        F.c("object.setProps", opacity=55)
        F.fx("distort.zigZag", size=1.2, relative=False, ridges=60, points="smooth")
    # sections of the outflow and the pycnocline on the cut faces
    poly([(0, 0, ZO + 10), (0.95 * L, 0, ZO + 4), (0.95 * L, 0, ZO - 2), (0, 0, ZO - 26)], color="#8f78c8", opacity=55)
    F.fx("distort.roughen", size=1.6, relative=False, detail=30, points="smooth", seed=41)
    # flows drawn on the cut face (section arrows): outflow at neutral buoyancy, warm inflow at depth
    ribbon([(0.24 * L, 0, ZO - 8), (0.66 * L, 0, ZO - 6)], 12, "z", "#5b3f8f", opacity=95)
    ribbon([(0.97 * L, 0, -0.80 * H), (0.55 * L, 0, -0.80 * H), (0.30 * L, 0, -0.74 * H)], 24, "z", "#b8322a", opacity=92)
    wav = [(x, 0, ZP + (7 * math.sin((x - L) / 34.0) * max(0.0, (x - 0.45 * L) / (0.55 * L)))) for x in
           [0.40 * L + i * 0.60 * L / 60 for i in range(61)]]
    F.path(d_of(wav, close=False), stroke="#4d4d4d", weight=1.3)
    F.path(d_of([(L, 0, ZP), (L, WD, ZP)], close=False), stroke="#4d4d4d", weight=1.3)
    # rock slab under the bed: near face and mouth face
    poly([(XG, 0, -H), (L, 0, -H), (L, 0, -H - SLAB), (XG, 0, -H - SLAB)], color=ROCK)
    F.fx("texture.texturizer", texture="sandstone", scaling=120, relief=4, lightDirection="top")
    F.fx("distort.roughen", size=2.0, relative=False, detail=20, points="corner", seed=7)
    poly([(L, 0, -H), (L, WD, -H), (L, WD, -H - SLAB), (L, 0, -H - SLAB)], color="#7b6c5b")
    F.fx("texture.texturizer", texture="sandstone", scaling=120, relief=4, lightDirection="top")
    # outline of the cut-away block
    F.path(d_of([(0, 0, 0), (L, 0, 0), (L, WD, 0)], close=False), stroke="#3b7fb3", weight=1.2)
    F.path(d_of([(L, 0, 0), (L, 0, -H - SLAB)], close=False), stroke="#555555", weight=0.8)


def iceberg(x, y, w, d, h, seed):
    rnd = random.Random(seed)
    n = 7
    base = [(x + w * (0.5 + 0.5 * math.cos(2 * math.pi * k / n)) * rnd.uniform(0.75, 1.0),
             y + d * (0.5 + 0.5 * math.sin(2 * math.pi * k / n)) * rnd.uniform(0.75, 1.0)) for k in range(n)]
    top = [(bx, by, h * rnd.uniform(0.75, 1.0)) for bx, by in base]
    bot = [(bx, by, 0) for bx, by in base]
    for k in range(n):                                       # visible side faces
        a, b = k, (k + 1) % n
        ex, ey = base[b][0] - base[a][0], base[b][1] - base[a][1]
        nx, ny = ey, -ex                                     # outward normal for a counter-clockwise base
        if nx * E[0] + ny * E[1] > 0:
            poly([bot[a], bot[b], top[b], top[a]], color="#d6e9f4", stroke="#8fb3c9", weight=0.6)
    poly(top, color="#f8fcfe", stroke="#8fb3c9", weight=0.6)
    F.fx("distort.roughen", size=1.0, relative=False, detail=18, points="corner", seed=seed)


def keels():
    """Submerged parts of the icebergs, drawn before the water so they are tinted by it."""
    for (x, y, w, d, h, s) in BERGS:
        poly([(x + 0.1 * w, y + 0.1 * d, 0), (x + 0.9 * w, y + 0.1 * d, 0), (x + 0.75 * w, y + 0.5 * d, -3.2 * h),
              (x + 0.25 * w, y + 0.5 * d, -3.0 * h)], color="#e3f1f8", opacity=70)


BERGS = [(14, 30, 36, 30, 12, 1), (60, 140, 28, 22, 9, 2), (24, 230, 40, 34, 14, 3), (110, 70, 22, 18, 7, 4),
         (150, 200, 30, 24, 10, 5), (230, 110, 20, 16, 6, 6), (90, 250, 18, 14, 6, 7), (320, 240, 24, 18, 8, 8)]


def surface_and_air():
    for (x, y, w, d, h, s) in sorted(BERGS, key=lambda b: -(b[1])):
        iceberg(x, y, w, d, h, s)
    rnd = random.Random(8)                                   # mélange: small bits near the face
    for k in range(40):
        x, y = rnd.uniform(4, 120), rnd.uniform(6, WD - 10)
        c0 = P(x, y, 0.5)
        F.c("shape.ellipse", x=c0[0] - 2.5, y=c0[1] - 1.2, width=5, height=2.4)
        F.fill(color="#ffffff", opacity=90)
    # wind (6) above the surface, and the surface water it pushes down-fjord
    for (x0, y0) in ((0.30 * L, 0.25 * WD), (0.48 * L, 0.70 * WD), (0.62 * L, 0.30 * WD)):
        ribbon([(x0, y0, 34), (x0 + 0.20 * L, y0, 34)], 14, "y", "#2b2b2b", opacity=85)
    ribbon([(0.52 * L, 0.5 * WD, -6), (0.86 * L, 0.5 * WD, -6)], 18, "y", "#2c6ea3", opacity=80)
    # shelf forcing (7): heave of the pycnocline at the mouth, and waves entering the fjord
    a, b = P(L + 40, 0.5 * WD, ZP + 46), P(L + 40, 0.5 * WD, ZP - 46)
    F.arrow(a[0], a[1], b[0], b[1], color=F.INK, weight=2.6, both=True)
    a, b = P(0.90 * L, 0, ZP - 22), P(0.62 * L, 0, ZP - 22)
    F.arrow(a[0], a[1], b[0], b[1], color="#4d4d4d", weight=1.6)


def annotations():
    # x0 between the portals, drawn at the foot of the ice face
    a, b = P(2, YP[0], -H + 4), P(2, YP[1], -H + 4)
    F.arrow(a[0], a[1], b[0], b[1], color="#ffffff", weight=1.8, both=True)
    m = P(2, 0.5 * WD, -H + 4)
    F.text(m[0] - 6, m[1] + 16, "x0", 11.5, color="#ffffff", style="Bold")
    # numbered processes
    marks = [(P(-150, 0, -150), 1), (P(30, YP[1] + 95, -H + 30), 2), (P(0, 0.07 * WD, -0.55 * H), 3),
             (P(0.33 * L, 0, ZO - 26), 4), (P(0.78 * L, 0, -0.66 * H), 5), (P(0.50 * L, 0.05 * WD, 40), 6),
             (P(L + 40, 0.5 * WD, ZP + 70), 7)]
    for (u, v), n in marks:
        F.marker(u, v, n)
    # region labels
    u, v = P(-180, 0.5 * WD, zs(-180) + 40)
    F.text(u, v, "Glacier", 14, style="Bold")
    u, v = P(0.62 * L, 0.85 * WD, 0)
    F.text(u, v, "Fjord", 14, style="Bold")
    u, v = P(L, 0.5 * WD, 0)
    F.text(u + 30, v - 30, "Shelf", 14, style="Bold")
    # labels written on the cut face, along the fjord axis
    a0, a1 = P(0, 0, 0), P(1, 0, 0)
    rot = math.degrees(math.atan2(a1[1] - a0[1], a1[0] - a0[0]))
    F.text(*P(0.42 * L, 0, -46), "Polar Water: cold, fresh", 11, color="#1f4e79", rotation=rot)
    F.text(*P(0.52 * L, 0, -0.93 * H), "Atlantic Water: warm, salty", 11, color="#7a2a14", rotation=rot)
    F.text(*P(0.74 * L, 0, ZP + 6), "pycnocline", 9.5, color="#4d4d4d", rotation=rot)
    q = P(0, YP[0], -H + 12)
    F.line(q[0], q[1], q[0] - 40, q[1] + 52, color=MW_D, weight=0.9)
    F.label(q[0] - 86, q[1] + 64, "portal", 9.5, color=MW_D, style="Bold")
    u, v = P(0, 0, -H - SLAB)
    F.text(u - 34, v + 16, "grounding line", 9.5, color=F.INK)
    F.line(*P(0, 0, -H), u, v + 4, color=F.INK, weight=0.8)
    for (pt, s_) in ((P(-250, 0.25 * WD, zs(-250)), "lake"), (P(-250, 0.72 * WD, -10), "moulin"),
                     (P(-200, 0.30 * WD, -H), "channel")):
        F.label(pt[0] + 8, pt[1] - 6, s_, 9.5, color=MW_D, style="Bold")
    u, v = P(0, WD, HC)
    F.line(u, v, u + 40, v - 26, color=F.INK, weight=0.9)
    F.label(u + 44, v - 28, "ice face (calving front)", 10.5, color="#2f5d78", style="Bold")
    u, v = P(140, 0.02 * WD, 14)
    F.text(u, v, "icebergs and mélange", 9.5, color="#466272")
    # down-fjord compass under the block
    a, b = P(0.12 * L, 0, -H - SLAB - 80), P(0.40 * L, 0, -H - SLAB - 80)
    F.arrow(a[0], a[1], b[0], b[1], color="#555555", weight=1.6)
    F.text(b[0] + 8, b[1] + 4, "down-fjord, toward the shelf", 9.5, color="#555555")
    # title and legend
    F.text(28, 34, "Processes at a tidewater glacier and in its fjord (3-D view)", 15, style="Bold")
    F.text(28, 54, "Schematic, not to scale: vertical exaggeration; plumes and icebergs enlarged. The near fjord "
                   "wall is cut away and the water drawn translucent; channels and moulins are shown through the ice.",
           9.5, color="#555555")
    items = ("subglacial discharge: lakes drain via moulins to channels at the bed and leave through portals",
             "buoyant plumes: rise along the ice; edge eddies engulf ambient water (entrainment)",
             "submarine melting: the plumes, and a thin meltwater layer, carry heat to the ice face",
             "neutral buoyancy: the plumes merge and spread down-fjord as an outflow layer",
             "warm Atlantic Water flows toward the glacier at depth, replacing the entrained water",
             "along-fjord wind drives a surface outflow and an exchange with the shelf",
             "shelf forcing: the pycnocline heaves over days and the motion enters the fjord")
    lx, ly = 40, CH - 108
    for k, s_ in enumerate(items):
        col, row = k % 2, k // 2
        x, y = lx + col * 680, ly + row * 24
        F.marker(x + 9, y, k + 1)
        F.text(x + 26, y + 4, s_, 10.5)


def build() -> list[dict]:
    F.cmds.clear()
    F.c("file.new", width=CW, height=CH, units="Points", backgroundContents="white", name="fig0_schematic3d")
    fit()
    background()
    back_wall()
    bed_top()
    glacier()
    hydrology_inside()
    hydrology_surface()
    interior_flows()
    keels()
    water()
    plumes()
    surface_and_air()
    annotations()
    return F.cmds


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cli", default=os.environ.get("VECTORCRAFT_CLI"))
    a = ap.parse_args()
    cmd_file = HERE / "fig0_schematic3d.commands.json"
    cmd_file.write_text(json.dumps(build(), indent=1))
    if not a.cli:
        print(f"wrote {cmd_file}; pass --cli to render")
        return
    outs = [FIG / f"fig0_schematic3d.{e}" for e in ("png", "svg", "pdf")]
    args = [a.cli, "run", "--commands", str(cmd_file), "--ppi", "150"]
    for o in outs:
        args += ["--export", str(o)]
    r = subprocess.run(args, capture_output=True, text=True, env={**os.environ, "VECTORCRAFT_NO_PREFS": "1"})
    bad = [l for l in r.stdout.splitlines() if '"error"' in l]
    if r.returncode or bad:
        raise SystemExit(f"vectorcraft-cli failed (rc {r.returncode}):\n" + "\n".join(bad) + r.stderr[-2000:])
    print("rendered", ", ".join(o.name for o in outs))


if __name__ == "__main__":
    main()
