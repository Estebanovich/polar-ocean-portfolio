"""Figure 0 of Experiment 4: Van Mijenfjorden and a plan-view process schematic, drawn with VectorCraft.

Writes the VectorCraft command list (``fig0_schematic.commands.json``), the facts the README caption quotes
(``figures/fig0_schematic.json``: image dates and sources) and, if the VectorCraft CLI is available, renders
``figures/fig0_schematic.pdf`` and a JPEG copy for the README (``.jpg``) with its headless batch mode:

    vectorcraft-cli run --commands fig0_schematic.commands.json --export fig0_schematic.png ...

Tool: VectorCraft (https://github.com/storytold/vectorcraft), CLI release v0.22.0, macOS universal.
Usage: fig0_schematic.py [--cli PATH_TO_vectorcraft-cli]   (or env VECTORCRAFT_CLI)

(a) is the Sentinel-2 image of fetch_sat.py in UTM 33N, with the real geometry of the experiment drawn on it:
the model domain and its open boundaries and sponge (cases.yaml configs.fjord_vmf), the strait sections that
close the inner-fjord box (figures/bathy_report.json, sounds.*.p0/p1), the inner-fjord polygon where the
late-winter runs put landfast ice (cases.yaml profiles.polygon_utm, the ice polygon of F2-F4), and the analysis
stations of V2a and V2b (figures/validation.json, cell indices -> UTM cell centres). Place names are positioned
from points of cases.yaml (the Akseloya island point) and by eye on the image. The inset is the MODIS image of
Svalbard (lat/lon), with the domain drawn from its UTM corners (pyproj).
(b) is a plan view, not to scale and drawn by hand, of the processes the experiment asks about. No model output
is drawn anywhere in this figure. Coordinates are points, y down.
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
EXP = FIG.parent
SAT = HERE / "sat"
W, H = 1100, 1395
FONT = "Helvetica"

INK, GREY = "#222222", "#555555"
WATER_D, WATER_L = "#2f6f9f", "#a7cde4"
ICE, ICE_LINE = "#f2f7fb", "#7d9db0"
LAND, LAND_D = "#d8d2c8", "#9c907f"
WIND, NI, KEL, RAD = "#3a3a3a", "#d1495b", "#5b3f8f", "#b8322a"
ACC = "#ffd400"

S2 = "vmf_s2_2024-03-25_domain.jpg"
MODIS = "svalbard_modis_2023-07-05.jpg"
BBOX_DOMAIN = (451000, 8590000, 559000, 8658000)        # as fetch_sat.py
SIZE_DOMAIN = (1404, 884)
BBOX_SVALBARD = (76.3, 9.0, 80.9, 30.0)                 # lat_min, lon_min, lat_max, lon_max
FACTS = {"sentinel2": {"file": S2, "date_text": "25 March 2024", "year": "2024",
                       "credit": "Contains modified Copernicus Sentinel data (2024), processed by ESA"},
         "modis": {"file": MODIS, "date_text": "5 July 2023", "credit": "NASA Worldview/GIBS, MODIS Terra"},
         "source": "exp4_arctic_bayiw/figures/src/fetch_sat.py, fig0_schematic.py"}

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
    return {"kind": "linear", "angle": angle, "stops": [{"offset": o, "color": col} for o, col in stops]}


def rect(x, y, w, h, **kw):
    c("shape.rectangle", x=x, y=y, width=w, height=h)
    fill(**kw)


def path(d, **kw):
    c("path.create", d=d)
    fill(**kw)


def poly(pts, close=True, **kw):
    path("M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in pts) + (" Z" if close else ""), **kw)


def line(x1, y1, x2, y2, color=INK, weight=1.0, dash=None, opacity=None):
    c("shape.line", x1=x1, y1=y1, x2=x2, y2=y2)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, cap="round", **({"dash": dash} if dash else {}))
    if opacity is not None:
        c("object.setProps", opacity=opacity)


def arrow(x1, y1, x2, y2, color=INK, weight=2.0, dash=None):
    c("shape.line", x1=x1, y1=y1, x2=x2, y2=y2)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, endArrow="Triangle", cap="round", **({"dash": dash} if dash else {}))


def curve_arrow(d, color, weight, dash=None, opacity=None):
    c("path.create", d=d)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, endArrow="Triangle", cap="round", join="round", **({"dash": dash} if dash else {}))
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
    """Text with simple markup: _{..} is a subscript, ^{..} a superscript."""
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
    w = width or 0.52 * size * len(s)
    rect(x - pad, y - size - pad + 2, w + 2 * pad, size + 2 * pad, color="#ffffff", opacity=opacity)
    text(x, y, s, size, color=color, style=style)


def marker(x, y, n, r=9):
    c("shape.ellipse", x=x - r, y=y - r, width=2 * r, height=2 * r)
    fill(color="#1d1d1d", stroke="#ffffff", weight=1.2)
    text(x, y + 4.2, str(n), 10.5, color="#ffffff", style="Bold", justify="center")


# ----------------------------------------------------------------------------------------------- (a)
def _inputs():
    import yaml
    cases = yaml.safe_load((EXP / "cases.yaml").read_text())["configs"]["fjord_vmf"]
    bathy = json.loads((FIG / "bathy_report.json").read_text())
    val = json.loads((FIG / "validation.json").read_text())
    return cases, bathy, val


def satellite(y0=56.0, ax0=40.0, aw=1020.0):
    from pyproj import Transformer
    cases, bathy, val = _inputs()
    hgt = aw * SIZE_DOMAIN[1] / SIZE_DOMAIN[0]
    X0, Y0, X1, Y1 = BBOX_DOMAIN

    def P(x, y):
        return (ax0 + (x - X0) / (X1 - X0) * aw, y0 + (Y1 - y) / (Y1 - Y0) * hgt)
    k = aw / ((X1 - X0) / 1e3)                                   # pt per km
    c("file.place", path=str(SAT / S2), rect=[ax0, y0, aw, hgt], link=False)
    # model domain, open boundaries (W, N, S) with the sponge, closed east side
    bx0, bx1, by0, by1 = cases["bbox_utm"]
    sp = cases["n_sponge"] * cases["dx_m"]
    for (xa, ya, xb, yb) in ((bx0, by0, bx0 + sp, by1), (bx0, by1 - sp, bx1, by1), (bx0, by0, bx1, by0 + sp)):
        (p, q), (r, s) = P(xa, yb), P(xb, ya)
        rect(p, q, r - p, s - q, color="#ffffff", opacity=38)
    (p, q), (r, s) = P(bx0, by1), P(bx1, by0)
    rect(p, q, r - p, s - q, stroke=INK, weight=3.6)
    rect(p, q, r - p, s - q, stroke=ACC, weight=2.0, dash=[7, 4])
    label(p + 8, q + 22, "model domain: open boundaries and sponge to the west, north and south", 9, width=330)
    # inner-fjord polygon: landfast ice in the late-winter runs (F2-F4)
    poly([P(x, y) for x, y in cases["profiles"]["polygon_utm"]], stroke="#38b6ff", weight=1.4, dash=[4, 3])
    # strait sections closing the inner-fjord box
    for nm, s_ in bathy["sounds"].items():
        a_, b_ = P(*s_["p0"]), P(*s_["p1"])
        line(*a_, *b_, color="#ff3b30", weight=3.2)
    # stations: V2a shelf station and the V2b south-coast stations (cell centres)
    def cell(j, i):
        return P(bx0 + (i + 0.5) * cases["dx_m"], by0 + (j + 0.5) * cases["dx_m"])
    pts = [cell(*val["V2a"]["station_ji"])] + [cell(*ji) for ji in val["V2b"]["stations"]["cells"]]
    for n, (x, y) in enumerate(pts):
        c("shape.ellipse", x=x - 3.4, y=y - 3.4, width=6.8, height=6.8)
        fill(color=ACC if n == 0 else "#ffffff", stroke=INK, weight=0.9)
    # place names
    ax_, ay_ = P(492000, 8627300)                                  # Akseloya (cases.yaml island point)
    line(ax_, ay_, ax_ + 48, ay_ + 64, color="#ffffff", weight=1.0)
    label(ax_ + 50, ay_ + 76, "Akseløya", 10, style="Bold")
    a_ = P(*bathy["sounds"]["Akselsundet"]["p1"])
    line(a_[0], a_[1], a_[0] - 30, a_[1] - 52, color="#ffffff", weight=1.0)
    label(a_[0] - 120, a_[1] - 56, "Akselsundet", 10, color="#b0201a", style="Bold")
    m_ = P(*bathy["sounds"]["Mariasundet"]["p1"])
    line(m_[0], m_[1], m_[0] + 40, m_[1] + 40, color="#ffffff", weight=1.0)
    label(m_[0] + 42, m_[1] + 52, "Mariasundet", 10, color="#b0201a", style="Bold")
    label(*P(515000, 8642800), "Van Mijenfjorden (inner fjord, landfast ice)", 11, style="Bold", width=270)
    text(*P(466000, 8628000), "Bellsund", 13, color="#ffffff", style="Bold")
    text(*P(466000, 8625600), "open water", 9.5, color="#ffffff")
    label(*P(505000, 8604500), "Van Keulenfjorden", 9.5)
    text(*P(454200, 8611000), "shelf", 10, color="#ffffff", style="Bold")
    sx_, sy_ = P(527000, 8598500)                                  # scale bar
    rect(sx_ - 6, sy_ - 18, 10 * k + 12, 26, color="#ffffff", opacity=80)
    line(sx_, sy_, sx_ + 10 * k, sy_, color=INK, weight=2.4)
    text(sx_, sy_ - 5, "10 km", 9, color=INK, style="Bold")
    nx_, ny_ = P(549000, 8597500)                                  # UTM grid north (convergence ~1 deg)
    arrow(nx_, ny_ + 14, nx_, ny_ - 14, color=INK, weight=2.0)
    text(nx_, ny_ - 19, "N", 9.5, color=INK, style="Bold", justify="center")
    # legend of the overlays
    lx, ly = ax0 + 10, y0 + hgt - 86
    rect(lx - 4, ly - 14, 300, 80, color="#ffffff", opacity=82)
    line(lx, ly + 44, lx + 22, ly + 44, color=INK, weight=3.6)
    line(lx, ly + 44, lx + 22, ly + 44, color=ACC, weight=2.0, dash=[7, 4])
    text(lx + 30, ly + 48, "model domain", 8.6)
    line(lx, ly - 4, lx + 22, ly - 4, color="#ff3b30", weight=3.2)
    text(lx + 30, ly, "strait sections: they close the inner-fjord box", 8.6)
    line(lx, ly + 12, lx + 22, ly + 12, color="#38b6ff", weight=1.8, dash=[4, 3])
    text(lx + 30, ly + 16, "inner-fjord polygon: landfast ice on its sea cells (late winter)", 8.6)
    for q_, (col, s_) in enumerate(((ACC, "shelf station (V2a)"), ("#ffffff", "south-coast stations (V2b)"))):
        c("shape.ellipse", x=lx + 8 + 120 * q_, y=ly + 26, width=6.8, height=6.8)
        fill(color=col, stroke=INK, weight=0.9)
        text(lx + 20 + 120 * q_, ly + 33, s_, 8.6)
    # inset: Svalbard (MODIS, lat/lon with the width scaled by cos 78.6)
    la0, lo0, la1, lo1 = BBOX_SVALBARD
    iw = 168.0
    ih = iw * 886 / 800
    ix, iy = ax0 + aw - iw - 10, y0 + 10
    c("file.place", path=str(SAT / MODIS), rect=[ix, iy, iw, ih], link=False)
    rect(ix, iy, iw, ih, stroke="#ffffff", weight=1.2)
    to_ll = Transformer.from_crs(32633, 4326, always_xy=True)
    corners = [(bx0, by0), (bx1, by0), (bx1, by1), (bx0, by1)]
    ll = [to_ll.transform(x, y) for x, y in corners]
    poly([(ix + (lo - lo0) / (lo1 - lo0) * iw, iy + (la1 - la) / (la1 - la0) * ih) for lo, la in ll],
         stroke="#ff3b30", weight=1.6)
    text(ix + iw / 2, iy + ih - 7, "Svalbard", 8.5, color="#ffffff", justify="center", style="Bold")
    text(ax0, y0 - 14, "(a)", 13, style="Bold")
    text(ax0 + 26, y0 - 14, "Van Mijenfjorden, Svalbard: Sentinel-2, 25 March 2024 (landfast ice), and the model "
         "domain", 10.5)
    text(ax0, y0 + hgt + 15, "Imagery: (a) contains modified Copernicus Sentinel data (2024), processed by ESA; "
         "inset: MODIS Terra true colour, 5 July 2023, NASA Worldview/GIBS. UTM 33N, grid north up.", 8.0,
         color=GREY)
    return y0 + hgt


# ----------------------------------------------------------------------------------------------- (b)
def plan_view(O):
    """Plan view, not to scale: the shelf and Bellsund at the left, the island across the mouth with two
    straits, the inner fjord to the right (its axis drawn horizontal), landfast ice over its eastern part."""
    L, R, B = 40, 1060, O + 440
    MX, HX = 300, 1000                           # mouth (island) and head of the fjord
    def yn(x):                                   # north coast of the inner fjord
        return O + 135 + 10 * math.sin((x - MX) / 70.0) + (x - MX) / (HX - MX) * 25
    def ys(x):                                   # south coast
        return O + 300 + 8 * math.sin((x - MX) / 55.0 + 1.0) - (x - MX) / (HX - MX) * 30
    # land: snow-covered mountains
    rect(L, O, R - L, B - O, gradient=grad((0, LAND_D), (1, LAND)))
    fx("texture.texturizer", texture="sandstone", scaling=110, relief=4, lightDirection="topLeft")
    # water: shelf and Bellsund (left) and the fjord
    xs = list(range(MX - 20, HX + 1, 10))
    north = [(x, yn(x)) for x in xs]
    south = [(x, ys(x)) for x in reversed(xs)]
    head = [(HX + 22 * math.sin(math.radians(a)), (yn(HX) + ys(HX)) / 2 + (ys(HX) - yn(HX)) / 2 * -math.cos(math.radians(a)))
            for a in range(10, 180, 15)]
    shelf_top, shelf_bot = O + 70, O + 380
    water = ([(L, shelf_top), (150, shelf_top - 8), (230, O + 95), (MX - 20, yn(MX - 20))] + north + head
             + south + [(230, O + 345), (150, shelf_bot + 6), (L, shelf_bot)])
    poly(water, gradient=grad((0, WATER_L), (0.45, "#6fa6cb"), (1, WATER_D), angle=0))
    fx("texture.grain", intensity=10, contrast=50, grainType="soft")
    fx("distort.roughen", size=1.4, relative=False, detail=22, points="smooth", seed=3)
    # the island across the mouth: two straits, Akselsundet (north, wider) and Mariasundet (south, narrow)
    iy0, iy1 = yn(MX) + 24, ys(MX) - 16
    poly([(MX - 8, iy0), (MX + 12, iy0 - 4), (MX + 18, (iy0 + iy1) / 2), (MX + 10, iy1 + 3), (MX - 10, iy1),
          (MX - 16, (iy0 + iy1) / 2)], gradient=grad((0, LAND_D), (1, LAND)), stroke="#857a6b", weight=0.8)
    fx("texture.texturizer", texture="sandstone", scaling=90, relief=3, lightDirection="topLeft")
    # (5) landfast ice over the eastern part of the fjord (it covers the whole inner fjord in F2-F4)
    IX = 760
    xi = list(range(IX, HX + 1, 10))
    ice = [(IX, yn(IX))] + [(x, yn(x)) for x in xi] + head + [(x, ys(x)) for x in reversed(xi)] + [(IX, ys(IX))]
    edge = [(IX - 6 + 7 * math.sin(k * 1.7), ys(IX) - (ys(IX) - yn(IX)) * k / 12) for k in range(13)]
    poly(ice + edge, gradient=grad((0, "#dce9f2"), (1, ICE)), stroke=ICE_LINE, weight=1.0)
    fx("texture.texturizer", texture="canvas", scaling=80, relief=2, lightDirection="topLeft")
    fx("distort.roughen", size=1.2, relative=False, detail=18, points="corner", seed=21)
    for x in (800, 860, 925):                                     # cracks / pressure ridges
        stroke_path(f"M {x} {yn(x) + 12} L {x + 14} {yn(x) + 40} L {x + 6} {ys(x) - 30} L {x + 20} {ys(x) - 10}",
                    "#c7d9e5", 1.0)
    # (1) the wind pulse, along the axis toward the mouth, over open water
    for yy, x1, x2 in ((O + 190, 720, 520), (O + 222, 700, 460), (O + 254, 720, 540)):
        arrow(x1, yy, x2, yy, color=WIND, weight=4.2)
    gx, gy = 520, O + 22                                          # the pulse glyph, tau(t)
    rect(gx - 8, gy - 8, 178, 72, color="#ffffff", opacity=85)
    line(gx, gy + 50, gx + 160, gy + 50, color=INK, weight=0.9)
    line(gx, gy + 50, gx, gy + 2, color=INK, weight=0.9)
    curve = [(gx + t, gy + 50 - 42 * math.exp(-0.5 * ((t - 30) / 10.0) ** 2) * (1 if t <= 30 else 1))
             for t in range(0, 161, 4)]
    curve = [(x, y if x - gx <= 30 else gy + 50 - 42 * math.exp(-0.5 * ((x - gx - 30) / 10.0) ** 2)) for x, y in curve]
    stroke_path("M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in curve), WIND, 2.0)
    text(gx + 4, gy + 6, "τ", 10, style="Bold")
    text(gx + 150, gy + 62, "t", 9.5, style="Bold")
    text(gx + 52, gy + 22, "one pulse of wind,", 8.6)
    text(gx + 52, gy + 34, "then free (no wind)", 8.6)
    # (2) near-inertial response: clockwise circles in open water
    for cx_, cy_ in ((432, O + 200), (432, O + 252)):
        r_ = 17
        pts = [(cx_ + r_ * math.cos(math.radians(a)), cy_ + r_ * math.sin(math.radians(a))) for a in range(200, 520, 12)]
        curve_arrow("M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in pts), NI, 2.2)
    # (3) internal Kelvin wave: trapped along the coasts, cyclonic (coast on the right)
    d = 13
    kel = ([(x, ys(x) - d) for x in range(390, HX - 10, 10)]
           + [(HX + (22 - d) * math.sin(math.radians(a)),
               (yn(HX) + ys(HX)) / 2 + ((ys(HX) - yn(HX)) / 2 - d) * math.cos(math.radians(a))) for a in range(20, 170, 15)]
           + [(x, yn(x) + d) for x in range(HX - 10, 390, -10)])
    curve_arrow("M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in kel), KEL, 2.6)
    for x in range(410, 960, 34):                                 # crests across the trapped band
        line(x, ys(x) - 2, x + 3, ys(x) - 24, color=KEL, weight=1.1, opacity=70)
        line(x + 12, yn(x + 12) + 2, x + 15, yn(x + 12) + 24, color=KEL, weight=1.1, opacity=70)
    bx_ = 690                                                    # trapping width ~ R1
    line(bx_, ys(bx_) - 1, bx_, ys(bx_) - 26, color=INK, weight=1.0)
    line(bx_ - 4, ys(bx_) - 26, bx_ + 4, ys(bx_) - 26, color=INK, weight=1.0)
    text(bx_ + 6, ys(bx_) - 9, "~R_{1}", 9, style="Bold")
    text(570, ys(570) + 28, "Kelvin wave: coast on its right", 9, color=KEL, style="Bold")
    arrow(470, ys(470) - d, 520, ys(520) - d, color=KEL, weight=2.6)
    # (4) radiation through the straits to the shelf
    for (y_, w_) in (((yn(MX) + iy0) / 2 + 2, 1.0), ((iy1 + ys(MX)) / 2, 0.7)):
        pts = [(MX + 30 - t, y_ + 4 * math.sin(t / 6.0)) for t in range(0, 150, 3)]
        curve_arrow("M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in pts), RAD, 2.6 * w_ + 0.6)
        for rr in (40, 70, 100):
            a0, a1 = 120, 240
            arc = [(MX - 20 + rr * math.cos(math.radians(a)), y_ + rr * 0.8 * math.sin(math.radians(a)))
                   for a in range(a0, a1 + 1, 6)]
            stroke_path("M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in arc), RAD, 0.9, opacity=45)
    # names
    text(70, O + 110, "shelf", 11, color="#ffffff", style="Bold")
    text(70, O + 124, "and Bellsund", 9.5, color="#ffffff")
    text(MX - 4, yn(MX) - 30, "Akselsundet", 9, color=INK, style="Bold", justify="center")
    text(MX - 4, ys(MX) + 32, "Mariasundet (narrow)", 9, color=INK, style="Bold", justify="center")
    text(MX + 24, (iy0 + iy1) / 2 + 4, "island", 8.5, color=INK)
    text(860, yn(860) - 12, "landfast ice (late winter)", 9.5, color="#1f4e79", style="Bold", justify="center")
    text(640, yn(640) - 14, "inner fjord", 11, color=INK, style="Bold")
    # numbered processes
    for (x, y, n) in ((735, O + 222, 1), (400, O + 226, 2), (540, ys(540) + 24, 3), (MX - 90, (yn(MX) + iy0) / 2 - 18, 4),
                      (900, (yn(900) + ys(900)) / 2, 5)):
        marker(x, y, n)
    text(20, O - 10, "(b)", 13, style="Bold")
    text(46, O - 10, "Processes after a wind pulse (plan view, not to scale; the fjord axis is drawn horizontal)", 11)
    # legend
    lx, ly = 40, B + 34
    items = (
        "a pulse of wind along the fjord axis, toward the mouth; then the fjord is left to ring",
        "near-inertial response: currents turning clockwise at the inertial period",
        "internal Kelvin wave, trapped within about one Rossby radius of the coast, travelling cyclonically",
        "radiation of internal-wave energy through the two straits to Bellsund and the shelf",
        "landfast ice: the wind does not reach the water; the fixed ice puts a drag on its surface",
    )
    for k_, s_ in enumerate(items):
        marker(lx + 9, ly + 24 * k_, k_ + 1)
        text(lx + 26, ly + 24 * k_ + 4, s_, 9.5)
    return B


def build() -> list[dict]:
    cmds.clear()
    c("file.new", width=W, height=H, units="Points", backgroundContents="white", name="fig0_schematic")
    bot = satellite()
    plan_view(bot + 70)
    return cmds


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cli", default=os.environ.get("VECTORCRAFT_CLI"))
    a = ap.parse_args()
    cmd_file = HERE / "fig0_schematic.commands.json"
    cmd_file.write_text(json.dumps(build(), indent=1))
    (FIG / "fig0_schematic.json").write_text(json.dumps(FACTS, indent=1) + "\n")
    if not a.cli:
        print(f"wrote {cmd_file}; pass --cli to render")
        return
    outs = [FIG / f"fig0_schematic.{e}" for e in ("png", "pdf")]
    args = [a.cli, "run", "--commands", str(cmd_file), "--ppi", "130"]
    for o in outs:
        args += ["--export", str(o)]
    r = subprocess.run(args, capture_output=True, text=True, env={**os.environ, "VECTORCRAFT_NO_PREFS": "1"})
    bad = [ln for ln in r.stdout.splitlines() if '"error"' in ln or '"warning"' in ln]
    if r.returncode or any('"error"' in ln for ln in bad):
        raise SystemExit(f"vectorcraft-cli failed (rc {r.returncode}):\n" + "\n".join(bad) + r.stderr[-2000:])
    for ln in bad:
        print("warning:", ln[:300])
    # PDF: the embedded images and textures are re-encoded as JPEG (200 ppi) with Ghostscript when it is
    # installed, as in Exp. 2 (only image streams change).
    import shutil
    if shutil.which("gs"):
        tmp = outs[1].with_suffix(".gs.pdf")
        subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7",
                        "-dDownsampleColorImages=true", "-dColorImageResolution=120", "-dColorImageDownsampleType=/Bicubic",
                        "-dAutoFilterColorImages=false", "-dPassThroughJPEGImages=false", "-dColorImageFilter=/DCTEncode", "-dJPEGQ=80",
                        f"-sOutputFile={tmp}", str(outs[1])], check=True)
        tmp.replace(outs[1])
    from PIL import Image
    jpg = FIG / "fig0_schematic.jpg"
    Image.open(outs[0]).convert("RGB").save(jpg, quality=82, optimize=True)
    outs[0].unlink()
    print("rendered", ", ".join(o.name for o in outs[1:] + [jpg]))


if __name__ == "__main__":
    main()
