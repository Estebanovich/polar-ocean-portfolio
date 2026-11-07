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
W, H = 1100, 1400
SAT = HERE / "sat"
O = 500                                  # top of the process schematic (c)
SL, PC, BD = O + 110, O + 250, O + 420   # sea level, pycnocline centre, sea floor (y, points)
FX = 230                                 # x of the ice face

ICE, ICE_LINE = "#eaf4fa", "#6b8fa3"
ICE_DEEP = "#b9d7ea"
PW, AW = "#a9cfe6", "#f2b48c"
PLUME = "#b9a7d8"
INK, WARM, COOL, BEDC, PURP = "#222222", "#c0392b", "#2c6ea3", "#a39282", "#5b3f8f"
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


def label(x, y, s, size=10, color=INK, style=None, pad=3):
    """Text on a translucent white box (for labels over imagery)."""
    w = 0.53 * size * len(s)
    rect(x - pad, y - size - pad + 2, w + 2 * pad, size + 2 * pad, color="#ffffff", opacity=78)
    text(x, y, s, size, color=color, style=style)


def satellite_row():
    # (a) MODIS overview, 900 x 1100 px placed at 300 x 367 pt
    ax0, ay0, sa = 40.0, 64.0, 300.0 / 900.0
    P = lambda px, py: (ax0 + sa * px, ay0 + sa * py)
    c("file.place", path=str(SAT / "sermilik_modis_2022-08-05.jpg"), rect=[ax0, ay0, 300, 366.7], link=False)
    x1, y1 = P(193, 90)
    x2, y2 = P(536, 270)
    c("shape.rectangle", x=x1, y=y1, width=x2 - x1, height=y2 - y1)
    fill(stroke="#ffd400", weight=1.6)
    label(*P(20, 470), "Greenland Ice Sheet", 9.5)
    label(*P(150, 60), "Helheim Glacier", 9.5)
    text(*P(560, 520), "Sermilik Fjord", 10, color="#ffffff", style="Bold")
    text(*P(470, 1075), "shelf (Irminger Sea)", 9.5, color="#ffffff")
    # scale bar 20 km (0.104 km per px across, 0.111 km per px along)
    bx, by = P(660, 1000)
    line(bx, by, bx + 20 / 0.104 * sa, by, color="#ffffff", weight=2.2)
    text(bx, by - 5, "20 km", 9, color="#ffffff")
    # model domain of (c) at the same scale: 40 km x 5.6 km, ice face at the top
    mx, my = P(70, 640)
    c("shape.rectangle", x=mx, y=my, width=5.6 / 0.104 * sa, height=40 / 0.111 * sa)
    fill(color="#f2b48c", stroke=INK, weight=1.2, opacity=85)
    line(mx, my, mx + 5.6 / 0.104 * sa, my, color=COOL, weight=3)
    label(mx + 26, my + 10, "fjord model domain,", 8.5)
    label(mx + 26, my + 24, "configuration (c),", 8.5)
    label(mx + 26, my + 38, "same scale, ice at top", 8.5)
    text(ax0, 52, "(a)", 13, style="Bold")
    text(ax0 + 26, 52, "Sermilik Fjord and Helheim Glacier, SE Greenland", 10.5)
    # (b) HLS zoom on the terminus, 1000 x 560 px placed at 700 x 392 pt
    bx0, by0, sb = 370.0, 64.0, 0.7
    Q = lambda px, py: (bx0 + sb * px, by0 + sb * py)
    c("file.place", path=str(SAT / "helheim_hls_2022-07-30.jpg"), rect=[bx0, by0, 700, 392], link=False)
    c("shape.rectangle", x=bx0, y=by0, width=700, height=392)
    fill(stroke="#ffd400", weight=1.6)
    arrow(*Q(120, 255), *Q(370, 300), color=COOL, weight=2.6)
    label(*Q(70, 235), "glacier ice flows to the sea", 10, color=COOL)
    line(*Q(415, 215), *Q(400, 395), color=WARM, weight=2.2, dash=[5, 3])
    label(*Q(330, 450), "calving front (terminus): the ice face in (c) and (d)", 10, color=WARM)
    label(*Q(600, 330), "ice mélange: icebergs and sea ice", 10)
    sx, sy = Q(40, 530)
    line(sx, sy, sx + 5 / 0.036 * sb, sy, color=INK, weight=2.2)
    text(sx, sy - 5, "5 km", 9, color=INK)
    text(bx0, 52, "(b)", 13, style="Bold")
    text(bx0 + 26, 52, "Helheim Glacier terminus, 30 m resolution", 10.5)
    text(40, 474, "Imagery: NASA Worldview / GIBS. (a) MODIS Terra true colour, 5 Aug 2022. (b) Harmonized Landsat "
         "Sentinel-2 (HLS S30), 30 Jul 2022; contains modified Copernicus Sentinel data 2022.", 8.5, color="#555555")


def fx(effect, **params):
    c("effect.apply", effect=effect, params=params)


def marker(x, y, n):
    c("shape.ellipse", x=x - 9, y=y - 9, width=18, height=18)
    fill(color="#1d1d1d", stroke="#ffffff", weight=1.2)
    text(x, y + 4.2, str(n), 10.5, color="#ffffff", style="Bold", justify="center")


def curl(x, y, r, color="#ffffff", weight=1.3, cw=True):
    """A small eddy: three-quarter circle with an arrow head."""
    s_ = 1 if cw else -1
    c("path.create", d=f"M {x + r} {y} C {x + r} {y + s_ * r * 0.55} {x + r * 0.55} {y + s_ * r} {x} {y + s_ * r} "
                       f"C {x - r * 0.55} {y + s_ * r} {x - r} {y + s_ * r * 0.55} {x - r} {y} "
                       f"C {x - r} {y - s_ * r * 0.55} {x - r * 0.55} {y - s_ * r} {x} {y - s_ * r}")
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, endArrow="Triangle", cap="round")


def curve_arrow(d, color, weight):
    c("path.create", d=d)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, endArrow="Triangle", cap="round", join="round")


def process_section():
    # sky and distant fjord walls
    rect(40, O, 1030, SL - O, gradient={"kind": "linear", "angle": 90,
                                         "stops": [{"offset": 0, "color": "#ffffff"}, {"offset": 1, "color": "#dfeaf3"}]})
    path(f"M {FX + 60} {SL} L {FX + 110} {O + 62} L {FX + 170} {O + 80} L {FX + 250} {O + 40} L {FX + 330} {O + 72} "
         f"L {FX + 420} {O + 50} L {FX + 520} {O + 84} L {FX + 610} {O + 58} L {FX + 700} {O + 90} L 1070 {O + 76} "
         f"L 1070 {SL} Z", gradient={"kind": "linear", "angle": 90, "stops": [{"offset": 0, "color": "#b9b4ac"},
                                                                          {"offset": 1, "color": "#e3dfd8"}]}, opacity=55)
    fx("distort.roughen", size=4, relative=False, detail=12, points="corner", seed=3)
    # water: Atlantic Water below, Polar Water above, a pycnocline between; fine grain
    rect(FX, SL, 1070 - FX, BD - SL, gradient={"kind": "linear", "angle": 90, "stops": [
        {"offset": 0, "color": "#e98f62"}, {"offset": 0.38, "color": "#f3b48e"}, {"offset": 0.52, "color": "#cdd6dc"},
        {"offset": 0.66, "color": "#9ec9e2"}, {"offset": 1, "color": "#c5e1f0"}]})
    fx("texture.grain", intensity=12, contrast=50, grainType="soft")
    # sea-surface ripples
    c("path.create", d=f"M {FX} {SL} L 1070 {SL}")
    c("paint.setFill", none=True)
    c("paint.setStroke", color="#3b7fb3")
    c("stroke.set", weight=1.3)
    fx("distort.zigZag", size=1.6, relative=False, ridges=40, points="smooth")
    # sea floor: rock with a rough top and a sandstone texture
    path(f"M 40 {BD} L 1070 {BD + 4} L 1070 {O + 470} L 40 {O + 470} Z", color="#8f7f6d")
    fx("texture.texturizer", texture="sandstone", scaling=120, relief=4, lightDirection="top")
    fx("distort.roughen", size=2.5, relative=False, detail=18, points="corner", seed=7)
    # glacier body, crevasses, supraglacial lake, moulin and subglacial channel
    path(f"M 40 {O + 26} C 110 {O + 30} 170 {O + 44} {FX} {O + 66} L {FX} {BD} L 40 {BD} Z",
         gradient={"kind": "linear", "angle": 90, "stops": [{"offset": 0, "color": ICE_DEEP}, {"offset": 1, "color": ICE}]},
         stroke=ICE_LINE, weight=1.4)
    fx("texture.texturizer", texture="canvas", scaling=80, relief=2, lightDirection="topLeft")
    for k in range(6):                                   # crevasses near the front
        x = 120 + k * 17
        line(x, O + 34 + k * 5.5, x + 6, O + 52 + k * 6, color="#7e9fb4", weight=1.1)
    for k in range(9):                                   # layering in the ice cliff
        y = O + 90 + k * 34
        line(FX - 50, y, FX - 2, y + 4, color="#c3dbea", weight=1.0)
    c("shape.ellipse", x=72, y=O + 24, width=46, height=8)
    fill(color="#4a9fd0")
    line(96, O + 32, 96, BD - 9, color=COOL, weight=1.8, dash=[5, 3])
    curve_arrow(f"M 96 {BD - 9} L {FX - 4} {BD - 9}", COOL, 3.2)
    # icebergs in front of the glacier (above and below the waterline)
    for (x0, w, h_up, h_dn) in ((FX + 22, 34, 12, 34), (FX + 78, 22, 8, 22), (FX + 122, 46, 15, 44), (FX + 196, 18, 6, 18)):
        path(f"M {x0} {SL} L {x0 + w * 0.2} {SL - h_up} L {x0 + w * 0.7} {SL - h_up * 0.8} L {x0 + w} {SL} "
             f"L {x0 + w * 0.85} {SL + h_dn} L {x0 + w * 0.25} {SL + h_dn * 0.9} Z", color="#f4fafd", stroke="#8fb3c9",
             weight=0.8, opacity=92)
        fx("distort.roughen", size=1.5, relative=False, detail=20, points="corner", seed=int(x0))
    # meltwater boundary layer on the ice face, and heat flux into the ice
    line(FX + 2, BD - 4, FX + 2, PC - 60, color="#e7f4fb", weight=2.6)
    for y in (BD - 50, BD - 110, BD - 170):
        curve_arrow(f"M {FX + 22} {y} C {FX + 16} {y - 4} {FX + 12} {y + 4} {FX + 5} {y}", WARM, 1.6)
    # the plume: turbulent half cone, overshoot, outflow tongue at neutral buoyancy
    path(f"M {FX} {BD} L {FX + 10} {BD} C {FX + 30} {BD - 80} {FX + 52} {PC + 60} {FX + 70} {PC + 6} "
         f"C {FX + 140} {PC + 16} {FX + 260} {PC + 12} {FX + 390} {PC - 4} C {FX + 260} {PC - 24} {FX + 140} {PC - 30} "
         f"{FX + 64} {PC - 30} C {FX + 50} {PC - 56} {FX + 24} {PC - 64} {FX} {PC - 52} Z",
         gradient={"kind": "linear", "angle": 90, "stops": [{"offset": 0, "color": "#7a63b8"}, {"offset": 1, "color": "#d5c9ef"}]},
         opacity=88)
    fx("distort.roughen", size=3.2, relative=False, detail=26, points="smooth", seed=11)
    for (x, y, r) in ((FX + 28, BD - 60, 7), (FX + 40, BD - 130, 9), (FX + 52, PC + 70, 10), (FX + 30, PC + 20, 8)):
        curl(x, y, r)
    curve_arrow(f"M {FX + 14} {BD - 14} C {FX + 18} {BD - 110} {FX + 24} {PC + 60} {FX + 30} {PC - 30}", "#ffffff", 3.0)
    curve_arrow(f"M {FX + 330} {PC - 4} L {FX + 520} {PC - 4}", PURP, 3.4)
    # entrainment of ambient water, and the warm inflow it draws from the shelf
    for (y, dx_) in ((BD - 40, 70), (BD - 100, 74), (BD - 160, 80)):
        curve_arrow(f"M {FX + dx_ + 40} {y + 10} C {FX + dx_ + 20} {y + 8} {FX + dx_} {y + 2} {FX + dx_ - 22} {y - 6}",
                    "#b8322a", 1.8)
    curve_arrow(f"M 1000 {BD - 52} C 860 {BD - 44} 700 {BD - 44} {FX + 190} {BD - 70}", "#b8322a", 5)
    # wind over the fjord and the surface water it pushes out
    for x0 in (FX + 240, FX + 420, FX + 600):
        curve_arrow(f"M {x0} {SL - 30} C {x0 + 50} {SL - 36} {x0 + 100} {SL - 26} {x0 + 140} {SL - 30}", INK, 2.4)
    curve_arrow(f"M {FX + 560} {SL + 20} L {FX + 700} {SL + 20}", "#2c6ea3", 2.0)
    # shelf: heaving pycnocline at the mouth, waves propagating into the fjord
    line(990, SL, 990, BD, color="#555555", weight=0.9, dash=[3, 3])
    arrow(1030, PC - 40, 1030, PC + 40, color=INK, weight=2.4, both=True)
    path(f"M {FX + 420} {PC} C {FX + 470} {PC - 12} {FX + 520} {PC + 12} {FX + 570} {PC} "
         f"C {FX + 620} {PC - 12} {FX + 670} {PC + 12} {FX + 720} {PC} C {FX + 740} {PC - 12} 960 {PC + 12} 990 {PC}",
         stroke="#4d4d4d", weight=1.4)
    curve_arrow(f"M 960 {PC + 26} L {FX + 700} {PC + 26}", "#4d4d4d", 1.6)
    # numbered processes
    for (x, y, n) in ((150, BD - 20, 1), (FX + 54, BD - 70, 2), (FX - 20, BD - 140, 3), (FX + 120, PC - 44, 4),
                      (760, BD - 74, 5), (FX + 420, SL - 52, 6), (1030, PC - 58, 7)):
        marker(x, y, n)
    # labels
    text(50, O + 18, "Glacier", 12, style="Bold")
    text(FX + 440, O + 18, "Fjord", 12, style="Bold", justify="center")
    text(1030, O + 18, "Shelf", 12, style="Bold", justify="center")
    text(FX + 560, SL + 46, "Polar Water: cold, fresh", 10.5, color="#1f4e79")
    text(FX + 560, BD - 14, "Atlantic Water: warm, salty", 10.5, color="#7a2a14")
    text(FX + 600, PC - 14, "pycnocline", 9.5, color="#4d4d4d")
    text(52, O + 470 + 18, "(c) is a schematic: not to scale; the plume and icebergs are enlarged.", 9, color="#555555")
    text(20, O - 8, "(c)", 13, style="Bold")
    text(46, O - 8, "Processes at a tidewater glacier and in its fjord", 11)
    # legend
    lx, ly = 40, O + 505
    items = ("subglacial discharge: surface meltwater drains to the bed and leaves at the grounding line",
             "buoyant plume: rises along the ice and entrains ambient water (turbulent eddies)",
             "submarine melting: the plume carries heat to the ice face",
             "neutral buoyancy: the plume detaches and flows out toward the sea",
             "warm Atlantic Water flows toward the glacier at depth, replacing the entrained water",
             "along-fjord wind drives an exchange with the shelf",
             "shelf forcing: the pycnocline heaves over days and the motion enters the fjord")
    for k, s_ in enumerate(items):
        col, row = k % 2, k // 2
        x, y = lx + col * 520, ly + row * 24
        marker(x + 9, y, k + 1)
        text(x + 26, y + 4, s_, 10)


def plume_geometry():
    """(d): line plume versus two point plumes on the ice face, drawn in perspective."""
    y0 = O + 650
    text(20, y0 - 8, "(d)", 13, style="Bold")
    text(46, y0 - 8, "Plume geometry at the ice face: configurations (a) and (b)", 11)
    for k, (x0, title) in enumerate(((60, "configuration (a): a line plume, discharge along the whole face"),
                                      (560, "configuration (b): two point plumes, x0 apart"))):
        top, bot, wid, sk = y0 + 40, y0 + 190, 380, 60
        # the ice face as a parallelogram seen obliquely, with a stippled texture
        path(f"M {x0} {top + 20} L {x0 + wid} {top + 20} L {x0 + wid + sk} {top - 10} L {x0 + sk} {top - 10} Z",
             color="#f6fbfe", stroke=ICE_LINE, weight=1.0)
        path(f"M {x0} {top + 20} L {x0 + wid} {top + 20} L {x0 + wid} {bot} L {x0} {bot} Z",
             gradient={"kind": "linear", "angle": 90, "stops": [{"offset": 0, "color": ICE_DEEP}, {"offset": 1, "color": ICE}]},
             stroke=ICE_LINE, weight=1.2)
        fx("texture.texturizer", texture="canvas", scaling=80, relief=2, lightDirection="topLeft")
        rect(x0 - 10, bot, wid + 20, 12, color="#8f7f6d")
        fx("texture.texturizer", texture="sandstone", scaling=100, relief=3, lightDirection="top")
        if k == 0:
            path(f"M {x0 + 10} {bot} L {x0 + wid - 10} {bot} L {x0 + wid - 10} {top + 60} L {x0 + 10} {top + 60} Z",
                 gradient={"kind": "linear", "angle": 90, "stops": [{"offset": 0, "color": "#7a63b8"},
                                                                    {"offset": 1, "color": "#d5c9ef"}]}, opacity=80)
            fx("distort.roughen", size=2.5, relative=False, detail=20, points="smooth", seed=5)
            for xx in range(int(x0) + 50, int(x0 + wid) - 20, 70):
                curve_arrow(f"M {xx} {bot - 8} L {xx} {top + 72}", "#ffffff", 2.2)
            text(x0 + wid / 2, bot + 28, "q0 per metre of ice face (2-D)", 10, justify="center")
        else:
            for xs in (x0 + wid / 2 - 45, x0 + wid / 2 + 45):
                path(f"M {xs - 4} {bot} L {xs + 4} {bot} L {xs + 60} {top + 60} L {xs - 60} {top + 60} Z",
                     gradient={"kind": "linear", "angle": 90, "stops": [{"offset": 0, "color": "#7a63b8"},
                                                                        {"offset": 1, "color": "#d5c9ef"}]}, opacity=72)
                fx("distort.roughen", size=2.5, relative=False, detail=20, points="smooth", seed=int(xs))
                curve_arrow(f"M {xs} {bot - 8} L {xs} {top + 74}", "#ffffff", 2.2)
                c("shape.ellipse", x=xs - 5, y=bot - 3, width=10, height=6)
                fill(color="#1f5f99")
            arrow(x0 + wid / 2 - 45, bot + 20, x0 + wid / 2 + 45, bot + 20, color=INK, weight=1.4, both=True)
            text(x0 + wid / 2, bot + 36, "x0", 11, justify="center", style="Bold")
            text(x0 + wid + 8, top + 100, "the cones widen", 9.5)
            text(x0 + wid + 8, top + 113, "as they rise and", 9.5)
            text(x0 + wid + 8, top + 126, "merge if x0 < x0c", 9.5)
            text(x0 + wid + 8, top + 139, "(Cenedese &", 9.5)
            text(x0 + wid + 8, top + 152, "Gatto 2016a)", 9.5)
        text(x0, top - 22, title, 10.5)
        text(x0 + wid + sk + 6, top - 4, "ice face", 9.5, color=ICE_LINE)


def build() -> list[dict]:
    c("file.new", width=W, height=H, units="Points", backgroundContents="white", name="fig0_schematic")
    satellite_row()
    process_section()
    plume_geometry()
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
