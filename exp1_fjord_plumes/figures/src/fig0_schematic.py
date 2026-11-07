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
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
PANELS = HERE / "panels"
SAT = HERE / "sat"
W, H = 1100, 1230
OY = 480                                # vertical offset of the section (c) below the satellite row
# panel (a): along-fjord section of configuration (c); x in km -> points, depth in m -> points
X0, KX = 120.0, 21.0                    # x of the ice face; points per km (0-40 km -> 120-960)
Y0, KZ = 110.0 + OY, 0.54                    # y of the sea surface; points per metre of depth


def X(km):
    return X0 + KX * km


def Y(m):
    return Y0 + KZ * m


ICE, ICE_LINE = "#e8f3f9", "#6b8fa3"
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
    label(mx + 26, my + 10, "model domain of (c)", 8.5)
    label(mx + 26, my + 24, "at the same scale", 8.5)
    label(mx + 26, my + 38, "(ice face at the top)", 8.5)
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
    label(*Q(330, 450), "calving front (terminus): the ice face in (c)-(e)", 10, color=WARM)
    label(*Q(600, 330), "ice mélange: icebergs and sea ice", 10)
    sx, sy = Q(40, 530)
    line(sx, sy, sx + 5 / 0.036 * sb, sy, color=INK, weight=2.2)
    text(sx, sy - 5, "5 km", 9, color=INK)
    text(bx0, 52, "(b)", 13, style="Bold")
    text(bx0 + 26, 52, "Helheim Glacier terminus, 30 m resolution", 10.5)
    text(40, 474, "Imagery: NASA Worldview / GIBS. (a) MODIS Terra true colour, 5 Aug 2022. (b) Harmonized Landsat "
         "Sentinel-2 (HLS S30), 30 Jul 2022; contains modified Copernicus Sentinel data 2022.", 8.5, color="#555555")


def build() -> list[dict]:
    c("file.new", width=W, height=H, units="Points", backgroundContents="white", name="fig0_schematic")
    satellite_row()
    o = OY
    # ---------------- (c) along-fjord section --------------------------------------------------
    rect(70, 40 + o, 1000, Y0 - 40 - o, color="#f6f9fb")
    rect(X0, Y0, 1060 - X0, Y(165) - Y0, color=PW)
    rect(X0, Y(235), 1060 - X0, Y(500) - Y(235), color=AW)
    rect(X0, Y(165), 1060 - X0, Y(235) - Y(165),
         gradient={"kind": "linear", "angle": 90, "stops": [{"offset": 0, "color": AW}, {"offset": 1, "color": PW}]})
    rect(70, Y(500), 1000, 14, color=BEDC)
    rect(X(38), Y0, X(40) - X(38), Y(500) - Y0, color="#7f7f7f", opacity=25)
    line(X(40), Y0, X(40), Y(500), color="#444444", weight=1.0, dash=[4, 3])
    line(X0, Y0, 1060, Y0, color=COOL, weight=1.4)
    path(f"M 70 {70 + o} L 95 {62 + o} L {X0} {56 + o} L {X0} {Y(500)} L 70 {Y(500)} Z", color=ICE, stroke=ICE_LINE, weight=1.4)
    zn = Y(181)
    path(f"M {X0} {Y(500)} L {X0 + 6} {Y(500)} C {X0 + 18} {Y(380)} {X0 + 30} {Y(260)} {X0 + 44} {zn + 8} "
         f"C {X(6)} {zn + 9} {X(11)} {zn + 6} {X(15)} {zn} C {X(11)} {zn - 6} {X(6)} {zn - 9} {X0} {zn - 8} Z",
         color=PLUME, opacity=92)
    arrow(X0 + 9, Y(480), X0 + 20, Y(260), color="#ffffff", weight=2.6)
    arrow(X(15) + 4, zn, X(19), zn, color=PURP, weight=2.6)
    line(X0, zn, X(3.2), zn, color=INK, weight=0.9, dash=[3, 2])
    for d in (440, 370, 300):
        arrow(X0 + 70, Y(d), X0 + 30 - (500 - d) * 0.025, Y(d) - 3, color=WARM, weight=1.4)
    arrow(X(31), Y(390), X(9), Y(390), color=WARM, weight=2.8)
    line(X(20), Y0 - 6, X(20), Y(500), color=INK, weight=1.2, dash=[6, 3])
    for k0 in (6.0, 14.5, 23.0):
        arrow(X(k0), 84 + o, X(k0) + 130, 84 + o, color=INK, weight=2.0)
    arrow(1010, Y(200) - 34, 1010, Y(200) + 34, color=INK, weight=2.0, both=True)
    path(f"M {X(13)} {Y(200)} C {X(18)} {Y(200) - 9} {X(23)} {Y(200) + 9} {X(28)} {Y(200)} "
         f"C {X(32)} {Y(200) - 9} {X(36)} {Y(200) + 9} {X(40)} {Y(200)}", stroke="#555555", weight=1.1)
    line(1068, Y0, 1068, Y(500), color=INK, weight=0.8)
    for d in range(0, 501, 100):
        line(1068, Y(d), 1073, Y(d), color=INK, weight=0.8)
        text(1077, Y(d) + 3.5, f"{d}", 9)
    text(1094, Y(250), "depth (m)", 9, rotation=-90, justify="center")
    line(X0, Y(500) + 18, X(40), Y(500) + 18, color=INK, weight=0.8)
    for k in range(0, 41, 10):
        line(X(k), Y(500) + 18, X(k), Y(500) + 23, color=INK, weight=0.8)
        text(X(k), Y(500) + 34, f"{k}", 9, justify="center")
    text(X(20), Y(500) + 47, "distance from the ice face (km)", 9, justify="center")
    text(20, 34 + o, "(c)", 13, style="Bold")
    text(48, 34 + o, "Along-fjord section of the idealised fjord, configuration (c); vertical scale exaggerated, plume enlarged", 11)
    text(76, 52 + o, "Glacier", 11, style="Bold")
    text(X(20), 70 + o, "along-fjord wind  τ = 0.2 N/m² for 12 h", 10.5, justify="center")
    text(X(20) + 5, Y(30), "Q_ex and heat toward the front, x = 20 km", 10)
    text(X(6), Y(120), "Polar Water  −1.5 °C, 32.9 psu", 10.5)
    text(X(6), Y(470), "Atlantic Water  4.0 °C, 34.6 psu", 10.5)
    text(X(3.4), zn + 3.5, "z_n ≈ 181 m", 9.5)
    text(X(4.2), zn - 12, "neutral buoyancy: plume detaches, outflow", 10, color=PURP)
    text(X0 + 74, Y(372), "entrainment", 10, color=WARM)
    text(X(9.5), Y(410), "warm inflow at depth", 10, color=WARM)
    text(X0 + 4, Y(500) + 11, "Q0 = 100 m³/s", 9.5, color="#ffffff")
    text(X(39), Y(470), "sponge", 9, rotation=-90, justify="center")
    text(1010, Y(200) + 50, "shelf", 10, justify="center", style="Bold")
    text(1010, Y(200) + 63, "z_c(t) = 200 m", 9, justify="center")
    text(1010, Y(200) + 75, "− 30 sin(2πt / 6 d)", 9, justify="center")
    c("shape.rectangle", x=X0 - 4, y=Y0 + 4, width=60, height=Y(500) - Y0 - 6)
    fill(stroke=INK, weight=1.0)
    c("stroke.set", dash=[4, 3])
    line(X0 - 4, Y(500) - 2, 155, 446 + o, color="#666666", weight=0.8, dash=[4, 3])
    line(X0 + 56, Y(500) - 2, 560, 446 + o, color="#666666", weight=0.8, dash=[4, 3])
    # ---------------- (d)-(g) data panels ---------------------------------------------------------
    py = 470 + o
    for tag, title, name, r in (
            ("(d)", "config (a): resolved line plume", "p_plume2d.png", [40, py, 245, 255]),
            ("(e)", "config (b): two plumes on the ice face", "p_face_w.png", [300, py, 280, 255]),
            ("(f)", "initial T and S profiles", "p_profiles.png", [600, py, 230, 255]),
            ("(g)", "forcing of config (c)", "p_forcing.png", [850, py, 235, 255])):
        c("file.place", path=str(PANELS / name), rect=r, link=False)
        text(r[0] + 4, py - 8, tag, 12, style="Bold")
        text(r[0] + 30, py - 8, title, 10.5)
    text(40, H - 10, "(d), (e): model output of this experiment (w, vertical velocity; black line in (d), the 1.25 °C "
         "isotherm between PW and AW). Schematic in (c) not to scale.", 9, color="#555555")
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
    subprocess.run([sys.executable, str(HERE / "fig0_panels.py")], check=True)
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
