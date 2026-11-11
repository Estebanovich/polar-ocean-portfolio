"""Figure 0 of Experiment 3: the emulator protocol (training set, held-out sets, 8 methods), drawn with VectorCraft.

Writes the VectorCraft command list (``fig0_protocol.commands.json``) and, if the VectorCraft CLI is available,
renders ``figures/fig0_protocol.png`` and ``.pdf`` with its headless batch mode, in the style and with the drawing
primitives of exp2_isomip_cavity/figures/src/fig0_schematic.py:

    vectorcraft-cli run --commands fig0_protocol.commands.json --export fig0_protocol.png ...

Tool: VectorCraft (https://github.com/storytold/vectorcraft), CLI release v0.22.0, macOS universal.
Usage: fig0_protocol.py [--cli PATH_TO_vectorcraft-cli]   (or env VECTORCRAFT_CLI)

Content (no results anywhere): (1) the 29 Exp. 2 runs in the (Tbot, z_c) plane, as designed (the 15-run factorial,
the 10 Latin-hypercube geometry runs with s = 0.8, 1.2, the 3 final-test runs at Tbot = 2.0 and the Ocean0
reference); (2) the held-out sets of pd0_preregistro.md; (3) the 8 methods with their sources; (4) the metrics and the
pre-registered hypotheses H3a-H3c. Coordinates are points, y down.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
W, H = 1100, 755
FONT = "Helvetica"

INK, WARM, COOL, PURP, GREY = "#222222", "#b8322a", "#2c6ea3", "#5b3f8f", "#555555"
PAPER, PANEL_LINE = "#f7f4ee", "#9a8f80"
TRAIN, GEOM, TEST, REF = "#2c6ea3", "#5aa0d8", "#b8322a", "#222222"
GROUPS = {"baseline": "#d9d9d9", "physics": "#cfe3f2", "stats": "#d6ecd9", "nn": "#f6dcc8"}

# The design of the runs (exp3 dataset; the case names encode Tbot, z_c and s): inputs, not results.
FACTORIAL = [(t, z) for z in (250, 400, 550) for t in (-0.5, 0.0, 0.5, 1.0, 1.5)]
GEOMETRY = [(0.75, 527, 0.8), (1.43, 295, 0.8), (-0.05, 328, 0.8), (-0.42, 402, 0.8), (0.34, 444, 0.8),
            (1.42, 481, 1.2), (0.02, 291, 1.2), (-0.39, 503, 1.2), (0.97, 386, 1.2), (0.66, 311, 1.2)]
FINAL = [(2.0, 250), (2.0, 400), (2.0, 550)]
OCEAN0 = (1.0, 360)          # ISOMIP+ Ocean0 linear profile (d = 720 m): profile test only

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
    return {"kind": "linear", "angle": angle, "stops": [{"offset": o, "color": col} for o, col in stops]}


def rect(x, y, w, h, radius=None, **kw):
    c("shape.rectangle", x=x, y=y, width=w, height=h, **({"radius": radius} if radius else {}))
    fill(**kw)


def line(x1, y1, x2, y2, color=INK, weight=1.0, dash=None, opacity=None):
    c("shape.line", x1=x1, y1=y1, x2=x2, y2=y2)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, cap="round", **({"dash": dash} if dash else {}))
    if opacity is not None:
        c("object.setProps", opacity=opacity)


def arrow(x1, y1, x2, y2, color=INK, weight=2.0):
    c("shape.line", x1=x1, y1=y1, x2=x2, y2=y2)
    c("paint.setFill", none=True)
    c("paint.setStroke", color=color)
    c("stroke.set", weight=weight, endArrow="Triangle", cap="round")


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


def panel(x, y, w, h, n, title):
    rect(x + 3, y + 4, w, h, radius=8, color="#000000", opacity=10)               # soft shadow
    rect(x, y, w, h, radius=8, color=PAPER, stroke=PANEL_LINE, weight=1.0)
    fx("texture.texturizer", texture="canvas", scaling=70, relief=2, lightDirection="topLeft")
    marker(x + 18, y + 20, n)
    text(x + 36, y + 25, title, 12, style="Bold")


# ----------------------------------------------------------------------------------- (1) the runs
def design_panel(x0=30, y0=60, w=360, h=430):
    panel(x0, y0, w, h, 1, "The runs: 29 MITgcm simulations (Exp. 2)")
    px, py, pw, ph = x0 + 62, y0 + 62, w - 92, h - 200           # plot box
    tmin, tmax, zmin, zmax = -0.8, 2.3, 220, 580
    X = lambda t: px + (t - tmin) / (tmax - tmin) * pw
    Y = lambda z: py + (z - zmin) / (zmax - zmin) * ph          # z_c downward, as a depth
    rect(px, py, pw, ph, color="#ffffff", stroke="#b5ab9c", weight=0.8)
    # warm_out band (Tbot >= 1.4) and the zc_out cuts (325, 475 m)
    rect(X(1.36), py, X(1.6) - X(1.36), ph, color="#f6d7cf", opacity=70)
    for zc in (325, 475):
        line(px, Y(zc), px + pw, Y(zc), color="#7d6f6a", weight=0.8, dash=[3, 3])
    text(px + pw - 4, Y(325) - 4, "zc_out cut 325 m", 7.5, color=GREY, justify="right")
    text(px + pw - 4, Y(475) - 4, "zc_out cut 475 m", 7.5, color=GREY, justify="right")
    text(X(1.48), py + 12, "warm_out", 7.5, color=WARM, justify="center")
    # axes
    for t in (-0.5, 0.0, 0.5, 1.0, 1.5, 2.0):
        line(X(t), py + ph, X(t), py + ph + 4, color=INK, weight=0.8)
        text(X(t), py + ph + 15, f"{t:g}", 8.5, justify="center")
    for z in (250, 400, 550):
        line(px - 4, Y(z), px, Y(z), color=INK, weight=0.8)
        text(px - 7, Y(z) + 3, f"{z}", 8.5, justify="right")
    text(px + pw / 2, py + ph + 31, "T_{bot}: CDW temperature below the thermocline (°C)", 9, justify="center")
    text(px - 40, py + ph / 2, "z_{c}: thermocline depth (m)", 9, justify="center", rotation=-90)
    # the runs as designed
    for t, z in FACTORIAL:
        dot(X(t), Y(z), 4.6, TRAIN)
    for t, z, s in GEOMETRY:
        dot(X(t), Y(z), 4.4, GEOM, shape="triangle")
    for t, z in FINAL:
        dot(X(t), Y(z), 4.6, TEST, shape="star")
    dot(X(OCEAN0[0]), Y(OCEAN0[1]), 4.4, REF, shape="diamond")
    # legend
    ly = py + ph + 52
    items = [("circle", TRAIN, "15 factorial runs, s = 1 (training)"),
             ("triangle", GEOM, "10 geometry runs, s = 0.8, 1.2 (training)"),
             ("star", TEST, "3 final-test runs, T_{bot} = 2.0 (held out)"),
             ("diamond", REF, "Ocean0 reference, linear profile (profile test)")]
    for k, (shp, col, lab) in enumerate(items):
        dot(x0 + 30, ly + k * 16, 4.4, col, shape=shp)
        text(x0 + 44, ly + k * 16 + 3.5, lab, 8.8)
    text(x0 + 18, y0 + h - 10, "s scales the ice draft; 4 km ISOMIP+ cavity; melt = mean of months 19–24", 7.8,
         color=GREY)


# ----------------------------------------------------------------------------------- (2) held-out sets
def folds_panel(x0=405, y0=60, w=300, h=430):
    panel(x0, y0, w, h, 2, "What each method is tested on")
    text(x0 + 18, y0 + 50, "Each fit sees only its training runs; standardisation,", 8.6, color=GREY)
    text(x0 + 18, y0 + 62, "EOFs, kernels and K, n are refitted inside each fold.", 8.6, color=GREY)
    rows = [
        ("loro", "25 folds: train 24 / test 1", "in range: leave one run out"),
        ("zc_out", "3 folds: thirds of z_{c}", "train 17/16/17, test 8/9/8"),
        ("warm_out", "1 fold: train T_{bot} ≤ 1.0 (20)", "test T_{bot} ≥ 1.4 (5)"),
        ("profile_test", "train 25 / test Ocean0", "different profile; per-cell methods"),
        ("final", "train 25 / test T_{bot} = 2.0 (3)", "evaluated once, frozen config"),
    ]
    y = y0 + 82
    for k, (name, a, b) in enumerate(rows):
        col = "#f6d7cf" if name == "final" else "#ffffff"
        rect(x0 + 16, y, w - 32, 50, radius=5, color=col, stroke="#b5ab9c", weight=0.8)
        # train / test bar
        bx, bw = x0 + 26, 70
        frac = {"loro": 24 / 25, "zc_out": 17 / 25, "warm_out": 20 / 25, "profile_test": 25 / 26,
                "final": 25 / 28}[name]
        rect(bx, y + 31, bw * frac, 9, color=TRAIN)
        rect(bx + bw * frac, y + 31, bw * (1 - frac), 9, color=TEST)
        text(x0 + 26, y + 18, name, 10, style="Bold", color=WARM if name == "final" else INK)
        text(x0 + 108, y + 19, a, 8.6)
        text(x0 + 108, y + 36, b, 8.6, color=GREY)
        y += 58
    rect(x0 + 26, y + 8, 10, 8, color=TRAIN)
    text(x0 + 40, y + 16, "training runs", 8.2)
    rect(x0 + 120, y + 8, 10, 8, color=TEST)
    text(x0 + 134, y + 16, "held-out runs", 8.2)
    text(x0 + 18, y + 36, "Sensitivity repeat: all folds but final without the", 8.2, color=GREY)
    text(x0 + 18, y + 48, "5 runs failing N1 (never used for the hypotheses).", 8.2, color=GREY)


# ----------------------------------------------------------------------------------- (3) methods
def methods_panel(x0=720, y0=60, w=350, h=430):
    panel(x0, y0, w, h, 3, "Eight ways to predict the melt map")
    groups = [
        ("baseline", "Baseline", [("M0", "mean of the training maps", "")]),
        ("physics", "Physics: one fitted constant", [
            ("M1", "γ c TF", "Favier et al. 2019, Eq. 2"),
            ("M2", "K c^{2} |TF| TF", "Burgard et al. 2022, Eq. 14"),
            ("M3", "K c^{2} TF |⟨TF⟩|", "Jourdain et al. 2020, Eq. 1")]),
        ("stats", "Statistics on (T_{bot}, z_{c}, s)", [
            ("M4", "per-cell regression, degree 1 or 2", ""),
            ("M5", "Gaussian process on EOFs", "Myren & Lawrence 2021")]),
        ("nn", "Neural networks, 5 members", [
            ("M6", "deep ensemble on 10 cell features", "Lakshminarayanan et al. 2017"),
            ("M7", "hybrid: K TF_{+}^{n} g_{θ}(geometry)", "g_{θ} sees no temperature")]),
    ]
    y = y0 + 46
    for key, title, items in groups:
        hgt = 22 + 30 * len(items)
        rect(x0 + 14, y, w - 28, hgt, radius=5, color=GROUPS[key], stroke="#b5ab9c", weight=0.7)
        text(x0 + 24, y + 15, title, 9.4, style="Bold")
        yy = y + 32
        for m, form, src in items:
            text(x0 + 26, yy + 3, m, 10, style="Bold", color=WARM if m == "M7" else INK)
            text(x0 + 58, yy + 3, form, 9.4)
            if src:
                text(x0 + w - 24, yy + 15, src, 7.8, color=GREY, justify="right")
            yy += 30
        y += hgt + 8
    text(x0 + 18, y0 + h - 10, "c = ρ_{sw}c_{p}/(ρ_{i}L_{i}); TF: thermal forcing at the ice base", 7.8, color=GREY)


# ----------------------------------------------------------------------------------- (4) metrics and hypotheses
def hypotheses_panel(x0=30, y0=510, w=1040, h=200):
    panel(x0, y0, w, h, 4, "How the predictions are scored (fixed before any fit on real data)")
    cx = x0 + 22
    text(cx, y0 + 54, "Metrics, per held-out run (Burgard et al. 2022):", 9.6, style="Bold")
    for k, s in enumerate(("RMSE of the melt map over the run's own ice shelf (Eq. 32)",
                           "relative error of the integrated melt (Eq. 33–34)",
                           "error of the mean melt near the grounding line (Eq. 35)",
                           "M6, M7: spread of the 5 network members; M5–M7: coverage, CRPS",
                           "noise floor: no difference below 0.157 m a^{-1} or 1.3 % is read")):
        text(cx + 10, y0 + 74 + 17 * k, "• " + s, 9)
    hx = x0 + 450
    text(hx, y0 + 54, "Pre-registered hypotheses:", 9.6, style="Bold")
    hyps = [("H3a", "On the 3 final-test runs, the hybrid M7 has a smaller integrated-melt error than M6,",
             "on average and in at least 2 of the 3 runs (a difference below 1.3 % is inconclusive)."),
            ("H3b", "Inside the training range (loro), the hybrid costs little: RMSE of M7 ≤ 1.10 × that of M6.", ""),
            ("H3c", "Out of range (zc_out, warm_out, final: 33 cases), the ensemble is never confidently wrong:",
             "no case with spread < 5 % and integrated error > 20 % (Ockenden et al. 2026 replication).")]
    y = y0 + 74
    for tag, a, b in hyps:
        rect(hx, y - 11, 34, 15, radius=3, color=PURP)
        text(hx + 17, y, tag, 8.6, color="#ffffff", style="Bold", justify="center")
        text(hx + 42, y, a, 9)
        if b:
            text(hx + 42, y + 14, b, 9)
        y += 34
    text(hx, y0 + h - 14, "The ranking of all 8 methods is reported as found, whichever wins.", 9, color=GREY)


def build() -> list[dict]:
    c("file.new", width=W, height=H, units="Points", backgroundContents="white", name="fig0_protocol")
    text(30, 36, "Experiment 3: can a melt emulator trained on 25 cavity runs extrapolate?", 15, style="Bold")
    design_panel()
    folds_panel()
    methods_panel()
    hypotheses_panel()
    for x in (393, 708):                                       # reading order (1) -> (2) -> (3)
        arrow(x - 1, 275, x + 10, 275, color=PANEL_LINE, weight=2.0)
    text(30, H - 16, "A schematic of the protocol: no model output or result is drawn in this figure.", 8.5,
         color=GREY)
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
