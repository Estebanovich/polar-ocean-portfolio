#!/opt/anaconda3/envs/pyocean/bin/python
r"""Render the DESCRIPTIVE Blender animation of Exp. 4 (Task 10) from a blender_export.py scene.

Two modes in one file:

1. 3-D panels, inside Blender (headless; numpy only, no other package)::

       blender -b -P blender_render.py -- render --scene DIR --frames K [K ...] | --stills
                                                 [--res W H] [--samples N]

   writes DIR/panels/<run>_<KKKK>.png (RGBA, transparent background), one per run of the shot.

2. Overlay and composition, in the pyocean env (matplotlib)::

       blender_render.py compose --scene DIR --frames K ... | --stills  --out-dir OUT [--name NAME]

   writes OUT/<name>_<HH>h.png (stills) or OUT/frames/<KKKK>.png: the panel(s) on white with the
   clock (hours and inertial periods), the run labels, the fixed colour bar, the vertical
   exaggeration, the isopycnal, the "descriptive" label, what the colour and the geometry show, and
   the wind-pulse indicator (peak |tau| of the forcing files, from meta.json).

   blender_render.py movie --scene DIR --out-dir OUT --name NAME  (ffmpeg H.264 MP4, Pillow GIF)

Scene (fixed for every frame and both shots): coordinates in km from the SW corner of the crop,
heights in km x VE; one oblique perspective camera from the southwest over Bellsund looking
up-fjord (ENE); one sun from the northwest (azimuth 315 deg, elevation 45 deg, the cartographic
hill-shade convention) and a uniform grey world light; Standard view transform (no tone curve, so
the colour-table colours are kept); no bloom, glare, depth of field or motion blur. Terrain: one
mesh on the cell centres, land = IBCAO elevation in a neutral warm grey, seafloor = - Depth in a
muted blue-grey ramp by depth. Isopycnal: mesh on the cell centres where the isopycnal exists
(quads with four valid corners), vertex colour = the filtered displacement on the fixed scale
(cmocean balance), alpha ISO_ALPHA. Ice (F2): flat slab over the exact AREA.bin >= 0.5 cells (cell
edges, not centres), thickness h_m x VE, white, alpha ICE_ALPHA.

Shot B (F3 | F2) is rendered as two separate scenes with the identical camera, lights and colour
scale, composed side by side by ``compose``: the two runs share the grid, so one scene with both
copies side by side would need a different (wider, more distant) camera than shot A, and two
viewports in one Blender frame need the compositor; two renders guarantee identical framing.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from exp4_arctic_bayiw.analysis.blender_export import lut_colour, quad_faces  # noqa: E402  numpy only

ISO_ALPHA = 0.85
ICE_ALPHA = 0.6
ICE_EMISSION = 0.35                                        # keeps the slab white under the sun
LAND_RGB = ((0.70, 0.68, 0.64), (0.86, 0.85, 0.83))       # coast .. highest land, neutral warm grey
SEA_RGB = ((0.62, 0.72, 0.79), (0.26, 0.35, 0.45))        # shallow .. deepest, muted blue-grey
NAN_RGB = (0.70, 0.70, 0.70)
SPONGE_RGB = (0.55, 0.55, 0.55)                            # masked open-boundary sponge                               # isopycnal without a filtered colour
CAM_LOC_KM = (-30.0, -38.0, 55.0)                          # southwest of Bellsund, above it
CAM_TARGET_KM = (50.0, 27.0, -3.0)                         # mid-fjord
CAM_LENS_MM = 30.0
RES = (1920, 1000)
CROP_TOP = 0.14                                            # fraction of the render above the domain


# ---- shared helpers (numpy) ---------------------------------------------------------------------

def load_scene(d: Path) -> tuple:
    meta = json.loads((Path(d) / "meta.json").read_text())
    ter = dict(np.load(Path(d) / "terrain.npz"))
    return meta, ter


def km_coords(ter: dict, meta: dict | None = None) -> tuple:
    """(X, Y) [km] of the terrain cell centres from the SW corner of the DATA crop (the frames), so
    the camera frames the same box whether or not the terrain extends over the sponge."""
    off = (meta or {}).get("view_offset", {"j": 0, "i": 0})
    x = (ter["xc_m"] - ter["xc_m"][off["i"]]) / 1e3
    y = (ter["yc_m"] - ter["yc_m"][off["j"]]) / 1e3
    return np.meshgrid(x, y)


def data_slice(meta: dict, shape) -> tuple:
    """Slices of the terrain grid that hold the frames (the data crop)."""
    off = meta.get("view_offset", {"j": 0, "i": 0})
    c = meta["crop"]
    return (slice(off["j"], off["j"] + c["j1"] - c["j0"]), slice(off["i"], off["i"] + c["i1"] - c["i0"]))


def terrain_height_rgb(ter: dict, ve: float, ve_land: float) -> tuple:
    """(height [km], rgb) of the terrain mesh on the cell centres: land elevation x ve_land,
    seafloor - Depth x ve."""
    land = ter["land"].astype(bool)
    z = np.where(land, ter["elev_m"] * ve_land, -ter["depth_m"] * ve) / 1e3
    em = max(float(ter["elev_m"].max()), 1.0)
    dm = max(float(ter["depth_m"].max()), 1.0)
    a = np.clip(ter["elev_m"] / em, 0, 1)[..., None]
    b = np.clip(ter["depth_m"] / dm, 0, 1)[..., None] ** 0.6
    lrgb = (1 - a) * np.array(LAND_RGB[0]) + a * np.array(LAND_RGB[1])
    srgb = (1 - b) * np.array(SEA_RGB[0]) + b * np.array(SEA_RGB[1])
    return z, np.where(land[..., None], lrgb, srgb)


def stills(meta: dict) -> list:
    return [s["index"] for s in meta["extra"]["stills"]]


# ---- Blender ------------------------------------------------------------------------------------

def _bl_mesh(bpy, name, verts, faces, rgb=None, smooth=True):
    me = bpy.data.meshes.new(name)
    me.from_pydata(verts.tolist(), [], faces.tolist())
    me.update()
    if rgb is not None:
        at = me.color_attributes.new("col", "FLOAT_COLOR", "POINT")
        rgba = np.concatenate([rgb.reshape(-1, 3), np.ones((rgb.reshape(-1, 3).shape[0], 1))], axis=1)
        at.data.foreach_set("color", rgba.astype(np.float32).ravel())
    if smooth:
        me.polygons.foreach_set("use_smooth", np.ones(len(me.polygons), bool))
    ob = bpy.data.objects.new(name, me)
    bpy.context.scene.collection.objects.link(ob)
    return ob


def _bl_material(bpy, name, alpha=1.0, rgb=None, roughness=0.8):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    p = nt.nodes["Principled BSDF"]
    if rgb is None:
        a = nt.nodes.new("ShaderNodeAttribute")
        a.attribute_name = "col"
        nt.links.new(a.outputs["Color"], p.inputs["Base Color"])
    else:
        p.inputs["Base Color"].default_value = (*rgb, 1.0)
    p.inputs["Roughness"].default_value = roughness
    p.inputs["Alpha"].default_value = alpha
    if "Specular IOR Level" in p.inputs:
        p.inputs["Specular IOR Level"].default_value = 0.2
    m.use_backface_culling = True
    if alpha < 1.0:
        try:
            m.surface_render_method = "BLENDED"
        except AttributeError:
            m.blend_method = "BLEND"
    return m


def _bl_setup(bpy, res, samples):
    import mathutils
    bpy.ops.wm.read_factory_settings(use_empty=True)
    sc = bpy.context.scene
    sc.render.engine = "BLENDER_EEVEE"
    sc.render.resolution_x, sc.render.resolution_y = res
    sc.render.resolution_percentage = 100
    sc.render.film_transparent = True
    sc.render.image_settings.file_format = "PNG"
    sc.render.image_settings.color_mode = "RGBA"
    sc.view_settings.view_transform = "Standard"
    sc.view_settings.look = "None"
    try:
        sc.eevee.taa_render_samples = samples
        sc.eevee.use_bloom = False
    except AttributeError:
        pass
    world = bpy.data.worlds.new("world")
    world.use_nodes = True
    world.node_tree.nodes["Background"].inputs["Color"].default_value = (0.55, 0.55, 0.55, 1.0)
    world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.6
    sc.world = world
    sun = bpy.data.lights.new("sun", "SUN")
    sun.energy = 2.2
    sun.angle = np.radians(2.0)
    so = bpy.data.objects.new("sun", sun)
    az, el = np.radians(315.0), np.radians(45.0)                  # from the NW, 45 deg up
    to_sun = mathutils.Vector((np.sin(az) * np.cos(el), np.cos(az) * np.cos(el), np.sin(el)))
    so.rotation_euler = (-to_sun).to_track_quat("-Z", "Y").to_euler()
    sc.collection.objects.link(so)
    cam = bpy.data.cameras.new("cam")
    cam.lens = CAM_LENS_MM
    cam.clip_start, cam.clip_end = 0.1, 1000.0
    co = bpy.data.objects.new("cam", cam)
    co.location = CAM_LOC_KM
    d = mathutils.Vector(CAM_TARGET_KM) - mathutils.Vector(CAM_LOC_KM)
    co.rotation_euler = d.to_track_quat("-Z", "Y").to_euler()
    sc.collection.objects.link(co)
    sc.camera = co
    return sc


def _north_arrow(bpy, X, Y, ve):
    """A dark-grey north arrow (+y) with an "N", just above the shelf water at the NW corner."""
    x0, y0 = float(X.max()) * 0.05, float(Y.max()) * 0.80
    z = 0.15
    w = 0.5
    verts = np.array([[x0 - w, y0, z], [x0 + w, y0, z], [x0 + w, y0 + 4, z], [x0 - w, y0 + 4, z],
                      [x0 - 1.4, y0 + 4, z], [x0 + 1.4, y0 + 4, z], [x0, y0 + 7, z]])
    m = _bl_material(bpy, "north_m", rgb=(0.12, 0.12, 0.12))
    for name, f in (("north", np.array([[0, 1, 2, 3]])), ("north_head", np.array([[4, 5, 6]]))):
        ob = _bl_mesh(bpy, name, verts, f, smooth=False)
        ob.data.materials.append(m)
    cu = bpy.data.curves.new("north_txt", "FONT")
    cu.body = "N"
    cu.size = 3.0
    cu.align_x = "CENTER"
    tx = bpy.data.objects.new("north_txt", cu)
    tx.location = (x0, y0 + 7.8, z)
    az = np.arctan2(CAM_TARGET_KM[0] - CAM_LOC_KM[0], CAM_TARGET_KM[1] - CAM_LOC_KM[1])
    tx.rotation_euler = (0.0, 0.0, -float(az))          # reads upright from the camera
    tx.data.materials.append(m)
    bpy.context.scene.collection.objects.link(tx)


def render_main(args) -> int:
    import bpy
    meta, ter = load_scene(args.scene)
    ve = float(meta["ve"])
    L = float(meta["colour_limit_m"])
    lut = np.asarray(meta["lut_diverging"])
    frames = stills(meta) if args.stills else (list(range(len(meta["frames"]))) if args.all else args.frames)
    X, Y = km_coords(ter, meta)
    ny, nx = X.shape
    dsl = data_slice(meta, X.shape)
    Xd, Yd = X[dsl], Y[dsl]
    out = Path(args.scene) / "panels"
    out.mkdir(exist_ok=True)
    timing = []
    for run in meta["runs"]:
        sc = _bl_setup(bpy, tuple(args.res), args.samples)
        zt, rgbt = terrain_height_rgb(ter, ve, float(meta["ve_land"]))
        terr = _bl_mesh(bpy, "terrain", np.stack([X, Y, zt], -1).reshape(-1, 3),
                        quad_faces(np.ones((ny, nx), bool)), rgbt)
        terr.data.materials.append(_bl_material(bpy, "terrain_m", roughness=0.9))
        _north_arrow(bpy, X, Y, ve)
        ice = ter.get(f"ice_{run}")
        if ice is not None:                      # slab on the exact ice cells (cell edges)
            h = float(meta["extra"]["ice_h_m"][run]) * ve / 1e3
            dx = float(X[0, 1] - X[0, 0])
            j, i = np.nonzero(ice >= 0.5)
            xc, yc = X[j, i], Y[j, i]
            c = np.stack([np.stack([xc - dx / 2, yc - dx / 2], -1), np.stack([xc + dx / 2, yc - dx / 2], -1),
                          np.stack([xc + dx / 2, yc + dx / 2], -1), np.stack([xc - dx / 2, yc + dx / 2], -1)], 1)
            # top face of the slab only, at 10 % of the ice thickness above sea level (90 % below):
            # the slab is h x VE thick, but a lower face seen through the upper one renders as
            # unlit grey, so only the white top is drawn
            vt = np.concatenate([c, np.full(c.shape[:2] + (1,), 0.1 * h)], -1).reshape(-1, 3)
            faces = np.arange(vt.shape[0]).reshape(-1, 4)
            slab = _bl_mesh(bpy, "ice", vt, faces, smooth=False)
            slab.visible_shadow = False          # no shadow-map noise on the isopycnal below
            im = _bl_material(bpy, "ice_m", alpha=ICE_ALPHA, rgb=(1.0, 1.0, 1.0), roughness=0.6)
            pr = im.node_tree.nodes["Principled BSDF"]          # flat white, independent of the sun
            pr.inputs["Emission Color"].default_value = (1.0, 1.0, 1.0, 1.0)
            pr.inputs["Emission Strength"].default_value = ICE_EMISSION
            slab.data.materials.append(im)
        sp = ter.get("sponge")
        if sp is not None:                       # open boundary + sponge (data.obcs) masked in grey
            dx = float(X[0, 1] - X[0, 0])
            j, i = np.nonzero(sp & ~ter["land"].astype(bool))
            xc, yc = X[j, i], Y[j, i]
            zs = float(meta["sigma"][run]["z_m"]) * ve / 1e3
            c = np.stack([np.stack([xc - dx / 2, yc - dx / 2], -1), np.stack([xc + dx / 2, yc - dx / 2], -1),
                          np.stack([xc + dx / 2, yc + dx / 2], -1), np.stack([xc - dx / 2, yc + dx / 2], -1)], 1)
            vs = np.concatenate([c, np.full(c.shape[:2] + (1,), zs)], -1).reshape(-1, 3)
            spo = _bl_mesh(bpy, "sponge", vs, np.arange(vs.shape[0]).reshape(-1, 4), smooth=False)
            spo.data.materials.append(_bl_material(bpy, "sponge_m", alpha=ISO_ALPHA, rgb=SPONGE_RGB,
                                                   roughness=0.8))
        iso_m = _bl_material(bpy, "iso_m", alpha=ISO_ALPHA, roughness=0.7)
        for k in frames:
            t0 = time.time()
            fr = np.load(Path(args.scene) / "frames" / run / f"{int(k):04d}.npz")
            z, col = fr["z_m"], fr["colour_m"]
            ok = np.isfinite(z)
            rgb = lut_colour(col, L, lut, nan_rgb=NAN_RGB)
            verts = np.stack([Xd, Yd, np.where(ok, z, 0.0) * ve / 1e3], -1).reshape(-1, 3)
            faces = quad_faces(ok)
            old = bpy.data.objects.get("iso")
            if old is not None:
                me = old.data
                bpy.data.objects.remove(old)
                bpy.data.meshes.remove(me)
            iso = _bl_mesh(bpy, "iso", verts, faces, rgb)
            iso.data.materials.append(iso_m)
            sc.render.filepath = str(out / f"{run}_{int(k):04d}.png")
            bpy.ops.render.render(write_still=True)
            timing.append({"run": run, "frame": int(k), "s": time.time() - t0})
            print(f"RENDERED {run} {k} {timing[-1]['s']:.1f} s", flush=True)
    (out / "timing.json").write_text(json.dumps(timing, indent=1))
    return 0


# ---- overlay (pyocean env) ----------------------------------------------------------------------

PANEL_TITLE = {"F1": "F1: summer profile, open water",
               "F3": "F3: winter profile, open water",
               "F2": "F2: winter profile, landfast ice (white slab = AREA.bin ice cells)"}


def compose_frame(meta: dict, scene: Path, k: int, out_png: Path, dpi: int = 100) -> None:
    """One overlay frame: header band (label, isopycnal, what colour and geometry show, VE, clock),
    the panel(s) side by side on white, footer band (wind pulse, colour bar, notes)."""
    import textwrap

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, Normalize
    from PIL import Image

    runs = meta["runs"]
    imgs = [Image.open(Path(scene) / "panels" / f"{r}_{k:04d}.png").convert("RGBA") for r in runs]
    w0, h0 = imgs[0].size
    imgs = [im.crop((0, int(CROP_TOP * h0), w0, h0)) for im in imgs]    # fixed crop: empty sky
    w, h = imgs[0].size
    HEAD, FOOT = 150, 150                                     # px
    W, H = w * len(imgs), h + HEAD + FOOT
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi, facecolor="white")
    for q, (r, im) in enumerate(zip(runs, imgs)):
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        bg.alpha_composite(im)
        ax = fig.add_axes([q * w / W, FOOT / H, w / W, h / H])
        ax.imshow(np.asarray(bg.convert("RGB")), interpolation="lanczos")
        ax.set_axis_off()
        ax.text(0.015, 0.985, PANEL_TITLE.get(r, r), transform=ax.transAxes, ha="left", va="top",
                fontsize=14, color="0.1")
    fr = meta["frames"][k]
    sig = meta["sigma"][runs[0]]
    L = float(meta["colour_limit_m"])
    tf = meta["extra"]["T_f_h"]
    y0 = 1 - 18 / H
    fig.text(12 / W, y0, "DESCRIPTIVE", fontsize=17, ha="left", va="top", color="0.05", weight="bold")
    fig.text(200 / W, y0, "Exp. 4, Van Mijenfjorden (MITgcm): one isopycnal after a single wind pulse", fontsize=16, ha="left", va="top", color="0.1")
    fig.text(1 - 14 / W, y0 - 44 / H, fr["label"], fontsize=17, ha="right", va="top",
             family="monospace", color="0.05")
    sub = (f"Isopycnal $\\rho$ = {sig['sigma']:.4f} kg m$^{{-3}}$ (depth {-sig['z_m']:.1f} m "
           f"at t = 0) over the model bathymetry, land from IBCAO. Geometry: raw isopycnal depth. "
           f"Colour: its displacement $\\zeta'$, with the running mean over one inertial period "
           f"(T$_f$ = {tf:.2f} h) removed. Vertical exaggeration {meta['ve']:g}x below sea level, "
           f"{meta['ve_land']:g}x on land. View from the southwest over Bellsund, looking up-fjord.")
    fig.text(12 / W, y0 - 40 / H, textwrap.fill(sub, int(0.66 * W / 7.6)), fontsize=11.5, ha="left",
             va="top", color="0.25", linespacing=1.4)
    # footer: wind pulse (left), colour bar (centre), notes (right)
    tt = (np.asarray(meta["extra"]["tau_t_s"]) - meta["onset_s"]) / 3600.0
    tau = np.asarray(meta["extra"]["tau_peak_Nm2"])
    axw = fig.add_axes([50 / W, 45 / H, min(0.2, 380 / W), 60 / H])
    axw.plot(tt, tau, color="0.35", lw=1.0)
    th = (fr["t_s"] - meta["onset_s"]) / 3600.0
    axw.axvline(th, color="0.1", lw=1.0)
    axw.set_xlim(tt[0], tt[-1])
    axw.set_ylim(0, tau.max() * 1.25)
    axw.set_yticks([0, round(float(tau.max()), 2)])
    axw.tick_params(labelsize=9)
    axw.set_xlabel("hours after the pulse onset", fontsize=9.5)
    axw.set_title(r"wind pulse: peak $|\tau|$ of the forcing (N m$^{-2}$)", fontsize=10, loc="left")
    for s_ in ("top", "right"):
        axw.spines[s_].set_visible(False)
    cax = fig.add_axes([0.5 - 0.17, 70 / H, 0.34, 22 / H])
    cm = ListedColormap(np.asarray(meta["lut_diverging"]))
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-L, L), cmap=cm), cax=cax,
                      orientation="horizontal", extend="both")
    cb.set_label(f"$\\zeta'$ (m, + up); fixed scale $\\pm${L:g} m for every frame", fontsize=11)
    cb.ax.tick_params(labelsize=10)
    notes = ["Descriptive only: no number from this animation enters the analysis.",
             "Grey on the isopycnal: no filtered value (isopycnal absent in part of the window); "
             "holes: outcropped or meeting the bed."]
    if "sponge_cells" in meta["extra"]:
        notes.append("Grey strip along the W, N and S edges: the open boundary and its "
                     f"{meta['extra']['sponge_cells']}-cell sponge (data.obcs), masked in both runs.")
    if "F2" in runs:
        notes.append("Under the F2 landfast ice the ocean receives no wind stress.")
    step = sig.get("step_kg_m3")
    if step is not None and step < 0.01:
        notes.append(f"The isopycnal sits on a {step:.4f} kg m$^{{-3}}$ step: metre-scale coastal values "
                     "may be mixing rather than waves.")
    fig.text(1 - 14 / W, 12 / H, "\n".join(textwrap.fill(n, 70) for n in notes), fontsize=9.5,
             ha="right", va="bottom", color="0.3", linespacing=1.35)
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, facecolor="white")
    plt.close(fig)


def compose_main(args) -> int:
    meta = json.loads((Path(args.scene) / "meta.json").read_text())
    if args.stills:
        for s_ in meta["extra"]["stills"]:
            compose_frame(meta, args.scene, s_["index"],
                          Path(args.out_dir) / f"{args.name}_{int(round(s_['requested_h'])):02d}h.png")
        return 0
    frames = list(range(len(meta["frames"]))) if args.all else args.frames
    for k in frames:
        compose_frame(meta, args.scene, int(k), Path(args.out_dir) / "frames" / f"{int(k):04d}.png")
    if args.all:                                   # saturation of the fixed scale, at render time
        from exp4_arctic_bayiw.analysis.blender_export import saturation_stats
        L = float(meta["colour_limit_m"])
        col = {r: np.stack([np.load(Path(args.scene) / "frames" / r / f"{k:04d}.npz")["colour_m"]
                            for k in frames]) for r in meta["runs"]}
        st = saturation_stats(col, L)
        for r in st:
            st[r]["worst_t_h"] = (meta["frames"][frames[st[r]["worst_index"]]]["t_s"] - meta["onset_s"]) / 3600
        st["limit_m"] = L
        st["frames"] = len(frames)
        (Path(args.out_dir) / "saturation.json").write_text(json.dumps(st, indent=1))
        print(json.dumps(st, indent=1))
    return 0


def movie_main(args) -> int:
    """H.264 MP4 (ffmpeg; crf raised until <= --max-mb-mp4) and a GIF preview (Pillow; the frame
    stride doubles until <= --max-mb-gif) of OUT/frames/*.png."""
    import shutil
    import subprocess

    from PIL import Image
    fr = Path(args.out_dir) / "frames"
    dest = Path(args.movie_out or args.out_dir)
    dest.mkdir(parents=True, exist_ok=True)
    mp4 = dest / f"{args.name}.mp4"
    pngs = sorted(fr.glob("*.png"))
    res = {"frames": len(pngs), "fps": args.fps}
    if shutil.which("ffmpeg") is None:
        raise SystemExit("ffmpeg not found")
    for crf in (23, 26, 29, 32, 35):
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(args.fps), "-i",
                        str(fr / "%04d.png"), "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", str(crf),
                        "-movflags", "+faststart", str(mp4)], check=True)
        mb = mp4.stat().st_size / 1e6
        res.update(mp4=str(mp4), mp4_mb=mb, crf=crf)
        if mb <= args.max_mb_mp4:
            break
    gif = dest / f"{args.name}.gif"
    stride = 1
    while True:
        sel = pngs[::stride]
        ims = []
        for p_ in sel:
            im = Image.open(p_).convert("RGB")
            w = args.gif_width
            im = im.resize((w, int(round(im.height * w / im.width))), Image.LANCZOS)
            ims.append(im.quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE))
        ims[0].save(gif, save_all=True, append_images=ims[1:], duration=int(1000 * stride / args.fps),
                    loop=0, optimize=True)
        mb = gif.stat().st_size / 1e6
        if mb <= args.max_mb_gif or stride >= 64:
            break
        stride *= 2
    res.update(gif=str(gif), gif_mb=mb, gif_stride=stride, gif_frames=len(sel), gif_width=args.gif_width)
    print(json.dumps(res, indent=1))
    (Path(args.out_dir) / "movie.json").write_text(json.dumps(res, indent=1))
    # the numbers the README quotes (render_readme.py), next to the MP4 and GIF
    meta = json.loads((Path(args.scene) / "meta.json").read_text())
    satp = Path(args.out_dir) / "saturation.json"
    r0 = meta["runs"][0]
    sig = meta["sigma"][r0]
    t = [f_["t_s"] for f_ in meta["frames"]]
    summary = {
        "shot": meta["shot"], "runs": meta["runs"],
        "sigma_kg_m3": sig["sigma"], "sigma_depth_t0_m": -sig["z_m"], "step_kg_m3": sig.get("step_kg_m3"),
        "sigma_rule": sig["rule"], "ve_sea": meta["ve"], "ve_land": meta["ve_land"],
        "colour_limit_m": meta["colour_limit_m"],
        "colour_limit_rule": meta["extra"].get("colour_limit_rule", "nice ceiling of the 99.5th percentile of |zeta'|"),
        "T_f_h": meta["extra"]["T_f_h"],
        "t_first_h": (t[0] - meta["onset_s"]) / 3600, "t_last_h": (t[-1] - meta["onset_s"]) / 3600,
        "frames": len(t), "sponge_cells": meta["extra"].get("sponge_cells"),
        "saturation": json.loads(satp.read_text()) if satp.exists() else None,
        "mp4": {"file": mp4.name, "mb": round(res["mp4_mb"], 2), "crf": res["crf"], "fps": args.fps},
        "gif": {"file": gif.name, "mb": round(res["gif_mb"], 2), "stride": res["gif_stride"],
                "width_px": args.gif_width},
    }
    (dest / f"{args.name}.json").write_text(json.dumps(summary, indent=1) + "\n")
    return 0 if res["mp4_mb"] <= args.max_mb_mp4 and mb <= args.max_mb_gif else 1


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[1:]
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("mode", choices=["render", "compose", "movie"])
    ap.add_argument("--scene", type=Path, required=True)
    ap.add_argument("--frames", type=int, nargs="+", default=[])
    ap.add_argument("--stills", action="store_true")
    ap.add_argument("--all", action="store_true", help="every exported frame")
    ap.add_argument("--max-mb-mp4", type=float, default=20.0)
    ap.add_argument("--max-mb-gif", type=float, default=5.0)
    ap.add_argument("--gif-width", type=int, default=720)
    ap.add_argument("--movie-out", type=Path, default=None, help="directory of the MP4 and GIF")
    ap.add_argument("--res", type=int, nargs=2, default=list(RES))
    ap.add_argument("--samples", type=int, default=32)
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--name", default="shot")
    ap.add_argument("--fps", type=int, default=12)
    a = ap.parse_args(argv)
    return {"render": render_main, "compose": compose_main, "movie": movie_main}[a.mode](a)


if __name__ == "__main__":
    rc = main()
    if "bpy" not in sys.modules:
        sys.exit(rc)
