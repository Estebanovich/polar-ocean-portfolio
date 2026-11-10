"""Figure 0d, step 2: a perspective render of the ISOMIP+ geometry with Blender (headless).

    blender -b --python fig0_geometry_blender.py -- [--work DIR] [--samples N] [--res W H]

Reads WORK/geom.npz (written by fig0_geometry_export.py from the official ISOMIP+ geometry file) and renders
WORK/geometry_render.png (transparent background), plus WORK/geometry_points.json: the image positions of the 3-D
points that fig0_geometry.py labels. Only the input geometry is drawn (bed, ice draft and surface, masks); the
arrows are a schematic of published circulation features, not model output.

Scene units: 1 unit = 1 km in x and y; z = elevation x VE / 1000, with vertical exaggeration VE = 100.
The model x = 320-800 km is drawn at X = x - 320. The geometry is symmetric about the centreline y = 40 km (to within
8 m in the file), so the near half (y < 40 km) is cut away and its outline drawn as thin lines: the cut face shows the
ice shelf, the cavity and the bed along the centre of the trough.
"""
import argparse
import json
import math
import sys
import tempfile
from pathlib import Path

import bmesh
import bpy
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector

VE = 100.0
VIEW = (0.30, -1.0, 0.50)       # direction from the scene centre to the camera
ZBASE = -0.72 * VE           # base of the block: the domain bottom, 720 m below sea level
ZCAP = 0.20 * VE             # the grounded-ice surface is drawn capped at 200 m above sea level (the file reaches
                             # 878 m, which at VE = 100 towers over the 720 m deep box); a drawing choice only


def args():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=str(Path(tempfile.gettempdir()) / "fig0_geometry"))
    ap.add_argument("--samples", type=int, default=96)
    ap.add_argument("--res", type=int, nargs=2, default=[2400, 1200])
    ap.add_argument("--lens", type=float, default=35.0)
    ap.add_argument("--margin", type=float, default=40.0)
    ap.add_argument("--points-only", action="store_true", help="write geometry_points.json without rendering")
    return ap.parse_args(argv)


# ------------------------------------------------------------------------------------------- mesh builders
def column_mesh(name, X, Y, ztop, zbot, mask):
    """A solid between two height fields on a vertex grid: top and bottom faces where all four corners are in
    `mask`, and vertical walls along the boundary of that region."""
    ny, nx = mask.shape
    bm = bmesh.new()
    vt, vb = {}, {}
    for j in range(ny):
        for i in range(nx):
            if mask[j, i]:
                vt[j, i] = bm.verts.new((X[i], Y[j], ztop[j, i]))
                vb[j, i] = bm.verts.new((X[i], Y[j], zbot[j, i]))
    quads = np.zeros((ny - 1, nx - 1), bool)
    for j in range(ny - 1):
        for i in range(nx - 1):
            if mask[j, i] and mask[j, i + 1] and mask[j + 1, i] and mask[j + 1, i + 1]:
                quads[j, i] = True
                c = ((j, i), (j, i + 1), (j + 1, i + 1), (j + 1, i))
                bm.faces.new([vt[k] for k in c])
                bm.faces.new([vb[k] for k in c[::-1]])
    edges = {}
    for j in range(ny - 1):
        for i in range(nx - 1):
            if quads[j, i]:
                c = ((j, i), (j, i + 1), (j + 1, i + 1), (j + 1, i))
                for a, b in zip(c, c[1:] + c[:1]):
                    key = tuple(sorted((a, b)))
                    edges[key] = edges.get(key, 0) + 1
    for (a, b), n in edges.items():
        if n == 1:
            try:
                bm.faces.new([vt[a], vt[b], vb[b], vb[a]])
            except ValueError:
                pass
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    me = bpy.data.meshes.new(name)
    bm.to_mesh(me)
    bm.free()
    ob = bpy.data.objects.new(name, me)
    bpy.context.scene.collection.objects.link(ob)
    for p in me.polygons:
        p.use_smooth = False
    return ob


def with_cut(ob, mat, cutmat, ycut):
    """Give the faces on the cut plane y = ycut (normal toward -y) their own material, so the section reads."""
    me = ob.data
    me.materials.append(mat)
    me.materials.append(cutmat)
    for p in me.polygons:
        if p.normal.y < -0.9 and abs(p.center.y - ycut) < 1e-3:
            p.material_index = 1
    return ob


def material(name, color, alpha=1.0, rough=0.6, noise_bump=0.0, emission=0.0):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    bsdf = nt.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*color, 1.0)
    bsdf.inputs["Roughness"].default_value = rough
    bsdf.inputs["Alpha"].default_value = alpha
    if emission:
        bsdf.inputs["Emission Color"].default_value = (*color, 1.0)
        bsdf.inputs["Emission Strength"].default_value = emission
    if noise_bump:
        tex = nt.nodes.new("ShaderNodeTexNoise")
        tex.inputs["Scale"].default_value = 0.15
        tex.inputs["Detail"].default_value = 8.0
        bump = nt.nodes.new("ShaderNodeBump")
        bump.inputs["Strength"].default_value = noise_bump
        nt.links.new(tex.outputs["Fac"], bump.inputs["Height"])
        nt.links.new(bump.outputs["Normal"], bsdf.inputs["Normal"])
        ramp = nt.nodes.new("ShaderNodeValToRGB")
        ramp.color_ramp.elements[0].color = (color[0] * 0.85, color[1] * 0.85, color[2] * 0.85, 1)
        ramp.color_ramp.elements[1].color = (min(1, color[0] * 1.1), min(1, color[1] * 1.1), min(1, color[2] * 1.1), 1)
        nt.links.new(tex.outputs["Fac"], ramp.inputs["Fac"])
        nt.links.new(ramp.outputs["Color"], bsdf.inputs["Base Color"])
    if alpha < 1:
        m.blend_method = "BLEND" if hasattr(m, "blend_method") else None
    return m


def depth_tint(m, z0, z1, c0, c1):
    """Multiply a material's base colour by a ramp in z (deep: c0, shallow: c1), so the bed trough reads."""
    nt = m.node_tree
    bsdf = nt.nodes["Principled BSDF"]
    src = bsdf.inputs["Base Color"].links[0].from_socket if bsdf.inputs["Base Color"].is_linked else None
    geo = nt.nodes.new("ShaderNodeNewGeometry")
    sep = nt.nodes.new("ShaderNodeSeparateXYZ")
    mr = nt.nodes.new("ShaderNodeMapRange")
    mr.inputs["From Min"].default_value, mr.inputs["From Max"].default_value = z0, z1
    ramp = nt.nodes.new("ShaderNodeValToRGB")
    ramp.color_ramp.elements[0].color = (*c0, 1)
    ramp.color_ramp.elements[1].color = (*c1, 1)
    mix = nt.nodes.new("ShaderNodeMix")
    mix.data_type = "RGBA"
    mix.blend_type = "MULTIPLY"
    mix.inputs["Factor"].default_value = 1.0
    nt.links.new(geo.outputs["Position"], sep.inputs["Vector"])
    nt.links.new(sep.outputs["Z"], mr.inputs["Value"])
    nt.links.new(mr.outputs["Result"], ramp.inputs["Fac"])
    if src is not None:
        nt.links.new(src, mix.inputs["A"])
    else:
        mix.inputs["A"].default_value = (1, 1, 1, 1)
    nt.links.new(ramp.outputs["Color"], mix.inputs["B"])
    nt.links.new(mix.outputs["Result"], bsdf.inputs["Base Color"])


def srgb(h):
    h = h.lstrip("#")
    c = [int(h[k:k + 2], 16) / 255 for k in (0, 2, 4)]
    return tuple(((v + 0.055) / 1.055) ** 2.4 if v > 0.04045 else v / 12.92 for v in c)


def tube(name, pts, radius, mat, head=True):
    """A curve with a round section through `pts`, with a cone at the end."""
    cu = bpy.data.curves.new(name, "CURVE")
    cu.dimensions = "3D"
    cu.bevel_depth = radius
    cu.bevel_resolution = 4
    sp = cu.splines.new("POLY")
    pts = [Vector(p) for p in pts]
    if head:
        d = (pts[-1] - pts[-2]).normalized()
        tip = pts[-1]
        pts[-1] = tip - d * radius * 3.2
    sp.points.add(len(pts) - 1)
    for k, p in enumerate(pts):
        sp.points[k].co = (*p, 1)
    ob = bpy.data.objects.new(name, cu)
    bpy.context.scene.collection.objects.link(ob)
    ob.data.materials.append(mat)
    if head:
        bpy.ops.mesh.primitive_cone_add(radius1=radius * 2.4, depth=radius * 4.2, location=tip - d * radius * 2.1)
        cone = bpy.context.active_object
        cone.rotation_mode = "QUATERNION"
        cone.rotation_quaternion = Vector((0, 0, 1)).rotation_difference(d)
        cone.data.materials.append(mat)
    return ob


def smooth_path(pts, n=12):
    """Catmull-Rom through the control points."""
    P = [Vector(p) for p in pts]
    P = [P[0]] + P + [P[-1]]
    out = []
    for k in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[k - 1], P[k], P[k + 1], P[k + 2]
        for t in [i / n for i in range(n)]:
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t
                              + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    out.append(P[-2])
    return [tuple(v) for v in out]


# ------------------------------------------------------------------------------------------- scene
def main():
    a = args()
    work = Path(a.work)
    g = np.load(work / "geom.npz")
    X = g["x"] / 1e3 - 320.0
    Y = g["y"] / 1e3
    bed = g["bed"] * VE / 1e3
    lower = g["lower"] * VE / 1e3
    upper = g["upper"] * VE / 1e3
    floating, grounded, openm = g["floating"] > 0, g["grounded"] > 0, g["open"] > 0
    ny, nx = bed.shape
    jc = ny // 2                                             # first row of the far half (y = 40.5 km)
    far = np.zeros_like(floating)
    far[jc:, :] = True
    D = 0.72 * VE                                            # domain depth, 720 m

    bpy.ops.wm.read_factory_settings(use_empty=True)
    sc = bpy.context.scene

    rock = material("rock", (1.0, 1.0, 1.0), rough=0.9, noise_bump=0.5)
    depth_tint(rock, -D, -0.1 * VE, srgb("#5f4f40"), srgb("#cbb9a1"))
    ice = material("ice", srgb("#f3f9fd"), alpha=0.75, rough=0.3)
    ice_g = material("ice_g", srgb("#f0f6fa"), rough=0.55)
    ice_cut = material("ice_cut", srgb("#cfe4f2"), rough=0.6)
    water = material("water", srgb("#3f86bf"), alpha=0.28, rough=0.1)
    water_cut = material("water_cut", srgb("#2f75b0"), alpha=0.5, rough=0.6)
    rock_cut = material("rock_cut", srgb("#6f5d4b"), rough=0.9, noise_bump=0.4)
    melt = material("melt", srgb("#8a6fd0"), alpha=0.9, rough=0.5, emission=0.3)
    warm = material("warm", srgb("#c0392b"), rough=0.5, emission=0.25)
    purp = material("purp", srgb("#5b3f8f"), rough=0.5, emission=0.25)
    restore = material("restore", srgb("#f39c4b"), alpha=0.5, rough=0.5, emission=0.2)

    # rock block: the bed over the whole domain, with walls down to ZBASE
    yc = float(Y[jc])
    with_cut(column_mesh("bed", X, Y, bed, np.full_like(bed, ZBASE), far), rock, rock_cut, yc)
    # ice, far half: grounded ice opaque (its base on the bed), the floating ice shelf translucent
    upper_g = np.minimum(upper, ZCAP)                        # capped grounded-ice surface (drawing only)
    with_cut(column_mesh("ice_grounded", X, Y, upper_g, bed + 0.05, grounded & far), ice_g, ice_cut, yc)
    with_cut(column_mesh("ice_shelf", X, Y, upper, lower, floating & far), ice, ice_cut, yc)
    # water, far half: under the ice shelf (up to the draft) and in the open ocean (up to sea level)
    wtop = np.where(floating, lower, 0.0)
    with_cut(column_mesh("water", X, Y, wtop - 0.03, bed + 0.03, (floating | openm) & far), water, water_cut, yc)
    # meltwater-tracer source: a thin sheet just under the ice base (floating ice, far half)
    column_mesh("tracer", X, Y, lower - 0.1, lower - 2.2, floating & far).data.materials.append(melt)
    # restoring zone x = 790-800 km, full width and depth
    bpy.ops.mesh.primitive_cube_add(location=(475.0, 60.25, -D / 2))
    rz = bpy.context.active_object
    rz.scale = (5.0, 19.75, D / 2)
    rz.data.materials.append(restore)

    # circulation sketch (process 7): the cyclonic (clockwise seen from above, Southern Hemisphere) cell that
    # Yung et al. (2026) find in most Ocean1 models, boundary-intensified in y > 50 km: inflow at depth along its
    # inner side, outflow along the y = 80 km wall under the ice base and out past the front.
    def zb(xk, yk, off):                                     # off below the ice base (or sea level)
        i, j = int(np.clip(round(xk), 0, nx - 1)), int(np.clip(round(yk), 0, ny - 1))
        return (lower[j, i] if floating[j, i] else 0.0) + off
    out = [(146, 74, zb(146, 74, -6)), (210, 74, zb(210, 74, -5)), (270, 74, zb(270, 74, -5)),
           (318, 74, zb(318, 74, -5)), (370, 72, -0.3 * VE), (430, 66, -0.3 * VE)]
    tube("outflow", smooth_path(out), 2.0, purp)
    def zbd(xk, yk, off):                                    # off above the bed
        return float(bed[int(round(yk)), int(round(xk))]) + off
    inflow = [(465, 53, zbd(465, 53, 9)), (380, 53, zbd(380, 53, 9)), (280, 53, zbd(280, 53, 9)),
              (190, 55, zbd(190, 55, 9)), (150, 62, zbd(150, 62, 12)), (140, 70, zb(140, 70, -8))]
    tube("inflow", smooth_path(inflow), 2.2, warm)
    # outline of the cut-away near half of the domain box (sea level to the bottom, y = 0-40 km)
    wire = material("wire", srgb("#33424d"), rough=0.6)
    for p0, p1 in (((0, 0, -D), (480, 0, -D)), ((480, 0, -D), (480, 0, 0)), ((480, 0, 0), (0, 0, 0)),
                   ((0, 0, 0), (0, 0, -D)), ((0, 0, -D), (0, 40, -D)), ((480, 0, -D), (480, 40, -D)),
                   ((480, 0, 0), (480, 40, 0)), ((0, 0, 0), (0, 40, 0))):
        tube("wire", [p0, p1], 0.22, wire, head=False)

    # camera: perspective along VIEW, fitted to the domain box
    cam_d = bpy.data.cameras.new("cam")
    cam_d.lens = a.lens
    cam_d.sensor_width = 36
    cam = bpy.data.objects.new("cam", cam_d)
    sc.collection.objects.link(cam)
    sc.camera = cam
    sc.render.resolution_x, sc.render.resolution_y = a.res
    sc.render.resolution_percentage = 100
    view = Vector(VIEW).normalized()
    cam.rotation_mode = "QUATERNION"
    cam.rotation_quaternion = (-view).to_track_quat("-Z", "Y")
    box = [Vector((x, y, z)) for x in (0, 480) for y in (0, 80) for z in (ZBASE, float(upper.max()))]   # framing as uncapped
    cam.location = Vector((240, 40, 0)) + view * 600
    bpy.context.view_layer.update()
    loc, _ = cam.camera_fit_coords(bpy.context.evaluated_depsgraph_get(), [c for v in box for c in v])
    cam.location = Vector(loc) + view * a.margin
    bpy.context.view_layer.update()                          # matrix_world must be current for world_to_camera_view
    sun = bpy.data.objects.new("sun", bpy.data.lights.new("sun", "SUN"))
    sun.data.energy = 2.2
    sun.data.angle = math.radians(10)
    sun.rotation_euler = (math.radians(48), math.radians(12), math.radians(-30))
    sc.collection.objects.link(sun)
    world = bpy.data.worlds.new("w")
    world.use_nodes = True
    world.node_tree.nodes["Background"].inputs["Color"].default_value = (0.80, 0.85, 0.90, 1)
    world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.65
    sc.world = world

    sc.render.engine = "CYCLES"
    sc.cycles.device = "CPU"
    sc.cycles.samples = a.samples
    sc.cycles.use_denoising = True
    sc.render.film_transparent = True
    sc.view_settings.view_transform = "Standard"
    sc.render.use_freestyle = True
    sc.render.line_thickness_mode = "ABSOLUTE"
    sc.render.line_thickness = 1.3 * a.res[0] / 2400
    ls = sc.view_layers[0].freestyle_settings.linesets[0]
    ls.select_by_visibility = True
    ls.select_crease = False
    ls.select_border = True
    ls.select_silhouette = True
    if ls.linestyle is None:
        ls.linestyle = bpy.data.linestyles.new("lines")
    ls.linestyle.color = (0.20, 0.25, 0.30)
    ls.linestyle.alpha = 0.75
    sc.render.filepath = str(work / "geometry_render.png")
    sc.render.image_settings.file_format = "PNG"
    sc.render.image_settings.color_mode = "RGBA"

    # image positions of the labelled points (fraction of the width from the left, of the height from the top)
    def P(x, y, z):
        v = world_to_camera_view(sc, cam, Vector((x, y, z)))
        return [round(v.x, 5), round(1 - v.y, 5)]
    glc = float(np.argmax(floating[jc]))                     # grounding line on the cut face
    pts = {
        "c_0_0_bot": P(0, 0, ZBASE), "c_480_0_bot": P(480, 0, ZBASE), "c_480_80_bot": P(480, 80, ZBASE),
        "c_480_0_sl": P(480, 0, 0), "c_480_0_720": P(480, 0, -D), "c_480_80_sl": P(480, 80, 0),
        "c_0_0_sl": P(0, 0, 0),
        "front_far_top": P(320, 80, float(upper[-1, 319])), "front_cut_top": P(319.5, Y[jc], float(upper[jc, 319])),
        "front_cut_base": P(319.5, Y[jc], float(lower[jc, 319])),
        "gl_cut": P(glc + 0.5, Y[jc], float(lower[jc, int(glc)])),
        "base_cut_mid": P(230, Y[jc], float(lower[jc, 230])),
        "cavity_cut": P(250, Y[jc], float(lower[jc, 250] + bed[jc, 250]) / 2),
        "restore_top": P(475, 80, 0), "restore_mid": P(475, 40.5, -D / 2),
        "bed_cut": P(150, Y[jc], float(bed[jc, 150])), "wall_far_ocean": P(400, 80, -0.25 * D),
        "near_box_mid": P(240, 0, -D / 2), "near_box_top": P(240, 0, 0), "c_0_40_sl": P(0, 40, 0),
        "c_0_0_bot2": P(0, 0, -D), "c_0_40_bot": P(0, 40, -D), "far_wall_mid": P(200, 80, float(bed[-1, 200])),
        "ice_grounded": P(50, 60, float(upper_g[60, 50])), "ice_shelf_top": P(240, 62, float(upper[62, 240])),
        "open_ocean": P(400, 75, 0), "outflow_end": P(430, 66, -0.3 * VE), "outflow_mid": P(270, 74, zb(270, 74, -5)),
        "inflow_mid": P(330, 53, zbd(330, 53, 9)), "inflow_end": P(150, 62, zbd(150, 62, 12)),
        "tracer_cut": P(200, Y[jc], float(lower[jc, 200])),
        "far_wall_ocean": P(430, 80, 0), "bed_grounded_cut": P(70, Y[jc], float(bed[jc, 70]) - 1.0),
        "outflow_open": P(370, 72, -0.3 * VE), "c_480_80_bot": P(480, 80, -D), "x480_mid_far": P(480, 80, -D / 2),
    }
    meta = {"VE": VE, "res": a.res, "points": pts, "gl_cut_x_km": glc + 320.0}
    (work / "geometry_points.json").write_text(json.dumps(meta, indent=1))
    if a.points_only:
        return
    bpy.ops.render.render(write_still=True)
    print("wrote", sc.render.filepath)


main()
