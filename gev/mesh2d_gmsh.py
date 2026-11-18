"""Body contour (geometry2d.Body) -> OpenFOAM 2D mesh via gmsh (2D, quad-dominant) + gmshToFoam.

Alternative to mesh2d_snappy (snappyHexMesh) with the same interface: mesh_case(case, body, profile) writes a
complete case from templates/2d with the same patch names (inlet, outlet, top, ground, body_<segment>,
front/back empty), forces and flow values (case2d.write_case_dicts). The 2D mesh is extruded to one cell
in z by write_extruded.

The profile gives absolute cell sizes instead of snappy refinement levels:
    size_far        largest cell (background)
    size_surface    cell size on the body contour, size_regions per segment regex overrides it
    size_corner     cell size at sharp corners (blunt trailing edges)
    size_ground     largest cell along the ground (caps the aspect ratio of the ground layers), None: size_far
    boxes           body -> [((x0, y0), (x1, y1), size), ...] refinement boxes
    envelopes       body -> [(segment regex, distance, size), ...]: size within the distance of those segments
    (size_regions and envelopes may be fixed or functions of the body)
    growth          cell size increase per unit distance away from contours and boxes (0.2 ~ ratio 1.2)
    layers          first cell, expansion ratio and count as in mesh2d_snappy (body_layers, ground_layers,
                    region_layers); quads, with fans at convex sharp corners
    concave_radius  sharp concave corners are rounded with this radius (gmsh's 2D boundary layer fails there)
    mesh_retry      rebuild a failed mesh with slightly scaled cell sizes (RETRY_SCALES)

Usage from Python:  mesh_case(case, body, scaled(TYRRELL, f))
"""

import re
import shutil
import time
from pathlib import Path

import gmsh
import numpy as np

from . import TEMPLATES
from . import mesh2d_snappy as snappy
from .case2d import write_case_dicts
from .foam import run, set_patch_types

# Tyrrell 026: the snappy profile mesh2d_snappy.TYRRELL as absolute sizes. Background cell B = Ahmed cell x c/L
# (16.6 mm); snappy level n gives B / 2^n: surface level 4 (1.04 mm), blunt trailing edge level 6 (0.26 mm),
# near box level 3, wake box level 2, gap box level 4 (gap ~ 1 mm cells). Nose box 0.52 mm (level 5) instead of
# snappy's level 6: 4 % fewer cells, C_l +0.4 % at h/c = 0.179 (check of 2026-10-06). The gap box reaches
# 10 mm above the trailing edge: it then also holds the separated region behind x/c ~ 0.75 near the ground.
# With the box ending below it (at 1.25 x the lowest wing point), C_l scattered by +-0.9 % between meshes
# whose cell sizes differed by only 2 % (h/c = 0.134 and 0.179; attached flow at 0.671: +-0.07 %); with the
# box, +-0.1 % (check of 2026-10-07).
# Same domain, layers and flow values as mesh2d_snappy.TYRRELL.
_C = snappy._C
_S = snappy._S
_B = 14 / 180 * _S


def tyrrell_boxes(body):
    P = body.points()
    (xn, yb), (xt, yt) = P.min(axis=0), P.max(axis=0)
    yle = P[P[:, 0].argmin(), 1]
    yte = body.segments["suction"][-1][1]      # suction-side trailing edge
    s, top = _S, 0.338
    return [((xn - 0.556 * s, 0.0), (xt + 0.2 * s, yt + (0.50 - top) * s), _B / 8),        # refinementBox
            ((xt - 0.35 * s, 0.0), (xt + 2.7 * s, yt + (1.0 - top) * s), _B / 4),          # wakeBox
            ((xn - 0.5 * _C, 0.0), (xt + 1.5 * _C, max(1.25 * yb, yte + 0.01)), _B / 16),  # gap + separation
            ((xn - 0.01, yle - 0.01), (xn + 0.01, yle + 0.01), _B / 32)]                   # nose


def tyrrell_double_boxes(body):
    """tyrrell_boxes for the double element wing (cases/tyrrell026/wing.py, double element): the same boxes over the whole
    wing, the gap box up to 10 mm above the flap trailing edge, nose boxes at both leading edges and a slot
    box around the main element trailing edge and the flap leading edge."""
    P = body.points()
    (xn, yb), (xt, yt) = P.min(axis=0), P.max(axis=0)
    yte = body.segments["flap_suction"][-1][1]
    s, top = _S, 0.338
    out = [((xn - 0.556 * s, 0.0), (xt + 0.2 * s, yt + (0.50 - top) * s), _B / 8),
           ((xt - 0.35 * s, 0.0), (xt + 2.7 * s, yt + (1.0 - top) * s), _B / 4),
           ((xn - 0.5 * _C, 0.0), (xt + 1.5 * _C, max(1.25 * yb, yte + 0.01)), _B / 16)]
    for el in ("main", "flap"):
        Q = np.vstack([body.segments[f"{el}_suction"], body.segments[f"{el}_pressure"]])
        le = Q[Q[:, 0].argmin()]
        out.append(((le[0] - 0.01, le[1] - 0.01), (le[0] + 0.01, le[1] + 0.01), _B / 32))
    te = body.segments["main_base"].mean(axis=0)
    out.append(((te[0] - 0.015, te[1] - 0.015), (te[0] + 0.01, te[1] + 0.025), _B / 32))   # slot
    return out


TYRRELL = dict(
    {k: snappy.TYRRELL[k] for k in ("x_min", "x_max", "y_max", "span", "n_procs", "mpi_args", "body_layers",
                                  "ground_layers", "region_layers", "layer_first", "layer_expansion",
                                  "U_inf", "k_inf", "omega_inf", "moving_ground")},
    size_far=_B, size_surface=_B / 16, size_regions=(), size_corner=_B / 64, size_ground=_B / 2,
    boxes=tyrrell_boxes, growth=0.2,
    corner_angle=30.0,     # turning angle [deg] above which a contour point is a corner (no spline through it)
    fan=7,                 # boundary-layer fan elements at convex corners
)


def scaled(p, f):
    """Profile with every cell size multiplied by f, for grid studies; the wall layers and the ground cap
    (which only limits their aspect ratio) stay."""
    boxes, regions, env = p["boxes"], p["size_regions"], p.get("envelopes")
    return dict(p, size_far=p["size_far"] * f, size_surface=p["size_surface"] * f, size_corner=p["size_corner"] * f,
                size_regions=lambda body: tuple((rx, s * f) for rx, s in _call(regions, body)),
                boxes=lambda body: [(a, b, s * f) for a, b, s in boxes(body)],
                envelopes=lambda body: [(rx, pad, s * f) for rx, pad, s in _call(env, body) or []])


def _call(v, body):
    """Profile entries may be fixed or functions of the body."""
    return v(body) if callable(v) else v


def per_segment(table, name, default):
    return next((v for rx, v in table if re.fullmatch(rx, name)), default)


def _turn(p, q, r):
    """Turning angle [deg] at q of the polyline p -> q -> r (positive = left, convex for a CCW contour)."""
    a, b = q - p, r - q
    return np.degrees(np.arctan2(a[0] * b[1] - a[1] * b[0], a @ b))


def _clean(seg, tol=5e-6):
    """Drop interior points closer than tol to the previous kept point or to the end point (near-duplicate
    vertices from the drawing give zero-length mesh edges, on which gmsh's 2D edge recovery fails)."""
    keep = [seg[0]]
    for q in seg[1:-1]:
        if np.linalg.norm(q - keep[-1]) > tol and np.linalg.norm(q - seg[-1]) > tol:
            keep.append(q)
    return np.array(keep + [seg[-1]])


def _densify(seg, ratio=2.0):
    """Split edges more than ratio times longer than a neighbouring edge into equal straight parts: the
    interpolating spline overshoots (and can loop) where the point spacing jumps."""
    d = np.linalg.norm(np.diff(seg, axis=0), axis=1)
    out = [seg[0]]
    for i, L in enumerate(d):
        nb = min(d[max(i - 1, 0)], d[min(i + 1, len(d) - 1)])
        k = max(1, int(np.ceil(L / (ratio * nb)))) if L > ratio * nb else 1
        out += [seg[i] + (seg[i + 1] - seg[i]) * j / k for j in range(1, k + 1)]
    return np.array(out)


def _round_concave(seg, corner_angle, radius, n_arc=6):
    """Replace sharp concave interior vertices (turning right on a CCW contour) by a circular arc of the given
    radius (cut back at most half of each adjacent edge): gmsh's 2D boundary layer fails at sharp concave
    corners. Convex corners stay sharp (they get a boundary-layer fan)."""
    out = [seg[0]]
    for i in range(1, len(seg) - 1):
        p, q, r = out[-1], seg[i], seg[i + 1]
        t = _turn(p, q, r)
        if t >= -corner_angle:
            out.append(q)
            continue
        u, v = (p - q) / np.linalg.norm(p - q), (r - q) / np.linalg.norm(r - q)
        half = np.radians(180 - abs(t)) / 2                     # half the interior angle
        d = min(radius / np.tan(half), 0.5 * np.linalg.norm(p - q), 0.5 * np.linalg.norm(r - q))
        a, b = q + d * u, q + d * v
        c = q + (u + v) / np.linalg.norm(u + v) * d / np.cos(half)  # arc centre
        ang = [np.arctan2(*(a - c)[::-1]), np.arctan2(*(b - c)[::-1])]
        if ang[1] - ang[0] > np.pi:
            ang[1] -= 2 * np.pi
        elif ang[0] - ang[1] > np.pi:
            ang[1] += 2 * np.pi
        rr = np.linalg.norm(a - c)
        out += [c + rr * np.array([np.cos(w), np.sin(w)]) for w in np.linspace(*ang, n_arc)]
    return np.array(out + [seg[-1]])


def _pieces(seg, corner_angle, concave_radius=0.0):
    """Split a segment polyline at its sharp interior vertices (concave ones rounded first)."""
    seg = _clean(seg)
    if concave_radius:
        seg = _round_concave(seg, corner_angle, concave_radius)
    seg = _densify(seg)
    cut = [0] + [i for i in range(1, len(seg) - 1) if abs(_turn(seg[i - 1], seg[i], seg[i + 1])) > corner_angle]
    return [seg[a:b + 1] for a, b in zip(cut, cut[1:] + [len(seg) - 1])]


def build(body, p, path):
    """Mesh the body in gmsh (2D) and write the mesh extruded to one cell as an msh 2.2 file."""
    gmsh.initialize(interruptible=False)    # no SIGINT handler: build() may run in a worker thread
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("section")
        geo = gmsh.model.geo
        lc = p["size_far"]

        # domain: the main surface sits on the ground layers, a structured strip of height T that write_extruded
        # builds from the nodes of the strip top (gmsh's own extrusion or boundary-layer field on the ground
        # fails at inlet/outlet or in the edge recovery)
        first, ratio, n_g = p["layer_first"], p["layer_expansion"], p["ground_layers"]
        ys = np.r_[0.0, np.cumsum(first * ratio ** np.arange(n_g))]
        T = ys[-1]
        corners = [geo.addPoint(x, y, 0, lc) for x, y in
                   ((p["x_min"], T), (p["x_max"], T), (p["x_max"], p["y_max"]), (p["x_min"], p["y_max"]))]
        strip_top = geo.addLine(corners[0], corners[1])
        dom = dict(outlet=[geo.addLine(corners[1], corners[2])], top=[geo.addLine(corners[2], corners[3])],
                   inlet=[geo.addLine(corners[3], corners[0])])
        loops = [geo.addCurveLoop([strip_top, dom["outlet"][0], dom["top"][0], dom["inlet"][0]])]

        # body contours: one spline per smooth piece, corners as shared points
        patch_curves = {}            # patch name -> [curve tags]
        names = {}                   # curve tag -> patch name
        sharp, convex = [], []       # point tags at corners; convex ones get a boundary-layer fan
        for loop in body.contours():
            segs = [body.segments[n] for n in loop]
            start = prev_pt = None
            curves = []
            for k, (name, seg) in enumerate(zip(loop, segs)):
                for piece in _pieces(seg, p["corner_angle"], p.get("concave_radius", 0.0)):
                    if prev_pt is None:
                        prev_pt = start = geo.addPoint(*piece[0], 0, lc)
                    last = k == len(segs) - 1 and np.allclose(piece[-1], segs[0][0])
                    end = start if last else geo.addPoint(*piece[-1], 0, lc)
                    inner = [geo.addPoint(x, y, 0, lc) for x, y in piece[1:-1]]
                    c = (geo.addSpline([prev_pt, *inner, end]) if inner else geo.addLine(prev_pt, end))
                    curves.append((c, piece))
                    patch_curves.setdefault(f"body_{name}", []).append(c)
                    names[c] = f"body_{name}"
                    prev_pt = end
            # corners between pieces (and segments)
            for (c0, a), (c1, b) in zip(curves, curves[1:] + curves[:1]):
                t = _turn(a[-2], a[-1], b[1])
                if abs(t) > p["corner_angle"]:
                    sharp.append((tuple(a[-1]), t > 0, names[c0]))
            loops.append(geo.addCurveLoop([c for c, _ in curves]))
        surf = geo.addPlaneSurface(loops)
        geo.synchronize()

        # point tags of the corners
        pts = {tuple(np.round(gmsh.model.getValue(0, t, [])[:2], 12)): t for _, t in gmsh.model.getEntities(0)}
        def tag(xy):
            return min(pts.items(), key=lambda kv: np.hypot(kv[0][0] - xy[0], kv[0][1] - xy[1]))[1]
        corner_tags = [tag(xy) for xy, _, _ in sharp]
        convex = [tag(xy) for xy, cv, _ in sharp if cv]

        # size field: min over contour distances, corners and boxes, each growing linearly with distance
        f = gmsh.model.mesh.field
        g = p["growth"]
        fields = []

        def threshold(dist_field, size):
            t = f.add("Threshold")
            f.setNumber(t, "InField", dist_field)
            f.setNumber(t, "SizeMin", size)
            f.setNumber(t, "SizeMax", lc)
            f.setNumber(t, "DistMin", 0.0)
            f.setNumber(t, "DistMax", (lc - size) / g)
            fields.append(t)

        regions = _call(p["size_regions"], body)
        surf_size = {name: per_segment(regions, name.removeprefix("body_"), p["size_surface"]) for name in patch_curves}
        by_size = {}
        for name, cs in patch_curves.items():
            by_size.setdefault(surf_size[name], []).extend(cs)
        for size, cs in by_size.items():
            d = f.add("Distance")
            f.setNumbers(d, "CurvesList", cs)
            f.setNumber(d, "Sampling", 400)
            threshold(d, size)
        by_corner = {}          # corners: size_corner, at most half the surface size next to them
        for (_, _, name), tg in zip(sharp, corner_tags):
            by_corner.setdefault(min(p["size_corner"], 0.5 * surf_size[name]), []).append(tg)
        for size, tags in by_corner.items():
            d = f.add("Distance")
            f.setNumbers(d, "PointsList", tags)
            threshold(d, size)
        for rx, dist, size in _call(p.get("envelopes"), body) or []:
            cs = [c for name, cl in patch_curves.items() if re.fullmatch(rx, name.removeprefix("body_")) for c in cl]
            if cs:
                d = f.add("Distance")
                f.setNumbers(d, "CurvesList", cs)
                f.setNumber(d, "Sampling", 400)
                t = f.add("Threshold")
                for k, v in dict(InField=d, SizeMin=size, SizeMax=lc, DistMin=dist, DistMax=dist + (lc - size) / g).items():
                    f.setNumber(t, k, v)
                fields.append(t)
        for (x0, y0), (x1, y1), size in p["boxes"](body):
            b = f.add("Box")
            for k, v in dict(VIn=size, VOut=lc, XMin=x0, XMax=x1, YMin=y0, YMax=y1, ZMin=-1, ZMax=1,
                             Thickness=(lc - size) / g).items():
                f.setNumber(b, k, v)
            fields.append(b)
        if p["size_ground"]:    # caps the aspect ratio of the ground layers in the far field
            b = f.add("Box")
            for k, v in dict(VIn=p["size_ground"], VOut=lc, XMin=p["x_min"] - 1, XMax=p["x_max"] + 1, YMin=-1,
                             YMax=T, ZMin=-1, ZMax=1, Thickness=(lc - p["size_ground"]) / g).items():
                f.setNumber(b, k, v)
            fields.append(b)
        m = f.add("Min")
        f.setNumbers(m, "FieldsList", fields)
        f.setAsBackgroundMesh(m)

        # body wall layers: one boundary-layer field per layer count
        by_n = {}
        for name, cs in patch_curves.items():
            by_n.setdefault(per_segment(p["region_layers"], name.removeprefix("body_"), p["body_layers"]), []).extend(cs)
        for n, cs in by_n.items():
            if not n:
                continue
            bl = f.add("BoundaryLayer")
            f.setNumbers(bl, "CurvesList", cs)
            f.setNumber(bl, "Size", first)
            f.setNumber(bl, "Ratio", ratio)
            f.setNumber(bl, "NbLayers", n)
            f.setNumber(bl, "Thickness", first * (ratio ** n - 1) / (ratio - 1))
            f.setNumber(bl, "Quads", 1)
            fan = [t for t in convex if any(t in gmsh.model.getAdjacencies(1, c)[1] for c in cs)]
            if fan:
                f.setNumbers(bl, "FanPointsList", fan)
            f.setAsBoundaryLayer(bl)

        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.BoundaryLayerFanElements", p["fan"])
        gmsh.option.setNumber("Mesh.Algorithm", 8)              # frontal-Delaunay for quads
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)  # blossom
        gmsh.option.setNumber("Mesh.Smoothing", 5)

        gmsh.model.mesh.generate(2)
        boundary = {name: cs for name, cs in list(dom.items()) + list(patch_curves.items())}
        write_extruded(path, surf, boundary, p["span"], (strip_top, ys))
    finally:
        gmsh.finalize()


def write_extruded(path, surf, boundary, dz, strip):
    """Add the ground layers below the 2D mesh of surf, extrude everything by dz to one cell layer and write it
    as msh 2.2 (ascii): quads -> hexahedra, triangles -> prisms, boundary edges of each named curve group ->
    quad faces, the 2D elements at z = 0 -> patch back, their copies at z = dz -> patch front.
    strip = (strip top curve, layer levels from y = 0 to its height): one column of quads below every node of
    the strip top, bottom edges -> ground, end columns -> inlet and outlet."""
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    pos = {int(t): x for t, x in zip(tags, xyz.reshape(-1, 3))}
    raw = []                                                    # node tags of the 2D elements at z = 0
    for s in [surf]:
        for etype, _, enodes in zip(*gmsh.model.mesh.getElements(2, s)):
            k = {2: 3, 3: 4}[int(etype)]                        # triangle, quad
            for row in enodes.reshape(-1, k):
                v = [int(a) for a in row]
                v = [a for i, a in enumerate(v) if a != v[i - 1]]   # boundary-layer fans: collapsed quads -> triangles
                raw.append(v)
    used = sorted({a for v in raw for a in v})                  # only nodes of elements (not geometry points)
    idx = {t: i + 1 for i, t in enumerate(used)}                # gmsh node tag -> 1-based node number
    xyz = np.array([pos[t] for t in used])
    n = len(used)
    cells = []
    for v in raw:                                               # counter-clockwise, so that the cells point to +z
        x, y = xyz[[idx[a] - 1 for a in v], 0], xyz[[idx[a] - 1 for a in v], 1]
        if np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
            v = v[::-1]
        cells.append((len(v), [idx[a] for a in v]))
    edges = {}
    for name, cs in boundary.items():
        for c in cs:
            for _, _, enodes in zip(*gmsh.model.mesh.getElements(1, c)):
                edges.setdefault(name, []).extend([idx[int(a)], idx[int(b)]] for a, b in enodes.reshape(-1, 2) if a != b)

    # ground layers: columns below the strip-top nodes, sorted by x
    top_tags, _, _ = gmsh.model.mesh.getNodes(1, strip[0], includeBoundary=True)
    top = sorted({idx[int(a)] for a in top_tags}, key=lambda i: xyz[i - 1, 0])
    ys = strip[1]
    new = np.array([[xyz[i - 1, 0], y, 0.0] for i in top for y in ys[:-1]])
    col = {i: [n + 1 + k * (len(ys) - 1) + j for j in range(len(ys) - 1)] + [i] for k, i in enumerate(top)}
    xyz = np.vstack([xyz, new])
    n = len(xyz)
    for a, b in zip(top, top[1:]):
        for j in range(len(ys) - 1):
            cells.append((4, [col[a][j], col[b][j], col[b][j + 1], col[a][j + 1]]))
        edges.setdefault("ground", []).append([col[a][0], col[b][0]])
    for name, i in (("inlet", top[0]), ("outlet", top[-1])):
        edges.setdefault(name, []).extend([col[i][j], col[i][j + 1]] for j in range(len(ys) - 1))

    # every edge used by exactly one 2D element must be a named boundary edge
    count = {}
    for _, v in cells:
        for a, b in zip(v, v[1:] + v[:1]):
            e = (min(a, b), max(a, b))
            count[e] = count.get(e, 0) + 1
    free = {e for e, m in count.items() if m == 1}
    named = {(min(a, b), max(a, b)) for es in edges.values() for a, b in es}
    if free != named:
        bad = sorted(free ^ named)[:5]
        raise RuntimeError(f"{len(free ^ named)} boundary edges not matched by the curve meshes, e.g. at "
                           + ", ".join(f"({xyz[a - 1][0]:.4f}, {xyz[a - 1][1]:.4f})" for a, _ in bad))

    names = ["internal", "back", "front", *edges]
    phys = {name: i + 1 for i, name in enumerate(names)}
    els = []
    for k, v in cells:
        top = [i + n for i in v]
        els.append((5 if k == 4 else 6, phys["internal"], v + top))       # hexahedron / prism
        els.append((3 if k == 4 else 2, phys["back"], v))
        els.append((3 if k == 4 else 2, phys["front"], top))
    for name, es in edges.items():
        els += [(3, phys[name], [a, b, b + n, a + n]) for a, b in es]
    with open(path, "w") as fh:
        fh.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n%d\n" % len(names))
        fh.write("".join(f'{3 if name == "internal" else 2} {phys[name]} "{name}"\n' for name in names))
        fh.write("$EndPhysicalNames\n$Nodes\n%d\n" % (2 * n))
        fh.write("".join(f"{i + 1} {x:.12g} {y:.12g} 0\n" for i, (x, y, _) in enumerate(xyz)))
        fh.write("".join(f"{i + 1 + n} {x:.12g} {y:.12g} {dz:.12g}\n" for i, (x, y, _) in enumerate(xyz)))
        fh.write("$EndNodes\n$Elements\n%d\n" % len(els))
        fh.write("".join(f"{i + 1} {t} 2 {g} {g} {' '.join(map(str, v))}\n" for i, (t, g, v) in enumerate(els)))
        fh.write("$EndElements\n")


BUILD_TIMEOUT = 600     # s; large multi-element sections took up to 210 s, the Tyrrell wing a few seconds


def _build_child(body, p, path, q):
    try:
        build(body, p, path)
        q.put(None)
    except Exception as e:  # noqa: BLE001
        q.put(f"{type(e).__name__}: {e}")


def build_isolated(body, p, path, timeout=BUILD_TIMEOUT):
    """build() in a forked child: gmsh can segfault (seen on aarch64) or not finish for some geometries, which
    would kill or block the calling process; both become an exception here."""
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    q = ctx.Queue()
    child = ctx.Process(target=_build_child, args=(body, p, path, q))
    child.start()
    child.join(timeout)
    if child.is_alive():
        child.kill()
        child.join()
        raise RuntimeError(f"gmsh did not finish within {timeout} s")
    if child.exitcode != 0:
        raise RuntimeError(f"gmsh crashed (exit code {child.exitcode})")
    try:
        err = q.get(timeout=10)
    except Exception:  # noqa: BLE001  queue.Empty
        err = "gmsh child returned no result"
    if err:
        raise RuntimeError(err)


def mesh_errors(case):
    """checkMesh failures that break the solver; the high aspect ratio of the ground layers is normal."""
    log = (Path(case) / "log.checkMesh").read_text()
    return [ln.strip() for ln in log.splitlines() if ln.strip().startswith("***") and "aspect ratio" not in ln]


# cell-size factors tried in turn when a profile has mesh_retry: about 6 % of a few hundred multi-element sections gave
# gmsh errors (edge not recovered, open cells, wrongly oriented faces, segfault on aarch64) that a 2 %
# change of every cell size removes; 2 % in cell size moved C_l by about 0.1 % in the Tyrrell checks
RETRY_SCALES = (1.0, 1.02, 0.98, 1.04)


def mesh_case(case, body, p=TYRRELL):
    """Same as mesh2d_snappy.mesh, with gmsh: complete case directory with a checked 2D mesh. With p["mesh_retry"]
    a failed mesh is rebuilt with slightly scaled cell sizes (RETRY_SCALES); the factor used is noted in
    log.gmsh. Returns that factor."""
    if not p.get("mesh_retry"):
        _mesh_once(case, body, p)
        return 1.0
    err = None
    for f in RETRY_SCALES:
        try:
            _mesh_once(case, body, p if f == 1.0 else scaled(p, f))
            bad = mesh_errors(case)
            if not bad:
                if f != 1.0:
                    with open(Path(case) / "log.gmsh", "a") as fh:
                        fh.write(f"cell sizes scaled by {f} (previous attempt: {err})\n")
                return f
            err = "checkMesh: " + "; ".join(bad)
        except RuntimeError as e:
            err = str(e)
    raise RuntimeError(f"no valid mesh at cell-size factors {RETRY_SCALES}: {err}")


def _mesh_once(case, body, p):
    case = Path(case)
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(TEMPLATES / "2d", case)
    write_case_dicts(case, body, p)
    t0 = time.time()
    build_isolated(body, p, case / "section.msh")
    (case / "log.gmsh").write_text(f"gmsh 2D + extrusion: {time.time() - t0:.1f} s\n")
    run(case, "gmshToFoam", "section.msh")
    types = {"front": "empty", "back": "empty", "ground": "wall", **{n: "wall" for n in
             re.findall(r"\n\s*(body_\w+)\s*\n", (case / "constant/polyMesh/boundary").read_text())}}
    set_patch_types(case / "constant/polyMesh/boundary", types)
    run(case, "renumberMesh", "-overwrite", "-constant")
    run(case, "checkMesh")
    (case / "section.msh").unlink()
    (case / f"{case.name}.foam").touch()


TYRRELL_DOUBLE = dict(TYRRELL, boxes=tyrrell_double_boxes)
