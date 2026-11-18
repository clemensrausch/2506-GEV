"""Diffuser-equipped bluff body of Senior (2002) and Ruhrmann & Zhang (2003), University of Southampton.

Dimensions from Senior's thesis (ch. 2.2, Fig. 4): length l = 1315 mm, width 2d = 314 mm, height H = 326 mm,
diffuser ramp 17 deg, height of the bluff base 159 mm, diffuser inlet faired at R 78.5 mm, nose profile
R 131.6 mm at the front blending into R 741.5 mm, upper long edges rounded R 61 mm (front view), lower long
edges sharp, side-plates along the diffuser down to the level of the flat underbody.

Construction: the solid is the intersection of three extrusions whose outlines are the views of Fig. 4:
side profile (x-y, extruded across the span), plan profile (x-z, extruded vertically) and front section
(y-z, extruded along x). The side-plates are added as plates of thickness PLATE_T, flush with the body sides.

Assumptions (not stated numerically in the sources):
- plan view of the nose: same radii as the side view (R 131.6 / R 741.5); Fig. 5(b) shows a similar shape;
- the ramp starts where 17 deg and the 159 mm base meet: 546 mm ahead of the base (Senior quotes 0.538 m as
  the diffuser length, 538 / cos 17 deg = 563 or 538 tan 17 deg = 164.5 mm do not close with 159 mm);
- side-plate thickness PLATE_T = 5 mm, no fairing; support struts (on top of the model) are not modelled.

CFD coordinates as for the Ahmed body: x = flow (nose at -l, base at 0), y = up (ground at y = 0, model
underbody at y = ride height), z = span (centre plane z = 0). Units m.

Usage (from the repository root):
    python cases/diffuser/diffuser_body.py OUTDIR [--hr 0.1]
writes OUTDIR/body.stl (solid "body", wall functions) and OUTDIR/underbody.stl (solid "underbody", y+ ~ 1),
two open surfaces that together close the model, after a check of the closed surface.
"""

import argparse
from pathlib import Path

import numpy as np
import shapely.geometry as sg
import trimesh

L = 1315.0          # model length, mm
D = 157.0           # half width d, mm
H = 326.0           # height, mm
RAMP_DEG = 17.0
BASE_H = 159.0      # height of the bluff base above the diffuser exit, mm
R_INLET = 78.5      # fairing radius at the diffuser inlet, mm
R_NOSE1, R_NOSE2 = 131.6, 741.5
R_TOP = 61.0        # upper long edges (front view)
PLATE_T = 5.0       # side-plate thickness, mm (assumed)
SEG = 2.0           # arc discretisation, mm

RAMP_RISE = H - BASE_H                                  # 167 mm
X_RAMP = L - RAMP_RISE / np.tan(np.radians(RAMP_DEG))   # sharp ramp corner, mm behind the nose tip


def _arc(c, r, a0, a1):
    n = max(4, int(np.ceil(abs(a1 - a0) * r / SEG)))
    a = np.linspace(a0, a1, n)
    return np.c_[c[0] + r * np.cos(a), c[1] + r * np.sin(a)]


def nose_curve(half):
    """Upper half of a nose profile of half-thickness `half` about axis 0: points from the tip (0, 0) to the
    tangent point with the flat side (x_t, half), and x_t."""
    r1, r2 = R_NOSE1, R_NOSE2
    c1 = np.array([r1, 0.0])
    yc2 = half - r2
    xc2 = r1 + np.sqrt((r2 - r1) ** 2 - yc2 ** 2)       # big arc tangent to the small one (inside) and the side
    c2 = np.array([xc2, yc2])
    u = (c1 - c2) / np.linalg.norm(c1 - c2)
    tp = c2 + r2 * u                                      # tangent point of the two arcs
    a_t = np.arctan2(tp[1] - c1[1], tp[0] - c1[0])
    p1 = _arc(c1, r1, np.pi, a_t)
    p2 = _arc(c2, r2, np.arctan2(u[1], u[0]), np.pi / 2)
    return np.vstack([p1, p2[1:]]), xc2


def side_profile():
    """Closed x-y outline (mm, nose tip at x = 0, underbody at y = 0), counter-clockwise."""
    up, _ = nose_curve(H / 2)
    top = up + [0, H / 2]
    # diffuser inlet fairing: convex corner between the flat bottom and the 17 deg ramp
    th = np.radians(RAMP_DEG)
    t = R_INLET * np.tan(th / 2)
    c = np.array([X_RAMP - t, R_INLET])
    fair = _arc(c, R_INLET, -np.pi / 2, -np.pi / 2 + th)
    lower = up * [1, -1] + [0, H / 2]                     # tip -> start of the flat bottom
    outline = np.vstack([lower, fair, [[L, RAMP_RISE], [L, H]], top[::-1][:-1]])
    return outline


def plan_profile():
    """Closed x-z outline (mm), counter-clockwise."""
    up, _ = nose_curve(D)
    right = up * [1, -1]                                   # tip -> z = -d
    left = up                                              # tip -> z = +d
    return np.vstack([right, [[L, -D], [L, D]], left[::-1][:-1]])


def front_section():
    """Closed y-z outline (mm): rectangle, upper corners rounded R_TOP, lower corners sharp. Points as (z, y)."""
    a = _arc((D - R_TOP, H - R_TOP), R_TOP, 0, np.pi / 2)
    b = _arc((-D + R_TOP, H - R_TOP), R_TOP, np.pi / 2, np.pi)
    return np.vstack([[[-D, 0.0], [D, 0.0]], a, b])


def _extrude(outline, lo, hi):
    poly = sg.Polygon(outline)
    if not poly.is_valid:
        poly = poly.buffer(0)
    m = trimesh.creation.extrude_polygon(poly, hi - lo)
    m.apply_translation([0, 0, lo])
    return m


def body(hr=100.0, plates=True):
    """Closed surface of the model with side-plates, mm -> m in CFD coordinates, ride height hr (mm)."""
    pad = 50.0
    # side view: outline in (x, y), extruded along z
    side = _extrude(side_profile(), -D - pad, D + pad)
    # plan view: outline in (x, z), extruded along y: (x, z, y) -> (x, y, z)
    plan = _extrude(plan_profile(), -pad, H + pad)
    plan.vertices = plan.vertices[:, [0, 2, 1]]
    # front view: outline in (z, y), extruded along x: (z, y, x) -> (x, y, z)
    front = _extrude(front_section(), -pad, L + pad)
    front.vertices = front.vertices[:, [2, 1, 0]]
    for m in (plan, front):
        if m.volume < 0:
            m.invert()
    solid = trimesh.boolean.intersection([side, plan, front], engine="manifold")
    if plates:
        x0 = X_RAMP - R_INLET * np.tan(np.radians(RAMP_DEG) / 2) - 1.0   # start of the fairing
        tri = np.array([[x0, -1.0], [L, -1.0], [L, H / 2], [x0, H / 4]])  # reaches into the body (union)
        tri[:, 1] = np.maximum(tri[:, 1], 0.0)
        for zs in (D - PLATE_T, -D):
            p = _extrude(tri, zs, zs + PLATE_T)
            solid = trimesh.boolean.union([solid, p], engine="manifold")
    solid.apply_translation([-L, hr, 0.0])
    solid.apply_scale(1e-3)
    return solid


def split(m, hr):
    """(upper, underbody) surfaces of the model built by body(hr) (m). Underbody: faces looking down (flat
    bottom, inlet fairing, ramp, lower nose; n_y < -0.3) and the side-plate faces below the ramp surface
    (inner faces and the lower part of the outer faces). Both are open; together they close the model."""
    hr_mm = hr * 1000.0
    n, c = m.face_normals, m.triangles_center * 1000.0 + [L, -hr_mm, 0.0]      # back to model mm coordinates
    vy = m.triangles[:, :, 1] * 1000.0 - hr_mm
    vx = m.triangles[:, :, 0] * 1000.0 + L
    ramp = np.maximum(vx - X_RAMP, 0.0) * np.tan(np.radians(RAMP_DEG))        # ramp height at each vertex
    x0 = X_RAMP - R_INLET * np.tan(np.radians(RAMP_DEG) / 2) - 1.0
    plate = (np.abs(n[:, 2]) > 0.9) & np.all(vy <= ramp + 0.5, axis=1) & (c[:, 0] > x0) & \
            (np.abs(c[:, 2]) > D - PLATE_T - 0.5)
    under = (n[:, 1] < -0.3) | plate
    return m.submesh([np.where(~under)[0]], append=True), m.submesh([np.where(under)[0]], append=True)


def frontal_area():
    """Projected area in the flow direction, m^2: the front section (side-plates lie inside the 2d width)."""
    return sg.Polygon(front_section()).area * 1e-6


def write_stl(mesh, path, solid):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = trimesh.exchange.stl.export_stl_ascii(mesh)
    path.write_text(text.replace("solid ", f"solid {solid}", 1).rsplit("endsolid", 1)[0] + f"endsolid {solid}\n")


def check(mesh, name):
    lo, hi = mesh.bounds
    ok = mesh.is_watertight and mesh.is_winding_consistent and mesh.volume > 0
    print(f"{name}: {len(mesh.faces)} triangles, watertight {mesh.is_watertight}, "
          f"winding {mesh.is_winding_consistent}, volume {mesh.volume * 1e3:.3f} dm3")
    print(f"  bounds x {lo[0]:.4f} .. {hi[0]:.4f}, y {lo[1]:.4f} .. {hi[1]:.4f}, z {lo[2]:.4f} .. {hi[2]:.4f} m")
    if not ok:
        raise RuntimeError(f"{name}: surface is not closed and outward oriented")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--hr", type=float, default=0.1, help="ride height of the flat underbody, m")
    a = ap.parse_args()
    m = body(a.hr * 1000.0)
    check(m, "body")
    print(f"  ramp corner {X_RAMP:.1f} mm behind the tip ({L - X_RAMP:.1f} mm ahead of the base), "
          f"frontal area {frontal_area():.5f} m2")
    upper, under = split(m, a.hr)
    print(f"  split: body {len(upper.faces)} faces {upper.area * 1e4:.1f} cm2, "
          f"underbody {len(under.faces)} faces {under.area * 1e4:.1f} cm2")
    write_stl(upper, a.outdir / "body.stl", "body")
    write_stl(under, a.outdir / "underbody.stl", "underbody")


if __name__ == "__main__":
    main()
