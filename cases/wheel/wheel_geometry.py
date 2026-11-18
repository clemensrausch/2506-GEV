"""Isolated rotating wheel of Patel, Garmory & Passmore (2022), half-scale wind-tunnel version.

Source: D. Patel, A. Garmory, M.A. Passmore, On the wake of an isolated rotating wheel: an experimental and
numerical investigation, J. Wind Eng. Ind. Aerodyn. (2022), CC BY 4.0; CAD and PIV data:
https://doi.org/10.17028/rd.lboro.16863232 (CC BY-NC 4.0). The CAD file is not needed: the geometry is built
from the dimensions below.

The dimensions below were read from the Parasolid file Wheel_WTVersion(withWT).x_t (analytic surfaces and edge
circles; no Parasolid reader is available here):
- tyre: tread R 165 mm (D 330), width 180 mm; shoulder of two R 60 mm arcs (centres r 105 / a 51.02 mm and
  r 120.31 / a 30 mm) meeting at r 160 / a 75 mm with a small crease; flat sidewall at a = 90 mm down to r 99.5 mm;
- rim: R 10 mm round into a lip at r 89.5 mm, R 2 mm fillet into the hub face, recessed to a = 76 mm;
- axle centre 164.5 mm above the floor: the tyre reaches 0.53 mm into the floor (contact patch);
- shaft R 12.5 mm to a = 350 mm; on each side a support plate (70 x 20 mm, a = 267.5 ... 287.5 mm, from the floor
  to 93.6 mm above the axle) and a bearing housing (box 62 x 111 mm to a = 299.5 mm, round R 33 mm to 315.5 mm,
  collar R 17 mm to 325.5 mm). Bolt heads, hub bolts and housing details are left out.
- contact patch: the tyre is cut by a plinth (Diasinos et al. 2015) of height --plinth (default 1 mm; Patel
  et al. used 0.1 mm with much finer cells), whose walls stand vertically on the floor;
- the wheel is driven by a drum (R 108.2 mm, a = +-127.5 mm) in a slot of the fixed floor, flush with the floor.

Coordinates as for the other 3D cases: x = flow, y = up (floor at y = 0), z = span (wheel centre plane z = 0, axle
along z), wheel axis at x = 0. Units m. Half model: only z >= 0.

Usage (from the repository root):
    python cases/wheel/wheel_geometry.py OUTDIR [--full]
writes OUTDIR/wheel.stl (tyre, rim and shaft: rotating wall) and OUTDIR/supports.stl (fixed).
"""

import argparse
from pathlib import Path

import numpy as np
import trimesh

R = 0.165                 # tread radius
AXLE_H = 0.16447          # axle height above the floor
HALF_W = 0.090            # half width of the tyre
U_TUNNEL = 30.0           # freestream and tread speed, m/s


def _arc(c, rad, t0, t1, n):
    """Points on a circle in the (r, a) plane, centre c, angles t0 -> t1 in degrees (excluding the end point)."""
    t = np.radians(np.linspace(t0, t1, n, endpoint=False))
    return np.c_[c[0] + rad * np.cos(t), c[1] + rad * np.sin(t)]


def _shoulder_crease():
    """Intersection of the two shoulder arcs (both R 60 mm), outer one."""
    c1, c2, rad = np.array([0.12031, 0.030]), np.array([0.105, 0.051021]), 0.060
    m, d = (c1 + c2) / 2, c2 - c1
    h = np.sqrt(rad ** 2 - (np.linalg.norm(d) / 2) ** 2)
    p = m + h * np.array([-d[1], d[0]]) / np.linalg.norm(d)
    return p if p[0] > m[0] else m - h * np.array([-d[1], d[0]]) / np.linalg.norm(d)


def profile(n=24):
    """Half profile (r, a) for a >= 0, from the axis on the hub face to the tread centre (a = 0)."""
    pc = _shoulder_crease()
    ang = lambda c, p: np.degrees(np.arctan2(p[1] - c[1], p[0] - c[0]))
    pts = [np.array([[0.0, 0.076]]),
           _arc((0.0875, 0.078), 0.002, -90, 0, 6),         # R 2 fillet, hub face -> lip
           _arc((0.0995, 0.080), 0.010, 180, 90, 10),       # R 10 round, lip -> sidewall
           np.array([[0.0995, 0.090]]),
           _arc((0.12031, 0.030), 0.060, 90, ang((0.12031, 0.030), pc), n),   # outer shoulder arc
           _arc((0.105, 0.051021), 0.060, ang((0.105, 0.051021), pc), 0, n),  # inner shoulder arc
           np.array([[R, 0.051021], [R, 0.0]])]
    return np.vstack(pts)


def r_outer(a):
    """Outer tyre radius at the axial position a (tread and the two shoulder arcs)."""
    a = np.abs(np.asarray(a, float))
    pc = _shoulder_crease()
    rA = 0.105 + np.sqrt(np.clip(0.060 ** 2 - (a - 0.051021) ** 2, 0, None))
    rB = 0.12031 + np.sqrt(np.clip(0.060 ** 2 - (a - 0.030) ** 2, 0, None))
    return np.where(a <= 0.051021, R, np.where(a <= pc[1], rA, rB))


def plinth(height, n=200):
    """Contact-patch plinth (Diasinos et al. 2015, as Patel et al.): the footprint of the tyre at y = height,
    extruded from below the floor up into the tyre, so that the tyre meets the floor with vertical walls instead
    of a sharp wedge. Returns the plinth mesh and its footprint (x, z) polygon."""
    import shapely.geometry as sg
    dy = AXLE_H - height
    a_max = brentq_r(dy)
    a = np.linspace(-a_max, a_max, n)
    x = np.sqrt(np.clip(r_outer(a) ** 2 - dy ** 2, 0, None))
    poly = sg.Polygon(np.r_[np.c_[x, a], np.c_[-x, a][::-1][1:-1]])   # x = 0 at both ends: no duplicates
    m = trimesh.creation.extrude_polygon(poly, height + 0.0035)   # local (x, a) plane, extruded along +local z
    # local (x, a, e) -> world (x, y, z) = (x, y0 + e, a)
    m.apply_transform(np.array([[1, 0, 0, 0], [0, 0, 1, -0.003], [0, 1, 0, 0], [0, 0, 0, 1.0]]))
    if m.volume < 0:                                              # the axis swap mirrors
        m.invert()
    return m, poly


def brentq_r(dy):
    """Axial half width of the footprint at distance dy below the axle (r_outer(a) = dy)."""
    from scipy.optimize import brentq
    return brentq(lambda a: float(r_outer(a)) - dy, 0.0, 0.0899)


def wheel(full=False, sections=180, plinth_h=0.001):
    """Tyre and rim (body of revolution about z) plus shaft and contact-patch plinth, as one closed mesh."""
    half = profile()
    prof = np.vstack([half, half[-2::-1] * [1, -1]])               # mirror to a < 0
    tyre = trimesh.creation.revolve(prof, sections=sections)
    if tyre.volume < 0:
        tyre.invert()                                              # profile runs clockwise in (r, a)
    shaft = trimesh.creation.cylinder(radius=0.0125, segment=[[0, 0, -0.35], [0, 0, 0.35]], sections=48)
    m = trimesh.boolean.union([tyre, shaft], engine="manifold")
    m.apply_transform(trimesh.transformations.translation_matrix([0, AXLE_H, 0]))
    if plinth_h:
        m = trimesh.boolean.union([m, plinth(plinth_h)[0]], engine="manifold")
    return m if full else _half(m)


def supports(full=False):
    """Support plate and bearing housing on each side (z > 0 side mirrored)."""
    parts = []
    for s in (1, -1):
        plate = trimesh.creation.box(bounds=[[-0.035, 0.0, 0.2675], [0.035, AXLE_H + 0.0936, 0.2875]])
        house = trimesh.creation.box(bounds=[[-0.031, AXLE_H - 0.0558, 0.2875], [0.031, AXLE_H + 0.0548, 0.2995]])
        ring = trimesh.creation.cylinder(radius=0.033, segment=[[0, AXLE_H, 0.2995], [0, AXLE_H, 0.3155]], sections=64)
        collar = trimesh.creation.cylinder(radius=0.017, segment=[[0, AXLE_H, 0.3155], [0, AXLE_H, 0.3255]], sections=48)
        m = trimesh.boolean.union([plate, house, ring, collar], engine="manifold")
        if s < 0:
            m.apply_transform(trimesh.transformations.reflection_matrix([0, 0, 0], [0, 0, 1]))
        parts.append(m)
    m = trimesh.util.concatenate(parts)
    return m if full else _half(m)


def _half(m):
    """Keep z >= 0, cut slightly below 0 so snappy sees the surface reach the symmetry plane."""
    return m.slice_plane([0, 0, -1e-4], [0, 0, 1], cap=False)


def write_stl(m, path, solid):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = trimesh.exchange.stl.export_stl_ascii(m).replace("solid ", f"solid {solid}", 1)
    text = text.rsplit("endsolid", 1)[0] + f"endsolid {solid}\n"
    path.write_text(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--full", action="store_true", help="whole wheel instead of the half model")
    ap.add_argument("--plinth", type=float, default=0.001, help="contact-patch plinth height, m (0: none)")
    a = ap.parse_args()
    w, s = wheel(a.full, plinth_h=a.plinth), supports(a.full)
    full_w = wheel(True, plinth_h=a.plinth)
    if a.plinth:
        _, poly = plinth(a.plinth)
        print(f"plinth {a.plinth * 1e3:g} mm: footprint {np.ptp(poly.bounds[::2]) * 1e3:.1f} x {np.ptp(poly.bounds[1::2]) * 1e3:.1f} mm")
    print(f"wheel: {len(w.faces)} faces, closed full wheel {full_w.is_watertight}, volume {full_w.volume * 1e3:.2f} dm3")
    print(f"  bounds {w.bounds.round(4).tolist()}")
    print(f"  lowest point y = {full_w.bounds[0, 1] * 1000:.2f} mm (floor at 0), omega = {U_TUNNEL / R:.1f} rad/s")
    print(f"supports: {len(s.faces)} faces, bounds {s.bounds.round(4).tolist()}")
    write_stl(w, a.outdir / "wheel.stl", "wheel")
    write_stl(s, a.outdir / "supports.stl", "supports")


if __name__ == "__main__":
    main()
