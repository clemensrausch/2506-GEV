"""Ahmed body (Ahmed et al. 1984; Lienhart & Becker 2003, ERCOFTAC case 82) as closed STL surfaces.

The nose is the original DLR surface (ERCOFTAC file ahmed-front-geo.dat: one quarter of the nose corner as a
structured 100 x 200 point grid, x = 0 ... 100 mm from the front). The rest of the nose follows from it: the
side and top/bottom strips are the grid's edge profiles extruded across the flat parts, the front face is flat.
Behind x = 100 mm the cross-section at x = 100 mm is swept to the rear; the slant (222 mm along the surface,
angle as parameter) lowers the roof and the rounded roof edges (r ~ 2.5 mm in the DLR data) together.
Stilts: 4 cylinders, diameter 30 mm, 202 mm and 672 mm behind the nose, 163.5 mm off the centre plane
(ERCOFTAC drawing); they reach below the ground and into the body, so snappyHexMesh cuts them at both ends.

CFD coordinates as in the 2D setup: x = flow (nose at -L, base at 0), y = up
(ground at y = 0), z = span (centre plane z = 0). Units m.

The nose file is not part of this repository (see SOURCES.md); put it into data/ahmed-front-geo.dat.

Usage (from the repository root):
    python cases/ahmed/ahmed_body.py OUTDIR [--slant 25] [--clearance 0.05] [--no-stilts]
writes OUTDIR/body.stl (solid "body") and OUTDIR/stilts.stl (solid "stilts") plus a check of the surface.
"""

import argparse
from pathlib import Path

import numpy as np
import trimesh

FRONT_DATA = Path(__file__).resolve().parent / "data/ahmed-front-geo.dat"
FRONT_URL = "http://cfd.mace.manchester.ac.uk/ercoftac/lib/exe/fetch.php?media=cdata:case082:ahmed-front-geo.dat"

L = 1044.0            # body length, mm
SLANT_LENGTH = 222.0  # length of the slanted surface, mm (same for all angles)
STILT_X = (202.0, 672.0)   # stilt centres behind the nose, mm
STILT_Z = 163.5            # stilt centres off the centre plane, mm
STILT_R = 15.0


def _grid_quads(n, m, offset=0):
    """Quads (as two triangles) of a structured n x m point grid stored row-major."""
    i, j = np.meshgrid(np.arange(n - 1), np.arange(m - 1), indexing="ij")
    a = (i * m + j).ravel() + offset
    b, c, d = a + 1, a + m + 1, a + m
    return np.vstack([np.c_[a, b, c], np.c_[a, c, d]])


class _Builder:
    def __init__(self):
        self.v, self.f = [], []
        self.n = 0

    def grid(self, P):
        """Add a structured grid P[n, m, 3]."""
        n, m, _ = P.shape
        self.v.append(P.reshape(-1, 3))
        self.f.append(_grid_quads(n, m, self.n))
        self.n += n * m

    def fan(self, loop, centre):
        """Add a planar polygon (ordered closed loop, star-shaped about centre) as a triangle fan."""
        k = len(loop)
        self.v.append(np.vstack([loop, centre]))
        idx = np.arange(k) + self.n
        self.f.append(np.c_[idx, np.roll(idx, -1), np.full(k, self.n + k)])
        self.n += k + 1

    def mesh(self):
        m = trimesh.Trimesh(np.vstack(self.v), np.vstack(self.f), process=False)
        m.merge_vertices(digits_vertex=6)          # patch borders share exact coordinates
        m.update_faces(m.nondegenerate_faces())    # the grid's first row is a single point
        m.remove_unreferenced_vertices()
        trimesh.repair.fix_normals(m)              # consistent winding, outward
        return m


def nose(path=FRONT_DATA, n_side=36, n_top=81):
    """Nose surface (data coordinates: x from the front, y span, z vertical about mid-height; mm), open at
    x = 100 mm. Returns the trimesh."""
    if not path.exists():
        raise FileNotFoundError(f"{path} missing; download it from {FRONT_URL}")
    Q = np.loadtxt(path, comments="#").reshape(100, 200, 3)   # rows: rings around (0, 80, 35); cols: z=35 .. y=80
    Q = Q[np.r_[0:99:2, 99]][:, np.r_[0:199:2, 199]]          # every 2nd point (<= 4 mm), edges kept
    nr = len(Q)
    side = Q[:, 0, [0, 1]]     # profile in the plane z = 35: (x, y)
    top = Q[:, -1, [0, 2]]     # profile in the plane y = 80: (x, z)
    y0, z0 = Q[0, 0, 1], Q[0, 0, 2]   # flat front face |y| <= 80, |z| <= 35

    b = _Builder()
    for sy in (1, -1):
        for sz in (1, -1):
            b.grid(Q * [1, sy, sz])
    zs = np.linspace(-z0, z0, n_side)
    for sy in (1, -1):         # side strips, extruded over |z| <= 35
        b.grid(np.stack([np.c_[side[:, 0], sy * side[:, 1], np.full(nr, z)] for z in zs], axis=1))
    ys = np.linspace(-y0, y0, n_top)
    for sz in (1, -1):         # top / bottom strips, extruded over |y| <= 80
        b.grid(np.stack([np.c_[top[:, 0], np.full(nr, y), sz * top[:, 1]] for y in ys], axis=1))
    face = np.vstack([np.c_[np.zeros(n_top), ys, np.full(n_top, -z0)],
                      np.c_[np.zeros(n_side), np.full(n_side, y0), zs][1:],
                      np.c_[np.zeros(n_top), ys[::-1], np.full(n_top, z0)][1:],
                      np.c_[np.zeros(n_side), np.full(n_side, -y0), zs[::-1]][1:-1]])
    b.fan(face, np.zeros(3))
    return b.mesh()


def _rear_loop(m, x=100.0):
    """Ordered boundary loop of the open nose at x (data coordinates)."""
    paths = m.outline().discrete
    loop = max(paths, key=len)[:-1]           # closed path repeats its first point
    assert np.allclose(loop[:, 0], x), "nose outline is not planar at x = 100 mm"
    return loop


def body(slant_deg=25.0, clearance=50.0, dx_prism=50.0, dx_slant=20.0, path=FRONT_DATA):
    """Closed body surface in CFD coordinates (m)."""
    m = nose(path)
    loop = _rear_loop(m)                      # (x=100, y, z)
    z_top = loop[:, 2].max()
    phi = np.radians(slant_deg)
    x_s = L - SLANT_LENGTH * np.cos(phi)      # slant starts here (sharp edge)

    # roof drop over the slant: points at or above z_hi (roof and its rounded edges) move down rigidly,
    # side-wall points between 0 and z_hi are compressed linearly (they stay in the side plane)
    z_hi = z_top - 6.0
    def ring(x):
        drop = max(0.0, (x - x_s) * np.tan(phi))
        y, z = loop[:, 1], loop[:, 2].copy()
        up = z >= z_hi
        mid = (z > 0) & ~up
        z[up] -= drop
        z[mid] *= (z_hi - drop) / z_hi
        r = np.c_[np.full(len(y), x), y, z]
        return np.vstack([r, r[:1]])          # closed ring

    xs = np.unique(np.r_[np.linspace(100.0, x_s, int(np.ceil((x_s - 100.0) / dx_prism)) + 1),
                         np.linspace(x_s, L, int(np.ceil((L - x_s) / dx_slant)) + 1)])
    b = _Builder()
    b.v.append(m.vertices); b.f.append(m.faces); b.n = len(m.vertices)
    b.grid(np.stack([ring(x) for x in xs]))
    base = ring(L)[:-1]
    b.fan(base, np.r_[L, 0.0, 0.5 * (base[:, 2].min() + base[:, 2].max())])
    s = b.mesh()

    # data -> CFD coordinates, m
    v = s.vertices
    h_half = -v[:, 2].min()
    s.vertices = np.c_[v[:, 0] - L, v[:, 2] + h_half + clearance, v[:, 1]] / 1000.0
    trimesh.repair.fix_normals(s)
    return s


def stilts(clearance=50.0, n=64, overlap=10.0):
    """Four closed cylinders (CFD coordinates, m) from overlap below the ground to overlap into the body."""
    parts = []
    for x in STILT_X:
        for z in (STILT_Z, -STILT_Z):
            c = trimesh.creation.cylinder(radius=STILT_R, height=clearance + 2 * overlap, sections=n)
            c.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0]))   # axis -> y
            c.apply_translation([x - L, 0.5 * clearance, z])
            parts.append(c)
    s = trimesh.util.concatenate(parts)
    s.vertices /= 1000.0
    return s


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
    ap.add_argument("--slant", type=float, default=25.0, help="slant angle, deg")
    ap.add_argument("--clearance", type=float, default=0.05, help="ground clearance, m")
    ap.add_argument("--no-stilts", action="store_true")
    a = ap.parse_args()
    b = body(a.slant, a.clearance * 1000.0)
    check(b, "body")
    write_stl(b, a.outdir / "body.stl", "body")
    if not a.no_stilts:
        s = stilts(a.clearance * 1000.0)
        check(s, "stilts")
        write_stl(s, a.outdir / "stilts.stl", "stilts")
    print(f"written to {a.outdir}")


if __name__ == "__main__":
    main()
