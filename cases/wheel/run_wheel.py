"""Isolated rotating wheel (Patel, Garmory & Passmore 2022) with the 3D template casetools.mesh3d: mesh, solve.

Flow as in the Loughborough Large Wind Tunnel test: 30 m/s (Re_D 0.67e6), turbulence intensity 0.2 %, fixed floor,
working section 1.92 m wide and 1.32 m high (blockage 3.3 %). The tunnel's side walls and roof are modelled as
slip walls at the same distance; the floor is a fixed no-slip wall that starts 3.65 m ahead of the axle, so that
a flat-plate boundary layer reaches about the measured displacement thickness of 7.2 mm at the wheel. The wheel,
rim and shaft are a rotating wall (tread speed = freestream); the supports are fixed. Patel et al. also left the
drum in the floor as a fixed wall in their validation simulation. Half model. Wall functions on the wheel, the
supports and the floor.

Mesh: the template profile scaled from the Ahmed body (cell 80 mm for L = 1.044 m) to the wheel diameter
(cell 25.4 mm, which also divides the tunnel height), plus near, wake, wheel-wake and contact boxes (boxes()).

Usage (from the repository root, OpenFOAM environment sourced):
    python cases/wheel/run_wheel.py mesh  [--name NAME] [--procs 8] [--dry-run] [--plinth 0.001]
    python cases/wheel/run_wheel.py solve [--name NAME] [--procs 8] [--max-iter 3000]
Cases go to $GEV_RUNS/wheel/wheel[_NAME].
"""

import argparse
import json

import numpy as np

import _path  # noqa: F401
import wheel_geometry as wg
from casetools import RUNS, mesh3d, run3d

CASES = RUNS / "wheel"
D = 2 * wg.R
TUNNEL_H, TUNNEL_HALF_W = 1.32, 0.96
FLOOR_AHEAD = 3.65            # floor length ahead of the axle, m (delta* ~ 7.2 mm at the wheel)
CELL = 0.0254                 # 0.08 m x D / 1.044 m, rounded so that 52 cells give the tunnel height

# frontal area of the whole wheel (tyre and rim, without shaft): 2 x integral over the height of the axial width
_r = np.linspace(0, wg.R, 4001)
_a = np.linspace(0, wg.HALF_W, 4001)
_ro = np.maximum(wg.r_outer(_a), 0)
FRONTAL_AREA = 2 * 2 * np.trapezoid([_a[_ro >= r].max() if (_ro >= r).any() else 0 for r in _r], _r)


def profile(cell=CELL, n_procs=8, half_span=0.35):
    """Template profile for the wheel; domain distances are in units of the wheel diameter (bbox length)."""
    return dict(mesh3d.DEFAULTS, U_inf=wg.U_TUNNEL, Tu=0.002, moving_ground=False, cell=cell, n_procs=n_procs,
                upstream=(FLOOR_AHEAD - wg.R) / D, downstream=9.6,
                roof=(TUNNEL_H - (wg.AXLE_H + wg.R)) / D, side=(TUNNEL_HALF_W - half_span) / D,
                ref_area=FRONTAL_AREA, ref_length=D, cofr=(0.0, wg.AXLE_H, 0.0),
                extra_boxes=boxes, layer_outer_iter=3)     # as validated (2026-10-19)


def boxes(bb, Lb):
    """Boxes on top of / instead of the template's (L = wheel diameter):
    - nearBox: level 3 around the wheel, up to 0.12 m (0.36 D) above the tyre for the top separation and its
      shear layer, but only 0.2 m out (the template box, scaled to the span with the shaft, reached 0.43 m and
      the first mesh, 6.2M cells, ran out of the 15 GB of the local machine during layer addition);
    - wakeBox: level 2 as the template, but from the nearBox start, so that no level-1 pocket is left above
      the wheel (first mesh: 13-25 mm cells there);
    - wheelWakeBox: level 3 from 0.25 m ahead of the axle (0.085 m ahead of the tyre) to 0.65 m behind it
      (PIV planes x/d = 1 and 1.5), up to 0.42 m;
    - contactBox: level 6 (0.4 mm) around the contact-patch plinth."""
    lo, hi = bb
    x0 = lo[0] - 0.533 * Lb
    return {"nearBox": ((x0, -0.015, -0.0008), (hi[0] + 0.192 * Lb, wg.AXLE_H + wg.R + 0.12, 0.20), 3),
            "wakeBox": ((x0, -0.015, -0.0008), (hi[0] + 2.586 * Lb, hi[1] + 0.634 * Lb, 0.526), 2),
            "wheelWakeBox": ((-0.25, -0.015, -0.0008), (0.65, 0.42, 0.25), 3),
            "contactBox": ((-0.045, -0.015, -0.0008), (0.045, 0.006, 0.075), 6)}


def setup(case, p, dry_run=False, plinth=0.001):
    geo = case.parent / f"{case.name}_geometry"
    wheel, sup = wg.wheel(not p["half"], plinth_h=plinth), wg.supports(not p["half"])
    wg.write_stl(wheel, geo / "wheel.stl", "wheel")
    wg.write_stl(sup, geo / "supports.stl", "supports")
    omega = wg.U_TUNNEL / wg.R
    rot = (f"type rotatingWallVelocity; origin (0 {wg.AXLE_H:g} 0); axis (0 0 1); omega {omega:.4g};"
           "   // tread speed = freestream, bottom moving downstream")
    parts = [mesh3d.Part("wheel", geo / "wheel.stl", level=(4, 5), feature_level=5, layers="wallfn", bc={"U": rot}),
             mesh3d.Part("supports", geo / "supports.stl", level=(3, 4), layers="wallfn")]
    return mesh3d.mesh(case, parts, p, dry_run=dry_run)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["mesh", "solve"])
    ap.add_argument("--name", default="")
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--max-iter", type=int, default=3000)
    ap.add_argument("--cell", type=float, default=CELL)
    ap.add_argument("--plinth", type=float, default=0.001, help="contact-patch plinth height, m")
    a = ap.parse_args()
    case = CASES / ("wheel" + (f"_{a.name}" if a.name else ""))
    p = profile(a.cell, a.procs)
    if a.step == "mesh":
        info = setup(case, p, a.dry_run, a.plinth)
        print(f"{case}: background {info['cells_background']} cells, frontal area {FRONTAL_AREA:.5f} m2")
        return
    r = run3d.solve(case, p, a.max_iter)
    w = run3d.part_forces(case, "wheel", r["window"], r["iterations"])
    q = 0.5 * 1.225 * wg.U_TUNNEL ** 2 * FRONTAL_AREA * 0.5          # half model
    r.update(Cd_wheel=w[0] / q, Cl_wheel=w[1] / q, plinth=a.plinth)
    (case / "result.json").write_text(json.dumps(r, indent=1, default=float))
    print(json.dumps(r, indent=1, default=float))


if __name__ == "__main__":
    main()
