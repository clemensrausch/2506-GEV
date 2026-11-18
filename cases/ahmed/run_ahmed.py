"""Ahmed body (ERCOFTAC case 82) with the 3D template gev.mesh3d: geometry, mesh, solve.

The validated runs are ahmed25_wf_lust and ahmed35_wf_lust: wall functions on the body, LUST for U, solver
preset "fast", background cell 0.08 m (grid study: --cell 0.104 / 0.08 / 0.0615).

Usage (from the repository root, OpenFOAM environment sourced):
    python cases/ahmed/run_ahmed.py mesh  --slant 25 --name wf_lust [--body-layers wallfn] [--cell 0.08] [--procs 8]
        [--dry-run] [--solver fast|strict] [--slant-box] [--layer-feature-angle 180] [--turb-scheme upwind]
    python cases/ahmed/run_ahmed.py solve --slant 25 --name wf_lust [--max-iter 3000] [--solver fast|strict]
        [--uniform-start] [--u-scheme LUST|linearUpwindV] [--turb-model kOmegaSST|realizableKE]
Cases go to $GEV_RUNS/ahmed/ahmed<slant>[_NAME].
"""

import argparse
import json
import re

import _path  # noqa: F401
import ahmed_body
from gev import RUNS, mesh3d, run3d

CASES = RUNS / "ahmed"


TURB_SCHEMES = {"upwind": "turbulence      bounded Gauss upwind;",
                "limitedLinear": "turbulence      bounded Gauss limitedLinear 1;"}


U_SCHEMES = {"linearUpwindV": "bounded Gauss linearUpwindV grad(U);",
             "LUST": "bounded Gauss LUST grad(U);"}       # 75 % linear + 25 % linearUpwind, less damping


def set_numerics(case, u_scheme, turb_model):
    """Before solving: convection scheme of U and turbulence model of an already meshed case.
    realizableKE gets an epsilon field (epsilon = Cmu k omega, wall functions) next to omega."""
    fs = case / "system/fvSchemes"
    t = re.sub(r"(div\(phi,U\)\s+).*;", lambda m: m.group(1) + U_SCHEMES[u_scheme], fs.read_text())
    if turb_model == "realizableKE" and "div(phi,epsilon)" not in t:
        t = t.replace("    div(phi,omega)  $turbulence;\n", "    div(phi,omega)  $turbulence;\n    div(phi,epsilon) $turbulence;\n")
    fs.write_text(t)
    tp = case / "constant/turbulenceProperties"
    tp.write_text(re.sub(r"RASModel\s+\w+;", f"RASModel        {turb_model};", tp.read_text()))
    if turb_model != "realizableKE":
        return
    fsol = case / "system/fvSolution"
    fsol.write_text(fsol.read_text().replace('"(k|omega)"', '"(k|omega|epsilon)"'))
    ic = case / "0/initialConditions"
    t = ic.read_text()
    if "turbulentEpsilon" not in t:
        k = float(re.search(r"turbulentKE\s+([^;]+);", t).group(1))
        om = float(re.search(r"turbulentOmega\s+([^;]+);", t).group(1))
        ic.write_text(t.replace("groundVelocity", f"turbulentEpsilon     {0.09 * k * om:.6g};\ngroundVelocity", 1))
    eps = (case / "0/omega").read_text()
    eps = eps.replace("object      omega;", "object      epsilon;").replace("[0 0 -1 0 0 0 0]", "[0 2 -3 0 0 0 0]")
    eps = eps.replace("$turbulentOmega", "$turbulentEpsilon")
    eps = re.sub(r"type\s+omegaWallFunction;\s*blending\s+\w+;\s*n\s+\d+;", "type            epsilonWallFunction;", eps)
    eps = re.sub(r"\n *// y\+ ~ 1 on some parts.*", "", eps)      # comment on the omega blending
    assert "omega" not in eps, "0/epsilon still mentions omega"
    (case / "0/epsilon").write_text(eps)


def slant_boxes(slant):
    """extra_boxes for mesh3d: level 4 (5 mm) around the slant, its side edges and the first 0.3 m of wake
    (shear layer of the roof edge and C-pillar vortices), level 3 (10 mm) in the near wake up to 1 L behind
    the base (the reference wake box is level 2 = 20 mm there)."""
    import numpy as np
    dx = 0.222 * np.cos(np.radians(slant))

    def f(bb, Lb):
        (lo, hi), za = bb, -0.0008
        zb = max(abs(lo[2]), abs(hi[2]))
        return {"slantBox": ((hi[0] - dx - 0.04, 0.10, za), (hi[0] + 0.30, hi[1] + 0.06, zb + 0.06), 4),
                "nearWakeBox": ((hi[0], -0.015, za), (hi[0] + 1.0 * Lb, hi[1] + 0.20, zb + 0.20), 3)}
    return f


def setup(case, slant, profile, dry_run=False, body_layers="wallfn"):
    geo = case.parent / f"{case.name}_geometry"
    b = ahmed_body.body(slant)
    ahmed_body.check(b, "body")
    ahmed_body.write_stl(b, geo / "body.stl", "body")
    ahmed_body.write_stl(ahmed_body.stilts(), geo / "stilts.stl", "stilts")
    parts = [mesh3d.Part("body", geo / "body.stl", layers=body_layers),
             mesh3d.Part("stilts", geo / "stilts.stl", layers="wallfn")]
    return mesh3d.mesh(case, parts, profile, dry_run=dry_run)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["mesh", "solve"])
    ap.add_argument("--slant", type=float, default=25.0)
    ap.add_argument("--name", default="")
    ap.add_argument("--solver", choices=list(mesh3d.SOLVERS), default="fast")
    ap.add_argument("--dry-run", action="store_true", help="write and check the dictionaries only")
    ap.add_argument("--max-iter", type=int, default=3000)
    ap.add_argument("--body-layers", choices=["resolved", "wallfn"], default="wallfn",
                    help="body wall: wall-function layers (y+ ~ 30-60, validated) or y+ ~ 1 layers")
    ap.add_argument("--procs", type=int, default=mesh3d.AHMED["n_procs"], help="MPI processes (physical cores)")
    ap.add_argument("--layer-feature-angle", type=float, default=mesh3d.AHMED["layer_feature_angle"],
                    help="snappy addLayers featureAngle: layers stop at edges sharper than this")
    ap.add_argument("--turb-scheme", choices=list(TURB_SCHEMES), default="upwind",
                    help="convection scheme of k and omega (template: first-order upwind as motorBike)")
    ap.add_argument("--cell", type=float, default=mesh3d.AHMED["cell"],
                    help="mesh: background cube size, m (all refinement levels scale with it; grid study 0.104/0.08/0.0615)")
    ap.add_argument("--slant-box", action="store_true",
                    help="mesh: refine slant, side edges and near wake (level 4 box, level 3 wake to 1 L)")
    ap.add_argument("--u-scheme", choices=list(U_SCHEMES), default="LUST",
                    help="solve: convection scheme of U (template: LUST; linearUpwindV as motorBike and 2D)")
    ap.add_argument("--turb-model", choices=["kOmegaSST", "realizableKE"], default="kOmegaSST",
                    help="solve: RAS model")
    ap.add_argument("--uniform-start", action="store_true", help="solve: start from uniform flow, no potentialFoam")
    a = ap.parse_args()
    case = CASES / (f"ahmed{a.slant:g}" + (f"_{a.name}" if a.name else ""))
    profile = dict(mesh3d.AHMED, solver=a.solver, n_procs=a.procs, layer_feature_angle=a.layer_feature_angle, cell=a.cell)
    if a.slant_box:
        profile["extra_boxes"] = slant_boxes(a.slant)
    if a.step == "mesh":
        info = setup(case, a.slant, profile, a.dry_run, a.body_layers)
        fs = case / "system/fvSchemes"
        fs.write_text(fs.read_text().replace(TURB_SCHEMES["upwind"], TURB_SCHEMES[a.turb_scheme]))
        print(f"{case}: background {info['cells_background']} cells, first layer {info['first_layer'] * 1e6:.1f} um")
        return
    if a.solver != "fast" or not (case / "system/fvSolution").read_text().count('preset "fast"'):
        mesh3d.write_fvsolution(case / "system", profile)       # solver preset of this run
    set_numerics(case, a.u_scheme, a.turb_model)
    r = run3d.solve(case, profile, a.max_iter, potential=not a.uniform_start)
    # body without stilts, coefficient on the full frontal area (half model: half force on half area)
    q_area = 0.5 * 1.225 * profile["U_inf"] ** 2 * profile["ref_area"] * 0.5
    for part in ("body", "stilts"):
        f = run3d.part_forces(case, part, r["window"], r["iterations"])
        r[f"Cd_{part}"], r[f"Cl_{part}"] = f[0] / q_area, f[1] / q_area
    r.update(solver=a.solver, u_scheme=a.u_scheme, turb_model=a.turb_model)
    (case / "result.json").write_text(json.dumps(r, indent=1, default=float))
    print(json.dumps(r, indent=1, default=float))


if __name__ == "__main__":
    main()
