"""Diffuser bluff body (Senior 2002, 17 deg) with the 3D template casetools.mesh3d: geometry, mesh, solve, compare.

One case per ride height h_r/d (d = 0.157 m). Flow as in Senior's force tests: 20 m/s (Re_l 1.75e6), moving
ground at freestream speed, half model (--full: whole model, for the asymmetric stall of type C). Underbody
(flat bottom, ramp, inner side-plates) with y+ ~ 1 layers, the rest with wall functions. Forces refer to the
frontal area 0.1008 m^2; C_L keeps the template sign (negative = downforce), as Senior publishes it.

Refinement on top of the template boxes:
- gap box under the model: cells of at most h_r / 8 (level 3 ... 6), so the gap has 8 or more cells;
- diffuser box: level 4 (5 mm) from 0.1 m ahead of the ramp to 0.4 m behind the base, up to the base top.

Usage (from the repository root, OpenFOAM environment sourced):
    python cases/diffuser/run_diffuser.py mesh  --hd 0.210 [--name NAME] [--procs 8] [--dry-run]
        [--underbody-layers resolved|wallfn] [--cell 0.08] [--full]
    python cases/diffuser/run_diffuser.py solve --hd 0.210 [--name NAME] [--procs 8] [--max-iter 3000] [--full]
    python cases/diffuser/run_diffuser.py plot
Cases go to $GEV_RUNS/diffuser/diffuser_hd<0210>[_NAME]. With data/senior2002_fig9.csv present (see
SOURCES.md) "plot" also draws the measured curves.
"""

import argparse
import json
from pathlib import Path

import numpy as np

import _path  # noqa: F401
import diffuser_body
from casetools import RUNS, mesh3d, run3d
from casetools.tables import read_csv

HERE = Path(__file__).resolve().parent
CASES = RUNS / "diffuser"
D = diffuser_body.D / 1000.0
DATA = HERE / "data/senior2002_fig9.csv"

# underbody surface cells: level 5 (2.5 mm) everywhere. With level (4, 5) the ramp had 5 mm cells and a 25-layer
# stack of 8.7 mm, which collapsed to 2-8 layers within 5 cm of the side-plates (test mesh 2026-10-19).
UNDERBODY_LEVEL = (5, 5)

# layers in one pass (nOuterIter 1): with the template's 3 passes the pass of the thinnest layers failed the
# quality checks everywhere and the first cell was 66 um (y+ ~ 3-4) instead of 18 um (test meshes 2026-10-19);
# in one pass flat bottom and ramp get all 21 layers, the stack thins out within 5 cm of the side-plates.
# Domain: cross-section of the Southampton 2.1 m x 1.5 m tunnel in which Senior measured the forces (blockage
# 2.86 %, not corrected in the data), side walls and roof as slip walls. The first sweep (2026-10-19) used the open
# template domain; its forces were 2-16 % low, partly from the missing blockage.
TUNNEL = (1.5, 1.05)
DIFFUSER = dict(mesh3d.DEFAULTS, U_inf=20.0, Tu=0.002, moving_ground=True, ref_area=diffuser_body.frontal_area(),
                layer_outer_iter=1, tunnel=TUNNEL)


def boxes(hr, cell, half=True):
    """extra_boxes for mesh3d: gap box and diffuser box (see module docstring)."""
    def f(bb, Lb):
        lo, hi = bb
        gap_level = int(np.clip(np.ceil(np.log2(8 * cell / hr)), 3, 6))
        zb = max(abs(lo[2]), abs(hi[2]))
        za = lambda zw: -0.0008 if half else -zw
        x_ramp = hi[0] - (diffuser_body.L - diffuser_body.X_RAMP) / 1000.0
        return {"gapBox": ((lo[0] - 0.02, -0.015, za(zb + 0.02)), (hi[0] + 0.05, hr + 0.004, zb + 0.02), gap_level),
                "diffuserBox": ((x_ramp - 0.1, -0.015, za(zb + 0.03)), (hi[0] + 0.4, hr + 0.18, zb + 0.03), 4)}
    return f


def case_dir(hd, name):
    return CASES / (f"diffuser_hd{round(hd * 1000):04d}" + (f"_{name}" if name else ""))


def setup(case, hd, profile, dry_run=False, underbody_layers="resolved"):
    hr = hd * D
    geo = case.parent / f"{case.name}_geometry"
    m = diffuser_body.body(hr * 1000.0)
    diffuser_body.check(m, "diffuser body")
    upper, under = diffuser_body.split(m, hr)
    diffuser_body.write_stl(upper, geo / "body.stl", "body")
    diffuser_body.write_stl(under, geo / "underbody.stl", "underbody")
    parts = [mesh3d.Part("body", geo / "body.stl", layers="wallfn"),
             mesh3d.Part("underbody", geo / "underbody.stl", layers=underbody_layers, level=UNDERBODY_LEVEL)]
    return mesh3d.mesh(case, parts, profile, dry_run=dry_run)


def plot():
    """C_L and C_D over h_r/d of all solved cases -> results/sweep.csv (CFD only) and results/sweep.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = []
    for r in sorted(CASES.glob("diffuser_hd[0-9][0-9][0-9][0-9]/result.json")):   # not the grid-study cases
        d = json.loads(r.read_text())
        rows.append((d["hd"], d["Cl"], d["Cd"], d["status"], r.parent.name))
    out = HERE / "results"
    out.mkdir(parents=True, exist_ok=True)
    (out / "sweep.csv").write_text("hr_d,CL,CD,status,case\n" + "".join(
        f"{h:.3f},{cl:.4f},{cd:.4f},{s},{n}\n" for h, cl, cd, s, n in rows))
    e = None
    if DATA.exists():
        e = read_csv(DATA)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for a, key, col in ((ax[0], "CL", 1), (ax[1], "CD", 2)):
        if e is not None:
            a.plot(e["hr_d"], e[key], "ks-", ms=4, lw=0.8, label="Senior (2002), Fig. 9")
        if rows:
            a.plot([r[0] for r in rows], [r[col] for r in rows], "o-", color="C0", label="CFD (SST, LUST)")
        a.set_xlabel("h_r / d"); a.set_ylabel(key); a.grid(alpha=0.3); a.legend()
    fig.tight_layout()
    fig.savefig(out / "sweep.png", dpi=110)
    print(f"{len(rows)} cases -> {out / 'sweep.png'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["mesh", "solve", "plot"])
    ap.add_argument("--hd", type=float, help="ride height h_r / d (d = 0.157 m)")
    ap.add_argument("--name", default="")
    ap.add_argument("--procs", type=int, default=DIFFUSER["n_procs"])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--max-iter", type=int, default=3000)
    ap.add_argument("--cell", type=float, default=DIFFUSER["cell"], help="background cube size, m")
    ap.add_argument("--layer-outer-iter", type=int, default=DIFFUSER["layer_outer_iter"],
                    help="mesh: snappy nOuterIter (layers added in this many passes)")
    ap.add_argument("--underbody-layers", choices=["resolved", "wallfn"], default="resolved")
    ap.add_argument("--full", action="store_true",
                    help="full model without symmetry plane (asymmetric diffuser stall, Senior type C)")
    a = ap.parse_args()
    if a.step == "plot":
        return plot()
    case = case_dir(a.hd, a.name)
    profile = dict(DIFFUSER, n_procs=a.procs, cell=a.cell, layer_outer_iter=a.layer_outer_iter, half=not a.full,
                   extra_boxes=boxes(a.hd * D, a.cell, half=not a.full))
    if a.step == "mesh":
        info = setup(case, a.hd, profile, a.dry_run, a.underbody_layers)
        print(f"{case}: background {info['cells_background']} cells, first layer {info['first_layer'] * 1e6:.1f} um")
        return
    r = run3d.solve(case, profile, a.max_iter)
    r.update(hd=a.hd, hr=a.hd * D, half=not a.full)
    (case / "result.json").write_text(json.dumps(r, indent=1, default=float))
    print(json.dumps(r, indent=1, default=float))


if __name__ == "__main__":
    main()
