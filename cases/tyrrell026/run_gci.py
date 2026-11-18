"""Grid convergence study (Celik et al. 2008) for the Tyrrell 026 single element wing.

Three systematically refined meshes per height, refinement ratio r = 1.4. gmsh: every cell size
(mesh2d_gmsh.scaled) is scaled by r; snappy: the background cells are scaled (all refinement levels and boxes
unchanged, so every non-layer cell scales with r). The wall layers (first cell 10 um, 20 layers, ratio 1.2) are
kept fixed so that the wall treatment (y+ ~ 1) is the same on all meshes. The medium mesh is the one of the
h/c sweep (run_sweep.py).

Usage (from the repository root, OpenFOAM environment sourced):
    python cases/tyrrell026/run_gci.py [--eval-only] [--mesher gmsh|snappy] [--heights 0.179 0.134 0.090]
Writes results/gci_cl<tag>.csv.
"""

import argparse
from pathlib import Path


import _path  # noqa: F401
import run_sweep
from gev import mesh2d_gmsh, mesh2d_snappy
from gev.foam import n_cells
from run_sweep import MESHERS

HERE = Path(__file__).resolve().parent
R = 1.4
HEIGHTS = [0.179, 0.134, 0.090]
LEVELS = {"_coarse": R, "": 1.0, "_fine": 1 / R}   # suffix -> cell size factor


def params(f, mesher="gmsh"):
    if mesher == "gmsh":
        return mesh2d_gmsh.scaled(mesh2d_gmsh.TYRRELL, f)
    t = mesh2d_snappy.TYRRELL
    return dict(t, nx=round(t["nx"] / f), ny=round(t["ny"] / f), span=t["span"] * f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eval-only", action="store_true")
    ap.add_argument("--mesher", choices=("gmsh", "snappy"), default="gmsh")
    ap.add_argument("--heights", type=float, nargs="+", default=HEIGHTS)
    a = ap.parse_args()
    fn, _, tag = MESHERS[a.mesher]
    levels = {tag + s: f for s, f in LEVELS.items()}
    for h in a.heights:
        for suf, f in levels.items():
            case = run_sweep.case_dir(h, suf)
            if a.eval_only or (suf == tag and (case / "postProcessing").exists()):
                continue
            print(f"h/c = {h}{suf}: cell size factor {f:.3f} ...", flush=True)
            run_sweep.run(h, params(f, a.mesher), suf, fn)
            d = run_sweep.result(h, suffix=suf)
            print(f"h/c = {h}{suf}: {n_cells(case)} cells, Cl = {d['cl']:.4f}, {d['it']} it, "
                  f"converged {d['converged']}", flush=True)

    rows = ["h_over_c,mesh,cells,cl,cl_min,cl_max,cd,iterations,converged"]
    for h in sorted(a.heights, reverse=True):
        for suf in levels:
            case = run_sweep.case_dir(h, suf)
            if not (case / "postProcessing").exists():
                continue
            d = run_sweep.result(h, suffix=suf)
            rows.append(f"{h},{suf.removeprefix(tag).strip('_') or 'medium'},{n_cells(case)},{d['cl']:.5f},"
                        f"{d['cl_min']:.5f},{d['cl_max']:.5f},{d['cd']:.5f},{d['it']},{d['converged']}")
    out = HERE / f"results/gci_cl{tag}.csv"
    out.write_text("\n".join(rows) + "\n")
    print("\n".join(rows))
    print("written:", out)


if __name__ == "__main__":
    main()
