"""Tyrrell 026 single element wing: ground-height sweep (Zerihan 2001, Table 15), mesh + simpleFoam + Cl(h/c).

Usage (from the repository root, OpenFOAM environment sourced):
    python cases/tyrrell026/run_sweep.py [h/c ...] [--mesher gmsh|snappy|gmsh_wf] [--procs 8]
    python cases/tyrrell026/run_sweep.py --plot-only [--mesher gmsh]

Solver: residualControl p 1e-5, U/k/omega 1e-6 (10x tighter than the Ahmed reference), at most 4000
iterations. Converged runs are averaged over the last 200 iterations, runs that end periodic over the
last 2000.
Writes results/sweep_cl<tag>.csv (CFD only). With data/zerihan_table15_cl.csv present (see SOURCES.md) it also
writes results/sweep_cl<tag>.png with the experiment and Zerihan's own 2D CFD.
"""

import argparse
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import _path  # noqa: F401
import wing
from gev import RUNS, mesh2d_gmsh, mesh2d_snappy
from gev.foam import sh
from gev.tables import read_csv

HERE = Path(__file__).resolve().parent
CASES = RUNS / "tyrrell026"
HEIGHTS = [0.671, 0.448, 0.313, 0.224, 0.179, 0.134, 0.090]
N_PROCS = mesh2d_snappy.TYRRELL["n_procs"]
# mesher -> (mesh function, profile, case and result suffix)
MESHERS = {"snappy": (mesh2d_snappy.mesh, mesh2d_snappy.TYRRELL, ""),
           "gmsh": (mesh2d_gmsh.mesh_case, mesh2d_gmsh.TYRRELL, "_gmsh"),
           # wall functions instead of y+ ~ 1: 3 layers from 0.6 mm (y+ ~ 25-50 at the cell centre) on wing and ground
           "gmsh_wf": (mesh2d_gmsh.mesh_case, dict(mesh2d_gmsh.TYRRELL, layer_first=6e-4, layer_expansion=1.2,
                                                   body_layers=3, ground_layers=3), "_gmsh_wf")}
EXPT = HERE / "data/zerihan_table15_cl.csv"


def case_dir(h, suffix=""):
    return CASES / (f"tyrrell_h{h:.3f}".replace(".", "") + suffix)


def solve(case, n_procs, max_iter=4000):
    """simpleFoam to the tightened residuals or max_iter, then reconstruct the latest time."""
    for f, v in [("p", "1e-5"), ("U", "1e-6"), ('"(k|omega)"', "1e-6")]:
        sh(case, f"foamDictionary system/fvSolution -entry 'SIMPLE/residualControl/{f}' -set {v}", "foamDictionary")
    sh(case, f"foamDictionary system/controlDict -entry endTime -set {max_iter}", "foamDictionary")
    sh(case, "decomposePar -force", "decomposePar.run")
    sh(case, f"mpirun -np {n_procs} simpleFoam -parallel", "simpleFoam")
    sh(case, "reconstructPar -latestTime", "reconstructPar")
    for d in case.glob("processor*"):
        shutil.rmtree(d)


def run(h, p=mesh2d_snappy.TYRRELL, suffix="", mesh_fn=mesh2d_snappy.mesh):
    case = case_dir(h, suffix)
    mesh_fn(case, wing.single(h_over_c=h, alpha=1.0), p)
    solve(case, p["n_procs"])


def result(h, avg=None, suffix=""):
    """Mean Cl (downforce positive) and Cd: converged runs over the last 200 iterations, runs that end periodic
    at the iteration limit over the last 2000; min and max of Cl over that window."""
    case = case_dir(h, suffix)
    rows = [l.split() for l in (case / "postProcessing/forceCoeffs_main/0/coefficient.dat").read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    f = np.array([list(map(float, r[:8])) for r in rows])
    it, cd, cl = f[:, 0], f[:, 1], -f[:, 4]
    converged = "SIMPLE solution converged" in (case / "log.simpleFoam").read_text()
    avg = avg or (200 if converged else 2000)
    last = it > it[-1] - avg
    drift = cl[-1] - cl[np.searchsorted(it, it[-1] - avg)]
    return dict(h=h, it=int(it[-1]), cl=cl[last].mean(), cl_min=cl[last].min(), cl_max=cl[last].max(),
                drift=drift, cd=cd[last].mean(), converged=converged, window=avg)


def write_results(res, tag=""):
    r = sorted(res, key=lambda d: d["h"])
    lines = ["h_over_c,Cl,Cl_min,Cl_max,Cd,iterations,converged,avg_window,Cl_drift_over_window"]
    for d in r:
        lines.append(f"{d['h']},{d['cl']:.5f},{d['cl_min']:.5f},{d['cl_max']:.5f},{d['cd']:.5f},{d['it']},"
                     f"{d['converged']},{d['window']},{d['drift']:.6f}")
    out = HERE / f"results/sweep_cl{tag}.csv"
    out.parent.mkdir(exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    print("written:", out)
    if not EXPT.exists():
        return
    ref = read_csv(EXPT)
    fin = np.isfinite(ref["h_over_c"])
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(ref["h_over_c"][fin], ref["Cl_expt"][fin], "ko-", label="Zerihan (2001) experiment")
    ax.plot(ref["h_over_c"][fin], ref["Cl_kw"][fin], "C3s--", mfc="none", label="Zerihan 2D CFD, k-ω")
    ax.plot(ref["h_over_c"][fin], ref["Cl_SA"][fin], "C2^:", mfc="none", label="Zerihan 2D CFD, S-A")
    ax.plot([d["h"] for d in r], [d["cl"] for d in r], "C0*-", ms=12, label="this setup (simpleFoam, kOmegaSST)")
    ax.set_xlabel("h/c"); ax.set_ylabel("sectional Cl (downforce +)")
    ax.set_title("Tyrrell 026, α = 1° (true 3.45°), Re = 4.62e5, moving ground")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / f"results/sweep_cl{tag}.png", dpi=140)
    print(f"{'h/c':>6} {'Cl':>7} {'expt':>6} {'dev':>7}")
    for d in r:
        e = ref["Cl_expt"][np.isclose(ref["h_over_c"], d["h"])][0]
        print(f"{d['h']:6.3f} {d['cl']:7.4f} {e:6.3f} {100 * (d['cl'] / e - 1):+6.1f}%")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("heights", type=float, nargs="*", default=HEIGHTS)
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--mesher", choices=list(MESHERS), default="gmsh")
    ap.add_argument("--procs", type=int, default=N_PROCS, help="MPI processes (physical cores)")
    a = ap.parse_args()
    fn, prof, suf = MESHERS[a.mesher]
    prof = dict(prof, n_procs=a.procs)
    if not a.plot_only:
        for h in a.heights:
            print(f"h/c = {h}: meshing + solving ...", flush=True)
            run(h, prof, suf, fn)
            d = result(h, suffix=suf)
            print(f"h/c = {h}: Cl = {d['cl']:.4f} after {d['it']} it (converged {d['converged']})", flush=True)
    write_results([result(h, suffix=suf) for h in HEIGHTS if (case_dir(h, suf) / "postProcessing").exists()], suf)


if __name__ == "__main__":
    main()
