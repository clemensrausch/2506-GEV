"""Tyrrell 026 double element wing (Zerihan 2001, chapter 8), datum flap deflection: ground-height sweep with the
gmsh mesh (mesh2d_gmsh.TYRRELL_DOUBLE), and an incidence sweep at the highest tested height.

Usage (from the repository root, OpenFOAM environment sourced):
    python cases/tyrrell026/run_double.py [h/c ...]          (default: HEIGHTS)
    python cases/tyrrell026/run_double.py --plot-only
    python cases/tyrrell026/run_double.py --alpha -5 -8 -11  (incidence sweep at h/c 0.592)

Same solver settings and averaging as run_sweep.py. C_L, C_D on the combined chord c = 380 mm; main element and
flap separately from the force groups. Writes results/double_cl_gmsh.csv and results/double_alpha_gmsh.csv
(CFD only); with data/zerihan_fig36b_cl_datum.csv present (see SOURCES.md) also results/double_cl_gmsh.png.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import _path  # noqa: F401
import wing
from gev import mesh2d_gmsh, mesh2d_snappy
from gev.foam import read_dat
from gev.tables import read_csv
from run_sweep import CASES, N_PROCS, solve

HERE = Path(__file__).resolve().parent
C = wing.CHORD_DOUBLE
HEIGHTS = [0.592, 0.394, 0.289, 0.237, 0.211, 0.159, 0.120, 0.098, 0.080, 0.067, 0.054]   # Fig. 36b markers
EXPT = HERE / "data/zerihan_fig36b_cl_datum.csv"
ALPHA_DATUM = 1.0           # incidence of the tested wing (Tables 3/4), deg
H_FREE = 0.592              # highest tested height


def case_dir(h, alpha=ALPHA_DATUM):
    a = "" if alpha == ALPHA_DATUM else f"_a{alpha:+.1f}".replace(".", "").replace("+", "p").replace("-", "m")
    return CASES / (f"tyrrell2_h{h:.3f}".replace(".", "") + a + "_gmsh")


def run(h, alpha=ALPHA_DATUM, n_procs=N_PROCS):
    case = case_dir(h, alpha=alpha)
    mesh2d_gmsh.mesh_case(case, wing.double(h_over_c=h, alpha=alpha), dict(mesh2d_gmsh.TYRRELL_DOUBLE, n_procs=n_procs))
    solve(case, n_procs)


def result(h, alpha=ALPHA_DATUM):
    """C_L, C_D (downforce positive) and the main / flap share; converged runs averaged over the last 200
    iterations, runs that reached 4000 over the last 2000."""
    case = case_dir(h, alpha)
    f = read_dat(case / "postProcessing/forceCoeffs_main/0/coefficient.dat")
    conv = "SIMPLE solution converged" in (case / "log.simpleFoam").read_text()
    win = 200 if conv else 2000
    sel = f[:, 0] > f[-1, 0] - win
    q = 0.5 * 1.225 * mesh2d_snappy.TYRRELL["U_inf"] ** 2 * C * mesh2d_snappy.TYRRELL["span"]
    out = dict(h=h, it=int(f[-1, 0]), converged=conv, cl=-f[sel, 4].mean(), cd=f[sel, 1].mean(),
               cl_min=-f[sel, 4].max(), cl_max=-f[sel, 4].min())
    for g in ("main", "flap"):
        F = read_dat(case / f"postProcessing/forces_{g}/0/force.dat")
        out[f"cl_{g}"] = -F[F[:, 0] > F[-1, 0] - win, 2].mean() / q
    return out


def plot():
    res = [result(h) for h in HEIGHTS if (case_dir(h) / "postProcessing").exists()]
    lines = ["h_over_c,cl,cl_min,cl_max,cl_main,cl_flap,cd,iterations,converged"]
    for r in res:
        lines.append(f"{r['h']},{r['cl']:.5f},{r['cl_min']:.5f},{r['cl_max']:.5f},{r['cl_main']:.5f},{r['cl_flap']:.5f},"
                     f"{r['cd']:.5f},{r['it']},{r['converged']}")
        print(lines[-1])
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results/double_cl_gmsh.csv").write_text("\n".join(lines) + "\n")
    if not EXPT.exists():
        return
    ex = read_csv(EXPT)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(ex["h_over_c"], ex["CL"], "s-", color="k", ms=4, label="Experiment (Zerihan Fig. 36b, datum flap)")
    if res:
        h = [r["h"] for r in res]
        ax.errorbar(h, [r["cl"] for r in res], yerr=[[r["cl"] - r["cl_min"] for r in res], [r["cl_max"] - r["cl"] for r in res]],
                    fmt="o-", color="C0", label="2D RANS (gmsh)")
        ax.plot(h, [r["cl_main"] for r in res], "--", color="C1", label="  main element")
        ax.plot(h, [r["cl_flap"] for r in res], ":", color="C2", label="  flap")
    ax.set_xlabel("h/c"); ax.set_ylabel("C_L (downforce)"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(HERE / "results/double_cl_gmsh.png", dpi=130)


def alpha_sweep(alphas):
    """Incidence sweep at H_FREE (2D lift curve of the wing far from the ground), including the datum run of the
    height sweep -> results/double_alpha_gmsh.csv."""
    for a in alphas:
        if not (case_dir(H_FREE, alpha=a) / "postProcessing").exists():
            run(H_FREE, a)
    lines = ["alpha,cl,cd,iterations,converged"]
    for a in sorted({*alphas, ALPHA_DATUM}):
        r = result(H_FREE, alpha=a)
        lines.append(f"{a},{r['cl']:.5f},{r['cd']:.5f},{r['it']},{r['converged']}")
        print(lines[-1], flush=True)
    (HERE / "results/double_alpha_gmsh.csv").write_text(f"# h/c {H_FREE}, gmsh, datum flap\n" + "\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("heights", nargs="*", type=float)
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--alpha", nargs="+", type=float, help=f"incidence sweep at h/c {H_FREE} (deg, datum {ALPHA_DATUM})")
    ap.add_argument("--procs", type=int, default=N_PROCS)
    a = ap.parse_args()
    if a.alpha:
        return alpha_sweep(a.alpha)
    if not a.plot_only:
        for h in a.heights or HEIGHTS:
            run(h, n_procs=a.procs)
            r = result(h)
            print(f"h/c {h}: C_L {r['cl']:.4f} (main {r['cl_main']:.3f}, flap {r['cl_flap']:.3f}), C_D {r['cd']:.4f}, "
                  f"{r['it']} it, converged {r['converged']}", flush=True)
    plot()


if __name__ == "__main__":
    main()
