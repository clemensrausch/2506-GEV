"""Mean wake of the solved wheel case vs. the PIV planes of Patel, Garmory & Passmore (2022).

1. Samples UMean of the decomposed case on the four PIV planes with OpenFOAM (postProcess, function object
   "surfaces", cutting planes, raw format) into CASE/postProcessing/pivPlanes/<time>/.
2. Maps the PIV data into the CFD frame and interpolates the CFD field linearly onto the PIV points.
3. Writes per plane the compared points (CSV) and an overview figure; prints the metrics.

PIV frames (from the data files; the centre planes have no documented origin, it was fitted to the masked wheel
outline, +-3 mm):
- centre planes "Y0mm" (z = 0) and "Y50mm" (z = 50 mm): in-plane x points upstream, y up; wheel centre at
  (292.1, -120.1) and (296.1, -123.1) mm; Vx = -U_x, Vy = U_y;
- cross planes "X330mm" / "X495mm" (x = 330 / 495 mm behind the axle, x/d = 1 / 1.5): in-plane x = span,
  y up with 0 at axle height; Vz = U_x, Vx = U_z, Vy = U_y. Both halves are compared to the half model
  (mirrored, z -> |z|).
The PIV mean is quoted with +-2 % uncertainty (99 % confidence).

The PIV files are not part of this repository (see SOURCES.md); put the *_Mean.dat files into data/.

Usage (from the repository root, OpenFOAM environment sourced):
    python cases/wheel/compare_piv.py [CASE] [--time 3000] [--procs 16]
writes CASE/piv_compare.png, CASE/piv_<plane>.csv, CASE/piv_metrics.csv.
"""

import argparse
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

import _path  # noqa: F401
from casetools import RUNS
from casetools.foam import write

HERE = Path(__file__).resolve().parent
import wheel_geometry as wg  # noqa: E402

PIV = HERE / "data"
U_INF = wg.U_TUNNEL
D = 2 * wg.R

# name: (PIV file, cutting plane point, normal)
PLANES = {
    "y0": ("Y0mm_Mean", (0, 0, 0.0005), (0, 0, 1)),        # just off the symmetry plane
    "y50": ("Y50mm_Mean", (0, 0, 0.050), (0, 0, 1)),
    "x330": ("X330mm_Mean", (0.330, 0, 0), (1, 0, 0)),
    "x495": ("X495mm_Mean", (0.495, 0, 0), (1, 0, 0)),
}
CENTRE = {"y0": (292.1, -120.1), "y50": (296.1, -123.1)}   # wheel centre in the PIV frame, mm


def sample(case, time, procs):
    """Cutting planes of UMean, written by OpenFOAM in raw format."""
    out = case / "postProcessing/pivPlanes" / str(time)
    if all((out / f"UMean_{n}.raw").exists() for n in PLANES):
        return out
    surfaces = "".join(f"""        {n}
        {{
            type        cuttingPlane;
            point       ({p[0]:g} {p[1]:g} {p[2]:g});
            normal      ({v[0]:g} {v[1]:g} {v[2]:g});
            interpolate true;
        }}
""" for n, (_, p, v) in PLANES.items())
    write(case / "system/pivPlanes", "pivPlanes", f"""\
functions
{{
readMean
{{
    type            readFields;
    libs            (fieldFunctionObjects);
    fields          (UMean);
}}
pivPlanes
{{
    type            surfaces;
    libs            (sampling);
    surfaceFormat   raw;
    fields          (UMean);
    interpolationScheme cellPoint;
    surfaces
    {{
{surfaces}    }}
}}
}}
""")
    cmd = ["mpirun", "--oversubscribe", "-np", str(procs), "postProcess", "-parallel", "-dict", "system/pivPlanes",
           "-time", str(time)]
    with open(case / "log.pivPlanes", "w") as fh:
        subprocess.run(cmd, cwd=case, stdout=fh, stderr=subprocess.STDOUT, check=True)
    return out


def load_cfd(out, name):
    f = out / f"UMean_{name}.raw"
    return np.loadtxt(f, comments="#")          # x y z Ux Uy Uz


def load_piv(name):
    lines = (PIV / f"{PLANES[name][0]}.dat").read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l[:1] in "-0123456789")
    return np.array([[float(x) for x in l.replace(",", " ").split()] for l in lines[i:] if l.strip()])


def piv_in_cfd(name):
    """PIV points and velocities in the CFD frame: dict(a, b plane coords in m, ux, uv (in-plane second
    velocity component), valid mask); a/b = (x, y) for the centre planes, (z, y) for the cross planes."""
    d = load_piv(name)
    if name.startswith("y"):
        xc, yc = CENTRE[name]
        a = (xc - d[:, 0]) / 1000.0                    # downstream of the axle
        b = (d[:, 1] - yc) / 1000.0 + wg.AXLE_H
        ux, uv = -d[:, 2], d[:, 3]
        valid = d[:, 4] > 0
    else:
        a = np.abs(d[:, 0]) / 1000.0
        b = d[:, 1] / 1000.0 + wg.AXLE_H
        ux, uv = d[:, 4], d[:, 3]
        valid = (np.abs(d[:, 2]) + np.abs(d[:, 3]) + np.abs(d[:, 4])) > 0
    return dict(a=a, b=b, ux=ux, uv=uv, valid=valid)


def compare(case, out):
    rows, fig, axs = [], None, None
    fig, axs = plt.subplots(len(PLANES), 3, figsize=(15, 3.6 * len(PLANES)))
    metrics = ["plane,n,mean_dU,rms_dU,share_within_5pct,wake_area_cfd,wake_area_piv"]
    for k, name in enumerate(PLANES):
        c = load_cfd(out, name)
        p = piv_in_cfd(name)
        if name.startswith("y"):
            ca, cb = c[:, 0], c[:, 1]
        else:
            ca, cb = c[:, 2], c[:, 1]
        m = p["valid"] & (p["b"] > 0.002)
        ucfd = griddata((ca, cb), c[:, 3], (p["a"][m], p["b"][m]), method="linear")
        ok = np.isfinite(ucfd)
        a, b, upiv, ucfd = p["a"][m][ok], p["b"][m][ok], p["ux"][m][ok], ucfd[ok]
        du = (ucfd - upiv) / U_INF
        da = np.median(np.diff(np.unique(np.round(a, 5)))) * np.median(np.diff(np.unique(np.round(b, 5))))
        w_cfd, w_piv = (ucfd < 0.7 * U_INF).sum() * da, (upiv < 0.7 * U_INF).sum() * da
        metrics.append(f"{name},{len(du)},{du.mean():.4f},{np.sqrt((du ** 2).mean()):.4f},"
                       f"{(np.abs(du) < 0.05).mean():.3f},{w_cfd:.5f},{w_piv:.5f}")
        np.savetxt(case / f"piv_{name}.csv", np.c_[a, b, upiv / U_INF, ucfd / U_INF], delimiter=",",
                   header="a_m,b_m,ux_piv,ux_cfd", comments="", fmt="%.5f")
        for j, (val, title) in enumerate(((upiv / U_INF, "PIV"), (ucfd / U_INF, "CFD (SST, LUST)"), (du, "CFD - PIV"))):
            ax = axs[k, j]
            sc = ax.scatter(a / D, b / D, c=val, s=1, cmap="RdBu_r" if j == 2 else "viridis",
                            vmin=-0.4 if j == 2 else -0.2, vmax=0.4 if j == 2 else 1.2, rasterized=True)
            ax.set_aspect("equal")
            ax.set_title(f"{name}: {title}", fontsize=9)
            ax.set_xlabel("x / D" if name.startswith("y") else "z / D")
            ax.set_ylabel("y / D")
            fig.colorbar(sc, ax=ax, shrink=0.8, label="U_x / U_inf" if j < 2 else "dU_x / U_inf")
            if name.startswith("y"):
                ax.add_patch(plt.Circle((0, wg.AXLE_H / D), wg.R / D, fill=False, lw=0.6, color="k"))
    fig.tight_layout()
    fig.savefig(case / "piv_compare.png", dpi=90)
    (case / "piv_metrics.csv").write_text("\n".join(metrics) + "\n")
    print("\n".join(metrics))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", type=Path, nargs="?", default=RUNS / "wheel/wheel")
    ap.add_argument("--time", default="3000")
    ap.add_argument("--procs", type=int, default=16, help="processor directories of the decomposed case")
    a = ap.parse_args()
    out = sample(a.case, a.time, a.procs)
    compare(a.case, out)


if __name__ == "__main__":
    main()
