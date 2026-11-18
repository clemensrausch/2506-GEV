"""Surface pressure of a solved Ahmed case vs. Lienhart & Becker (ERCOFTAC case 82, ahmed-<slant>-press.dat).

Reads the reconstructed case (constant/polyMesh and <time>/p, binary) and takes, for every pressure tap, the
nearest face of the "body" patch; the cell value of p is used (zeroGradient wall). Cp = (p - p_inf) / (U^2/2)
with p_inf the mean pressure on the inlet patch. ERCOFTAC coordinates: x = 0 at the base, y span (taps on
the y < 0 side, mirrored to the half model), z up from the ground.

The tap data are not part of this repository (see SOURCES.md); put ahmed-25-press.dat and ahmed-35-press.dat
into data/.

Usage (from the repository root, after reconstructPar -latestTime -fields '(pMean)'):
    python cases/ahmed/compare_pressure.py CASE [--slant 25] [--field p|pMean] [--time T]
writes CASE/cp_compare.png and CASE/cp_compare.csv.
"""

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

DATA = Path(__file__).resolve().parent / "data"


def _list(b, dtype, start=0, ncomp=1):
    """Binary OpenFOAM list after position start (ncomp components per entry): (array, end position)."""
    m = re.compile(rb"\n(\d+)\n\(").search(b, start)
    n, s = int(m.group(1)) * ncomp, m.end()
    size = np.dtype(dtype).itemsize
    return np.frombuffer(b[s:s + size * n], dtype=dtype), s + size * n


def read_mesh(case):
    pm = case / "constant/polyMesh"
    pts = _list((pm / "points").read_bytes(), "<f8", ncomp=3)[0].reshape(-1, 3)
    b = (pm / "faces").read_bytes()
    off, end = _list(b, "<i4")
    lab = _list(b, "<i4", end)[0]
    owner = _list((pm / "owner").read_bytes(), "<i4")[0]
    bd = {}
    for m in re.finditer(r"\n\s+(\w+)\n\s+\{(.*?)\}", (pm / "boundary").read_text(errors="ignore"), re.S):
        d = dict(re.findall(r"(\w+)\s+([^;]+);", m.group(2)))
        if "startFace" in d:
            bd[m.group(1)] = (int(d["nFaces"]), int(d["startFace"]))
    return pts, off, lab, owner, bd


def patch_centres(pts, off, lab, n, s):
    return np.array([pts[lab[off[f]:off[f + 1]]].mean(axis=0) for f in range(s, s + n)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", type=Path)
    ap.add_argument("--slant", type=int, default=25)
    ap.add_argument("--field", default="p")
    ap.add_argument("--time", default=None)
    ap.add_argument("--U", type=float, default=40.0)
    a = ap.parse_args()
    case = a.case.expanduser()
    t = a.time or max((d.name for d in case.iterdir() if re.fullmatch(r"[0-9.e+-]+", d.name) and d.name != "0"),
                      key=float)
    pts, off, lab, owner, bd = read_mesh(case)
    p = _list((case / t / a.field).read_bytes(), "<f8")[0]

    n, s = bd["body"]
    xc = patch_centres(pts, off, lab, n, s)
    p_body = p[owner[s:s + n]]
    ni, si = bd["inlet"]
    p_inf = p[owner[si:si + ni]].mean()
    cp = (p_body - p_inf) / (0.5 * a.U ** 2)

    e = np.loadtxt(DATA / f"ahmed-{a.slant}-press.dat", comments="#")
    tap = np.c_[e[:, 0], e[:, 2], np.abs(e[:, 1])] / 1000.0      # -> CFD (x, y up, z span >= 0)
    d, i = cKDTree(xc).query(tap)
    cp_cfd = cp[i]

    # groups of taps: roof ahead of the slant, slant, base, side wall
    x, z_up = e[:, 0], e[:, 2]
    xs = -222.0 * np.cos(np.radians(a.slant)) - 2.0     # start of the slant (222 mm along the surface)
    groups = {
        "roof ahead of the slant": (x < xs) & (z_up > 300),
        "slant": (x > xs) & (x < -1) & (np.abs(e[:, 1]) < 194) & (z_up > 200),
        "base": np.abs(x) < 1,
        "side wall": np.abs(np.abs(e[:, 1]) - 194.5) < 0.6,
    }
    lines = ["x_mm,y_mm,z_mm,Cp_expt,Cp_cfd,dist_mm,group"]
    gname = np.full(len(e), "", object)
    for g, sel in groups.items():
        gname[sel & (gname == "")] = g
    for k in range(len(e)):
        lines.append(f"{e[k, 0]},{e[k, 1]},{e[k, 2]},{e[k, 3]},{cp_cfd[k]:.4f},{d[k] * 1000:.2f},{gname[k]}")
    (case / "cp_compare.csv").write_text("\n".join(lines) + "\n")

    fig, ax = plt.subplots(2, 2, figsize=(13, 10))
    print(f"time {t}, field {a.field}, p_inf {p_inf:.2f} m2/s2, taps {len(e)}, max tap-face distance {d.max() * 1000:.1f} mm")
    for axi, (g, sel) in zip(ax.ravel(), groups.items()):
        sel = gname == g
        axi.scatter(e[sel, 3], cp_cfd[sel], s=12)
        lim = [min(e[sel, 3].min(), cp_cfd[sel].min()) - 0.05, max(e[sel, 3].max(), cp_cfd[sel].max()) + 0.05]
        axi.plot(lim, lim, "k--", lw=0.8)
        rms = np.sqrt(np.mean((cp_cfd[sel] - e[sel, 3]) ** 2))
        bias = np.mean(cp_cfd[sel] - e[sel, 3])
        axi.set_title(f"{g}: {sel.sum()} taps, mean deviation {bias:+.3f}, RMS {rms:.3f}")
        axi.set_xlabel("Cp Experiment"); axi.set_ylabel("Cp CFD"); axi.grid(alpha=0.3)
        print(f"  {g:22s} n={sel.sum():3d}  Cp expt mean {e[sel, 3].mean():+.3f}  CFD mean {cp_cfd[sel].mean():+.3f}  "
              f"bias {bias:+.3f}  rms {rms:.3f}")
    fig.suptitle(f"Ahmed {a.slant}°: wall pressure CFD ({case.name}, {a.field}, t = {t}) vs. Lienhart & Becker")
    fig.tight_layout()
    fig.savefig(case / "cp_compare.png", dpi=90)


if __name__ == "__main__":
    main()
