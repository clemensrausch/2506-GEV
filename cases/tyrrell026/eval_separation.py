"""Trailing-edge separation on the suction side, compared with Zerihan (2001) Table 9 (oil flow).

Evaluates the existing sweep and grid-study cases with gev.separation (no new runs): skin friction along the
suction side and the start of reversed flow x_sep/c, measured along the chord line from the leading edge. For
runs that end periodic (h/c = 0.090) the values are a snapshot of the last iteration.

Usage (from the repository root, OpenFOAM environment sourced):
    python cases/tyrrell026/eval_separation.py [--mesher gmsh|snappy]
Writes results/separation_xsep<tag>.csv and results/separation_cf<tag>.csv (CFD only); with
data/zerihan_table9_xsep.csv present (see SOURCES.md) also results/separation<tag>.png.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import _path  # noqa: F401
import wing
from gev import separation
from gev.tables import read_csv
from run_sweep import HEIGHTS, case_dir

HERE = Path(__file__).resolve().parent
U_INF = 31.0
MESHES = {"coarse": "_coarse", "medium": "", "fine": "_fine"}


def chord_coordinate(h):
    """(x, y) -> x/c along the chord line, leading edge 0, trailing edge (base centre) 1."""
    b = wing.single(h_over_c=h, alpha=1.0)
    le = b.segments["suction"][0]
    d = 0.5 * (b.segments["suction"][-1] + b.segments["pressure"][0]) - le
    return lambda x, y: ((np.c_[x, y] - le) @ d) / (d @ d)


def evaluate(h, suffix=""):
    s = separation.patch_separation(case_dir(h, suffix), "body_suction", U_INF)
    xc = chord_coordinate(h)(s.x, s.y)
    xs = None
    if s.x_sep is not None:
        xs = chord_coordinate(h)(np.array([s.x_sep]), np.array([np.interp(s.x_sep, s.x, s.y)]))[0]
    return s, xc, xs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesher", choices=("snappy", "gmsh"), default="gmsh")
    tag = "_gmsh" if ap.parse_args().mesher == "gmsh" else ""
    meshes = {m: tag + s for m, s in MESHES.items()}
    t9_file = HERE / "data/zerihan_table9_xsep.csv"
    t9 = read_csv(t9_file) if t9_file.exists() else None
    expt = dict(zip(np.round(t9["h_over_c"], 3), t9["x_sep_over_c"])) if t9 is not None else {}
    out = HERE / "results"

    rows = ["h_over_c,mesh,x_sep_over_c,separated_fraction,n_regions"]
    cf_rows = ["h_over_c,x_over_c,cf"]
    curves = {}
    for h in HEIGHTS:
        for mesh_name, suf in meshes.items():
            if not (case_dir(h, suf) / "postProcessing").exists():
                continue
            s, xc, xs = evaluate(h, suf)
            e = expt.get(round(h, 3), np.nan)
            rows.append(f"{h},{mesh_name},{'' if xs is None else f'{xs:.4f}'},{s.frac:.4f},{len(s.regions)}")
            print(f"h/c = {h:5.3f} {mesh_name:6s}: x_sep/c = {'attached' if xs is None else f'{xs:.3f}':>8}"
                  f"  reversed {100 * s.frac:4.1f} % of surface  expt {e}")
            if mesh_name == "medium":
                curves[h] = (xc, s.cf)
                cf_rows += [f"{h},{a:.5f},{b:.6f}" for a, b in zip(xc, s.cf)]
    (out / f"separation_xsep{tag}.csv").write_text("\n".join(rows) + "\n")
    (out / f"separation_cf{tag}.csv").write_text("\n".join(cf_rows) + "\n")
    print("written:", out / f"separation_xsep{tag}.csv")
    if t9 is None:
        return

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(12, 4.5))
    for i, (h, (xc, cf)) in enumerate(curves.items()):
        a0.plot(xc, cf, color=f"C{i}", label=f"h/c = {h}")
        if round(h, 3) in expt:
            a0.axvline(expt[round(h, 3)], color=f"C{i}", ls=":", lw=1)
    a0.set_xlim(0.5, 1.0); a0.set_ylim(-0.002, 0.008); a0.axhline(0, color="k", lw=0.8)
    a0.set_xlabel("x/c"); a0.set_ylabel("cf (attached +)"); a0.legend(fontsize=8); a0.grid(alpha=0.3)
    r = np.genfromtxt(out / f"separation_xsep{tag}.csv", delimiter=",", names=True, dtype=None, encoding=None)
    med = r[(r["mesh"] == "medium") & np.isfinite(r["x_sep_over_c"])]
    a1.plot(med["h_over_c"], med["x_sep_over_c"], "C3o-", label="this setup (medium)")
    for m, mk in (("fine", "D"), ("coarse", "s")):
        q = r[(r["mesh"] == m) & np.isfinite(r["x_sep_over_c"])]
        a1.plot(q["h_over_c"], q["x_sep_over_c"], "C3" + mk, mfc="none", label=f"this setup ({m})")
    a1.plot(t9["h_over_c"], t9["x_sep_over_c"], "ko", mfc="none", label="Zerihan oil flow (Table 9)")
    a1.set_xlabel("h/c"); a1.set_ylabel("x_sep/c"); a1.legend(); a1.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / f"separation{tag}.png", dpi=130)


if __name__ == "__main__":
    main()
