"""Plot a 2D OpenFOAM mesh (ascii polyMesh) from the faces of one empty patch.

Usage:  python -m casetools.plot_mesh <case> out.png [--patch back] [--window x0 x1 y0 y1 ...]
"""

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection


def _body(text):
    """Contents of the outermost list of an ascii OpenFOAM list file."""
    text = text[text.index("// * *"):] if "// * *" in text else text
    start = re.search(r"\n\s*\d+\s*\n\s*\(", text)
    return text[start.end(): text.rindex(")")]


def read_patch(case, patch):
    m = Path(case) / "constant/polyMesh"
    pts = np.array(re.findall(r"\(([^()]+)\)", _body((m / "points").read_text())), dtype=object)
    pts = np.array([list(map(float, p.split())) for p in pts])
    faces = [list(map(int, f.split())) for f in re.findall(r"\d+\(([^()]+)\)", _body((m / "faces").read_text()))]
    bnd = (m / "boundary").read_text()
    b = re.search(rf"\n\s*{patch}\s*\n\s*\{{(.*?)\}}", bnd, re.S).group(1)
    n = int(re.search(r"nFaces\s+(\d+)", b).group(1))
    s = int(re.search(r"startFace\s+(\d+)", b).group(1))
    walls = {}
    for name, blk in re.findall(r"\n\s*(\w+)\s*\n\s*\{(.*?)\}", bnd, re.S):
        if re.search(r"type\s+wall;", blk):
            wn = int(re.search(r"nFaces\s+(\d+)", blk).group(1))
            ws = int(re.search(r"startFace\s+(\d+)", blk).group(1))
            walls[name] = faces[ws: ws + wn]
    return pts, faces[s: s + n], walls


def plot(case, out, patch="back", windows=None, title=None, overlay=None):
    pts, faces, walls = read_patch(case, patch)
    polys = [pts[f][:, :2] for f in faces]
    windows = windows or [None]
    fig, axes = plt.subplots(len(windows), 1, figsize=(12, 4.2 * len(windows)), squeeze=False)
    for ax, w in zip(axes[:, 0], windows):
        ax.add_collection(PolyCollection(polys, facecolors="none", edgecolors="0.35", linewidths=0.25))
        for i, (name, fs) in enumerate(walls.items()):
            segs = [pts[f][:, :2] for f in fs]  # wall faces are quads spanning z; project to edges
            ax.add_collection(LineCollection(segs, colors=f"C{i}", linewidths=1.6, label=name))
        if overlay is not None:
            overlay(ax)
        if w:
            ax.set_xlim(w[0], w[1]); ax.set_ylim(w[2], w[3])
        else:
            ax.autoscale()
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=7)
    if title:
        axes[0, 0].set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--patch", default="back")
    ap.add_argument("--window", type=float, nargs=4, action="append")
    a = ap.parse_args()
    plot(a.case, a.out, a.patch, a.window)
