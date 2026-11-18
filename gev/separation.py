"""Flow separation on wall patches of a 2D case, from the wall shear stress.

On each patch the faces are ordered along the direction d (x; for a patch that does not run in x, such
as a steep flap, its main direction, the principal axis of the face centres oriented downstream) and
the tangential skin friction

    cf = -(tau_w . t) / (0.5 U_inf^2)

is evaluated, t being the unit surface tangent pointing downstream (along d). OpenFOAM's
wallShearStress is the stress the wall exerts on the fluid (kinematic, m^2/s^2), so it
points against an attached flow: cf > 0 is attached, cf < 0 is reversed (separated) flow.
The tangent orientation is only defined for surfaces that run along one direction; patches that turn
around (blunt bases, rectangle front/rear faces, a whole closed profile) are not evaluated.
Between a rounded nose and its stagnation point the attached flow runs upstream and shows up
as a short reversed region at the upstream end of the patch (e.g. the wing's pressure side).
Evaluate the patches where separation matters (suction side, diffuser roof), not those.

Usage (OpenFOAM environment sourced):
    python -m gev.separation <case> [--patches body_suction ...] [--U 31]
"""

import argparse
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def latest_time(case):
    times = [d.name for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.eE+-]+", d.name)]
    return max(times, key=float)


def ensure_fields(case, time):
    """Write wallShearStress and the cell/face centres C for `time` if they are missing."""
    tdir = case / time
    if not (tdir / "wallShearStress").exists():
        app = re.search(r"application\s+(\w+);", (case / "system/controlDict").read_text()).group(1)
        _run(case, [app, "-postProcess", "-func", "wallShearStress", "-time", time], "wallShearStress")
    if not (tdir / "C").exists():
        _run(case, ["postProcess", "-func", "writeCellCentres", "-time", time], "writeCellCentres")


def _run(case, cmd, log):
    with open(case / f"log.{log}", "w") as fh:
        r = subprocess.run(cmd, cwd=case, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError(f"{' '.join(cmd)} failed, see {case / f'log.{log}'}")


def patch_vectors(field_file, patch):
    """Face values of a vector field on one patch (ascii, nonuniform List<vector>)."""
    text = Path(field_file).read_text()
    m = re.search(rf"\n\s*{re.escape(patch)}\s*\n\s*\{{(.*?)\n\s*\}}", text.split("boundaryField", 1)[1], re.S)
    if not m:
        raise KeyError(f"patch {patch} not in {field_file}")
    v = re.search(r"value\s+nonuniform\s+List<vector>\s*(\d+)\s*\((.*)\)\s*;", m.group(1), re.S)
    if not v:
        raise ValueError(f"{field_file}: patch {patch} has no nonuniform vector value")
    n = int(v.group(1))
    vals = np.array(re.findall(r"\(([^()]+)\)", v.group(2))[:n], dtype=object)
    return np.array([list(map(float, s.split())) for s in vals])


@dataclass
class PatchSeparation:
    patch: str
    x: np.ndarray          # face centres, ordered by x
    y: np.ndarray
    cf: np.ndarray         # tangential skin friction, + attached
    regions: list          # [(x_start, x_end), ...] of reversed flow, x ordered
    frac: float            # reversed-flow fraction of the patch's surface length

    @property
    def x_sep(self):
        """Start of the first reversed-flow region, None if the patch is attached."""
        return self.regions[0][0] if self.regions else None


def patch_separation(case, patch, U_inf, time=None, min_faces=2):
    """Skin friction and reversed-flow regions on one wall patch.

    Reversed-flow runs shorter than `min_faces` faces are ignored (single-face sign flips
    at mesh irregularities).
    """
    case = Path(case)
    time = time or latest_time(case)
    ensure_fields(case, time)
    c = patch_vectors(case / time / "C", patch)[:, :2]
    tau = patch_vectors(case / time / "wallShearStress", patch)[:, :2]
    d = np.array([1.0, 0.0])                        # x for patches that run in x (as validated) ...
    tx = np.gradient(c[np.argsort(c[:, 0])], axis=0)
    if np.mean(np.abs(tx[:, 0]) / np.linalg.norm(tx, axis=1)) < 0.7:   # ... the main direction for steep ones
        d = np.linalg.svd(c - c.mean(axis=0), full_matrices=False)[2][0]
        d *= np.sign(d[0]) or 1.0                   # downstream
    order = np.argsort(c @ d)
    c, tau = c[order], tau[order]

    t = np.gradient(c, axis=0)                      # surface tangent from neighbouring faces
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    if np.mean(np.abs(t @ d)) < 0.7:
        raise ValueError(f"patch {patch} does not run along one direction; tangent orientation undefined")
    t *= np.sign(t @ d)[:, None]                    # downstream orientation
    cf = -np.sum(tau * t, axis=1) / (0.5 * U_inf ** 2)

    ds = np.linalg.norm(np.gradient(c, axis=0), axis=1)   # surface length per face
    rev = cf < 0
    regions, length = [], 0.0
    i = 0
    while i < len(rev):
        if rev[i]:
            j = i
            while j + 1 < len(rev) and rev[j + 1]:
                j += 1
            if j - i + 1 >= min_faces:
                # interpolated zero crossings of cf at both ends of the run
                xa = c[i, 0] if i == 0 else np.interp(0, [cf[i], cf[i - 1]], [c[i, 0], c[i - 1, 0]])
                xb = c[j, 0] if j == len(rev) - 1 else np.interp(0, [cf[j], cf[j + 1]], [c[j, 0], c[j + 1, 0]])
                regions.append((xa, xb))
                length += ds[i:j + 1].sum()
            i = j + 1
        else:
            i += 1
    return PatchSeparation(patch, c[:, 0], c[:, 1], cf, regions, length / ds.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", type=Path)
    ap.add_argument("--patches", nargs="+", default=["body_suction"])
    ap.add_argument("--U", type=float, default=31.0, help="free-stream velocity [m/s]")
    ap.add_argument("--time")
    a = ap.parse_args()
    for p in a.patches:
        s = patch_separation(a.case, p, a.U, a.time)
        reg = ", ".join(f"{x0:.5f}..{x1:.5f}" for x0, x1 in s.regions) or "none"
        print(f"{p}: reversed flow {100 * s.frac:.1f} % of surface, regions (x): {reg}")


if __name__ == "__main__":
    main()
