"""Smoke test without third-party data: NACA 4412, inverted, at h/c = 0.2 above a moving ground, meshed with the
validated gmsh profile of the Tyrrell wing at doubled cell sizes, 100 iterations of simpleFoam.

Usage (from the repository root, OpenFOAM environment sourced):  python tests/smoke_naca.py [--procs 2]
Checks that the mesh has no errors and that the forces are written and finite. Takes a few seconds on 4 cores.
"""

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from casetools import RUNS, geometry2d, mesh2d_gmsh  # noqa: E402
from casetools.foam import n_cells, read_dat, sh  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=2)
    a = ap.parse_args()
    case = RUNS / "tests/naca4412_h0200"
    c = 0.2234
    suction, pressure = geometry2d.naca4("4412")
    body = geometry2d.single_element(suction, pressure, h_over_c=0.2, alpha=0.0, chord=c, x_le=-c)
    p = dict(mesh2d_gmsh.scaled(mesh2d_gmsh.TYRRELL, 2.0), n_procs=a.procs)
    mesh2d_gmsh.mesh_case(case, body, p)
    errors = mesh2d_gmsh.mesh_errors(case)
    assert not errors, errors
    sh(case, "foamDictionary system/controlDict -entry endTime -set 100", "foamDictionary")
    sh(case, "decomposePar -force", "decomposePar.run")
    sh(case, f"mpirun -np {a.procs} simpleFoam -parallel", "simpleFoam")
    f = read_dat(case / "postProcessing/forceCoeffs_main/0/coefficient.dat")
    assert len(f) and np.isfinite(f[-1]).all(), "forces missing or not finite"
    print(f"ok: {n_cells(case)} cells, after {int(f[-1, 0])} iterations Cl = {-f[-1, 4]:.3f} (downforce +), "
          f"Cd = {f[-1, 1]:.4f}")
    for d in case.glob("processor*"):
        shutil.rmtree(d)


if __name__ == "__main__":
    main()
