"""Solve a meshed 3D case (mesh3d.mesh): rebalance, potentialFoam start, simpleFoam with force-based stop rules.

Stop rules (casetools.convergence): the run stops when the forces stand (plateau: Cl(f), Cl(r), Cd within 0.1 % over
150 iterations) or when the means of two consecutive windows (whole shedding periods, at least 400 iterations)
differ by less than 0.5 % (lift) / 1 % (drag) after iteration 1200. Forces are averaged over the last 200
iterations (plateau) or the last window (settled, or unsettled at the iteration limit = "unsteady").

Rebalancing: snappy adds the layer cells on the processors that own the walls; after layers the largest
processor had ~1.9x the mean cell count (Ahmed). The mesh is reconstructed and decomposed again.

Usage from Python: solve(case, profile, max_iter=3000) -> dict with the force means.
"""

import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from . import convergence
from .foam import mpirun, run, write

COEFFS = "forceCoeffs_car"
AVG_START = 400     # fieldAverage of p and U from this iteration (Ahmed 25: force means stood at ~300)


def coeffs(case):
    """Rows (Time Cd Cd(f) Cd(r) Cl Cl(f) Cl(r)) of the whole-body coefficients written so far, or None."""
    return convergence.coeffs(case, COEFFS)


def rebalance(case, p):
    """Reconstruct the snappy mesh (unless that is done already) and decompose it again, fields from 0/.
    constant/polyMesh holds the blockMesh background until reconstructParMesh, so compare file times."""
    whole, part = case / "constant/polyMesh/owner", case / "processor0/constant/polyMesh/owner"
    if part.exists() and (not whole.exists() or part.stat().st_mtime > whole.stat().st_mtime):
        run(case, "reconstructParMesh", "-constant")
    if "car" not in (case / "constant/polyMesh/boundary").read_text():
        raise RuntimeError(f"{case}: constant/polyMesh has no car patches (snappy mesh missing)")
    for proc in case.glob("processor*"):
        shutil.rmtree(proc)
    run(case, "foamDictionary", "system/decomposeParDict", "-entry", "numberOfSubdomains", "-set",
               str(p["n_procs"]), log="foamDictionary")   # the mesh may come from a machine with other core count
    run(case, "decomposePar", "-force", log="decomposePar.solve")


def write_field_average(case, start=AVG_START):
    """Mean p and U from iteration start (system/fieldAverage, included by controlDict)."""
    write(case / "system/fieldAverage", "fieldAverage", f"""\
fieldAverage1
{{
    type            fieldAverage;
    libs            (fieldFunctionObjects);
    timeStart       {start};
    writeControl    writeTime;
    fields
    (
        p {{ mean on; prime2Mean off; base iteration; }}
        U {{ mean on; prime2Mean off; base iteration; }}
    );
}}
""")
    cd = case / "system/controlDict"     # foamDictionary -set expands includes, so hook in at "functions {"
    text = cd.read_text()
    if '#include "fieldAverage"' not in text and "fieldAverage1" not in text:
        cd.write_text(text.replace("\nfunctions\n{\n", '\nfunctions\n{\n    #include "fieldAverage"\n', 1))


def solve(case, p, max_iter=3000, poll=10.0, potential=True):
    """potential: start from potentialFoam (attached flow) instead of the uniform 0/ fields."""
    case = Path(case)
    t0 = time.time()
    rebalance(case, p)
    run(case, "foamDictionary", "system/controlDict", "-entry", "endTime", "-set", str(max_iter),
               log="foamDictionary")
    run(case, "foamDictionary", "system/controlDict", "-entry", "stopAt", "-set", "endTime",
               log="foamDictionary")       # a stopped run leaves writeNow behind
    write_field_average(case)
    if potential:
        mpirun(case, p, "potentialFoam", "-writephi", log="potentialFoam")
    t1 = time.time()

    cmd = ["mpirun", *p["mpi_args"], "-np", str(p["n_procs"]), "simpleFoam", "-parallel"]
    stop = (None, None)
    with open(case / "log.simpleFoam", "w") as fh:
        proc = subprocess.Popen(cmd, cwd=case, stdout=fh, stderr=subprocess.STDOUT)
        while proc.poll() is None:
            time.sleep(poll)
            if stop[0] is None:
                a = coeffs(case)
                if (it := convergence.plateau(a)) is not None:
                    stop = ("plateau", it)
                elif (st := convergence.settled(a)) is not None:
                    stop = ("settled", st)
                if stop[0]:
                    run(case, "foamDictionary", "system/controlDict", "-entry", "stopAt", "-set", "writeNow",
                               log="foamDictionary")
        if proc.returncode:
            raise RuntimeError(f"simpleFoam failed, see {case / 'log.simpleFoam'}")
    t2 = time.time()

    a = coeffs(case)
    if a is None or not np.abs(a[:, 1]).any():
        raise RuntimeError(f"no forces on the car patches, check the mesh in {case}")
    t = a[:, 0]
    log = (case / "log.simpleFoam").read_text()
    if stop[0] == "settled":
        win, status = stop[1][1], "ok"
    elif stop[0] == "plateau" or "SIMPLE solution converged" in log:
        win, status = 200, "ok"
    else:
        win, status = convergence.window(a), "unsteady"
    m = a[t > t[-1] - win].mean(axis=0)
    return dict(status=status, stop=stop[0], iterations=int(t[-1]), window=int(win),
                Cd=m[1], Cl=m[4], Cl_front=m[5], Cl_rear=m[6],
                time_prepare=t1 - t0, time_solve=t2 - t1, s_per_iter=(t2 - t1) / max(t[-1], 1))


def part_forces(case, part, window, it_end):
    """Mean force vector (N, modelled geometry) of one part over the last window iterations."""
    f = case / f"postProcessing/forces_{part}/0/force.dat"
    a = np.array([[float(x) for x in l.replace("(", " ").replace(")", " ").split()[:4]]
                  for l in f.read_text().splitlines() if l.strip() and not l.startswith("#")])
    sel = (a[:, 0] > it_end - window) & (a[:, 0] <= it_end)
    return a[sel, 1:4].mean(axis=0)
