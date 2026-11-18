# Ground-Effect Validation for OpenFOAM v2506

Steady RANS setups for OpenFOAM v2506 (openfoam.com) for bodies close to the ground: wings in ground
effect, bluff bodies, underbody diffusers and rotating wheels. Each setup was checked against a published
wind-tunnel experiment. The repository contains the settings, the Python scripts that turn a geometry into
a complete case, the results of the validation runs and a write-up of each comparison, including where the
setup fails.

The cases share one set of physics: k-ω SST, SIMPLEC, the same wall functions and the same mesh logic.
The 2D setup started from the Ahmed reference mesh scaled to the wing. The 3D template then took over the
2D settings and was tested on three bodies of increasing difficulty.

| Case | Tests | Result | Not covered |
|---|---|---|---|
| [Tyrrell 026 wing, 2D](cases/tyrrell026/) | wing over a moving ground, y⁺ ≈ 1, 7 heights | lift 1 to 12 % above the measurement, maximum at the measured height; separation 0.09 to 0.15c too early | flow below the lift maximum, transition |
| [Ahmed body, 3D](cases/ahmed/) | separation at sharp rear edges, 25° and 35° slant | drag within 8 %, wall pressure within ΔC_p = 0.05 except on the 25° slant (0.13) | lift, moving ground |
| [Diffuser body, 3D](cases/diffuser/) | underbody diffuser over a moving ground, 10 heights | downforce within 6.7 % (RMS) at and above the force maximum | the asymmetric stall below the maximum |
| [Rotating wheel, 3D](cases/wheel/) | rotating wall, contact patch, wake | wake below axle height within 10 % | upper wake, wheel forces |

Each case folder has a README with experiment, setup, verification, comparison and limits, and a
SOURCES.md for the third-party data.

## Requirements

- OpenFOAM v2506 from openfoam.com, sourced in the shell: `source /usr/lib/openfoam/openfoam2506/etc/bashrc`
- Python 3.12 with the packages in `requirements.txt` (`pip install -r requirements.txt`)
- For the 3D cases: 8 to 16 cores and 32 GB of memory; the wheel mesh needs about 2.8 GB per million cells

Check the installation with the smoke test. It meshes a NACA 4412 in ground effect with the validated 2D
settings at doubled cell sizes and runs 100 iterations, in a few seconds:

```
python tests/smoke_naca.py
```

Cases are written to `$GEV_RUNS` (default `~/gev-runs`). All scripts run from the repository root.

## Layout

```
casetools/                Python tools that build, mesh, run and evaluate the cases
  geometry2d.py           2D contours: Body, wings from tables, NACA 4-digit sections
  mesh2d_gmsh.py          2D mesh with gmsh (quads, prism layers, fans at sharp corners), profile TYRRELL
  mesh2d_snappy.py        2D mesh with snappyHexMesh on a one-cell slab, profile TYRRELL
  case2d.py               forces and flow values of a 2D case
  mesh3d.py               3D case from STL parts with snappyHexMesh, profiles DEFAULTS and AHMED
  run3d.py                3D solve: rebalance, potentialFoam start, simpleFoam with force-based stop rules
  convergence.py          the stop rules: force plateau and settled window means
  separation.py           separation point on 2D walls from the wall shear stress
  templates/2d, 3d        OpenFOAM dictionaries (0, constant, system)
cases/<case>/             geometry and run scripts, README, SOURCES.md, results/validated, figures
tests/smoke_naca.py
```

## Using the settings for another geometry

In 2D, build a `geometry2d.Body` (closed, counter-clockwise, named segments become patches) and mesh it
with a gmsh profile. The profile `TYRRELL` is set up for a chord of 223.4 mm with the leading edge at
x = −c: domain, boxes and cell sizes are absolute lengths. A wing of that size can use it directly; for
another size copy the profile and scale these lengths with the chord.

```python
from casetools import geometry2d, mesh2d_gmsh
suction, pressure = geometry2d.naca4("4412")
wing = geometry2d.single_element(suction, pressure, h_over_c=0.2, alpha=2.0, chord=0.2234, x_le=-0.2234)
mesh2d_gmsh.mesh_case(case, wing, mesh2d_gmsh.TYRRELL)
```

In 3D, describe the body as closed STL parts and pick the wall treatment per part:

```python
from casetools import mesh3d, run3d
parts = [mesh3d.Part("body", "body.stl", layers="wallfn"),
         mesh3d.Part("underbody", "underbody.stl", layers="resolved", level=(5, 5))]
profile = dict(mesh3d.DEFAULTS, U_inf=20.0, moving_ground=True, ref_area=0.1)
mesh3d.mesh(case, parts, profile)
result = run3d.solve(case, profile, max_iter=3000)
```

The case scripts show complete examples, including extra refinement boxes, tunnel walls and rotating
parts. Every profile entry is commented in `mesh3d.py` and `mesh2d_gmsh.py`.

## Settings that decided results

These are the settings where the validation showed that the obvious choice gives a wrong answer. The case
READMEs give the evidence.

- Convection of U in 3D: LUST. With linearUpwindV, the scheme of the motorBike tutorial, the 35° flow
  on the Ahmed slant stayed attached and the drag came out 50 % too high. In 2D linearUpwindV is kept;
  the wing results were validated with it.
- Wall functions that work at any y⁺. nutUSpaldingWallFunction and omegaWallFunction with binomial
  blending, because a body can carry y⁺ ≈ 1 layers on some parts and wall-function layers on others, and
  where layers break off y⁺ lies in between.
- Stop on forces, not residuals. Residuals of separated cases creep long after the forces stand. A
  residual limit of 1e-4 stopped the diffuser runs while C_L still drifted by up to 1.4 %. The rules in
  `convergence.py` stop on a force plateau or on two consecutive window means that agree, the windows
  being whole shedding periods.
- Refinement boxes that hold the separated region. On the wing a gap box that ended below the
  trailing edge left the separation in coarse cells; meshes with 2 % different cell sizes then scattered
  by ±0.9 % in C_l instead of ±0.1 %.
- Tunnel walls where the experiment had them. Modelling the cross-section of the tunnel changed the
  diffuser downforce by 2.5 to 4.6 %.
- Layers in one pass on thin y⁺ ≈ 1 stacks (`nOuterIter 1`). In three passes the thinnest layers
  failed the quality checks and the first cell grew from 18 to 66 µm.

## Hardware

The 2D cases run on a desktop: 100 000 to 180 000 cells, 2 to 5 minutes per case on 8 cores. The 3D
cases have 1.5 to 6.7 million cells and need about 2 GB of memory per million cells while solving;
snappyHexMesh takes more for a short time, up to about 2.8 GB per million cells on the wheel. Meshes above
about 6 million cells need more than 16 GB.

The 3D runs were made on rented cloud servers. A comparison on the coarse Ahmed mesh (0.75 million
cells, 8 MPI processes):

| Server | s per iteration | Cost per 1000 iterations and million cells |
|---|---|---|
| x86, 8 shared cores (AMD EPYC Genoa) | 0.25 | 1.2 ct |
| x86, 4 shared cores | 0.66 | 1.7 ct |
| ARM, 8 shared cores | 0.72 | 1.1 ct |
| x86, 4 dedicated cores with 2 threads each | 0.50 | 4.9 ct |

Shared x86 servers with 8 cores were the best trade-off: three times faster than ARM at almost the same
cost per iteration. Doubling the cores from 4 to 8 made the run 2.65 times faster. On the 1.5 million cell
Ahmed mesh a desktop with a Ryzen 7 5800X (8 cores) needed 2.3 s per iteration, the 8-core x86 server
0.51 s. Spreading one case over several servers does not pay off; one case per server does. Each value is a
single measurement on shared hardware.

## Reproducibility

The scripts here reproduce the validation runs: for the 3D cases they write the same geometry
(byte-identical STL files) and the same dictionaries as the code the runs were made with, and the
evaluation scripts give the archived results from the finished cases. Two exceptions are documented in the
case READMEs. gmsh's quad mesher does not give identical meshes from run to run; a rerun differs by a few
tens of cells and should agree with the archived lift to about 0.1 %.

## What is not in the repository

No meshes, no fields and no third-party data. The wing coordinates, the nose surface of the Ahmed body
and all measurements belong to their authors and are not redistributed here. Each `SOURCES.md` says
where to get them, how to name them and which script reads them; they go into the `data/` folder of the
case, which git ignores. The diffuser body and the wheel are built from published dimensions and run
without any third-party file.

The results in `results/validated/` contain only values of these computations. The figures show the
measurements as cited comparison points.

## Licence

GPL-3.0, see `LICENSE`. The third-party data listed in the `SOURCES.md` files keep their own terms.

Clemens Rausch, [clemensrausch.de](https://clemensrausch.de)
