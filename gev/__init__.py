"""gev: meshing, solver settings and evaluation of validated OpenFOAM cases for ground-effect aerodynamics.

Requires OpenFOAM v2506 (openfoam.com) sourced in the shell that runs the scripts.
"""

import os
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent
TEMPLATES = PACKAGE / "templates"
# case directories live outside the repository; override with $GEV_RUNS
RUNS = Path(os.environ.get("GEV_RUNS", Path.home() / "gev-runs"))
