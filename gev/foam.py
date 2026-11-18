"""Small helpers around OpenFOAM: dictionary files, running utilities, reading function-object output."""

import re
import subprocess
from pathlib import Path

import numpy as np

HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2506                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
FOOTER = "\n// " + "*" * 73 + " //\n"


def write(path, obj, body):
    """Write an OpenFOAM dictionary with the standard header."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER.format(obj=obj) + body + FOOTER)


def run(case, *cmd, log=None):
    """Run one command in the case directory, output to log.<name>; raise on failure."""
    case = Path(case)
    log = case / f"log.{log or Path(str(cmd[0])).name}"
    with open(log, "w") as fh:
        r = subprocess.run([str(c) for c in cmd], cwd=case, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed, see {log}")


def mpirun(case, p, *cmd, log=None):
    """Run an OpenFOAM application in parallel with the profile's process count and MPI arguments."""
    run(case, "mpirun", *p["mpi_args"], "-np", str(p["n_procs"]), *cmd, "-parallel", log=log or cmd[0])


def sh(case, cmd, log):
    """Run a shell command line in the case directory, output to log.<log>; raise on failure."""
    case = Path(case)
    with open(case / f"log.{log}", "w") as fh:
        r = subprocess.run(cmd, cwd=case, stdout=fh, stderr=subprocess.STDOUT, shell=True)
    if r.returncode:
        raise RuntimeError(f"{cmd} failed, see {case / f'log.{log}'}")


def set_entry(case, file, entry, value):
    """foamDictionary -entry <entry> -set <value> on a case file."""
    run(case, "foamDictionary", file, "-entry", entry, "-set", str(value), log="foamDictionary")


def set_patch_types(boundary_file, types):
    """Rewrite 'type' of patches in an (ascii) polyMesh/boundary file."""
    boundary_file = Path(boundary_file)
    text = boundary_file.read_text()
    for name, typ in types.items():
        m = re.search(rf"(\n\s*{name}\s*\n\s*\{{)(.*?)(\}})", text, re.S)
        if not m:
            raise KeyError(f"patch {name} not in {boundary_file}")
        body = re.sub(r"\n\s*physicalType\s+\w+;", "", m.group(2))
        body = re.sub(r"\n\s*inGroups\s+[^;]*;", "", body)
        body = re.sub(r"type\s+\w+;", f"type            {typ};", body)
        text = text[:m.start(2)] + body + text[m.end(2):]
    boundary_file.write_text(text)


def read_dat(path, ncol=None):
    """Numeric rows of a function-object .dat file (comment lines and vector brackets ignored)."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        v = line.replace("(", " ").replace(")", " ").split()
        try:
            rows.append([float(x) for x in (v[:ncol] if ncol else v)])
        except ValueError:          # line being written
            pass
    return np.array(rows)


def latest_time(case):
    """Name of the latest time directory of a (reconstructed) case."""
    case = Path(case)
    times = [d.name for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.eE+-]+", d.name)]
    return max(times, key=float)


def n_cells(case, log="log.checkMesh"):
    """Cell count from a checkMesh log."""
    for line in (Path(case) / log).read_text().splitlines():
        if line.strip().startswith("cells:"):
            return int(line.split()[1])
    return None
