"""Body contour (geometry2d.Body) -> OpenFOAM 2D mesh via blockMesh + snappyHexMesh + extrudeMesh.

Meshing settings follow the 3D Ahmed reference case: same background cell size, refinement levels, boxes,
layer and quality settings. 2D is obtained by meshing a slab one background cell thick and extruding its
back patch (z = 0) to a single cell.

Coordinates: x = flow, y = vertical (ground at y = 0), z = span. Requires a sourced OpenFOAM v2506
environment.

Usage from Python:  mesh(case, body, TYRRELL)    (body from geometry2d, profile = dict of settings below)
Command line (demo):  python -m gev.mesh2d_snappy CASE   (rectangle with the Ahmed body's proportions)
"""

import argparse
import re
import shutil
from pathlib import Path

import numpy as np

from . import TEMPLATES, geometry2d
from .case2d import write_case_dicts
from .foam import run, set_patch_types, write

# Domain and background mesh as the 3D Ahmed reference case: x -4..10 (180), y 0..4 (80, grading 0.8)
DEFAULTS = dict(
    x_min=-4.0, x_max=10.0, y_max=4.0,
    nx=180, ny=80, y_grading=0.8,
    span=0.1,             # slab / 2D thickness (= Ahmed background cell size in z)
    surface_level=(4, 5), feature_level=4,
    near_level=3, wake_level=2,
    body_layers=3, ground_layers=2,
    layer_first=None, layer_expansion=None,   # None: relative sizes as Ahmed
    box_scale=1.0,        # Ahmed box offsets are scaled by this (body size / Ahmed length)
    le_box=None,          # (half size, refinement) box around the nose
    gap_box=None,         # (upstream, downstream, [refinement, ...]) box between body and ground
    n_procs=8,
    # flow as the Ahmed reference case: 60 m/s, fixed ground
    U_inf=60.0, k_inf=0.135, omega_inf=180.1, moving_ground=False,
    layer_relaxed=None,   # e.g. "maxNonOrtho 75;": quality limits used by the layer addition after nRelaxedIter
    mpi_args=(),          # extra mpirun arguments, e.g. ("--bind-to", "none") when several cases run at once
    x_grading="1",        # blockMesh grading in x (OpenFOAM syntax, multi-grading allowed)
    region_levels=(),     # ((segment regex, (min, max)), ...): surface levels per patch, others surface_level
    region_layers=(),     # ((segment regex, n), ...): prism layers per patch, others body_layers
                          # (written after the "body.*" entry: snappy applies the entries in order)
    extra_boxes=None,     # body -> {name: ((x0, y0), (x1, y1), refinement spec)}
)

AHMED_LENGTH = 1.044

# Tyrrell 026 (c = 0.2234, LE at x = -c, TE at x = 0): the Ahmed mesh scaled by c/L, i.e. background
# cells and boxes x c/L with the Ahmed refinement levels unchanged. The domain follows Zerihan's
# 2D CFD (thesis ch. 9.2): 8c upstream of the LE, 8c downstream of the TE, roof at 7c.
# Extra: nose and gap refinement, wall-resolved layers (y+ ~ 1, first layer 10 um).
# Nose and gap use directional (x, y only) refinement: isotropic refinement would also split the
# slab in z, which only inflates the intermediate 3D mesh (factor 2 per level).
_C = 0.2234
_S = _C / AHMED_LENGTH
TYRRELL = dict(
    DEFAULTS,
    x_min=-9 * _C, x_max=8 * _C, y_max=7 * _C,
    nx=round(17 * _C / (14 / 180 * _S)), ny=round(7 * _C / (4 / 80 * _S)),   # Ahmed cell size x c/L
    span=0.1 * _S,
    box_scale=_S,
    le_box=(0.01, "levelIncrement (3 5 (1 1 0));"),              # nose: x2 in x, y
    gap_box=(0.5 * _C, 1.5 * _C,                                 # gap: every cell to ~level 4
             ["levelIncrement (2 2 (2 2 0));", "levelIncrement (3 3 (1 1 0));"]),
    body_layers=20, ground_layers=20, layer_first=1e-5, layer_expansion=1.2,
    # Zerihan: Re_c = 4.62e5 (his CFD, upper end of the tests) -> U = Re nu / c = 31.0 m/s with
    # the Ahmed nu = 1.5e-5; tunnel turbulence intensity 0.2 % (thesis p. 17), turbulent/laminar
    # viscosity ratio 10; moving belt at U_inf
    U_inf=31.0, k_inf=1.5 * (0.002 * 31.0) ** 2, omega_inf=1.5 * (0.002 * 31.0) ** 2 / (10 * 1.5e-5),
    moving_ground=True,
)


def multi_grading(sections):
    """blockMesh multi-grading for sections [(length, first cell, last cell), ...] -> (cells, grading string).
    Geometric expansion within a section; cell count from the continuous estimate L ln(r) / (d1 - d0)."""
    total_L = sum(L for L, _, _ in sections)
    ns, parts = [], []
    for L, d0, d1 in sections:
        n = round(L / d0) if abs(d1 - d0) < 1e-12 else round(L * np.log(d1 / d0) / (d1 - d0))
        ns.append(max(n, 1))
    for (L, d0, d1), n in zip(sections, ns):
        parts.append(f"({L / total_L:.6g} {n / sum(ns):.6g} {d1 / d0:.6g})")
    return sum(ns), "(" + " ".join(parts) + ")"


def write_stl(body, path, z0, z1):
    """Contour swept from z0 to z1 as an ASCII STL, one solid per named segment."""
    out = []
    for name, seg in body.segments.items():
        out.append(f"solid {name}")
        for a, b in zip(seg[:-1], seg[1:]):
            n = np.array([b[1] - a[1], a[0] - b[0], 0.0])  # outward for a CCW contour
            n /= np.linalg.norm(n)
            quad = [(*a, z0), (*b, z0), (*b, z1), (*a, z1)]
            for tri in ((0, 1, 2), (0, 2, 3)):
                out.append(f"  facet normal {n[0]:.9g} {n[1]:.9g} {n[2]:.9g}\n    outer loop")
                out += [f"      vertex {quad[i][0]:.9g} {quad[i][1]:.9g} {quad[i][2]:.9g}" for i in tri]
                out.append("    endloop\n  endfacet")
        out.append(f"endsolid {name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out) + "\n")


def write_dicts(case, body, p):
    sysd = case / "system"
    pts = body.points()
    xn, xt = pts[:, 0].min(), pts[:, 0].max()     # nose / tail
    yb, yt = pts[:, 1].min(), pts[:, 1].max()     # bottom / top
    s = p["box_scale"]
    ahmed_top = 0.338                             # Ahmed body top; box heights are relative to it
    z0, z1, dz = -0.5 * p["span"], 1.5 * p["span"], p["span"]
    lo, hi = f"{-0.01 * dz:g}", f"{1.01 * dz:g}"  # boxes span the whole slab

    write(sysd / "blockMeshDict", "blockMeshDict", f"""\
scale   1;

xMin {p['x_min']}; xMax {p['x_max']};
yMin 0;          yMax {p['y_max']};
zMin 0;          zMax {dz};

vertices
(
    ($xMin $yMin $zMin) ($xMax $yMin $zMin) ($xMax $yMax $zMin) ($xMin $yMax $zMin)
    ($xMin $yMin $zMax) ($xMax $yMin $zMax) ($xMax $yMax $zMax) ($xMin $yMax $zMax)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({p['nx']} {p['ny']} 1) simpleGrading ({p['x_grading']} {p['y_grading']} 1)
);

boundary
(
    inlet  {{ type patch; faces ((0 4 7 3)); }}
    outlet {{ type patch; faces ((1 2 6 5)); }}
    top    {{ type patch; faces ((3 7 6 2)); }}
    ground {{ type wall;  faces ((0 1 5 4)); }}
    back   {{ type patch; faces ((0 3 2 1)); }}
    front  {{ type patch; faces ((4 5 6 7)); }}
);
""")

    write(sysd / "surfaceFeatureExtractDict", "surfaceFeatureExtractDict", """\
body.stl
{
    extractionMethod    extractFromSurface;
    includedAngle       150;
    geometricTestOnly   yes;    // no feature lines at segment (patch) boundaries on a smooth contour
    subsetFeatures
    {
        nonManifoldEdges    no;
        openEdges           no;
    }
    writeObj            yes;
}
""")

    sl, fl = p["surface_level"], p["feature_level"]
    def per_segment(table, n, default):
        return next((v for rx, v in table if re.fullmatch(rx, n)), default)

    regions = "\n".join(
        f"                {n} {{ name {n}; level ({lv[0]} {lv[1]}); patchInfo {{ type wall; }} }}"
        for n in body.segments for lv in [per_segment(p["region_levels"], n, sl)])
    extra_geo, extra_reg = "", ""
    boxes = {}   # name -> (min, max, refinement spec)
    if p["le_box"]:
        d, spec = p["le_box"]
        yle = pts[pts[:, 0].argmin(), 1]
        boxes["leBox"] = ((xn - d, yle - d), (xn + d, yle + d), spec)
    if p["gap_box"]:
        up, down, specs = p["gap_box"]
        for i, spec in enumerate(specs):
            boxes[f"gapBox{i}"] = ((xn - up, -0.015), (xt + down, 1.25 * yb), spec)
    if p["extra_boxes"]:
        boxes.update(p["extra_boxes"](body))
    for name, (a, b, spec) in boxes.items():
        extra_geo += (f"    {name}\n    {{\n        type searchableBox;\n"
                      f"        min  ({a[0]:g} {a[1]:g} {lo});\n        max  ({b[0]:g} {b[1]:g} {hi});\n    }}\n")
        extra_reg += (f"        {name}\n        {{\n            mode inside;\n            levels ((1E15 0));\n"
                      f"            {spec}\n        }}\n")

    if p["layer_first"]:   # absolute sizes: first layer + expansion (wall-resolved)
        layer_size = (f"relativeSizes false;\n    firstLayerThickness {p['layer_first']:g};\n"
                      f"    expansionRatio {p['layer_expansion']:g};\n    minThickness {0.1 * p['layer_first']:g};")
    else:                  # Ahmed: relative sizes
        layer_size = ("relativeSizes true;\n    expansionRatio 1.0;\n"
                      "    finalLayerThickness 0.3;\n    minThickness 0.1;")

    # thin wall-resolved layers on the coarse far-field ground cells have a tiny determinant
    # (aspect ratio ~1e3-1e4); relaxed as in the v2506 tutorial mesh/snappyHexMesh/airfoilWithLayers
    det = "    minDeterminant 1e-8;" if p["layer_first"] else ""
    if p["layer_relaxed"]:
        det += f"\n    relaxed\n    {{\n        {p['layer_relaxed']}\n    }}"

    layer_regions = "".join(
        f"        body_{n}\n        {{\n            nSurfaceLayers {per_segment(p['region_layers'], n, p['body_layers'])};\n        }}\n"
        for n in body.segments) if p["region_layers"] else ""

    write(sysd / "snappyHexMeshDict", "snappyHexMeshDict", f"""\
castellatedMesh true;
snap            true;
addLayers       true;

geometry
{{
    body.stl
    {{
        type triSurfaceMesh;
        name body;
    }}
    refinementBox
    {{
        type searchableBox;
        min  ({xn - 0.556 * s:g} -0.015 {lo});
        max  ({xt + 0.2 * s:g} {yt + (0.50 - ahmed_top) * s:g} {hi});
    }}
    wakeBox
    {{
        type searchableBox;
        min  ({xt - 0.35 * s:g} -0.015 {lo});
        max  ({xt + 2.7 * s:g} {yt + (1.0 - ahmed_top) * s:g} {hi});
    }}
{extra_geo}}}

castellatedMeshControls
{{
    maxLocalCells 2000000;
    maxGlobalCells 20000000;
    minRefinementCells 10;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 3;

    features
    (
        {{
            file "body.eMesh";
            level {fl};
        }}
    );

    refinementSurfaces
    {{
        body
        {{
            level ({sl[0]} {sl[1]});
            patchInfo
            {{
                type wall;
            }}
            regions
            {{
{regions}
            }}
        }}
    }}

    resolveFeatureAngle 30;

    refinementRegions
    {{
        refinementBox
        {{
            mode inside;
            levels ((1E15 {p['near_level']}));
        }}
        wakeBox
        {{
            mode inside;
            levels ((1E15 {p['wake_level']}));
        }}
{extra_reg}
    }}

    locationInMesh ({p['x_max'] - 0.05 * (p['x_max'] - p['x_min']):g} {0.5625 * p['y_max']:g} {0.5 * dz:g});
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 30;
    nRelaxIter 5;
    nFeatureSnapIter 10;
    implicitFeatureSnap false;
    explicitFeatureSnap true;
    multiRegionFeatureSnap false;
}}

addLayersControls
{{
    {layer_size}
    layers
    {{
        "body.*"
        {{
            nSurfaceLayers {p['body_layers']};
        }}
{layer_regions}
        "ground.*"
        {{
            nSurfaceLayers {p['ground_layers']};
        }}
    }}
    nGrow 0;
    featureAngle 60;
    slipFeatureAngle 30;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}}

meshQualityControls
{{
    #include "meshQualityDict"
{det}
    nSmoothScale 4;
    errorReduction 0.75;
}}

writeFlags
(
    scalarLevels
    layerSets
    layerFields
);

mergeTolerance 1e-6;
""")

    write(sysd / "meshQualityDict", "meshQualityDict", """\
#includeEtc "caseDicts/meshQualityDict"

maxNonOrtho 65;
maxBoundarySkewness 20;
maxInternalSkewness 4;
maxConcave 80;
minVol -1e30;
minTetQuality 1e-15;
minTwist 0.02;
minDeterminant 0.001;
minFaceWeight 0.05;
minVolRatio 0.01;
""")

    write(sysd / "extrudeMeshDict", "extrudeMeshDict", f"""\
constructFrom    patch;
sourceCase       "<case>";
sourcePatches    (back);
exposedPatchName front;

extrudeModel     linearDirection;
direction        (0 0 -1);   // along the outward normal of back (z = 0)
nLayers          1;
thickness        {dz};

flipNormals      false;
mergeFaces       false;
""")

    write_stl(body, case / "constant/triSurface/body.stl", z0, z1)
    write_case_dicts(case, body, p)


def mesh(case, body, p=DEFAULTS):
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(TEMPLATES / "2d", case)
    write_dicts(case, body, p)

    run(case, "blockMesh")
    run(case, "surfaceFeatureExtract")
    run(case, "decomposePar", "-force", "-no-fields")
    run(case, "mpirun", *p["mpi_args"], "-np", str(p["n_procs"]), "snappyHexMesh", "-parallel", "-overwrite",
        log="snappyHexMesh")
    run(case, "reconstructParMesh", "-constant")
    for d in case.glob("processor*"):
        shutil.rmtree(d)
    run(case, "checkMesh", log="checkMesh.3D")

    # slab -> single-cell 2D mesh from the z = 0 patch
    run(case, "extrudeMesh")
    set_patch_types(case / "constant/polyMesh/boundary", {"front": "empty", "back": "empty"})
    run(case, "checkMesh")
    (case / f"{case.name}.foam").touch()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", type=Path)
    args = ap.parse_args()
    # bodies end at x ~ 0, like the Ahmed body (nose at x = -1.044, tail at x = 0)
    mesh(args.case.resolve(), geometry2d.rectangle(x0=-1.044), DEFAULTS)
    print(f"mesh ready: {args.case}")


if __name__ == "__main__":
    main()
