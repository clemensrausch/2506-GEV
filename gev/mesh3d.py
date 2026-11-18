"""Parts (closed STL surfaces) -> 3D OpenFOAM case: blockMesh + snappyHexMesh, decomposed, ready for run3d.

Physics and meshing follow the 2D setup (mesh2d_snappy, validated against the Tyrrell 026 wing) and the Ahmed
reference case: kOmegaSST, same schemes, same snappy snap and quality settings, Ahmed refinement boxes scaled
with the body. In 3D:
- half model with a symmetry plane at z = 0 (default; full model with half=False);
- background cells are cubes;
- prism layers per part: "resolved" (y+ ~ 1, absolute first layer from a flat-plate estimate) or "wallfn"
  (wall functions, layers relative to the surface cells). The wall BCs (nutUSpalding, omega with binomial
  blending) are valid for both;
- one STL per part; the part name becomes the patch name, all parts are in the patch group "car";
- the mesh stays decomposed (no reconstructParMesh): the solver runs on the processor directories.

Coordinates as in 2D: x = flow, y = up (ground at y = 0), z = span (symmetry plane z = 0). Units m.

Usage from Python:
    parts = [Part("body", stl, level=(4, 5), layers="resolved"), Part("stilts", stl2, layers="wallfn")]
    mesh(case, parts, AHMED)
Requires a sourced OpenFOAM v2506 environment.
"""

import shutil
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import trimesh

from . import TEMPLATES
from .foam import mpirun, run, write


@dataclass
class Part:
    """One closed (or, together with other parts, closing) STL surface = one wall patch."""
    name: str
    stl: Path
    level: tuple = None          # snappy surface level (min, max); None: profile surface_level
    feature_level: int = None    # level on the feature edges; None: profile feature_level
    layers: str = "wallfn"       # "resolved" (y+ ~ 1), "wallfn" or None (no layers)
    bc: dict = field(default_factory=dict)   # optional per-part BC overrides, {"U": "type ...; ..."}


# Profile: everything that is not geometry (the settings of the validation cases, see their READMEs).
# Box and domain sizes are in units of the body length.
DEFAULTS = dict(
    half=True,
    # domain around the body bounding box, as the Ahmed reference case (x -4 ... 10 m, roof and side 4 m:
    # 2.8 L ahead of the nose, 9.6 L behind the base, 3.5 L above, 3.6 L beside the body; blockage 0.35 %)
    upstream=3.0, downstream=9.6, roof=3.5, side=3.6,
    tunnel=None,                  # (height, half width) in m: wind-tunnel cross-section with slip walls instead of
                                  # roof/side (blockage as in the experiment); cells then only nearly cubic
    cell=0.08,                    # background cell (cube), m: Ahmed reference cell size in x
    surface_level=(4, 5), feature_level=4,
    near_level=3, wake_level=2,   # Ahmed refinementBox / wakeBox
    extra_boxes=None,             # f(bbox, L) -> {name: ((x0, y0, z0), (x1, y1, z1), level)}
    # layers
    y_plus=1.0, layer_expansion=1.2, layer_final=0.3,     # resolved: last layer = layer_final x surface cell
    wallfn_layers=4, wallfn_final=0.4, wallfn_expansion=1.2,
    ground_layers="wallfn",
    layer_relaxed="maxNonOrtho 75;\n        minTetQuality -1e30;",   # quality limits after nRelaxedIter (20)
                                  # layer iterations; first Ahmed meshes lost layers at edges on these two
    layer_medial_ratio=1.0,       # maxThicknessToMedialRatio; 0.3 (Ahmed reference) thinned the layers on the
                                  # slant, the base and around edges of the first Ahmed meshes
    layer_outer_iter=1,           # nOuterIter: layers in one pass. With 3 passes (outer layers first, then the
                                  # thin ones) the pass of the thinnest layers failed the quality checks on the
                                  # diffuser body and the first cell was 66 instead of 18 um (2026-10-19). Ahmed and
                                  # the wheel were validated with 3 and keep it in their profiles.
    layer_feature_angle=180,      # layers wrap around convex edges; 60 (Ahmed reference) stops them at every
                                  # 90 deg edge and thins them over ~0.1 m (first Ahmed mesh, 2026-10-17)
    # flow
    U_inf=40.0, nu=1.5e-5, Tu=0.0025, nut_ratio=10.0, moving_ground=False,
    ref_area=None,                # full-car frontal area, m^2 (coefficients refer to the whole car)
    ref_length=None,              # default: body length
    cofr=None,                    # moment centre (x, y, z); default: body centre on the ground
    # numerics
    solver="fast",                # "fast" (motorBike tutorial) or "strict" (settings of the validated 2D cases)
    n_procs=8,
    mpi_args=(),
)

# Ahmed body (ERCOFTAC case 82, Lienhart & Becker 2003): 40 m/s, fixed floor, Tu < 0.25 %.
AHMED = dict(DEFAULTS, U_inf=40.0, Tu=0.0025, moving_ground=False, ref_area=0.389 * 0.288,
             layer_outer_iter=3)          # as validated (wall functions; the pass count hardly matters there)


def first_layer(p, length):
    """First-cell height for the target y+ from the turbulent flat-plate skin friction at x = length
    (Cf = 0.026 Re^-1/7, White)."""
    re = p["U_inf"] * length / p["nu"]
    u_tau = p["U_inf"] * np.sqrt(0.5 * 0.026 * re ** (-1 / 7))
    return p["y_plus"] * p["nu"] / u_tau


def bbox(parts):
    lo, hi = np.full(3, np.inf), np.full(3, -np.inf)
    for part in parts:
        m = trimesh.load(part.stl, force="mesh")
        lo, hi = np.minimum(lo, m.bounds[0]), np.maximum(hi, m.bounds[1])
    return lo, hi


def layer_spec(kind, p, level, first):
    """snappy per-patch layer entries."""
    if kind is None:
        return "nSurfaceLayers 0;"
    cell = p["cell"] / 2 ** level
    if kind == "resolved":
        final = p["layer_final"] * cell
        n = int(np.ceil(np.log(final / first) / np.log(p["layer_expansion"])))
        return (f"nSurfaceLayers {n};\n            relativeSizes false;\n"
                f"            thicknessModel firstAndExpansion;\n"
                f"            firstLayerThickness {first:.3g};\n            expansionRatio {p['layer_expansion']:g};\n"
                f"            minThickness {first:.3g};      // = one first layer: keep thin partial stacks")
    return (f"nSurfaceLayers {p['wallfn_layers']};\n            relativeSizes true;\n"
            f"            thicknessModel finalAndExpansion;\n"
            f"            finalLayerThickness {p['wallfn_final']:g};\n            expansionRatio {p['wallfn_expansion']:g};\n"
            f"            minThickness 0.1;")


def domain(lo, hi, p):
    """Domain around the body bounding box, snapped to whole background cells:
    (x0, x1, y1, z0, z1) and background cell counts (nx, ny, nz)."""
    Lb = hi[0] - lo[0]
    h = p["cell"]
    x0 = lo[0] - p["upstream"] * Lb
    nx = int(np.ceil((p["downstream"] * Lb + Lb + p["upstream"] * Lb) / h))
    ny = int(np.ceil((hi[1] + p["roof"] * Lb) / h))
    nz = int(np.ceil((max(abs(lo[2]), abs(hi[2])) + p["side"] * Lb) / h))
    x1, y1, z1 = x0 + nx * h, ny * h, nz * h
    if p.get("tunnel"):
        y1, z1 = p["tunnel"]
        ny, nz = max(1, round(y1 / h)), max(1, round(z1 / h))
    z0 = 0.0 if p["half"] else -z1
    return (x0, x1, y1, z0, z1), (nx, ny, nz if p["half"] else 2 * nz)


def refinement_boxes(lo, hi, p):
    """{name: (min, max, level)}: the Ahmed reference boxes relative to the body (L = 1.044, top 0.338,
    half width 0.1945), plus the profile's extra boxes."""
    Lb = hi[0] - lo[0]
    zb = max(abs(lo[2]), abs(hi[2]))
    def box(xa, xb, ytop, zw):
        za = -0.01 * p["cell"] if p["half"] else -zw
        return (xa, -0.015, za), (xb, ytop, zw)
    boxes = {
        "nearBox": (*box(lo[0] - 0.533 * Lb, hi[0] + 0.192 * Lb, hi[1] + 0.155 * Lb, zb + 0.245 * Lb),
                    p["near_level"]),
        "wakeBox": (*box(hi[0] - 0.335 * Lb, hi[0] + 2.586 * Lb, hi[1] + 0.634 * Lb, zb + 0.532 * Lb),
                    p["wake_level"]),
    }
    if p["extra_boxes"]:
        boxes.update(p["extra_boxes"]((lo, hi), Lb))
    return boxes


def write_dicts(case, parts, p):
    sysd = case / "system"
    lo, hi = bbox(parts)
    Lb = hi[0] - lo[0]
    L = p["ref_length"] or Lb
    h = p["cell"]
    (x0, x1, y1, z0, z1), (nx, ny, nzt) = domain(lo, hi, p)
    sym = ("    symmetry { type symmetry; faces ((0 3 2 1)); }\n" if p["half"]
           else "")
    side_faces = "(4 5 6 7)" if p["half"] else "(4 5 6 7) (0 3 2 1)"

    write(sysd / "blockMeshDict", "blockMeshDict", f"""\
scale   1;

vertices
(
    ({x0:g} 0 {z0:g}) ({x1:g} 0 {z0:g}) ({x1:g} {y1:g} {z0:g}) ({x0:g} {y1:g} {z0:g})
    ({x0:g} 0 {z1:g}) ({x1:g} 0 {z1:g}) ({x1:g} {y1:g} {z1:g}) ({x0:g} {y1:g} {z1:g})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nzt}) simpleGrading (1 1 1)
);

boundary
(
    inlet  {{ type patch; faces ((0 4 7 3)); }}
    outlet {{ type patch; faces ((1 2 6 5)); }}
    top    {{ type patch; faces ((3 7 6 2)); }}
    side   {{ type patch; faces ({side_faces}); }}
    ground {{ type wall;  faces ((0 1 5 4)); }}
{sym});
""")

    tri = case / "constant/triSurface"
    tri.mkdir(parents=True, exist_ok=True)
    for part in parts:
        shutil.copy(part.stl, tri / f"{part.name}.stl")

    write(sysd / "surfaceFeatureExtractDict", "surfaceFeatureExtractDict", "".join(f"""\
{part.name}.stl
{{
    extractionMethod    extractFromSurface;
    includedAngle       150;
    geometricTestOnly   yes;    // 2D lesson: no feature lines at region borders on a smooth surface
    subsetFeatures
    {{
        nonManifoldEdges    no;
        openEdges           yes;
    }}
    writeObj            yes;
}}

""" for part in parts))

    boxes = refinement_boxes(lo, hi, p)
    geo_boxes = "".join(f"""    {n}
    {{
        type searchableBox;
        min  ({a[0]:g} {a[1]:g} {a[2]:g});
        max  ({b[0]:g} {b[1]:g} {b[2]:g});
    }}
""" for n, (a, b, _) in boxes.items())
    reg_boxes = "".join(f"""        {n}
        {{
            mode inside;
            levels ((1E15 {lv}));
        }}
""" for n, (_, _, lv) in boxes.items())

    geo_parts = "".join(f"""    {part.name}.stl
    {{
        type triSurfaceMesh;
        name {part.name};
    }}
""" for part in parts)
    features = "".join(f"""        {{
            file "{part.name}.eMesh";
            level {part.feature_level if part.feature_level is not None else p['feature_level']};
        }}
""" for part in parts)
    surfaces = ""
    for part in parts:
        a, b = part.level or p["surface_level"]
        surfaces += f"""        {part.name}
        {{
            level ({a} {b});
            patchInfo
            {{
                type wall;
                inGroups (car);
            }}
        }}
"""

    first = first_layer(p, L)
    layers = ""
    for part in parts:
        lv = (part.level or p["surface_level"])[0]
        layers += f"""        {part.name}
        {{
            {layer_spec(part.layers, p, lv, first)}
        }}
"""
    layers += f"""        ground
        {{
            {layer_spec(p['ground_layers'], p, p['near_level'], first)}
        }}
"""
    resolved = any(part.layers == "resolved" for part in parts) or p["ground_layers"] == "resolved"
    det = "    minDeterminant 1e-8;     // 2D lesson: thin absolute layers on coarse cells" if resolved else ""
    if p["layer_relaxed"]:
        det += f"\n    relaxed\n    {{\n        {p['layer_relaxed']}\n    }}"

    write(sysd / "snappyHexMeshDict", "snappyHexMeshDict", f"""\
castellatedMesh true;
snap            true;
addLayers       true;

geometry
{{
{geo_parts}{geo_boxes}}}

castellatedMeshControls
{{
    maxLocalCells 4000000;
    maxGlobalCells 40000000;
    minRefinementCells 10;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 3;

    features
    (
{features}    );

    refinementSurfaces
    {{
{surfaces}    }}

    resolveFeatureAngle 30;

    refinementRegions
    {{
{reg_boxes}    }}

    locationInMesh ({x1 - 0.05 * (x1 - x0) + 0.137 * h:g} {0.5625 * y1 + 0.371 * h:g} {(z0 + z1) / 2 + 0.4321 * h if not p['half'] else 0.5 * z1 + 0.293 * h:g});   // off the background cell faces
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
    // defaults; every patch with layers has its own entry below
    relativeSizes true;
    expansionRatio 1.0;
    finalLayerThickness 0.3;
    minThickness 0.1;

    layers
    {{
{layers}    }}

    nGrow 0;
    featureAngle {p['layer_feature_angle']:g};
    slipFeatureAngle 30;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio {p['layer_medial_ratio']:g};
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
    nOuterIter {p['layer_outer_iter']};
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

    write(sysd / "decomposeParDict", "decomposeParDict", f"""\
numberOfSubdomains {p['n_procs']};
method          scotch;
""")

    write_fvsolution(sysd, p)
    write_forces(sysd, parts, p, (lo, hi), L)
    write_initial(case, p)
    write_part_bcs(case, parts)
    return dict(bbox=(lo, hi), first_layer=first, cells_background=nx * ny * nzt)


SOLVERS = {
    # settings of the validated 2D cases and the Ahmed reference: p solved to 1e-7 every iteration, two
    # non-orthogonal correctors, SIMPLEC with p relaxed to 0.3
    "strict": dict(p_relTol=0, U_relTol=0, kw_relTol=0.01, nonOrth=2, p_relax=0.3, U_relax=0.8, kw_relax=0.8),
    # tutorial incompressible/simpleFoam/motorBike (v2506): SIMPLEC without p relaxation, inner solvers to relTol
    "fast": dict(p_relTol=0.01, U_relTol=0.1, kw_relTol=0.1, nonOrth=0, p_relax=None, U_relax=0.9, kw_relax=0.7),
}


def write_fvsolution(sysd, p):
    s = SOLVERS[p["solver"]]
    p_relax = f"    fields\n    {{\n        p               {s['p_relax']};\n    }}\n" if s["p_relax"] else ""
    write(sysd / "fvSolution", "fvSolution", f"""\
// solver preset "{p['solver']}" (gev.mesh3d.SOLVERS)
solvers
{{
    p
    {{
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          {s['p_relTol']};
    }}

    Phi
    {{
        $p;
    }}

    U
    {{
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          {s['U_relTol']};
        nSweeps         1;
    }}

    "(k|omega)"
    {{
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          {s['kw_relTol']};
        nSweeps         1;
    }}
}}

SIMPLE
{{
    nNonOrthogonalCorrectors {s['nonOrth']};
    consistent      yes;

    // no residual stop: the force-based rules of gev.convergence decide (plateau, or settled windows after 1200
    // iterations). With p 1e-4 the diffuser runs stopped at 640-840 iterations while C_L still drifted
    // by 0.5-1.4 % over the last window (2026-10-19).
    residualControl
    {{
        p               0;
        U               0;
        "(k|omega)"     0;
    }}
}}

potentialFlow
{{
    nNonOrthogonalCorrectors 10;
}}

relaxationFactors
{{
{p_relax}    equations
    {{
        U               {s['U_relax']};
        "(k|omega)"     {s['kw_relax']};
    }}
}}

cache
{{
    grad(U);
}}
""")


def write_forces(sysd, parts, p, bb, L):
    """Forces per part (dimensional, of the modelled half in a half model) and coefficients of the whole car."""
    lo, hi = bb
    cofr = p["cofr"] or (0.5 * (lo[0] + hi[0]), 0.0, 0.0)
    area = p["ref_area"] * (0.5 if p["half"] else 1.0)   # half model: half the force on half the area
    blocks = f"""forceCoeffs_car
{{
    type            forceCoeffs;
    libs            (forces);
    patches         (car);
    rho             rhoInf;
    rhoInf          1.225;
    dragDir         (1 0 0);
    liftDir         (0 1 0);        // downforce comes out negative
    pitchAxis       (0 0 1);
    magUInf         {p['U_inf']:g};
    lRef            {L:g};
    Aref            {area:g};
    CofR            ({cofr[0]:g} {cofr[1]:g} {cofr[2]:g});
    log             yes;
    writeControl    timeStep;
    writeInterval   1;
}}
"""
    for part in parts:
        blocks += f"""
forces_{part.name}
{{
    type            forces;
    libs            (forces);
    patches         ({part.name});
    rho             rhoInf;
    rhoInf          1.225;
    CofR            ({cofr[0]:g} {cofr[1]:g} {cofr[2]:g});
    writeControl    timeStep;
    writeInterval   1;
}}
"""
    write(sysd / "forces", "forces", blocks)


def write_initial(case, p):
    k = 1.5 * (p["Tu"] * p["U_inf"]) ** 2
    omega = k / (p["nut_ratio"] * p["nu"])
    ug = p["U_inf"] if p["moving_ground"] else 0.0
    write(case / "0/initialConditions", "initialConditions", f"""\
flowVelocity         ({p['U_inf']:g} 0 0);
pressure             0;
turbulentKE          {k:.6g};
turbulentOmega       {omega:.6g};
groundVelocity       ({ug:g} 0 0);
""")
    (case / "constant/transportProperties").write_text(
        (case / "constant/transportProperties").read_text().replace("1.5e-05", f"{p['nu']:g}"))


def write_part_bcs(case, parts):
    text = "// per-part velocity overrides, written by mesh3d.py\n"
    for part in parts:
        if "U" in part.bc:
            text += f"{part.name}\n{{\n    {part.bc['U']}\n}}\n"
    (case / "0/partBCs.U").write_text(text)


def mesh(case, parts, p=DEFAULTS, dry_run=False):
    """Write the case and mesh it in parallel. The mesh stays decomposed. dry_run: only write the
    dictionaries and check them (blockMesh, surfaceFeatureExtract, snappyHexMesh -dry-run)."""
    case = Path(case)
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(TEMPLATES / "3d", case)
    info = write_dicts(case, parts, p)

    run(case, "blockMesh")
    run(case, "surfaceFeatureExtract")
    if dry_run:
        run(case, "snappyHexMesh", "-dry-run", log="snappyHexMesh.dryRun")
        return info
    run(case, "decomposePar", "-force", "-no-fields")
    mpirun(case, p, "snappyHexMesh", "-overwrite", log="snappyHexMesh")
    mpirun(case, p, "checkMesh", "-constant", log="checkMesh")
    for proc in case.glob("processor*"):           # fields onto the decomposed mesh
        shutil.copytree(case / "0", proc / "0", dirs_exist_ok=True)
    (case / f"{case.name}.foam").touch()
    return info
