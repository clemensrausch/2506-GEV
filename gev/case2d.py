"""Mesher-independent dictionaries of a 2D case: decomposition, forces, flow values."""

from .foam import FOOTER, write


def write_case_dicts(case, body, p):
    sysd = case / "system"
    dz = p["span"]
    write(sysd / "decomposeParDict", "decomposeParDict", f"""\
numberOfSubdomains {p['n_procs']};
method          scotch;
""")

    # force coefficients; reference sizes come from the body (2D: per slab thickness dz)
    ref_length, ref_area = body.params["ref_length"], body.params["ref_area_per_span"] * dz
    cofr = body.params.get("quarter_chord", (0.0, 0.0))   # moments about c/4 for a wing
    write(sysd / "forceCoeffs", "forceCoeffs", f"""\
forceCoeffs_main
{{
    type            forceCoeffs;
    libs            (forces);
    patches         ("body_.*");
    rho             rhoInf;
    rhoInf          1.225;
    dragDir         (1 0 0);
    liftDir         (0 1 0);
    pitchAxis       (0 0 1);
    magUInf         {p['U_inf']:g};
    lRef            {ref_length:g};
    Aref            {ref_area:g};
    CofR            ({cofr[0]:g} {cofr[1]:g} 0);
    log             yes;
    writeControl    timeStep;
    writeInterval   2;
}}
""")
    # dimensional forces and moments per component (force_groups: name -> patch regexes), moments about
    # body.params["cofr_forces"]; per unit span = value / slab thickness dz
    groups = body.params.get("force_groups", {})
    if groups:
        cf = body.params.get("cofr_forces", (0.0, 0.0))
        blocks = "".join(f"""
forces_{name}
{{
    type            forces;
    libs            (forces);
    patches         ({" ".join(f'"{r}"' for r in regexes)});
    rho             rhoInf;
    rhoInf          1.225;
    CofR            ({cf[0]:g} {cf[1]:g} 0);
    writeControl    timeStep;
    writeInterval   1;
}}
""" for name, regexes in groups.items())
        f = sysd / "forceCoeffs"
        f.write_text(f.read_text().replace(FOOTER, blocks + FOOTER))

    ug = p["U_inf"] if p["moving_ground"] else 0.0
    write(case / "0/initialConditions", "initialConditions", f"""\
flowVelocity         ({p['U_inf']:g} 0 0);
pressure             0;
turbulentKE          {p['k_inf']:.6g};
turbulentOmega       {p['omega_inf']:.6g};
groundVelocity       ({ug:g} 0 0);
""")
