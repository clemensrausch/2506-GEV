"""Tyrrell 026 wings of Zerihan (2001), inverted and positioned as in the experiment.

The profile coordinates are not part of this repository (see SOURCES.md). Put them into data/:
    data/tyrrell026_table2.csv     single element, Table 2 of the thesis
    data/zerihan_table3_main.csv   double element, main element, Table 3
    data/zerihan_table4_flap.csv   double element, flap, Table 4
"""

from pathlib import Path

import numpy as np

import _path  # noqa: F401  (repository root on sys.path)
from gev import geometry2d
from gev.tables import read_csv

DATA = Path(__file__).resolve().parent / "data"
SINGLE = DATA / "tyrrell026_table2.csv"
DOUBLE_MAIN = DATA / "zerihan_table3_main.csv"
DOUBLE_FLAP = DATA / "zerihan_table4_flap.csv"
CHORD = 0.2234          # single element, m
CHORD_DOUBLE = 0.380    # double element, combined chord, m


def _table(path):
    """(suction, pressure) as (N, 2) arrays from a file with columns x_suction, y_suction, x_pressure, y_pressure."""
    if not path.exists():
        raise FileNotFoundError(f"{path} missing: transcribe it from the thesis as described in SOURCES.md")
    t = read_csv(path)
    return np.c_[t["x_suction"], t["y_suction"]], np.c_[t["x_pressure"], t["y_pressure"]]


def single(h_over_c=0.179, alpha=1.0, chord=CHORD, x_le=-CHORD, n_side=200):
    """Single element wing. Table 2 gives the profile at the reference incidence alpha = 1 deg (true incidence
    3.45 deg; quoted alpha = true - 2.45 deg). h is the distance from the lowest point of the suction surface
    to the ground at alpha = 0; the incidence is set about the quarter chord. Blunt trailing edge 0.007c."""
    suction, pressure = _table(SINGLE)
    return geometry2d.single_element(suction, pressure, h_over_c, alpha, chord, x_le, n_side, alpha_table=1.0,
                                     params=dict(alpha_true=alpha + 2.45))


def double(h_over_c=0.263, alpha=1.0, chord=CHORD_DOUBLE, x_le=-CHORD_DOUBLE, n_side=200, n_flap=160):
    """Double element wing (main element = Table 2 scaled by 0.588, flap), datum flap deflection. Tables 3 and 4
    give both elements assembled as tested at alpha = 1 deg, normalised by the combined chord of 380 mm."""
    main, flap = _table(DOUBLE_MAIN), _table(DOUBLE_FLAP)
    return geometry2d.multi_element({"main": (*main, n_side), "flap": (*flap, n_flap)}, h_over_c, alpha, chord, x_le,
                                    alpha_table=1.0)
