"""Reading the third-party CSV files of the cases: any number of '#' comment lines, one header line with column
names, numeric rows (non-numeric entries become NaN)."""

from pathlib import Path

import numpy as np


def read_csv(path):
    """{column name: array} of a CSV file as described in the module doc."""
    lines = [l for l in Path(path).read_text().splitlines() if l.strip() and not l.lstrip().startswith("#")]
    names = [n.strip() for n in lines[0].split(",")]
    cols = {n: [] for n in names}
    for line in lines[1:]:
        vals = [v.strip() for v in line.split(",")]
        for n, v in zip(names, vals + [""] * (len(names) - len(vals))):
            try:
                cols[n].append(float(v))
            except ValueError:
                cols[n].append(np.nan)
    return {n: np.array(v) for n, v in cols.items()}
