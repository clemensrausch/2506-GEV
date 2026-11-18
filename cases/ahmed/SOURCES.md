# Sources: Ahmed body

The body is built by `ahmed_body.py` from the dimensions of the ERCOFTAC drawing. Only the nose is taken
from a measured surface, and that file, like the measured wall pressures, is not part of this repository.
Both are free to download from the ERCOFTAC database.

**Data set:** H. Lienhart, S. Becker and C. Stoots, *Flow Around a Simplified Car Body (Ahmed Body)*,
ERCOFTAC Classic Collection, Case 82.
http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case082

**Experiment:** H. Lienhart and S. Becker, Flow and turbulence structure in the wake of a simplified car
model, SAE Technical Paper 2003-01-0656, 2003. https://doi.org/10.4271/2003-01-0656

**Drag values:**
- S. R. Ahmed, G. Ramm and G. Faltin, Some salient features of the time-averaged ground vehicle wake,
  SAE Technical Paper 840300, 1984. https://doi.org/10.4271/840300
- S. Möller, D. Suzzi and W. Meile, Investigation of the flow around the Ahmed body using RANS and URANS
  with various turbulence models, 3rd OpenFOAM Workshop, Politecnico di Milano, 2008 (slides with the
  drag measured by W. Meile at TU Graz, 40 m/s).
- W. Meile, G. Brenn, A. Reppenhagen, B. Lechner and A. Fuchs, Experiments and numerical simulations on
  the aerodynamics of the Ahmed body, *CFD Letters* 3(1), 2011, 32-39.

## Files

All files go into `data/`, unchanged, with the names they have in the ERCOFTAC archive.

| File | Needed for |
|---|---|
| `ahmed-front-geo.dat` | geometry of the nose (required) |
| `ahmed-25-press.dat`, `ahmed-35-press.dat` | wall pressure comparison of `compare_pressure.py` (optional) |

The nose file is linked on the case page. The pressure files are in the archive with all files of the
case (`ahby-allfiles.zip`).

## What the repository contains instead

- `ahmed_body.py` builds the closed body and the stilts from the drawing and the nose file.
- `results/validated/` holds the results of this setup only, without measured values.
- The figures in `figures/` show the measurements as cited comparison points.
