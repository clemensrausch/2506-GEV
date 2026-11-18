# Sources: diffuser bluff body

The body is built by `diffuser_body.py` from the dimensions given in the thesis below, so the case runs
without any third-party file. Where the thesis gives no number, the script states the assumption it
makes. Only the measured forces, used for the comparison plot, are not part of this repository.

**Reference:** A. E. Senior, *The Aerodynamics of a Diffuser Equipped Bluff Body in Ground Effect*,
PhD thesis, University of Southampton, 2002. https://eprints.soton.ac.uk/47113/

Supporting paper on the hysteresis near the force maximum: A. Ruhrmann and X. Zhang, Influence of
diffuser angle on a bluff body in ground effect, *Journal of Fluids Engineering* 125(2), 2003, 332-338. https://doi.org/10.1115/1.1537252

## Files

| File | Needed for | Source in the thesis | Columns |
|---|---|---|---|
| `senior2002_fig9.csv` | comparison plot of `run_diffuser.py plot` (optional) | Fig. 9, p. 45, digitised | `hr_d,hr_m,CL,CD` (a further text column `regime` is ignored) |

Format: comment lines starting with `#`, one header line, comma-separated rows. `hr_d` is the ride height
over the half width `d` = 157 mm, `hr_m` the ride height in metres. `CL` keeps Senior's sign (negative is
downforce).

## What the repository contains instead

- `diffuser_body.py` builds the model from the dimensions of Senior's Fig. 4.
- `results/validated/` holds the results of this setup only, without measured values.
- The figures in `figures/` show the measurements as cited comparison points.
