# Sources: Tyrrell 026 wing in ground effect

The wing geometry and all measured values of this case come from one PhD thesis. They may not be
redistributed, so this repository contains neither the coordinates nor the measurements. The thesis
itself is freely available; the files below have to be transcribed from it once.

**Reference:** J. D. C. Zerihan, *An Investigation into the Aerodynamics of Wings in Ground Effect*,
PhD thesis, University of Southampton, 2001. https://eprints.soton.ac.uk/426058/

Supporting paper: J. Zerihan and X. Zhang, Aerodynamics of a single element wing in ground effect,
*Journal of Aircraft* 37(6), 2000, 1058-1064. https://doi.org/10.2514/2.2711

## Files

All files go into `data/`. Format: any number of comment lines starting with `#`, then one header line
with the column names below, then one row per point, comma-separated. Page numbers are those printed in
the thesis; the PDF page is 18 higher.

| File | Needed for | Source in the thesis | Columns |
|---|---|---|---|
| `tyrrell026_table2.csv` | geometry of the single element wing (required) | Table 2, p. 25 | `x_suction,y_suction,x_pressure,y_pressure` |
| `zerihan_table3_main.csv` | double element wing, main element (required for `run_double.py`) | Table 3, p. 26 | `x_suction,y_suction,x_pressure,y_pressure` |
| `zerihan_table4_flap.csv` | double element wing, flap (required for `run_double.py`) | Table 4, p. 27 | `x_suction,y_suction,x_pressure,y_pressure` |
| `zerihan_table15_cl.csv` | comparison plot of `run_sweep.py` (optional) | Table 15, p. 179 | `h_over_c,Cl_expt,Cl_SA,Cl_kw` |
| `zerihan_table9_xsep.csv` | comparison plot of `eval_separation.py` (optional) | Table 9, p. 56 | `h_over_c,x_sep_over_c` |
| `zerihan_fig36b_cl_datum.csv` | comparison plot of `run_double.py` (optional) | Fig. 36b (PDF p. 135), digitised | `h_over_c,CL` |

Notes on the transcription:

- Tables 2 to 4 give both surfaces from the leading to the trailing edge, normalised by the chord, at
  the reference incidence of 1° and with the wing inverted (suction surface `y < 0`). Enter them as
  printed; `wing.py` handles incidence, ground height and scale. Table 2 uses the chord of 223.4 mm,
  Tables 3 and 4 the combined chord of 380 mm of the assembled wing at the datum flap deflection.
- The scan has OCR errors. After transcribing, plot the profile and check that both surfaces are smooth;
  a single wrong digit shows up as a kink.
- In Table 9 the column with `h/c` is missing from the OCR text; read it from the scan.
- Fig. 36b has no table. Digitise the datum-flap curve, for example with WebPlotDigitizer.

## What the repository contains instead

- `wing.py` builds the wing from the transcribed tables.
- `results/validated/` holds the results of this setup only, without measured values.
- The figures in `figures/` show the measurements as cited comparison points.
