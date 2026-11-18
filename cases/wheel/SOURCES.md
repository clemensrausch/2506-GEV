# Sources: isolated rotating wheel

The wheel, its supports and the contact-patch plinth are built by `wheel_geometry.py` from dimensions read
from the authors' CAD file, so the case runs without any third-party file. The PIV measurements used for
the comparison are published under a non-commercial licence and are not part of this repository.

**Experiment:** D. Patel, A. Garmory and M. Passmore, On the wake of an isolated rotating wheel: an
experimental and numerical investigation, *Journal of Wind Engineering and Industrial Aerodynamics* 227,
2022, 105049. https://doi.org/10.1016/j.jweia.2022.105049 (CC BY 4.0)

**Data:** D. Patel, A. Garmory and M. Passmore, CAD and experiment PIV data, Loughborough University,
https://doi.org/10.17028/rd.lboro.16863232 (CC BY-NC 4.0)

**Contact patch:** S. Diasinos, T. J. Barber and G. Doig, The effects of simplifications on isolated wheel
aerodynamics, *Journal of Wind Engineering and Industrial Aerodynamics* 146, 2015, 90-101.
https://doi.org/10.1016/j.jweia.2015.08.004

**Force comparison:** J. E. Fackrell, *The Aerodynamics of an Isolated Wheel Rotating in Contact with the
Ground*, PhD thesis, Imperial College London, 1974, quoted from A. P. Mears, *The Aerodynamic
Characteristics of an Exposed Racing Car Wheel*, PhD thesis, Durham University, 2004.

## Files

All files go into `data/`, unchanged, with the names of the data set.

| File | Needed for |
|---|---|
| `Y0mm_Mean.dat`, `Y50mm_Mean.dat` | PIV comparison in the centre plane and 50 mm beside it (`compare_piv.py`) |
| `X330mm_Mean.dat`, `X495mm_Mean.dat` | PIV comparison in the cross planes 1 and 1.5 diameters behind the axle |

## What the repository contains instead

- `wheel_geometry.py` builds the geometry from the dimensions in its header.
- `results/validated/` holds the results of this setup only, and the comparison metrics derived from the
  PIV data, without the PIV data themselves.
- The figures in `figures/` show the measurements as cited comparison points.
