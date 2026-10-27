# abl_stable_kosovic

## Purpose
Flat-ground reference for the terrain surface-stability options. GABLS1 stable
boundary layer with the Kosovic model: the surface cools at 0.25 K/hr from 265 K
and the `ABLWallFunction` (`ABL.wall_shear_stress_type = local`) gives the
surface momentum and heat fluxes. `abl_stable_kosovic_offset` runs the same case
on ground raised by 50 m with `TerrainDrag`.

## Files
- `abl_stable_kosovic.inp`: input deck (400 m x 400 m x 400 m, 32^3 cells).

## Data and acknowledgements
No external data are used. The set-up follows the GABLS1 intercomparison:
R. J. Beare et al., "An intercomparison of large-eddy simulations of the stable
boundary layer", Boundary-Layer Meteorology 118, 247-272 (2006),
https://doi.org/10.1007/s10546-004-2820-6. We thank the GABLS participants for
the case definition.

## How to run
```
mpirun -np 4 kynema_sgf abl_stable_kosovic.inp
```
The deck stops after 10 steps for the regression test. Run the full 9 hours with
`time.max_step=-1`; plane averages are written to `post_processing/` every 100
steps.

## Expected Results
After 9 hours the flat and offset runs give the same profiles of wind speed,
direction, potential temperature and resolved heat flux; see
`abl_stable_kosovic_offset/README.md`.

## CTest
`abl_stable_kosovic` (`add_test_r`), 6 steps, compared with the gold plotfile.
