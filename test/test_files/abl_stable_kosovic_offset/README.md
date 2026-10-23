# abl_stable_kosovic_offset

## Purpose
Checks that the terrain surface layer with `ABL.surface_stability_type =
heating_rate` reproduces the flat-wall `ABLWallFunction`. It is the GABLS1
stable boundary layer of `abl_stable_kosovic` (Kosovic model, surface cooling of
0.25 K/hr from 265 K) on flat ground raised by 50 m with `TerrainDrag`.
`DragForcing` and `DragTempForcing` solve the Monin-Obukhov surface layer at each
terrain cell from the cell above it. Heights in the deck are the flat-case
heights plus 50 m.

## Files
- `abl_stable_kosovic_offset.inp`: input deck (400 m x 400 m x 500 m, 32 x 32 x 40
  cells, dz = 12.5 m; the terrain top is a cell face).
- `terrain.amrwind`: flat terrain at 50 m.
- `compare_gabls_offset.py`: averages the ABLStats plane averages of the flat and
  offset runs over a time window, shifts the offset heights by 50 m, prints the
  largest differences and plots wind speed, direction, potential temperature and
  resolved heat flux.

## Data and acknowledgements
No external data are used. The set-up follows GABLS1: R. J. Beare et al., "An
intercomparison of large-eddy simulations of the stable boundary layer",
Boundary-Layer Meteorology 118, 247-272 (2006),
https://doi.org/10.1007/s10546-004-2820-6. We thank the GABLS participants for
the case definition.

## How to run
```
cd flat   && mpirun -np 4 kynema_sgf abl_stable_kosovic.inp time.max_step=-1
cd offset && mpirun -np 4 kynema_sgf abl_stable_kosovic_offset.inp time.max_step=-1
python3 compare_gabls_offset.py flat offset --offset 50 --tstart 28800 --tend 32400
```
Each 9-hour run takes about 45 minutes on 4 CPU ranks.

## Expected Results
Mean over 8-9 h (32 x 32 x 32 cells for the flat run, the same cell size for the
offset run):

| Quantity                                  | Flat (ABLWallFunction) | Offset 50 m (heating_rate) |
|-------------------------------------------|------------------------|----------------------------|
| Wind speed at the first cell (6.25 m)     | 3.01 m/s               | 3.40 m/s                   |
| Wind direction at the first cell          | 33.4 deg               | 29.3 deg                   |
| Potential temperature at the first cell   | 263.42 K               | 263.26 K                   |
| Low-level jet                             | 9.40 m/s at 194 m      | 9.22 m/s at 169 m          |
| Height where the resolved heat flux is 5% of its peak | 194 m      | 169 m                      |

Largest differences over the profile: 0.79 m/s in wind speed and 0.71 K in
potential temperature, both at the top of the boundary layer (206 m); 4.1 deg
in direction near the surface; 1e-3 K m/s in resolved heat flux.

The wind-speed difference at 206 m comes from the height of the jet nose
(169 m offset, 194 m flat); the jet speeds differ by 0.18 m/s. The offset
boundary layer is one to two cells shallower from the first hour on, and the
square root of the total (resolved plus SGS) stress at 18.8 m is 4-5% lower
(0.231 against 0.241 m/s). This is the known effect of the binary (blanked)
terrain representation, found earlier for the form7 cases in kynema-sgf; it is
not caused by the surface-stability option. The flat run also uses the
plane-averaged surface state of `ABLWallFunction` and the offset run the
per-cell state, so the two are not expected to agree exactly.

## CTest
`abl_stable_kosovic_offset` (`add_test_r`), 6 steps, compared with the gold
plotfile.
