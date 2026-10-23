# abl_stable_kosovic_offset_surface_temp

## Purpose
Shows one of the four `ABL.surface_stability_type` options on terrain. It is the
GABLS1 set-up of `abl_stable_kosovic_offset` (Kosovic model, flat ground raised
by 50 m with `TerrainDrag`) with the terrain surface layer set by
`ABL.surface_stability_type = surface_temperature` with `ABL.surface_temp_init = 263` K. Each terrain cell solves the Monin-Obukhov surface layer for the friction velocity, Obukhov length and heat flux.

## Files
- `abl_stable_kosovic_offset_surface_temp.inp`: input deck (400 m x 400 m x 500 m, 32 x 32 x 40 cells).
- `terrain.amrwind`: flat terrain at 50 m.

## Data and acknowledgements
No external data are used. The set-up follows GABLS1 (Beare et al. 2006,
https://doi.org/10.1007/s10546-004-2820-6); see `abl_stable_kosovic/README.md`.

## How to run
```
mpirun -np 4 kynema_sgf abl_stable_kosovic_offset_surface_temp.inp
```

## Expected Results
The run completes without non-finite values. The option and its inputs are
listed in `docs/sphinx/user/inputs_ABL.rst` under `ABL.surface_stability_type`.

## CTest
`abl_stable_kosovic_offset_surface_temp` (`add_test_r`), 6 steps, compared with the gold plotfile.
