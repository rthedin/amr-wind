# klaxell_stability_abort

## Purpose
`KLAxell` takes the friction velocity from `ABL.monin_obukhov_length` and the TKE
boundary condition from `ABL.surface_temp_flux`, both from the same 1-D RANS
solution. It therefore supports only `ABL.surface_stability_type =
monin_obukhov_length`. This deck sets `heat_flux` and must abort.

## Files
- `klaxell_stability_abort.inp`: KLAxell on flat terrain raised to 20 m.
- `terrain.amrwind`, `rans_1d.info`: copied from `forest_drag_roughness`.

## Data and acknowledgements
No external data are used.

## How to run
```
kynema_sgf klaxell_stability_abort.inp
```

## Expected Results
The run aborts with
`KLAxell-ConstTransport uses ABL.monin_obukhov_length ... ABL.surface_stability_type = heat_flux is not supported`.

## CTest
`klaxell_stability_abort` passes when the output matches that message
(`PASS_REGULAR_EXPRESSION`). The message names `KLAxell`, so the test fails if
only the `KransAxell` check is left.
