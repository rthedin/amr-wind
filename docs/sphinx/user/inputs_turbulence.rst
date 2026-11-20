.. _inputs_turbulence:

Section: turbulence
~~~~~~~~~~~~~~~~~~~

This section is for setting turbulence model parameters

.. input_param:: turbulence.model

   **type:** String, optional, default = Laminar

   Specifies which turbulence model to use, by default "Laminar" is
   chosen (effectively no turbulence model).  Currently the supported
   turbulence models are "Smagorinsky", "AMD", "Kosovic", 
   "OneEqKsgsM84", "KOmegaSST", "KOmegaSSTIDDES" or "KLAxell".

   
.. input_param:: Smagorinsky_coeffs.Cs

   **type:** Real, optional, default = 0.135

   Specifies the coefficient used in the `Smagorinsky` turbulence model.

The `Kosovic` model adds the nonlinear backscatter stress of Kosovic (1997) to
a Smagorinsky eddy viscosity. The nonlinear stress acts through the
`NonLinearSGSTerm` momentum source, which must be listed in
`ICNS.source_terms`. With `C_1 = C_2` the stress is

.. math::

   \tau^{NL}_{ij} = -(C_s \Delta)^2 C_1
   \left( S_{ik} S_{kj} + S_{ik} R_{kj} - R_{ik} S_{kj} \right),

where :math:`S` and :math:`R` are the symmetric and antisymmetric parts of
the velocity gradient :math:`\partial u_i / \partial x_j`, and the isotropic
part is absorbed in the pressure. The model also reads
:input_param:`ABL.surface_roughness_z0`, :input_param:`ABL.kappa`,
`ABL.mo_beta_m`, `ABL.mo_gamma_m` and :input_param:`ABL.surface_stability_type`.

.. input_param:: Kosovic.Cb

   **type:** Real, optional, default = 0.36

   Backscatter coefficient. It sets :math:`C_s = \sqrt{8 (1 + C_b) / (27 \pi^2)}`
   and :math:`C_1 = C_2 = \sqrt{960}\, C_b / (7 (1 + C_b) S_k)` with
   :math:`S_k = 0.5`. Must be >= 0; 0 turns the nonlinear stress off.

.. input_param:: Kosovic_coeffs.Cs

   **type:** Real, optional, default = value set by `Kosovic.Cb`

   Replaces the Smagorinsky coefficient computed from `Kosovic.Cb`. It does
   not change :math:`C_1`. Must be > 0.

.. input_param:: Kosovic.surfaceRANS

   **type:** Boolean, optional, default = false

   Blends the LES length scale with the near-wall RANS length scale
   :math:`(\kappa z / \phi_m(z/L))^2`, where :math:`z` is the height of the
   top face of the cell above the ground (the lower domain face, or the
   terrain with `TerrainDrag`).

.. input_param:: Kosovic.switchLoc

   **type:** Real, optional, default = 24.0

   Height (m) of the RANS-LES blend, :math:`f = \exp(-z/z_{switch})`. Read
   only with `Kosovic.surfaceRANS = true`. Must be > 0.

.. input_param:: Kosovic.surfaceRANSExp

   **type:** Real, optional, default = 2.0

   Exponent :math:`n` of the blend, :math:`(1-f)^n` for LES and :math:`f^n`
   for RANS. Read only with `Kosovic.surfaceRANS = true`. Must be > 0.

.. input_param:: Kosovic.LESTurbOff_transition_zaglstart

   **type:** Real, optional, default = no limit

   Height above the ground (m) where the model starts to be tapered off with
   a cosine-squared function. Must be <= `Kosovic.LESTurbOff_transition_zaglend`.

.. input_param:: Kosovic.LESTurbOff_transition_zaglend

   **type:** Real, optional, default = no limit

   Height above the ground (m) above which the eddy viscosity and the
   nonlinear stress are zero. Equal start and end heights give a step.

.. input_param:: Kosovic.muCoeff

   **type:** Real, optional, default = 3.0

   Ratio of the turbulent thermal diffusivity to the eddy viscosity (the
   inverse turbulent Prandtl number). Must be >= 0. The stability correction
   of the eddy viscosity, :math:`\sqrt{|S|^2 - 3 N^2}`, keeps the factor 3
   whatever this value is.

.. input_param:: Kosovic.writeTerms

   **type:** Boolean, optional, default = false

   Adds the fields `Nij` (the tensor
   :math:`4 (S_{ik} S_{kj} + S_{ik} R_{kj} - R_{ik} S_{kj})`, unscaled) and
   `divNij` (the momentum source) to the plot files.


   

.. input_param:: KLAxell_coeffs.Cmu

   **type:** Real, optional, default = 0.556

   Neutral stability constant :math:`C_\mu` of the `KLAxell` model. It sets the
   eddy viscosity, and the `KransAxell` source uses the same value for the
   dissipation and for the wall value of the turbulent kinetic energy. Must be
   > 0; the run aborts otherwise. `ABL.Cmu`, if given, must equal it.

.. input_param:: KLAxell_coeffs.Cb_stable

   **type:** Real, optional, default = 0.25

   Coefficient of the buoyancy length scale :math:`C_b \sqrt{k}/N` in stable
   stratification. Must be > 0; the run aborts otherwise.

.. input_param:: KLAxell_coeffs.Cb_unstable

   **type:** Real, optional, default = 0.35

   Coefficient of the length scale in unstable stratification. Must be > 0;
   the run aborts otherwise.

.. input_param:: KLAxell_coeffs.legacy_alphaeff

   **type:** Boolean, optional, default = false

   The turbulent heat diffusivity is :math:`\mu_t / Pr_t(R_t)`. By default
   :math:`R_t` is the value used for the eddy viscosity in the same step,
   including the neutral switch (:math:`R_t = 0` when the magnitude of
   `ABL.surface_temp_flux` is below :math:`10^{-5}` and the height is at or
   below :input_param:`ABL.meso_sponge_start`). The heat flux then
   matches the buoyancy production, which uses :math:`C_\mu'(R_t)`. With
   `true`, :math:`R_t` is recomputed from the updated length scale without
   the neutral switch, as in earlier versions.

`KLAxell_coeffs.Cmu_prime` and `KLAxell_coeffs.prandtl` were read but never used: the
stability function :math:`C_\mu'` and the turbulent Prandtl number follow from `Cmu` and
the turbulent Richardson number. Setting either key now aborts the run. The `KLAxell`
length scale uses :input_param:`ABL.kappa`, like the `KransAxell` wall boundary condition.

The following inputs change how the ``KLAxell`` model treats the cells next to
the terrain with ``TerrainDrag`` (binary blanking). They are all off by
default, and with every one of them off the model is unchanged;
:input_param:`TerrainDrag.wall_treatment` ``= improved`` turns them all on.
They have no effect on flat ground; turned on without ``TerrainDrag`` in
:input_param:`incflo.physics` they abort. Use them with
``incflo.diffusion_type = 2`` (see :input_param:`TerrainDrag.wall_treatment`).

.. input_param:: KLAxell.terrain_wall_stencil

   **type:** Boolean, optional, default = false (true with
   :input_param:`TerrainDrag.wall_treatment` ``= improved``)

   Computes the strain rate of the fluid cells next to the blanked cells with
   the stencil flat ground uses at its bottom wall, the wall being the face of
   the blanked cell. A central difference into a blanked cell reads the
   near-zero velocity that the drag holds there, which for a horizontal wall
   gives :math:`\partial U/\partial z = U_{k+1}/(2\Delta z)`,
   :math:`\ln(1.5\Delta z/z_0)/4` times the log-law shear (1.4 at
   :math:`\Delta z` = 2 m, :math:`z_0` = 0.01 m), and a wall-cell TKE above the value flat ground reaches.
   With this option, in each direction that has a blanked neighbor the
   tangential velocity takes the one-sided derivative
   :math:`(-3u_0 + 4u_1 - u_2)/(2\Delta x)` over the wall cell and the next
   two fluid cells, and the wall-normal velocity the derivative through its
   zero wall value, :math:`(u_0 + u_1/3)/\Delta x`. Between two blanked
   neighbors the derivative is zero; at a domain face the usual boundary
   stencil is kept.

.. input_param:: KLAxell.terrain_blanked_face_length

   **type:** Boolean, optional, default = false (true with
   :input_param:`TerrainDrag.wall_treatment` ``= improved``)

   Measures the mixing length from the top face of the blanked column instead
   of the terrain height. A cell is blanked when its center is at or below
   the terrain height :math:`h`, so the resolved wall is the face at
   :math:`z_f = z_{lo} + \max(0, \Delta z\,\lfloor (h - z_{lo})/\Delta z + 1/2 \rfloor)`,
   and the height entering the mixing length is
   :math:`\max(z_c - z_f, \Delta z/2)`, as flat ground measures it from its
   wall. The count of blanked cells is checked against the same test that
   ``TerrainDrag`` uses, so that it also matches when :math:`h` is at a cell
   center. Not available in the ocean-wave mode of ``TerrainDrag`` (aborts
   if set there). Measured from :math:`h` it is up to :math:`\Delta z/2` smaller in the
   first fluid cells, depending on where :math:`h` falls within its cell.

.. input_param:: KLAxell.terrain_face_stress

   **type:** Boolean, optional, default = false (true with
   :input_param:`TerrainDrag.wall_treatment` ``= improved``)

   Sizes the turbulent viscosity of the drag cells (the first fluid cell above
   a blanked column) so that the face above each one carries the wall stress
   of the ``DragForcing`` wall law. ``DragForcing`` relaxes the drag cell
   toward the log-law velocity of the cell above it, which fixes the velocity
   difference across that face; with the model viscosity the resolved flux
   through it falls short of :math:`u_*^2`. With this option

   .. math::

      u_* = \frac{\kappa |U_{k+1}|}{\ln(1.5\Delta z/z_0) - \psi_m(1.5\Delta z/L)},
      \qquad
      \mu_f = \frac{\rho u_*^2 \Delta z}{\max(|U_{k+1} - U_k|, 10^{-2})},
      \qquad
      \mu_k = \max(2\mu_f - \mu_{k+1}, 0),

   where :math:`U` is the horizontal velocity, :math:`\mu_f` the face value
   (the mean of the two cells) and :math:`z_0` is floored at
   :math:`10^{-4}` m, as in ``DragForcing``. :math:`\psi_m` is used only with
   ``ABL.wall_het_model = mol``, and :math:`\kappa`, :math:`L` and the
   stability constants (``ABL.kappa``, ``ABL.monin_obukhov_length``,
   ``ABL.mo_beta_m``, ``ABL.mo_gamma_m``) are the ``ABL`` inputs the model
   already reads for its stability functions. The shear production of the
   drag cell keeps the model viscosity.

.. input_param:: KLAxell.terrain_face_heat_flux

   **type:** Boolean, optional, default = false (true with
   :input_param:`TerrainDrag.wall_treatment` ``= improved``)

   With ``ABL.wall_het_model = mol``, sizes the heat diffusivity of the drag
   cells so that the face above each one carries the surface heat flux given
   by the Obukhov length, as :input_param:`KLAxell.terrain_face_stress` does
   for the stress. ``DragTempForcing`` relaxes the drag cell toward the
   surface-layer temperature from the cell above it, which fixes the
   temperature difference across that face. With :math:`u_*` as above,

   .. math::

      \theta_* = \frac{\theta_k u_*^2}{\kappa g L}, \qquad
      q = -u_*\theta_*, \qquad
      \alpha_f = \frac{\rho q \Delta z}{\theta_k - \theta_{k+1}}, \qquad
      \alpha_k = \max(2\alpha_f - \alpha_{k+1}, \alpha_{lam}).

   The diffusivity is left unchanged where the temperature difference is below
   :math:`10^{-4}` K or runs against the flux. Near neutral (large
   :math:`|L|`) :math:`q` is small and the floor :math:`\alpha_{lam}` is
   reached. Without ``mol`` the option is turned off when it comes from
   :input_param:`TerrainDrag.wall_treatment`, and aborts when set on its own.
