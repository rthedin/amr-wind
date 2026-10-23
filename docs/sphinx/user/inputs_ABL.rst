.. _inputs_abl:

Section: ABL
~~~~~~~~~~~~

This section is for setting atmospheric boundary layer parameters.

.. input_param:: ABL.kappa

   **type:** Real, optional, default = 0.41

   Wall model coefficient.

.. input_param:: ABL.surface_roughness_z0

   **type:** Real, optional, default = 0.1

   Wall model surface roughness length (in meters). When specified, this sets both the aerodynamic 
   (pertaining to velocity) and thermal (pertaining to temperature) roughness lengths.
   These roughness lengths can be set separately using the two parameters listed below.
   Specifying these separately will override the combine parameter; however,
   if this parameter (which sets both) is used with only one of the roughness parameters below,
   the code will abort due to the ambiguity.

.. input_param:: ABL.aerodynamic_roughness_length

   **type:** Real, optional, default = 0.1

   Wall model surface roughness length (in meters) for aerodynamic processes, i.e., the wall momentum flux/shear stress model.

.. input_param:: ABL.thermal_roughness_length

   **type:** Real, optional, default = 0.1

   Wall model surface roughness length (in meters) for thermal processes, i.e., the wall temperature flux model.

.. input_param:: ABL.mo_alpha_h

   **type:** Real, optional, default = 1.0

   Scaling of the Monin-Obukhov temperature profile, which plays the role of a
   turbulent Prandtl number in the surface layer:
   :math:`\theta(z) = \theta_s - \alpha_h \, q \left[ \ln(z/z_{0t}) - \psi_h(z/L) \right] / (\kappa u_*)`.
   It is used both in the Monin-Obukhov iteration and in the wall-model heat flux,
   so the flux the wall model applies is the one the iteration computes.
   The default leaves the results unchanged.

.. input_param:: ABL.normal_direction

   **type:** Integer, optional, default = 2

   Wall model normal direction. x-direction = 0, y-direction = 1, z-direction = 2.

.. input_param:: ABL.wall_position

   **type:** Real, optional, default = lower domain boundary in :input_param:`ABL.normal_direction` (typically zlo)

   Physical location of the wall-modeled boundary along the wall-normal direction.
   When not specified, this is taken from the lower problem-domain boundary in the
   selected normal direction. This value is used together with
   :input_param:`ABL.log_law_height` to determine the sampling location for the
   planar averages used in the wall model.

.. input_param:: ABL.log_law_height

   **type:** Real, optional

   Height to evaluate the log law for the wall model.
   If this parameter is not specified in the input file, the first half cell height is set
   as the log law height. For meshes with multiple levels, the first half cell height
   on the coarsest level is used, but the :input_param:`ABL.stats_max_level` parameter
   overrides this behavior, making the code use the first half cell height on the specified level.
   Therefore the log law height depends on the domain size and number of elements in the normal direction.
   If this parameter is set to a number the log law is evaluated at a fixed height.

.. input_param:: ABL.stats_max_level

   **type:** Integer, optional, default = 0

   Maximum level to consider in planar averages for ABLStats. This also controls what is
   considered as the first cell height, which is the default log law height in the wall function.
   The default value of 0 means the coarsest level is used, which has been the legacy behavior of
   the code. Setting a value of -1 means there is no limit on the maximum level, allowing the code
   to select the finest level available on the mesh.

.. input_param:: ABL.surface_temp_flux

   **type:** Real, optional

   Surface temperature flux setting for the ABL wall function. Specifies a constant temperature flux
   at the wall-modeled boundary. This is not a required argument because there are other options for
   setting up the surface temperature condition. With :input_param:`ABL.surface_stability_type` =
   heat_flux it also sets the heat flux of the terrain surface, in K m/s.

.. input_param:: ABL.surface_temp_flux_timetable

   **type:** String, optional

   File name of a surface temperature flux time table, allowing the heat flux
   (in K m/s) to change with time. The file starts with a one-line header,
   followed by one line per entry with the time and the flux; the flux is
   interpolated linearly in time. This option cannot be combined with another
   surface temperature boundary condition.

.. input_param:: ABL.surface_temp_timetable

   **type:** String, optional

   File name of surface temperature time table, allowing the surface temperature
   to change with time without specifying a surface temperature rate.

.. input_param:: ABL.near_surface_temp_timetable

   **type:** String, optional

   File name of a time table of the temperature measured near the surface, at
   :input_param:`ABL.near_surface_height` (for example a 2 m station
   temperature). The file starts with a one-line header, followed by one line per
   entry with the time and the temperature. At every timestep, the
   Monin-Obukhov profile is required to pass through this temperature and the
   planar-averaged temperature at the log-law height, which gives the surface
   temperature and the surface heat flux together. This option cannot be combined
   with another surface temperature boundary condition.

.. input_param:: ABL.near_surface_height

   **type:** Real, mandatory with :input_param:`ABL.near_surface_temp_timetable`

   Height above the surface (in meters) of the temperatures in
   :input_param:`ABL.near_surface_temp_timetable`. It must be above the thermal
   roughness length and differ from the log-law height.

.. input_param:: ABL.surface_temp_rate

   **type:** Real, optional

   Constant rate at which the surface temperature changes, in K/hr. With
   :input_param:`ABL.surface_stability_type` = heating_rate it also sets the rate for the
   terrain surface.

.. input_param:: ABL.surface_temp_init

   **type:** Real, optional

   Initial temperature of the wall-modeled surface. This parameter is only active
   when the surface temperature rate is specified. If this parameter is active but not
   specified, the initial temperature will be set to the reference temperature of the simulation.
   With :input_param:`ABL.surface_stability_type` = heating_rate or surface_temperature it also
   sets the terrain surface temperature.

.. input_param:: ABL.surface_temp_rate_tstart

   **type:** Real, optional

   Start time of the surface temperature rate of change. Prior to this time,
   the surface temperature remains at the initial value. This parameter is only
   active when the surface temperature rate is specified. The default start time is 0.


.. input_param:: ABL.temperature_heights

   **type:** List of Reals, mandatory

   Height(s) in meters at which temperature values are prescribed.

.. input_param:: ABL.temperature_values

   **type:** List of Reals (has to be same length as :input_param:`ABL.temperature_heights`), mandatory

   Temperature values in Kelvin at the corresponding :input_param:`ABL.temperature_heights`.
   The temperature below the first height is assumed to be constant and equal to the
   first temperature value.
   The temperature between values is initialized to have linear variation.
   The final temperature is constant above the last specified height.


.. input_param:: ABL.perturb_velocity

   **type:** Boolean, optional, default = true

   If true this flag turns on perturbations to the freestream flow.

.. input_param:: ABL.pertub_ref_height

   **type:** Real, optional, default = 50.0

   Reference height for velocity perturbations,
   perturbations exist below this height and decay above this height.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.Uperiods

   **type:** Real, optional, default = 4.0

   Number of sinusoidal waves in x-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.Vperiods

   **type:** Real, optional, default = 4.0

   Number of sinusoidal waves in y-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.deltaU

   **type:** Real, optional, default = 1.0

   Amplitude of fluctuations in x-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.deltaV

   **type:** Real, optional, default = 1.0

   Amplitude of fluctuations in y-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.perturb_temperature

   **type:** Boolean, optional, default = false

   Perturb temperature field with random fluctuations.

.. input_param:: ABL.theta_amplitude

   **type:** Real, optional, default = 0.8 K

   Amplitude of the temperature perturbations added to the initial field. Only
   active when :input_param:`ABL.perturb_temperature` is true.

.. input_param:: ABL.cutoff_height

   **type:** Real, optional, default = domain height

   Height below which temperature perturbations are added

.. input_param:: ABL.random_gauss_mean

   **type:** Real, optional, default = 0.0

   Mean for the Gaussian random number generator

.. input_param:: ABL.random_gauss_var

   **type:** Real, optional, default = 1.0

   Variance for the Gaussian random number generator


.. input_param:: ABL.bndry_file

   **type:** String, optional, default = ""

   NetCDF-4 file name for ABL inflow

.. input_param:: ABL.bndry_io_mode

   **type:** Int, optional, default = -1

   IO mode (0=output, 1=input)

.. input_param:: ABL.bndry_planes

   **type:** String, optional, default = ""

   IO planes for ABL inflow

.. input_param:: ABL.bndry_output_start_time

   **type:** Real, optional, default = 0.0

   Time at which to start ABL inflow output

.. input_param:: ABL.bndry_var_names

   **type:** String, optional, default = ""

   Variables for IO for ABL inflow

.. input_param:: ABL.wall_shear_stress_type

   **type:** String, optional, default = "Moeng"

   Wall shear stress model: options include
   "constant", "local", "Schumann", and "Moeng"

.. input_param:: ABL.bndry_output_format

   **type:** String, optional, default = "native"

   Output of boundary plane files. Valid values are ``netcdf`` and ``native``.

.. input_param:: ABL.initial_condition_input_file

   **type:** String, optional, default= ""

   File that contains initial conditions for the
   velocity field in netcdf file format.
   This file is expected to have the same dimensions as the simulation.
   Values are passed directly from the file to the velocity field inside the code.
   Only spanwise velocity components are supported.

.. input_param:: ABL.anelastic

   **type:** Boolean, optional, default= false

   Activate anelastic behavior. This adds `reference_density` and
   `reference_pressure` fields.

.. input_param:: ABL.bottom_reference_pressure

   **type:** Real, optional, default = 1.01325e5

   Reference pressure at the bottom of the domain. Used for anelastic ABL.

.. input_param:: ABL.initial_wind_profile 

   **type:** Boolean, optional, default= false

   Activates the reading of wind speed profile from a file. Recommended for 
   RANS models and also for wind conditions input from climate model.

.. input_param:: ABL.rans_1dprofile_file 

   **type:** String, optional, default = ""

   This input is required when the ABL.initial_wind_profile is set to True. 

.. input_param:: ABL.minimum_vertical_position 

   **type:** Real, optional, default = 0.1

   When using an initial wind profile, this parameter limits how small the vertical
   position can be. It affects calculations in the initialization and in DragForcing.

.. input_param:: ABL.meso_sponge_start 

   **type:** Real, optional, default = 600

   Approximate height of the planetary boundary layer height to enable the forcing 
   in the free atmosphere. Recommended for use with RANS model and optionally to run 
   LES with non canonical flow conditions. If not specified in the input file, the 
   default value is active for the TKE source term but not for velocity or temperature.
   For this value to be used, the source terms must be included as input arguments 
   as follows, for temperature, velocity, or TKE, respectively.
   
   `Temperature.source_terms  = TemperatureFreeAtmosphereForcing`

   `ICNS.source_terms  = VelocityFreeAtmosphereForcing`

   `TKE.source_terms = KransAxell`

.. input_param:: ABL.velocity_sponge_start

   **type:** Real, optional

   Instead of specifying the same sponge start location for velocity, temperature, and TKE,
   this input argument allows the velocity sponge start (involved in `VelocityFreeAtmosphereForcing`)
   to be differentiated from the value specified by :input_param:`ABL.meso_sponge_start`.

.. input_param:: ABL.temp_sponge_start

   **type:** Real, optional

   Instead of specifying the same sponge start location for velocity, temperature, and TKE,
   this input argument allows the temperature sponge start (involved in `TemperatureFreeAtmosphereForcing`)
   to be differentiated from the value specified by :input_param:`ABL.meso_sponge_start`.

.. input_param:: ABL.wall_het_model

   **type:** String, optional, default = "none"

   Allows the use of different surface model options for the Monin-Obukhov length. Currently supports two options:
   (i) "none" - original model in the code and (ii) "mol" - Monin-Obukhov length is constant while heat-flux varies 

.. input_param:: ABL.monin_obukhov_length

   **type:** Real, optional, default = 1e-10 times the largest Real (neutral)

   Used in conjunction with `ABL.wall_het_model`. The default value runs a neutral boundary layer.
   With `ABL.wall_het_model = "mol"` it must satisfy \|L\| > machine epsilon; the run aborts otherwise.

.. input_param:: ABL.surface_stability_type

   **type:** String, optional, default = "monin_obukhov_length"

   Sets the stability of the terrain surface layer used by the `DragForcing` and
   `DragTempForcing` source terms and by the `Kosovic` model on terrain. Only these and
   `KLAxell`/`KransAxell` read this key; the flat-wall `ABLWallFunction` does not. The options
   follow the amrTerrain 1-D RANS solver and reuse the `ABLWallFunction` keys: the flat wall has
   the heat-flux and heating-rate inputs, and a constant surface temperature there is a zero
   :input_param:`ABL.surface_temp_rate`.
   Options other than the default require `ABL.wall_het_model = "mol"`; the run aborts otherwise.
   The run also aborts if a surface temperature key that belongs to another option
   (:input_param:`ABL.surface_temp_flux`, :input_param:`ABL.surface_temp_rate`,
   :input_param:`ABL.surface_temp_timetable`, :input_param:`ABL.surface_temp_flux_timetable`,
   :input_param:`ABL.near_surface_temp_timetable`) is also set, so that the flat wall and the
   terrain cannot use different inputs. :input_param:`ABL.kappa` must be > 0 and
   `ABL.mo_gamma_m`, `ABL.mo_beta_m`, `ABL.mo_gamma_h` and `ABL.mo_beta_h` must be >= 0.

   ======================== ============================================================ ===========================================================================
   Value                    Inputs read                                                  Effect
   ======================== ============================================================ ===========================================================================
   monin_obukhov_length     :input_param:`ABL.monin_obukhov_length`                      The same Obukhov length in every terrain cell (the behaviour before this option).
   heat_flux                :input_param:`ABL.surface_temp_flux` (required, \|q\| <= 10) Each terrain cell solves for the friction velocity, Obukhov length and surface temperature.
   heating_rate             :input_param:`ABL.surface_temp_rate` (required),             The surface temperature changes at the given rate (K/hr) after the start time;
                            :input_param:`ABL.surface_temp_init` (> 0, default =         each terrain cell solves for the heat flux and Obukhov length.
                            reference temperature),
                            :input_param:`ABL.surface_temp_rate_tstart` (>= 0)
   surface_temperature      :input_param:`ABL.surface_temp_init` (required, > 0)         A constant surface temperature; each terrain cell solves for the heat flux and
                                                                                         Obukhov length.
   ======================== ============================================================ ===========================================================================

   The per-cell solve uses the wind speed and potential temperature of the cell above the
   terrain cell, at 1.5 times the cell height, and the terrain roughness for both momentum
   and heat. It finds :math:`\zeta = z/L` by bisection in :math:`[\zeta_{lo}, 2]`, where
   :math:`\zeta_{lo}` is -5 or, if larger, the value at which
   :math:`\ln(z/z_0) - \psi_h` falls to half of :math:`\ln(z/z_0)`. When there is no root
   in this range (a stable surface below the critical wind speed, or a very unstable surface
   at low wind speed) :math:`\zeta` is held at the bound: with `heat_flux` the heat flux is
   then reduced to the flux the bound can carry, and with the surface temperature inputs the
   heat flux follows from the surface temperature. In the heat-flux option the requested
   flux is delivered in full for wind speeds of 3 m/s and above at 1.5 times the cell
   height; at 1 m/s and below a large unstable flux is reduced strongly. With `Kosovic` the near-wall RANS length scale uses the
   mean of 1/L over the terrain cells of each level (so it can differ between levels), and
   options other than the default need `TerrainDrag`. The roughness length used in the solve is
   limited to 1/e of the reference height. The `KLAxell` model and the `KransAxell` source term take the friction
   velocity from :input_param:`ABL.monin_obukhov_length` and the TKE boundary condition from
   :input_param:`ABL.surface_temp_flux`, so they abort for any option other than the default.

.. input_param:: ABL.terrain_aligned_profile 

   **type:** Boolean, optional, default= false

   Used in conjunction with immersed forcing for terrain. This option allows the user to align the wind, temperature and turbulence profiles to be aligned with the terrain.
