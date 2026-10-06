.. _inputs_momentum_sources:
   
Section: Momentum Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. input_param:: ICNS.source_terms

   **type:** String(s), optional
   
   Activates source terms for the incompressible Navier-Stokes momentum
   equations. These strings can be entered in any order with a space between
   each. Please consult the :doc:`../doxygen/html/index` for a
   comprehensive list of all momentum source terms available. Note that the
   following input arguments specific to each source term will only be active
   if the corresponding source term (the root name) is listed in 
   :input_param:`ICNS.source_terms`.

.. input_param:: CoriolisForcing.latitude 

   **type:** Real, mandatory
   
   Latitude in degrees where the Coriolis forcing is computed. Positive values
   indicate northern hemisphere.
   
.. input_param:: CoriolisForcing.rotational_time_period 

   **type:** Real, optional, default = 86400.0
   
   Rotational time period of a day in seconds.
   
.. input_param:: CoriolisForcing.east_vector

   **type:** List of 3 reals, optional, default = 1.0 0.0 0.0
   
   East vector that gives the orientation of the grid w.r.t. to planetary coordinate system.
   This vector is automatically normalized within Kynema-SGF.
   
.. input_param:: CoriolisForcing.north_vector

   **type:** List of 3 reals, optional, default = 0.0 1.0 0.0
   
   North vector that gives the orientation of the grid w.r.t. to planetary coordinate system.
   This vector is automatically normalized within Kynema-SGF.

.. input_param:: GeostrophicForcing.geostrophic_wind

   **type:** List of 3 reals, optional

   The user has to choose between GeostrophicForcing and ABLForcing. 
   CoriolisForcing input must be present when using GeostrophicForcing.
   These checks are not enforced for now.

.. input_param:: GeostrophicForcing.geostrophic_wind_timetable

   **type:** String, optional
   
   Input file name for table that lists time in seconds, wind speed 
   in meters per second, and horizontal wind direction in degrees of the Geostrophic 
   forcing velocity. Each line in the file should be a sequence of 
   three floats specifying the inputs in that order (e.g., 0.0 8.0 -5.0). If this
   argument is present, the :input_param:`GeostrophicForcing.geostrophic_wind`
   will be ignored. Note that the code expects there to be a single-line header
   at the beginning of the geostrophic wind timetable file; if no header exists, 
   the first line of data will be ignored.
   
.. input_param:: ABLForcing.abl_forcing_height

   **type:** Real, mandatory
   
   Height in meters at which the flow is forced to maintain the freestream
   inflow velocities specified through :input_param:`incflo.velocity`.

.. input_param:: ABLForcing.velocity_timetable

   **type:** String, optional
   
   Input file name for table that lists time in seconds, wind speed 
   in meters per second, and horizontal wind direction in degrees of the ABL 
   forcing velocity. Each line in the file should be a sequence of 
   three floats specifying the inputs in that order (e.g., 0.0 8.0 -5.0). If this
   argument is present, the :input_param:`incflo.velocity` argument
   will be ignored. Note that the code expects there to be a single-line header
   at the beginning of the velocity timetable file; if no header exists, the first 
   line of data will be ignored.

.. input_param:: ABLForcing.forcing_timetable_output_file

   **type:** String, optional
   
   Output file name for writing the ABL forcing vector to a text file over the course
   of a simulation. This output is primarily intended for replicating the ABL forcing
   from a precursor simulation in a subsequent inflow-outflow simulation by providing 
   an input file for Body Forcing. The output file will contain the time and three vector
   components of the force.

   With :input_param:`ABLForcing.free_atmosphere_damping`, the file also
   contains the geostrophic wind that balances the forcing and the
   free-atmosphere height.

.. input_param:: ABLForcing.free_atmosphere_damping

   **type:** Boolean, optional, default = false

   Relax the planar-averaged velocity above a free-atmosphere height :math:`h`
   toward the geostrophic wind that balances the ABL forcing :math:`\mathbf{F}`,
   :math:`U_g = F_y / f` and :math:`V_g = -F_x / f`, with the source term
   :math:`(U_g - \langle u \rangle(z)) / \tau` and :math:`(V_g - \langle v \rangle(z)) / \tau`
   for :math:`z \geq h`. This removes the inertial oscillations of the free
   atmosphere that the time-varying forcing otherwise excites. The Coriolis
   parameter :math:`f = 2 \Omega \sin \phi` is computed from the
   ``CoriolisForcing`` inputs, and ``CoriolisForcing`` must be one of the
   ``ICNS.source_terms``. As in ``GeostrophicForcing``, :math:`x` is east and
   :math:`y` is north. The forcing of a timestep corrects the drift of the
   previous one, so the damping starts at the second timestep after start-up.

.. input_param:: ABLForcing.free_atmosphere_height

   **type:** Real, mandatory with :input_param:`ABLForcing.free_atmosphere_damping`
   unless :input_param:`ABLForcing.detect_free_atmosphere_height` is true

   Height (in the domain coordinates) above which the free atmosphere is damped.

.. input_param:: ABLForcing.detect_free_atmosphere_height

   **type:** Boolean, optional, default = false

   Set the free-atmosphere height to the capping inversion height every
   timestep. The capping inversion height is the planar average of the height of
   the largest vertical potential temperature gradient on level 0, as in the ABL
   statistics output. It is meant for boundary layers under a capping inversion.

.. input_param:: ABLForcing.free_atmosphere_damping_time_scale

   **type:** Real, optional, default = 100.0

   Relaxation time scale :math:`\tau` in seconds.

.. input_param:: ABLForcing.free_atmosphere_damping_start_time

   **type:** Real, optional, default = 0.0

   Time at which the free-atmosphere damping starts.

.. input_param:: ABLForcing.free_atmosphere_damping_end_time

   **type:** Real, optional, default = no end

   Time at which the free-atmosphere damping ends.

.. input_param:: ABLForcing.forcing_timetable_frequency

   **type:** Int, optional
   
   The interval of timesteps for writing to the forcing timetable output file. The default
   is 1, i.e., writing every step, which is also the default of the boundary plane output feature.

.. input_param:: ABLForcing.forcing_timetable_start_time

   **type:** Real, optional
   
   The start time for writing to the forcing timetable output file. The default is 0.

.. input_param:: ABLForcing.abl_forcing_off_height

   **type:** Real, required for multiphase simulations with ABL
   
   This parameter indicates the vertical distance above the water level that the ABL
   forcing term should be turned off. This tuning parameter is used to avoid applying 
   the ABL forcing to ocean waves. This is not used when the volume fraction field (vof)
   is not present in the simulation.

.. input_param:: ABLForcing.abl_forcing_ramp_height

   **type:** Real, required for multiphase simulations with ABL
   
   This parameter indicates the vertical distance above the water level and the "off height"
   that the ABL forcing term should ramp up from zero to full strength. This is not used
   when the volume fraction field (vof) is not present in the simulation.

.. input_param:: ABLForcing.abl_forcing_band

   **type:** Real, optional for multiphase simulations with ABL
   
   This parameter is an additional safeguard against applying ABL forcing within the waves.
   This specifies the number of computational cells in a band around the air-water interface
   that the ABL forcing should be deactivated. While the other arguments relate to the height coordinate
   within the domain, this argument is relative to the actual position of water in the simulation.
   The default value is 2.

.. input_param:: BodyForce.type

   **type:** String, optional
   
   The type of body force being used. The default is uniform_constant, which applies a single constant
   force vector over the entire domain. Other available types are height_varying, oscillatory, and
   uniform_constant.

.. input_param:: BodyForce.magnitude

   **type:** List of 3 reals, conditionally mandatory
   
   The force vector to be applied as a body force. This argument is mandatory for uniform_constant 
   (default) and oscillatory body force types.

.. input_param:: BodyForce.angular_frequency

   **type:** Real, conditionally mandatory
   
   The angular frequency to be used for applying sinusoidal time variation to the body force. This 
   argument is mandatory for the oscillatory body force type and is only active for the oscillatory type.

.. input_param:: BodyForce.bodyforce_file

   **type:** String, conditionally mandatory
   
   The text file for specifying the body force vector as a function of height. This text file must contain
   heights (z coordinate values), force components in x, and force components in y. This argument is mandatory for
   the height_varying body force type and is only active for the height_varying type.

.. input_param:: BodyForce.uniform_timetable_file

   **type:** String, conditionally mandatory
   
   The text file for specifying the body force vector as a function of time. This text file must contain
   times, force components in x, force components in y, and force components in z. This argument is mandatory for
   the uniform_timetable body force type and is only active for the uniform_timetable type.  Note that the code 
   expects there to be a single-line header at the beginning of the uniform timetable file; if no header exists, 
   the first line of data will be ignored.

.. input_param:: DragForcing.drag_coefficient

   **type:** Real, optional

   This value specifies the coefficient for the forcing term in the immersed boundary forcing method. It is currently
   recommended to use the default value to avoid initial numerical stability. 

.. input_param:: DragForcing.sponge_strength

   **type:** Real, optional

   The value of the sponge layer coefficient. It is recommended to use the default value of 1.0.  

.. input_param:: DragForcing.sponge_density

   **type:** Real, optional

   The value of the sponge layer density. It is recommended to use the default value of 1.0.  

.. input_param:: DragForcing.bc_forcing_time_factor

   **type:** Real, optional, default = 5.0

   This value modifies the time scale of the BC forcing component of DragForcing relative to
   the time step size.

.. input_param:: DragForcing.sponge_west

   **type:** Boolean, optional, default = false

   This term turns on the sponge layer in the west (-x) boundary.

.. input_param:: DragForcing.sponge_east

   **type:** Boolean, optional, default = false

   This term turns on the sponge layer in the east (+x) boundary.

.. input_param:: DragForcing.sponge_south

   **type:** Boolean, optional, default = false

   This term turns on the sponge layer in the south (-y) boundary.

.. input_param:: DragForcing.sponge_north

   **type:** Boolean, optional, default = false

   This term turns on the sponge layer in the north (+y) boundary.

.. input_param:: DragForcing.sponge_distance_west

   **type:** Real, mandatory if west sponge is active

   This value is specified as a negative value when the inflow x-velocity is <=0. 

.. input_param:: DragForcing.sponge_distance_east

   **type:** Real, mandatory if east sponge is active

   This value is specified as a positive value when the inflow x-velocity is >=0. 

.. input_param:: DragForcing.sponge_distance_south

   **type:** Real, mandatory if south sponge is active

   This value is specified as a negative value when the inflow y-velocity is <=0. 

.. input_param:: DragForcing.sponge_distance_north

   **type:** Real, mandatory if north sponge is active

   This value is specified as a positive value when the inflow y-velocity is >=0.

.. input_param:: DragForcing.is_laminar

   **type:** int, optional

   This term turns off the sponge layer. This term is required for terrain simulations with periodic 
   boundary conditions. The default value is 0. 

.. input_param:: DragForcing.wave_model_inviscid_form_drag

   **type:** Boolean, optional, default = false

   This input file option turns on or off an inviscid model for the form drag of waves in the domain. 
   The formulation of this model is adapted from the Moving Surface Drag (MOSD) model developed by
   `Ayala et al (2024) <https://doi.org/10.1007/s10546-024-00884-8>`_.
   
   When the OceanWaves physics module is active, and the volume fraction variable ("vof") is not in the simulation,
   DragForcing will represent ocean waves as moving terrain. This is automatic and independent of DragForcing
   input arguments. When waves are represented as moving terrain and
   there is sufficient mesh resolution to resolve the shape of the wave, the blanking of cells performed
   by the DragForcing routine will naturally introduce the form drag of the waves into the flow. However,
   when the waves are not sufficiently resolved, such as when the wave amplitude is less than the cell height,
   the analytical model for the form drag, activated by setting this option to true, can be used to compensate
   for the lack of resolution. Therefore, this option should remain set to false except in scenarios
   when the form drag is known to be under-resolved.

The following arguments are influential when ``MetMastForcing`` is included in
:input_param:`ICNS.source_terms`. The source term relaxes the velocity towards
met-mast and lidar measurements. Each measurement station is a horizontal
location and a profile of one or more heights. Station :math:`i` has the weight

.. math::

   w_i = \exp\left(-\frac{1}{4}\left[\frac{r_i^2}{R_h^2} + \frac{d_i^2}{R_z^2}\right]\right)

where :math:`r_i` is the horizontal distance to the station and :math:`d_i` is
the vertical distance outside the measured height range, which is zero inside
it. The target :math:`t_i` is the measured velocity :math:`\bar{u}_i`,
interpolated linearly in height, widened to the band
:math:`\bar{u}_i \pm \alpha \sigma_i` so that only the part of the velocity
outside the measured variability is forced. The stations are combined into one
relaxation

.. math::

   S = -\min\left(\frac{W}{\tau}, \frac{1}{\Delta t}\right) (u - \bar{t}), \quad
   W = \min\left(1, \sum_i w_i\right), \quad
   \bar{t} = \frac{\sum_i w_i t_i}{\sum_i w_i}

When ``TerrainDrag`` is active, the heights are above the local terrain and the
cells inside the terrain are not forced.

.. input_param:: ABL.metmast_1dprofile_file

   **type:** String, optional

   File with met-mast points, one per line: ``x y z u v w T``, where ``z`` is
   the height above the terrain. The temperature ``T`` is not used. At least one
   of this file and :input_param:`ABL.metmast_profile_files` is required.

.. input_param:: ABL.metmast_profile_files

   **type:** List of strings, optional

   Lidar profile files, one per lidar. The first line is the location ``x y``,
   followed by one line per height ``z u v w su sv sw`` in increasing order of
   ``z``, the height above the terrain. ``su``, ``sv`` and ``sw`` are the
   standard deviations of the velocity components. Above and below the measured
   heights, the end values are used and the weight is tapered with
   :input_param:`ABL.metmast_vertical_radius`.

.. input_param:: ABL.metmast_timescale

   **type:** Real, optional, default = ``ABL.meso_timescale`` or 30.0

   Relaxation time scale :math:`\tau` in seconds. The relaxation rate is
   limited to :math:`1/\Delta t`.

.. input_param:: ABL.metmast_horizontal_radius

   **type:** Real, optional, default = 500.0

   Horizontal radius :math:`R_h` of the forcing around each station.

.. input_param:: ABL.metmast_vertical_radius

   **type:** Real, optional, default = 25.0

   Vertical radius :math:`R_z` of the forcing around a met-mast point, and of
   the taper above and below a lidar profile.

.. input_param:: ABL.metmast_damping_radius

   **type:** Real, optional, default = 1400.0

   Cutoff on the normalized squared distance
   :math:`r_i^2/R_h^2 + d_i^2/R_z^2` beyond which a station does not force.

.. input_param:: ABL.metmast_sigma_factor

   **type:** Real, optional, default = 1.0

   Half-width :math:`\alpha` of the tolerance band in standard deviations. A
   value of 0 relaxes the velocity to the measured mean. Met-mast points have
   no standard deviation and are always relaxed to the mean. Only used by the
   ``relaxation`` forcing.

.. input_param:: ABL.metmast_forcing_type

   **type:** String, optional, default = relaxation (or body_force, monitor)

   ``relaxation`` pulls each cell's own velocity towards the measurements as
   described above. In LES this
   also damps the resolved turbulence at the rate :math:`2k/\tau`, which is
   larger than the shear production for the time scales that hold the mean.
   ``body_force`` compares the measurements with the model's footprint average
   instead, like a virtual lidar, and applies the difference as a body force
   that does not depend on the local velocity. For every level :math:`l` of
   every station :math:`i`, the velocity is averaged over the station weights
   (with the horizontal radius :input_param:`ABL.metmast_averaging_radius`) and
   filtered in time,

   .. math::

      \bar{U}_{il} \leftarrow \bar{U}_{il} + \frac{\Delta t}{T_{avg}}
      \left(\frac{\sum w_{il} u \Delta V}{\sum w_{il} \Delta V} -
      \bar{U}_{il}\right)

   and a proportional-integral controller gives the force

   .. math::

      F_{il} = \frac{e_{il}}{\tau} + \frac{1}{\tau_I}\int e_{il}\,dt, \qquad
      e_{il} = \bar{u}_{il} - \bar{U}_{il}

   which is spread to the cells with the station weights,
   :math:`S = \sum_i w_i F_i / \max(1, \sum_i w_i)`. Cells covered by a finer
   AMR level and cells inside the terrain are left out of the average. The
   averaging and the controller are updated once per time step.

   A body force cannot add mass flux: the projection keeps the flow
   divergence-free, so a faster stream through the footprint is balanced by
   slower flow around it, and the accelerated stream continues downstream of
   the footprint. Large differences between the measurements and the
   unforced flow therefore give unrealistic flow around the station; such
   differences should be corrected through the inflow or the large-scale
   forcing instead.

   ``monitor`` computes and writes the footprint averages like ``body_force``
   but applies no force. It gives the model's virtual-lidar values, for
   example to compare a run with lidar data or to build a target with the same
   averaging as the controller.

   The controller is a loop with lags: the time filter :math:`T_{avg}` and the
   time :math:`\tau_d` for the forced air to reach the lidar. Air stays in the
   footprint for about :math:`T_p \approx R_h / U`, which is the gain of the
   force on the footprint average. For a stable loop, keep the integral slow
   compared with these lags, :math:`\tau_I \gtrsim 2 T_p (T_{avg} + \tau_d)`.
   In LES, an averaging time shorter than the large-eddy time scale (several
   minutes) makes the controller follow individual large eddies and add slow
   variability; average over about the lidar averaging period instead.

.. input_param:: ABL.metmast_averaging_time

   **type:** Real, optional, default = 120.0

   Time filter :math:`T_{avg}` of the footprint averages in seconds, for
   example the lidar averaging period or the time to cross the footprint.

.. input_param:: ABL.metmast_averaging_radius

   **type:** Real, optional, default = :input_param:`ABL.metmast_horizontal_radius`

   Horizontal radius of the footprint average. A radius matching the lidar
   measurement volume makes the controller hold the velocity at the lidar
   while the force is still spread over
   :input_param:`ABL.metmast_horizontal_radius`.

.. input_param:: ABL.metmast_gate_length

   **type:** Real, optional, default = 0.0

   Range-gate length of the footprint average in meters. When positive, each
   gate averages only the cells within half this length of its height, like a
   lidar range gate. When zero, the average uses the forcing weights, which
   interpolate between gates and taper above and below the end gates; that
   biases the end gates towards the flow outside the measured range. Use a
   length of at least the cell height, so that every gate contains cells.

.. input_param:: ABL.metmast_start_time

   **type:** Real, optional, default = 0.0

   Simulation time in seconds from which the footprint averages and the force
   start. Before it the body force is zero and nothing is averaged; at the
   first update after it the time filter starts from the current footprint
   average and the integral from zero. Use it to start the controller after
   the flow has spun up, for example after a restart from a precursor, so that
   the controller does not integrate the spin-up transient.

.. input_param:: ABL.metmast_gain_schedule

   **type:** Boolean, optional, default = false

   Set the gains of each gate from the local flow instead of
   :input_param:`ABL.metmast_timescale` and
   :input_param:`ABL.metmast_integral_timescale`. With the filtered
   footprint speed :math:`U` (at least :input_param:`ABL.metmast_min_speed`),
   the air spends :math:`T_p = R_h / U` under the force, which is also taken
   as the delay, and

   .. math::

      \tau_I = k_I \, T_p \, (T_{avg} + T_p), \qquad \tau = \tau_I / r

   with :math:`k_I` = :input_param:`ABL.metmast_integral_factor` and
   :math:`r` = :input_param:`ABL.metmast_integral_ratio`. Slow air, which a
   force accelerates the most, gets the gentlest gains. The integral is kept as
   a force, so changing gains do not rescale its past. The gains in use are
   written as the last two columns of the output file.

.. input_param:: ABL.metmast_integral_factor

   **type:** Real, optional, default = 2.0

   :math:`k_I` of the gain schedule. About 2 gives a phase margin near
   60 degrees; larger is slower and more robust.

.. input_param:: ABL.metmast_integral_ratio

   **type:** Real, optional, default = 10.0

   Ratio :math:`\tau_I / \tau` of the gain schedule.

.. input_param:: ABL.metmast_min_speed

   **type:** Real, optional, default = 0.5

   Smallest footprint speed in m/s used by the gain schedule, which bounds the
   residence time in nearly stagnant air.

.. input_param:: ABL.metmast_integral_timescale

   **type:** Real, optional, default = 4 :input_param:`ABL.metmast_timescale`

   Integral time scale :math:`\tau_I` in seconds. Zero or a negative value
   turns the integral term off.

.. input_param:: ABL.metmast_max_force

   **type:** Real, optional, default = 0.0

   Largest force per component in :math:`m/s^2`. The integral stops
   accumulating while the force is capped. Zero or a negative value means no
   cap.

.. input_param:: ABL.metmast_force_vertical

   **type:** Boolean, optional, default = false

   Also force the vertical velocity with ``body_force``.

.. input_param:: ABL.metmast_output_frequency

   **type:** Integer, optional, default = 10

   Steps between outputs of the ``body_force`` forcing to
   ``post_processing/metmast_body_force.txt``: per station level the measured
   velocity and standard deviation, the footprint average and standard
   deviation, and the force. The footprint standard deviation includes the
   spatial variation within the footprint. Zero turns the output off.

.. input_param:: ABL.metmast_restart_state

   **type:** String, optional

   Controller state to continue from on a restart. The state is written to
   ``post_processing/metmast_state<step>.txt`` at every checkpoint step; give
   the file of the checkpoint being restarted.


The following arguments are influential when ``GravityForcing`` is included in :input_param:`ICNS.source_terms`.

   .. input_param:: ICNS.use_perturb_pressure

   **type:** Boolean, optional, default = false
   
   When this option is off, the GravityForcing term is simply :math:`g`, which becomes
   :math:`\rho g` when included in the momentum equation. By activating this option,
   the momentum term applied by GravityForcing will become :math:`(\rho - \rho_0) g`,
   where :math:`rho_0` is some constant reference density profile. The reference density field
   can be created by either MultiPhase physics or anelastic ABL physics. By using the
   reference density, the pressure field seen by the solver is represented as a
   perturbation from a reference pressure field, enabling pressure_outflow boundary
   conditions to better handle certain flows, e.g., those with equilibrium pressure gradients
   parallel to the outflow plane.

   .. input_param:: ICNS.reconstruct_true_pressure

   **type:** Boolean, optional, default = false
   
   This option is only valid when the perturbational pressure form is being used, i.e.,
   :input_param:`ICNS.use_perturb_pressure` = true. Reconstructing the true pressure
   adds back the reference pressure profile to obtain the full pressure after the
   pressure solve has been performed. This makes no difference to the flow evolution,
   but it changes the field available for post-processing or coupling to overset solvers.
