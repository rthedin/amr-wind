.. _inputs_forestdrag:

Section: ForestDrag
~~~~~~~~~~~~~~~~~~~

These parameters are active when ``ForestDrag`` is included in
:input_param:`incflo.physics`. Two forest representations are supported:

1. A legacy cylinder-based forest file, controlled by ``ForestDrag.forest_file``.
2. A point-cloud forest description, controlled by ``ForestDrag.point_cloud_files``.

The two approaches are mutually exclusive and cannot be used together in the same run.

When ``TerrainDrag`` is also included in :input_param:`incflo.physics`, the
forests stand on the terrain: tree heights in the forest file and the ``z``
coordinates of the point-cloud files are heights above the local ground, and no
forest drag is applied inside the terrain. The ground height is taken
column by column from the ``terrain_height`` field, so a forest patch on a
slope follows the slope. ``TerrainDrag`` must be listed before ``ForestDrag``
in :input_param:`incflo.physics`.

.. input_param:: ForestDrag.forest_file

   **type:** String, optional, default = ``forest.amrwind``

   Input file for the legacy cylinder-based forest representation. If
   ``ForestDrag`` is enabled and
   :input_param:`ForestDrag.point_cloud_files` is not provided,
   Kynema-SGF reads this file.

   Each non-comment line of the file must contain eight values,

   ``TreeType xc yc height diameter cd lai laimax``

   where ``TreeType`` selects the vertical LAD model, ``xc`` and ``yc`` are
   the forest-center coordinates, ``height`` and ``diameter`` define the
   cylindrical footprint, ``cd`` is the drag coefficient, ``lai`` is the leaf
   area index, and ``laimax`` gives the normalized height of the LAD maximum
   for the type-2 profile.

   The supported legacy LAD models are the same as those described in the
   theory documentation: a uniform profile and an analytical vertically varying
   profile.

.. input_param:: ForestDrag.point_cloud_files

   **type:** List of strings, optional

   List of point-cloud files, with one file per forest patch. When this
   parameter is present, the legacy cylinder-based input
   :input_param:`ForestDrag.forest_file` must not be specified.

   Each non-comment line of a point-cloud file must contain four values,

   ``x y z lad``

   where ``x``, ``y``, and ``z`` are the point coordinates and ``lad`` is the
   leaf area density sampled at that point.

   Blank lines are ignored. Text following ``#`` on a line is treated as a
   comment.

.. input_param:: ForestDrag.coefficients_of_drag

   **type:** List of reals, mandatory when
   :input_param:`ForestDrag.point_cloud_files` is used

   Drag coefficient for each point-cloud forest. The number of entries in this
   list must match the number of files listed in
   :input_param:`ForestDrag.point_cloud_files`.

.. input_param:: ForestDrag.point_neighbors

   **type:** Integer, optional, default = 4

   Number of nearest point-cloud samples used to reconstruct the local LAD at a
   cell center. Values are internally clamped to the range 1 through 8.

.. input_param:: ForestDrag.point_interp_eps

   **type:** Real, optional, default = 1.0e-12

   Regularization parameter used in the inverse-distance weighting for the
   point-cloud forest model. This value avoids singular weights when a cell
   center is extremely close to a point-cloud sample.

.. input_param:: ForestDrag.terrain_aware

   **type:** Boolean, optional, default = true

   Place the forests on the ``TerrainDrag`` terrain. Only used when
   ``TerrainDrag`` is active. Set it to false when the point-cloud ``z``
   coordinates are absolute heights rather than heights above the ground.

.. input_param:: ForestDrag.model

   **type:** String, optional, default = ``canopy``

   Forest representation. ``canopy`` resolves the forest as a drag
   :math:`C_d \, L` in the momentum equation (``forest_drag`` field, applied by
   the ``ForestForcing`` source term). ``roughness`` applies no canopy drag:
   each forest footprint writes its roughness length into the ``terrainz0``
   field of ``TerrainDrag``, which the terrain wall model (``DragForcing``,
   ``KransAxell``, ``DragTempForcing`` and the Kosovic model) uses in the
   first cell above the terrain. The ``forest_drag`` field stays zero, so
   ``ForestForcing`` and the canopy turbulence terms have no effect, and
   ``forest_id`` marks the cells between the ground and the tree top.

   The footprint is the cylinder cross section of a legacy forest, or the
   :math:`x-y` convex hull of a point-cloud forest. No displacement height is
   applied.

   The ``roughness`` model requires ``TerrainDrag`` listed before
   ``ForestDrag`` in :input_param:`incflo.physics`, otherwise the run aborts.
   The roughness is then set in this order:

   1. ``TerrainDrag`` fills ``terrainz0`` from
      ``TerrainDrag.roughness_file``, or with 0.1 m where no file is given.
   2. ``ForestDrag`` overwrites it inside every forest footprint, whole
      columns, in the order of the forest file or of
      :input_param:`ForestDrag.point_cloud_files`; where footprints overlap,
      the later forest wins.

   ``terrainz0`` only acts in cells above blanked terrain cells. On a flat
   domain bottom without terrain cells, the ABL wall function uses the uniform
   ``ABL.surface_roughness_z0`` and the forest roughness has no effect; raise
   the ground by a whole number of cells in the terrain file to use it.

.. input_param:: ForestDrag.roughness_z0

   **type:** List of reals, optional

   Roughness length (m) of the forests with
   :input_param:`ForestDrag.model` = ``roughness``. One value applies to all
   forests; otherwise give one value per forest, in the order of the forest
   file or of :input_param:`ForestDrag.point_cloud_files`. Exactly one of
   this parameter and :input_param:`ForestDrag.roughness_height_fraction`
   must be given in roughness mode.

.. input_param:: ForestDrag.roughness_height_fraction

   **type:** List of reals, optional

   Roughness length as a fraction :math:`c` of the tree height,
   :math:`z_0 = c \, h`, with :input_param:`ForestDrag.model` =
   ``roughness``. One value for all forests or one per forest. The height of
   a point-cloud forest is the highest sample. A common choice is
   :math:`c = 0.1`.

