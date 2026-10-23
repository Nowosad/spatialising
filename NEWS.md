spatialising 0.7.0
==================

New features
------------

* Added spatially varying external pressure through the new `external_field` argument to `kinetic_ising()`.
* Added support for spatial masks: `NA` cells are excluded from the simulation and no longer contribute neighbor pressure.

Bug fixes
---------

* Cell locations are now sampled uniformly during simulation. Previously, cells
  on the outer rows and columns were selected less often than interior cells.
* `composition_index()` now ignores masked cells consistently for matrices and rasters.
