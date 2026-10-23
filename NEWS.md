spatialising 0.7.0
==================

New features
------------

* Added spatially varying external pressure through the new `external_field` argument to `kinetic_ising()`.
* Added support for spatial masks: `NA` cells are excluded from the simulation and no longer contribute neighbor pressure.
