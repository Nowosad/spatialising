# devtools::load_all()
data(r_start, package = "spatialising")
r_startr = terra::rast(r_start)

set.seed(2022-04-04)
ts10v1 = spatial_ising(r_startr, B = -0.3, J = 0.7, updates = 10, version = 1)

set.seed(2022-04-04)
ts10v2 = spatial_ising(r_startr, B = -0.3, J = 0.7, updates = 10, version = 2)

test_that("both versions give the same results", {
  expect_equal(terra::values(ts10v1), terra::values(ts10v2))
})

#
# devtools::load_all()
# bench::mark({set.seed(1);spatial_ising(r_startr, B = -0.3, J = 0.7, updates = 10, version = 1)->a},
#             {set.seed(1);spatial_ising_version1b(r_startr, B = -0.3, J = 0.7, updates = 10)->b},
#             {set.seed(1);spatial_ising(r_startr, B = -0.3, J = 0.7, updates = 10, version = 2)->ab},
#             check = FALSE)
