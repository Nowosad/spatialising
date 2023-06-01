# devtools::load_all()
data(r_start, package = "spatialising")
r_startr = terra::rast(r_start)

set.seed(2022-04-04)
ts10v1 = kinetic_ising(r_startr, B = -0.3, J = 0.7, updates = 1, version = 1)
# set.seed(2022-04-04)
# ts10v1b = kinetic_ising(r_startr, B = -0.3, J = 0.7, updates = 1, version = 1)

# plot(ts10v1)
# plot(ts10v1b)

set.seed(2022-04-04)
ts10v2 = kinetic_ising(r_startr, B = -0.3, J = 0.7, updates = 1, version = 2)
# plot(ts10v2)
test_that("both versions give the same results", {
  expect_equal(terra::values(ts10v1), terra::values(ts10v2))
})

