# devtools::load_all()
set.seed(2022-04-04)
data(r_start, package = "spatialising")
data(r_end, package = "spatialising")
ts1 = kinetic_ising(r_start, B = -0.3, J = 0.7, rule = "metropolis")

# r_start2 = r_start
# data(r_start, package = "spatialising")
# terra::plot(c(terra::rast(r_start2), terra::rast(r_start)))

set.seed(2022-04-04)
# data(r_start, package = "spatialising")
ts1v2 = kinetic_ising(terra::rast(r_start), B = -0.3, J = 0.7, rule = "metropolis", version = 2)
ts1v2 = terra::as.matrix(ts1v2, wide = TRUE)

# r_start2 = r_start
# data(r_start, package = "spatialising")
# terra::plot(c(terra::rast(r_start2), terra::rast(r_start)))

ts10 = kinetic_ising(r_start, B = -0.3, J = 0.7, updates = 10)

# terra::plot(terra::rast(ts1))

test_that("calculations on matrices works", {
  expect_equal(dim(ts1), c(50, 50))
  expect_equal(dim(ts10), c(50, 50, 10))
  expect_equal(ts1, r_end)
  expect_equal(ts1v2, r_end)
  expect_equal(ts1v2, ts1)
})

# terra::plot(c(terra::rast(r_end), terra::rast(ts1v2)))
# terra::plot(c(terra::rast(ts1), terra::rast(r_end), terra::rast(ts1v2), ts10v1, ts10v2))
# devtools::load_all()
# all.equal(rns1, as.matrix(rns2, wide = T))
# terra::plot(c(terra::rast(rns1), rns2))
#
# bench::mark(a = {set.seed(2022-04-04); kinetic_ising(r_start, B = -0.3, J = 0.7, rule = "metropolis")},
#             b = {set.seed(2022-04-04); ts1v2 = kinetic_ising(terra::rast(r_start), B = -0.3, J = 0.7, rule = "metropolis", version = 2); terra::as.matrix(ts1v2, wide = TRUE)})

