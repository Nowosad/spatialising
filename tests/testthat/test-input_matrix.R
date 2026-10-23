set.seed(2022-04-04)
data(r_start, package = "spatialising")
data(r_end, package = "spatialising")
ts1 = kinetic_ising(r_start, B = -0.3, J = 0.7, rule = "metropolis")

set.seed(2022-04-04)
ts1v2 = kinetic_ising(terra::rast(r_start), B = -0.3, J = 0.7, rule = "metropolis", version = 2)
ts1v2 = terra::as.matrix(ts1v2, wide = TRUE)

ts10 = kinetic_ising(r_start, B = -0.3, J = 0.7, updates = 10)

test_that("calculations on matrices works", {
  expect_equal(dim(ts1), c(50, 50))
  expect_equal(dim(ts10), c(50, 50, 10))
  expect_equal(ts1, r_end)
  expect_equal(ts1v2, r_end)
  expect_equal(ts1v2, ts1)
})

test_that("invalid simulation inputs fail clearly", {
  x = matrix(c(-1, 1, 1, -1), nrow = 2)

  expect_error(kinetic_ising(x, 0, 1, rule = "invalid"), "one of")
  expect_error(kinetic_ising(x, 0, 1, updates = 0), "positive integer")
  expect_error(kinetic_ising(x, 0, 1, iter = 0), "positive integer")
  expect_error(kinetic_ising(x, 0, 1, version = 3), "must be 1 or 2")
  expect_error(kinetic_ising(x, 0, 1, version = 2), "requires a SpatRaster")

  r = c(terra::rast(x), terra::rast(x))
  expect_error(kinetic_ising(r, 0, 1), "single-layer")
})

test_that("binary validation does not change the random-number state", {
  set.seed(1)
  state = .Random.seed
  spatialising:::check_if_proper_binary(matrix(c(-1, 1), nrow = 1))
  expect_identical(.Random.seed, state)
})
