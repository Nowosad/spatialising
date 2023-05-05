set.seed(2022-04-04)
data(r_start, package = "spatialising")
data(r_end, package = "spatialising")
ts1 = spatial_ising(r_start, B = -0.3, J = 0.7, rule = "metropolis")
ts1v2 = spatial_ising(r_start, B = -0.3, J = 0.7, rule = "metropolis", version = 2)
# ts1v2 = terra::as.matrix(ts1v2, wide = TRUE)
# attributes(ts1v2) = NULL

ts10 = spatial_ising(r_start, B = -0.3, J = 0.7, updates = 10)

test_that("calculations on matrices works", {
  expect_equal(dim(ts1), c(50, 50))
  expect_equal(dim(ts10), c(50, 50, 10))
  expect_equal(ts1, r_end)
  # expect_equal(ts1v2, r_end)
})
