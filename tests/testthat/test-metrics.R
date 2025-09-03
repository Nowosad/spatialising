set.seed(2022-04-04)
data(r_start, package = "spatialising")
data(r_end, package = "spatialising")
ts1 = kinetic_ising(r_start, B = -0.3, J = 0.7)
ts10 = kinetic_ising(r_start, B = -0.3, J = 0.7, updates = 10)
ts_oneclass = matrix(rep(1, 9), nrow = 3)

m1 = composition_index(ts1)
m2 = composition_index(ts10)
m3 = composition_index(ts_oneclass)

ti1 = texture_index(ts1)
ti2 = texture_index(ts10)
ti3 = texture_index(ts_oneclass)

test_that("calculations of metrics works", {
  expect_equal(length(m1), 1)
  expect_equal(length(m2), 10)
  expect_equal(length(m3), 1)
  expect_true(ti3 == 1)
})
