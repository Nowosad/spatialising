test_that("a constant external field is added to B", {
  x = matrix(rep(c(-1, 1), 8), nrow = 4)
  field = matrix(0.125, nrow = 4, ncol = 4)

  for (rule in c("glauber", "metropolis")) {
    set.seed(1)
    with_field = kinetic_ising(x, B = -0.25, J = 0.7, updates = 2,
                               rule = rule, external_field = field)
    set.seed(1)
    with_adjusted_B = kinetic_ising(x, B = -0.125, J = 0.7, updates = 2,
                                    rule = rule)
    expect_equal(with_field, with_adjusted_B)
  }
})

test_that("both versions use the same spatial field", {
  x = terra::rast(matrix(rep(c(-1, 1), 8), nrow = 4))
  field = terra::rast(matrix(seq(-0.4, 0.4, length.out = 16), nrow = 4))

  set.seed(2)
  v1 = kinetic_ising(x, B = -0.3, J = 0.7, external_field = field)
  set.seed(2)
  v2 = kinetic_ising(x, B = -0.3, J = 0.7, version = 2,
                     external_field = field)

  expect_equal(terra::values(v1), terra::values(v2))
})

test_that("external field dimensions must match x", {
  x = matrix(rep(c(-1, 1), 8), nrow = 4)
  expect_error(
    kinetic_ising(x, B = -0.3, J = 0.7, external_field = matrix(0, 2, 2)),
    "same dimensions"
  )

  r = terra::rast(x)
  field = terra::rast(nrows = 2, ncols = 2)
  expect_error(
    kinetic_ising(r, B = -0.3, J = 0.7, external_field = field),
    "same geometry"
  )
})
