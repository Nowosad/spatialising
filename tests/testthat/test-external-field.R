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
  field_matrix = matrix(c(-20, 20, -20, 20,
                           20, 20, -20, -20,
                          -20, -20, 20, 20),
                        nrow = 3, byrow = TRUE)
  x_matrix = matrix(-1, nrow = 3, ncol = 4)
  x_matrix[2, 2] = NA
  field_matrix[2, 2] = NA
  x = terra::rast(x_matrix)
  field = terra::rast(field_matrix)
  expected = terra::values(terra::rast(ifelse(field_matrix > 0, 1, -1)))

  for (rule in c("glauber", "metropolis")) {
    set.seed(2)
    v1 = kinetic_ising(x, B = 0, J = 0, iter = 200, rule = rule,
                       external_field = field)
    set.seed(2)
    v2 = kinetic_ising(x, B = 0, J = 0, iter = 200, rule = rule,
                       version = 2, external_field = field)

    expect_equal(terra::values(v1), terra::values(v2))
    expect_equal(unname(terra::values(v1)), unname(expected))
  }
})

test_that("external field dimensions must match x", {
  x = matrix(rep(c(-1, 1), 8), nrow = 4)
  expect_error(
    kinetic_ising(x, B = -0.3, J = 0.7, external_field = matrix(0, 2, 2)),
    "type and geometry"
  )

  r = terra::rast(x)
  field = terra::rast(nrows = 2, ncols = 2)
  expect_error(
    kinetic_ising(r, B = -0.3, J = 0.7, external_field = field),
    "type and geometry"
  )

  field = matrix(0, nrow = 4, ncol = 4)
  field[1, 1] = NA
  expect_error(
    kinetic_ising(x, B = -0.3, J = 0.7, external_field = field),
    "cannot be NA"
  )
})
