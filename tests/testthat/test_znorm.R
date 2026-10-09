test_that("znorm rejects non-numeric input", {
  expect_error(znorm("abc", limits = c(2, 8)), "'x' has to be numeric")
})

test_that("znorm calculates values for a small dataset", {
  x <- 2
  expect_equal(znorm(x, limits = c(2, 8), lambda = 1), -1.959964, tolerance = 1e-6)
})
