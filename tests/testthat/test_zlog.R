test_that("zlog rejects non-numeric input", {
  expect_error(zlog("abc", limits = c(2, 8)), "x must be numeric")
})

test_that("zlog calculates values for a small dataset", {
  x <- 2
  expect_equal(zlog(x, limits = c(2, 8)), -1.959964 ,tolerance = 1e-6)
})
