test_that("max-functionality works", {
  x <- rnorm(10000)
  y <- x
  M <- 100
  expect_gt(ot_indices_1d(data.frame(x), y, M)$indices, 0.95)
})

test_that("zero-indipendence works", {
  x <- data.frame(rnorm(10000))
  y <- rnorm(10000)
  M <- 100
  expect_lt(ot_indices_1d(x, y, M)$indices, 0.05)
})
