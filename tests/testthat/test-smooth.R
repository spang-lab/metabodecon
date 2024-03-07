library(testthat)

test_that("1. smooth(fast = FALSE) works", {
  y <- c(1, 4, 2, 4, 3)
  exp1 <- sapply(list(y[1:2], y[1:3], y[2:4], y[3:5], y[4:5]), mean)
  exp2 <- sapply(list(exp1[1:2], exp1[1:3], exp1[2:4], exp1[3:5], exp1[4:5]), mean)
  exp3 <- sapply(list(exp2[1:2], exp2[1:3], exp2[2:4], exp2[3:5], exp2[4:5]), mean)
  z1 <- smooth(y = y, k = 3, reps = 1, fast = FALSE)
  z2 <- smooth(y = y, k = 3, reps = 2, fast = FALSE)
  z3 <- smooth(y = y, k = 3, reps = 3, fast = FALSE)
  expect_equal(z1, exp1)
  expect_equal(z2, exp2)
  expect_equal(z3, exp3)
})

test_that("2. smooth(fast = TRUE) works", {
  y <- c(1, 4, 2, 4, 3)
  exp1 <- sapply(list(y[1:2], y[1:3], y[2:4], y[3:5], y[4:5]), mean)
  exp2 <- sapply(list(exp1[1:2], exp1[1:3], exp1[2:4], exp1[3:5], exp1[4:5]), mean)
  exp3 <- sapply(list(exp2[1:2], exp2[1:3], exp2[2:4], exp2[3:5], exp2[4:5]), mean)
  z1 <- smooth(y = y, k = 3, reps = 1, fast = TRUE)
  z2 <- smooth(y = y, k = 3, reps = 2, fast = TRUE)
  z3 <- smooth(y = y, k = 3, reps = 3, fast = TRUE)
  expect_equal(z1, exp1)
  expect_equal(z2, exp2)
  expect_equal(z3, exp3)
})

test_that("3. smooth(fast = TRUE) is faster than smooth(fast = FALSE)", {
    y <- rnorm(5000)
    runtime_v1 <- system.time(z_v1 <- smooth(y = y, k = 5, reps = 10, fast = FALSE))
    runtime_v2 <- system.time(z_v2 <- smooth(y = y, k = 5, reps = 10, fast = TRUE))
    expect_true(runtime_v1[3] > runtime_v2[3] * 10) # at least 10 times faster
    expect_equal(z_v1, z_v2)
    # In manual testting the speedup was factor 27 for y <- rnorm(130000). There we had runtime_v1[3] == 1.09 and runtime_v2[3] == 0.04.
})
