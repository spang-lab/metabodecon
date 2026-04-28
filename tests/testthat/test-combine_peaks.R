library(testthat)

test_that("combine_peaks merges neighbouring compatible columns", {
    M <- rbind(
        c(0.13, 0.00, 0.00, 0.11, 0.00),
        c(0.00, 0.88, 0.00, 0.12, 0.00),
        c(0.07, 0.56, 0.30, 0.00, 0.00),
        c(0.08, 0.00, 0.07, 0.00, 0.07),
        c(0.04, 0.00, 0.00, 0.04, 0.00)
    )

    obj <- combine_peaks(M, range = 1)
    expect_equal(obj$long[, 3], c(0, 0, 0, 0, 0))
    expect_equal(obj$long[, 4], c(0.11, 0.12, 0.30, 0.07, 0.04))
    expect_equal(ncol(obj$short), 4)
    expect_equal(obj$short[, 3], obj$long[, 4])
})

test_that("combine_peaks with range zero keeps columns and normalizes NAs", {
    M <- rbind(
        c(1, NA, 0),
        c(0, 0, 0),
        c(NA, 2, 0)
    )

    obj <- combine_peaks(M, range = 0)
    exp_long <- rbind(
        c(1, 0, 0),
        c(0, 0, 0),
        c(0, 2, 0)
    )
    exp_short <- exp_long[, 1:2, drop = FALSE]
    expect_equal(obj$long, exp_long)
    expect_equal(obj$short, exp_short)
})

test_that("combine_peaks snaps to reference positions and drops far peaks", {
    M <- rbind(
        c(0, 2, 0, 0, 3, 0, 0),
        c(0, 0, 4, 0, 0, 5, 0)
    )

    obj <- combine_peaks(M, range = 1, ref = c(2, 5))
    exp <- rbind(
        c(0, 2, 0, 0, 3, 0, 0),
        c(0, 4, 0, 0, 5, 0, 0)
    )
    expect_equal(obj$long, exp)
    expect_equal(obj$short, exp[, c(2, 5), drop = FALSE])
})

test_that("combine_peaks snaps midpoint ties to the left reference", {
    M <- matrix(c(0, 0, 0, 0, 7, 0, 0, 0, 0), nrow = 1)

    obj <- combine_peaks(M, range = 2, ref = c(3, 7))
    exp <- matrix(c(0, 0, 7, 0, 0, 0, 0, 0, 0), nrow = 1)
    expect_equal(obj$long, exp)

    obj2 <- combine_peaks(M, range = 1, ref = c(3, 7))
    expect_equal(obj2$long, matrix(0, nrow = 1, ncol = 9))
    expect_length(obj2$short, 0)
})