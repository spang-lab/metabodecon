library(testthat)

test_that("combine_peaks merges neighbouring compatible columns", {
    M <- rbind(
        c(0.13, 0.00, 0.00, 0.11, 0.00),
        c(0.00, 0.88, 0.00, 0.12, 0.00),
        c(0.07, 0.56, 0.30, 0.00, 0.00),
        c(0.08, 0.00, 0.07, 0.00, 0.07),
        c(0.04, 0.00, 0.00, 0.04, 0.00)
    )

    obj <- metabodecon:::combine_peaks(M, range = 1)
    expect_equal(obj[, 3], c(0, 0, 0, 0, 0))
    expect_equal(obj[, 4], c(0.11, 0.12, 0.30, 0.07, 0.04))
    expect_equal(sum(colSums(obj != 0) > 0), 4)
    expect_equal(obj[, which(colSums(obj != 0) > 0)[3]], obj[, 4])
})

test_that("combine_peaks with range zero keeps columns and normalizes NAs", {
    M <- rbind(
        c(1, NA, 0),
        c(0, 0, 0),
        c(NA, 2, 0)
    )

    obj <- metabodecon:::combine_peaks(M, range = 0)
    exp_long <- rbind(
        c(1, 0, 0),
        c(0, 0, 0),
        c(0, 2, 0)
    )
    expect_equal(obj, exp_long)
})

test_that("combine_peaks snaps to reference positions and drops far peaks", {
    M <- rbind(
        c(0, 2, 0, 0, 3, 0, 0),
        c(0, 0, 4, 0, 0, 5, 0)
    )

    obj <- metabodecon:::combine_peaks(M, range = 1, ref = c(2, 5))
    exp <- rbind(
        c(0, 2, 0, 0, 3, 0, 0),
        c(0, 4, 0, 0, 5, 0, 0)
    )
    expect_equal(obj, exp)
})

test_that("combine_peaks snaps midpoint ties to the left reference", {
    M <- matrix(c(0, 0, 0, 0, 7, 0, 0, 0, 0), nrow = 1)

    obj <- metabodecon:::combine_peaks(M, range = 2, ref = c(3, 7))
    exp <- matrix(c(0, 0, 7, 0, 0, 0, 0, 0, 0), nrow = 1)
    expect_equal(obj, exp)

    obj2 <- metabodecon:::combine_peaks(M, range = 1, ref = c(3, 7))
    expect_equal(obj2, matrix(0, nrow = 1, ncol = 9))
    expect_equal(sum(colSums(obj2 != 0) > 0), 0)
})