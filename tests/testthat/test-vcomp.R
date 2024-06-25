library(testthat)

test_that("vcomp correctly identifies identical vectors", {
    capture.output(x <- vcomp(c(1, 2, 3), c(1, 2, 3)))
    expect_identical(x, 0)
})

test_that("vcomp correctly identifies all.equal vectors", {
    capture.output(x <- vcomp(c(1, 2, 3), c(1, 2, 3) + 1e-8))
    expect_identical(x, 1)
})

test_that("vcomp correctly identifies different vectors", {
    capture.output(x <- vcomp(c(1, 2, 3), c(4, 5, 6)))
    expect_identical(x, 2)
})
