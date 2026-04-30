library(testthat)

# Stashed: tests for decon1/decons1 collection methods.
# Moved here when depr.R (which defines generate_lorentz_curves_sim returning
# decon1 objects) was removed from the project.

test_that("c methods work for decon1/decons1 collections", {
    d1 <- generate_lorentz_curves_sim(sim[[1]])

    cc1 <- c(d1, d1)
    expect_true(is_decons1(cc1))
    expect_equal(length(cc1), 2)
})
