library(testthat)

test_that("deconvolute works", {
    
    gspec <- as_gspec(metabodecon_file("sim_subset/sim_01"))
    gdecon <- deconvolute_gspec(gspec, verbose = FALSE)
    expect_true(inherits(gdecon, "gdecon"))
})
