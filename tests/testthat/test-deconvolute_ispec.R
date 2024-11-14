library(testthat)

test_that("deconvolute works", {
    
    ispec <- as_ispec(metabodecon_file("sim_subset/sim_01"))
    idecon <- deconvolute_ispec(ispec, verbose = FALSE)
    expect_true(inherits(idecon, "idecon"))
})
