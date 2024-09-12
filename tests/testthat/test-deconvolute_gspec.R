library(testthat)

test_that("deconvolute works", {

    # Define Inputs
    gspec <- as_gspec(metabodecon_file("sim_subset/sim_01"))
    nfit <- 3
    smopts <- c(1, 5)
    delta <- 0.1
    sfr <- c(3.58, 3.42)
    wshw <- 0
    force <- FALSE
    bwc <- 1

    # Call Function
    obj <- deconvolute_gspec(gspec, nfit, smopts, delta, sfr, wshw, force, bwc)

    # Test Outputs
    expect_true(inherits(obj, "gdecon"))

    # TODO: obtain true parameters used to simulate spectrum and calculate PRARP
})
