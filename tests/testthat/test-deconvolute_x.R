library(testthat)

test_that("deconvolute works", {

    # Define Inputs
    x <- read_spectrum(metabodecon_file("sim_subset/sim_01"))
    nfit <- 3
    smopts <- c(1, 5)
    delta <- 0.1
    sfr <- c(3.58, 3.42)
    wsr <- c(3.50, 3.50)
    rtyp <- "decons3"
    rmwsv <- 1
    force <- FALSE

    # Call Function
    obj <- deconvolute_x(x, nfit, smopts, delta, sfr, wsr, rtyp, rmwsv, force)

    # Test Outputs
    expect_true(inherits(obj, rtyp))
    expect_true(length(obj, length(x)))
})
