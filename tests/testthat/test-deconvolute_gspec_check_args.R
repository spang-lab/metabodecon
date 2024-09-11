library(testthat)

test_that("deconvolute works", {

    # Initialize Helpers Variables
    p <- metabodecon_file("sim_subset/sim_01")
    s <- read_spectrum(p)

    # Define Inputs
    x <- as_gspec(s)
    nfit <- 3
    smopts <- c(1, 5)
    delta <- 0.1
    sfr <- c(3.58, 3.42)
    wsr <- c(3.50, 3.50)
    rtyp <- "decons3"
    rm_ws_version <- 1
    force <- FALSE

    # Call Function
    args <- named(x, nfit, smopts, delta, sfr, wsr, rtyp, rm_ws_version, force)
    runf <- function(...) do.call(check_args, set(args, ...))
    # styler: off
    expect_no_error(runf())
    expect_error(runf(x      = 0  ), "x must be a spectrum or gspec, not .*")
    expect_error(runf(nfit   = "a"), "nfit must be a integer.1., not .*")
    expect_error(runf(smopts = "a"), "smopts must be a integer.2., not .*")
    expect_error(runf(delta  = log), "delta must be a numeric.1., not .*")
    expect_error(runf(sfr    = 0  ), "sfr must be a numeric.2., not .*")
    expect_error(runf(wsr    = NA ), "wsr must be a numeric.2., not .*")
    expect_error(runf(rtyp   = 0  ), "rtyp must match 'decons.1-3.', but is 0")
    expect_error(runf(rm_ws_version  = NA ), "rm_ws_version must be a integer.1., not .*")
    expect_error(runf(force  = "a"), "force must be a logical.1., not .*")
    # styler: on

    # Test Outputs
    expect_true(inherits(x, "gspec"))
    expect_true(inherits(obj, rtyp))
    expect_true(length(obj, length(x)))
})
