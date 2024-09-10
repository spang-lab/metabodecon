library(testthat)

test_that("deconvolute works", {
    x <- read_spectra(metabodecon_file("sim_subset"))
    rtyp <- "decons3";
    sfrs <- list(sim_01 = c(3.58, 3.42), sim_02 = c(3.58, 3.42));
    wsrs <- list(sim_01 = c(3.50, 3.50), sim_02 = c(3.50, 3.50));
    nfit <- 3; smopts <- c(1, 5); delta <- 0.1;
    ask <- FALSE; nworkers <- 1; verbose <- TRUE; force <- FALSE;
    obj <- deconvolute_xx(x, rtyp, sfrs, wsrs, nfit, smopts, delta, ask, nworkers, verbose, force)
    expect_true(inherits(x, "spectra"))
    expect_true(inherits(obj, rtyp))
    expect_true(length(obj, length(x)))
})
