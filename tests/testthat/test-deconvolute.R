library(testthat)

test_that("deconvolute works", {
    x <- read_spectra(metabodecon_file("sim_subset"))
    rtyp <- "decons2"; sfr <- c(3.58, 3.42); wshw <- 0;
    nfit <- 3; smopts <- c(1, 5); delta <- 0.1;
    ask <- FALSE; nworkers <- 1; verbose <- TRUE; force <- FALSE;
    obj <- deconvolute_gspecs(x, rtyp, sfr, wshw, nfit, smopts, delta, ask, nworkers, verbose, force)
    expect_true(inherits(x, "spectra"))
    expect_true(inherits(obj, rtyp))
    expect_true(length(obj, length(x)))
})
