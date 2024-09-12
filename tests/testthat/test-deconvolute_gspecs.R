library(testthat)

test_that("deconvolute works", {
    gspecs <- as_gspecs(read_spectra(metabodecon_file("bruker/urine")))
    nfit <- 3
    smopts <- c(1, 5)
    delta <- 0.1
    sfr <- c(3.58, 3.42)
    wshw <- 0
    ask <- FALSE
    nworkers <- 2
    verbose <- TRUE
    force <- FALSE
    bwc <- 1
    rt(obj <- deconvolute_gspecs(
        gspecs, nfit, smopts, delta, sfr, wshw,
        ask, force, verbose, bwc,
        nworkers
    ))
    expect_true(inherits(x, "spectra"))
    expect_true(inherits(obj, rtyp))
    expect_true(length(obj, length(x)))
})
