library(testthat)

test_that("deconvolute_gspecs works", {
    spectra <- sim[1:2]
    truepars <- lapply(spectra, function(x) x$meta$simpar)
    gspecs <- as_gspecs(spectra)
    decons <- deconvolute_gspecs(
        gspecs = gspecs,
        nfit = 3, sfr = c(3.55, 3.35), wshw = 0,
        ask = FALSE, force = FALSE, verbose = TRUE, rtyp = "decons1"
    )
    prarps <- mapply(calc_prarp, decons, truepars)
    expect_true(all(prarps > 0.8))
    calc_prarp(decons[[1]], truepars[[1]], show = TRUE)
    expect_equal(class(decons), "decons1")
    expect_equal(length(decons), 2)
})

test_that("deconvolute_gspecs is backwards compatible to MetaboDecon1D", {
    # CONTINUE HERE
})
