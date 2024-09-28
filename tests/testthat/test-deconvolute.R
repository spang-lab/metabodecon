library(testthat)

test_that("deconvolute_gspecs works", {
    obj <- deconvolute_gspecs(
        gspecs = as_gspecs(read_spectra(metabodecon_file("sim_subset"))),
        nfit = 3,
        smopts = c(1, 5),
        delta = 0.1,
        sfr = c(3.58, 3.42),
        wshw = 0,
        ask = FALSE,
        force = FALSE,
        verbose = TRUE,
        rtyp = "decons1"
    )
    expect_equal(class(obj), "decons1")
    expect_equal(length(obj), length(x))
})
