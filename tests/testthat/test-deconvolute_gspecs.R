library(testthat)

sim_subset <- as_gspecs(read_spectra(metabodecon_file("sim_subset")))
sim_one <- sim_subset[[1]]

test_that("deconvolute_gspecs for: 1 spec, 1 core", {
    obj <- deconvolute_gspecs(sim_one)
    expect_true(inherits(obj, "gdecons"))
    expect_equal(length(obj), length(sim_subset))
})

test_that("deconvolute_gspecs for: 1 spec, 2 cores", {
    obj <- deconvolute_gspecs(sim_one, nworkers = 2)
    expect_true(inherits(obj, "gspecs"))
    expect_equal(length(obj), length(sim_subset))
})

test_that("deconvolute_gspecs for: 2 specs, 1 core", {
    nworkers <- 1
    obj <- deconvolute_gspecs(sim_subset, nfit, smopts, delta, sfr, wshw, ask, force, verbose, bwc, nworkers)
    expect_true(inherits(obj, "gspecs"))
    expect_equal(length(obj), length(sim_subset))
})
