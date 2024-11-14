library(testthat)

test_that("deconvolute_ispecs works: ", {
    sim_subset <- as_ispecs(read_spectra(metabodecon_file("sim_subset")))
    sim_one <- sim_subset[[1]]
    d1 <- deconvolute_ispecs(sim_one, nworkers = 1, verbose = FALSE)
    d2 <- deconvolute_ispecs(sim_one, nworkers = 2, verbose = FALSE)
    d3 <- deconvolute_ispecs(sim_subset, nworkers = 1, verbose = FALSE)
    d4 <- deconvolute_ispecs(sim_subset, nworkers = 2, verbose = FALSE)
    expect_true(inherits(d1, "idecons"))
    expect_true(inherits(d3, "idecons"))
    expect_equal(length(d1), 1)
    expect_equal(length(d3), length(sim_subset))
    d2$sim_01$dcp$nworkers <- 1
    d4$sim_01$dcp$nworkers <- 1
    d4$sim_02$dcp$nworkers <- 1
    expect_identical(d1, d2)
    expect_identical(d3, d4)
})
