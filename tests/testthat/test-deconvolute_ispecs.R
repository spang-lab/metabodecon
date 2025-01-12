library(testthat)

check_prarp <- test_that("deconvolute_ispecs works", {
    # Call deconvolute_ispecs
    ispecs <- as_ispecs(sim[1:2])
    idecons <- deconvolute_ispecs(ispecs, nfit = 3, wshw = 0.00, sfr = c(3.55, 3.35), verbose = FALSE)

    # Test for correct types
    expect_identical(class(idecons), "idecons")
    expect_identical(class(idecons[[1]]), "idecon")
    expect_identical(class(idecons[[2]]), "idecon")
    expect_identical(names(idecons), c("sim_01", "sim_02"))
    expect_identical(names(idecons[[1]]), idecon_members)
    expect_identical(names(idecons[[2]]), idecon_members)

    # Test for good PRARP
    truepars <- lapply(sim[1:2], function(x) x$meta$simpar)
    obj1 <- calc_prarp(idecons[[1]], truepar = sim[[1]]$meta$simpar)
    obj2 <- calc_prarp(idecons[[2]], truepar = sim[[2]]$meta$simpar)
    expect_true(obj1$prarpx > 0.777)
    expect_true(obj2$prarpx > 0.750)

    # Plot if testing interactively
    if (environment() %==% .GlobalEnv) {
        plot_spectrum(idecons[[1]])
        plot_spectrum(idecons[[2]])
    }
})

skip_if_slow_tests_disabled()

check_bwc <- test_that("deconvolute_ispecs is backwards compatible to MetaboDecon1D", {

    # Inputs
    sim_subset <- metabodecon_file("sim_subset")
    spectra <- read_spectra(sim_subset)
    ispecs <- as_ispecs(spectra)

    # Decovolute with old function
    decons0_old <- MetaboDecon1D_silent( # takes 0.26 sec on tux15
        filepath = sim_subset,
        number_iterations = 3,
        range_water_signal_ppm = 0.00,
        signal_free_region = c(3.55, 3.35),
        debug = FALSE
    )
    decons1_old <- as_decons1(x = decons0_old, spectra = spectra)

    # Decovolute with new function
    idecons_new <- deconvolute_ispecs( # takes 0.14 sec on tux15
        ispecs = ispecs,
        nfit = 3,
        wshw = 0.00,
        sfr = c(3.55, 3.35),
        verbose = FALSE,
        bwc = 0 # backwards compatibility to version 0.x
    )
    decons0_new <- as_decons0(x = idecons_new, spectra = spectra)
    decons1_new <- as_decons1(x = idecons_new, spectra = spectra)

    # Checks
    expect_equal(decons0_old[[1]], decons0_new[[1]])
    expect_equal(decons0_old[[2]], decons0_new[[2]])
    expect_equal(decons1_old[[1]], decons1_new[[1]])
    expect_equal(decons1_old[[2]], decons1_new[[2]])
})
