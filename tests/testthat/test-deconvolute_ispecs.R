library(testthat)

check_prarp <- test_that("deconvolute_ispecs works", {
    sim_subset <- sim[1:2]
    ispecs <- as_ispecs(sim_subset)
    idecons <- deconvolute_ispecs(
        ispecs = ispecs,
        nfit = 3,
        wshw = 0.00,
        sfr = c(3.55, 3.35),
        verbose = FALSE,
        bwc = 2
    )
    truepars <- lapply(sim_subset, function(x) x$meta$simpar)
    prarps <- mapply(calc_prarp, idecons, truepars)
    expect_true(all(prarps > 0.8))
    expect_equal(class(idecons), "idecons")
    expect_equal(length(idecons), 2)
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
