library(testthat)

test1 <- test_that("deconvolute_ispecs works", {
    sim_subset <- sim[1:2]
    ispecs <- as_ispecs(sim_subset)
    idecons <- deconvolute_ispecs( # takes 0.14 sec on tux15
        ispecs = ispecs,
        nfit = 3,
        wshw = 0.00, # (1)
        sfr = c(3.55, 3.35),
        verbose = FALSE,
        bwc = 2 # (2)
    )
    truepars <- lapply(spectra, function(x) x$meta$simpar)
    prarps <- mapply(calc_prarp, decons1_new, truepars)
    expect_true(all(prarps > 0.8))
    expect_equal(class(idecons), "idecons")
    expect_equal(length(idecons), 2)
})

test2 <- test_that("deconvolute_ispecs is backwards compatible to MetaboDecon1D", {
    sim_subset <- metabodecon_file("sim_subset")
    spectra <- read_spectra(sim_subset)
    ispecs <- as_ispecs(spectra)
    decons0_old <- MetaboDecon1D_silent( # takes 0.26 sec on tux15
        filepath = sim_subset,
        number_iterations = 3,
        range_water_signal_ppm = 0.00, # (1)
        signal_free_region = c(3.55, 3.35),
        debug = TRUE
    )
    decons1_old <- as_decons1(x = decons0_old, spectra = spectra)
    idecons_new <- deconvolute_ispecs( # takes 0.14 sec on tux15
        ispecs = ispecs,
        nfit = 3,
        wshw = 0.00, # (1)
        sfr = c(3.55, 3.35),
        verbose = FALSE,
        bwc = 0 # (2)
    )
    # (1) The old MetaboDecon1D function has a bug in the calculation of the Water
    # Signal Range, that causes the borders to be shifted by 1-2 datapoints. I.e.
    # even for WSHW = 0, one datapoint gets set to zero. This introduces a gap in
    # the spectrum and introduces additional peaks in the deconvoluted spectrum.
    # Therefore, to make the comparison fair, we need to set a WSHW for the new
    # function as well.
    # (2) Set bwc to 0 to achieve backwards compatibility with MetaboDecon1D.
    decons0_new <- as_decons0(x = idecons_new) # CONTINUE HERE: implement!
    decons1_new <- as_decons1(x = idecons_new)
    plot_spectrum(decons1_old[[1]], foc_rgn=c(0.55, 0.45))
    plot_spectrum(decons1_new[[1]], foc_rgn=c(0.55, 0.45))
    expect_equal(decons0_old[[1]], decons_new[[1]])
    expect_equal(decons0_old[[2]], decons_new[[2]])
})
