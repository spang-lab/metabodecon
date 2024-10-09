library(testthat)

test_deconvolute_gspecs

sim_subset <- metabodecon_file("sim_subset")
spectra <- read_spectra(sim_subset)
gspecs <- as_gspecs(spectra)
decons0 <- MetaboDecon1D_silent( # takes 0.26 sec on tux15
    filepath = sim_subset,
    number_iterations = 3,
    range_water_signal_ppm = 0.01, # (1)
    signal_free_region = c(3.55, 3.35)
)
decons1 <- deconvolute_gspecs( # takes 0.14 sec on tux15
    gspecs = gspecs,
    nfit = 3,
    wshw = 0.01, # (1)
    sfr = c(3.55, 3.35),
    rtyp = "decons1",
    verbose = FALSE
)
# (1) The old MetaboDecon1D function has a bug in the calculation of the Water Signal Range, that causes the borders to be shifted by 1-2 datapoints. I.e. even for WSHW = 0, one datapoint gets set to zero. This introduces a gap in the spectrum and introduces additional peaks in the deconvoluted spectrum. Therefore, to make the comparison fair, we need to set a WSHW for the new function as well.

truepars <- lapply(spectra, function(x) x$meta$simpar)
prarps0 <- mapply(calc_prarp, decons0, truepars)
prarps1 <- mapply(calc_prarp, decons1, truepars)

test_that("deconvolute_gspecs works", {
    expect_true(all(prarps > 0.8))
    expect_equal(class(decons), "decons1")
    expect_equal(length(decons), 2)
})

test_that("deconvolute_gspecs is backwards compatible to MetaboDecon1D", {
    for (name in names(decons0[[1]])) {
        expect_equal()
    }
})
