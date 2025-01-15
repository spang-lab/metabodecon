# Inputs #####
args0 <- args1 <- args2 <- list(
    ispec   = as_ispec(sap[[1]]),
    nfit    = 3,
    sfr     = c(3.2, -3.2),
    smopts  = c(1, 3),
    delta   = 3,
    bwc     = 0,
    verbose = FALSE
)
args1$bwc <- 1
args2$bwc <- 2

# Call #####
idecon0 <- do.call(deconvolute_ispec, args0)
idecon1 <- do.call(deconvolute_ispec, args1)
idecon2 <- do.call(deconvolute_ispec, args2)

# Checks #####
types <- test_that("returned objects have correct type", {
    expect_identical(object = names(idecon0), expected = idecon_members)
    expect_identical(object = names(idecon1), expected = idecon_members)
    expect_identical(object = names(idecon2), expected = idecon_members)
    expect_identical(object = class(idecon0), expected = "idecon")
    expect_identical(object = class(idecon1), expected = "idecon")
    expect_identical(object = class(idecon2), expected = "idecon")
})

prarps <- test_that("PRARPs are good", {
    obj0 <- calc_prarp(x = idecon0, truepar = sap[[1]]$meta$simpar)
    obj1 <- calc_prarp(x = idecon1, truepar = sap[[1]]$meta$simpar)
    obj2 <- calc_prarp(x = idecon2, truepar = sap[[1]]$meta$simpar)
    expect_true(obj0$prarpx >= 0.507) # MetaboDecon1D has a PRARPX of 0.507. See test-MetaboDecon1d.R.
    expect_true(obj1$prarpx >= 0.961)
    expect_true(obj2$prarpx >= 0.961)
    expect_true(obj2$prarpx >= obj1$prarpx)
    expect_true(obj1$prarpx >= obj0$prarpx)
})

testthat::skip_on_cran()
testthat::skip_on_ci()

bwc <- test_that("Deconvolute with bwc=0 returns the same as MetaboDecon1D", {
    decon0_deconvolute <- as_decon0(idecon0, optional = FALSE)
    decon0_MetaboDecon1D <- MetaboDecon1D_silent(
        filepath = metabodecon_file("bruker/sap"),
        filename = "sap_01",
        number_iterations = 3,
        range_water_signal_ppm = 0,
        signal_free_region = c(3.2, -3.2),
        smoothing_param = c(1, 3),
        delta = 3
    )
    spec <- read_spectrum(metabodecon_file("bruker/sap/sap_01"))
    expect_equal(decon0_deconvolute, decon0_MetaboDecon1D)
})
