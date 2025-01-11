# IMPORTANT: we dont test for jcampdx files because after calling
# `read_spectrum()`, the data is the same as for bruker, which is tested in
# `test-read_spectrum.R`. Also, the calculations in the old `MetaboDecon1D()`
# function are slightly different for jcampdx and bruker files (in the jcampdx
# case it calculates with n-1 instead of n) and our `compare_spectra` function
# currently only accounts for bruker-type errors of MetaboDecon1D, but not
# jcampdx errors.

single_spectrum <- test_that("GLC works for a single spectrum", {
    x <- generate_lorentz_curves(
        data_path = sim[[1]],
        nfit = 3,
        sfr = c(3.55, 3.35),
        wshw = 0,
        ask = FALSE,
        verbose = FALSE
    )
    expect_identical(object = names(x), expected = decon1_members)
    expect_identical(object = class(x), expected = "decon1")
    obj <- calc_prarp(x, truepar = sim[[1]]$meta$simpar)
    expect_true(obj$prarp  >= 0.836)
    expect_true(obj$prarpx >= 0.779)
})

bruker_folder <- test_that("GLC works for a folder with bruker spectra", {
    data_path <- metabodecon_file("bruker/sim_subset")
    x <- generate_lorentz_curves(
        data_path,
        nfit = 3,
        sfr = c(3.55, 3.35),
        wshw = 0,
        ask = FALSE,
        verbose = FALSE
    )
    expect_identical(object = names(x), expected = c("sim_01", "sim_02"))
    expect_identical(object = class(x), expected = "decons1")
    expect_identical(object = class(x$sim_01), expected = "decon1")
    expect_identical(object = class(x$sim_02), expected = "decon1")

    obj1 <- calc_prarp(x$sim_01, truepar = sim[[1]]$meta$simpar)
    expect_true(obj1$prarp  >= 0.836)
    expect_true(obj1$prarpx >= 0.779)

    obj2 <- calc_prarp(x$sim_02, truepar = sim[[2]]$meta$simpar)
    expect_true(obj2$prarp  >= 0.908)
    expect_true(obj2$prarpx >= 0.709)
})

wrong_sfr <- test_that("GLC works when no peaks are filtered out", {
    x <- simulate_spectrum(ndp = 256, npk = 3)
    expect_error(
        generate_lorentz_curves(
            data_path = x,
            sfr = c(Inf, -Inf),
            wshw = 0,
            smopts = c(0, 3),
            ask = FALSE,
            verbose = FALSE
        )
    )
    decon <- generate_lorentz_curves(
        data_path = x,
        sfr = c(Inf, -Inf),
        wshw = 0,
        smopts = c(0, 3),
        ask = FALSE,
        force = TRUE,
        verbose = FALSE
    )
    expect_identical(length(decon), 32L)
})

broken_A <- simulate_spectrum(
    # TODO: deconvolute this spectrum and check whether we get A's with
    # different signs (I think so.) If yes, fix this. Mixed signs in A (and
    # lambda) should not be possible!
    "sap",
    cs     = seq(0.5, -0.5, by = -0.1),
    x0     = c(0.2),
    A      = 2000,
    lambda = 0.08,
    noise  = c(0, 0, 0, 0, 0, 5000, 0, 1000, 0, 0, 0),
    fqref  = 6e8
)
