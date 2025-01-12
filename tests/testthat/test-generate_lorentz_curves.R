# IMPORTANT: we dont test for jcampdx files because after calling
# `read_spectrum()`, the data is the same as for bruker, which is tested in
# `test-read_spectrum.R`. Also, the calculations in the old `MetaboDecon1D()`
# function are slightly different for jcampdx and bruker files (in the jcampdx
# case it calculates with n-1 instead of n) and our `compare_spectra` function
# currently only accounts for bruker-type errors of MetaboDecon1D, but not
# jcampdx errors.

sap2 <- test_that("GLC works for single spectrum", {
    decon1 <- generate_lorentz_curves(
        data_path = sap2,
        nfit = 3,
        sfr = c(3.2, -3.2),
        wshw = 0,
        smopts = c(1, 3),
        delta = 3,
        ask = FALSE,
        verbose = FALSE
    )
    expect_identical(object = names(decon1), expected = decon1_members)
    expect_identical(object = class(decon1), expected = "decon1")
    decon2 <- as_decon2(decon1, spectrum = sap2, sfr = c(3.2, -3.2), wshw = 0)
    obj2 <- calc_prarp(x = decon2, truepar = sap2$meta$simpar)
    expect_true(obj2$prarpx >= 0.961) # MetaboDecon1D has a PRARPX of 0.507. See test-MetaboDecon1d.R.
})

sim_subset <- test_that("MetaboDecon1D works for multiple spectra", {
    decons1 <- generate_lorentz_curves(
        data_path = sim[1:2],
        nfit = 3,
        sfr = c(3.55, 3.35),
        wshw = 0,
        ask = FALSE,
        verbose = FALSE
    )
    expect_identical(names(decons1), c("sim_01", "sim_02"))
    expect_identical(class(decons1), "decons1")
    expect_identical(names(decons1[[1]]), decon1_members)
    expect_identical(class(decons1[[1]]), "decon1")
    expect_identical(names(decons1[[2]]), decon1_members)
    expect_identical(class(decons1[[2]]), "decon1")
    decons2 <- as_decons2(decons1, spectra = sim[1:2])
    obj1 <- calc_prarp(decons2[[1]], truepar = sim[[1]]$meta$simpar)
    obj2 <- calc_prarp(decons2[[2]], truepar = sim[[2]]$meta$simpar)
    expect_true(obj1$prarpx >= 0.777) # MetaboDecon1D has a PRARPX of 0.732. See test-MetaboDecon1d.R.
    expect_true(obj2$prarpx >= 0.750) # MetaboDecon1D has a PRARPX of 0.710. See test-MetaboDecon1d.R.
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
