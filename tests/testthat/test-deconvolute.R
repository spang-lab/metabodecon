sap <- test_that("deconvolute works for single spectrum", {
    decon2 <- deconvolute(sap[[1]], nfit = 3, sfr = c(3.2, -3.2), smopts = c(1, 3), delta = 3)
    expect_identical(object = names(decon2), expected = decon2_members)
    expect_identical(object = class(decon2), expected = "decon2")
    obj2 <- calc_prarp(x = decon2, truepar = sap[[1]]$meta$simpar)
    expect_true(obj2$prarpx >= 0.961) # MetaboDecon1D has a PRARPX of 0.507. See test-MetaboDecon1d.R.
})

sim_subset <- test_that("deconvolute works for multiple spectra", {
    decons2 <- deconvolute(x = sim[1:2], sfr = c(3.55, 3.35))
    expect_identical(class(decons2), "decons2")
    expect_identical(class(decons2[[1]]), "decon2")
    expect_identical(class(decons2[[2]]), "decon2")
    expect_identical(names(decons2), c("sim_01", "sim_02"))
    expect_identical(names(decons2[[1]]), decon2_members)
    expect_identical(names(decons2[[2]]), decon2_members)
    obj1 <- calc_prarp(decons2[[1]], truepar = sim[[1]]$meta$simpar)
    obj2 <- calc_prarp(decons2[[2]], truepar = sim[[2]]$meta$simpar)
    expect_true(obj1$prarpx >= 0.777) # MetaboDecon1D has a PRARPX of 0.732. See test-MetaboDecon1d.R.
    expect_true(obj2$prarpx >= 0.750) # MetaboDecon1D has a PRARPX of 0.710. See test-MetaboDecon1d.R.
})

wrong_sfr <- test_that("deconvolute works when no peaks are filtered out", {
    x <- simulate_spectrum(ndp = 256, npk = 3)
    expect_error(deconvolute(x, sfr = c(Inf, -Inf), smopts = c(0, 3)))
    deconForc <- deconvolute(x, sfr = c(Inf, -Inf), smopts = c(0, 3), force = TRUE)
})
