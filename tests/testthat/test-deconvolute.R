single_spectrum <- test_that("deconvolute works for 1 spectrum", {
    x <- deconvolute(x = sim[[1]])
    x1 <- as_decon1(x)
    prarp <- calc_prarp(decon = x, truepar = sim[[1]]$meta$simpar)
    expect_identical(object = names(x), expected = decon1_members)
    expect_identical(object = class(x), expected = "decon1")
    expect_true(prarp$peak_ratio > 0.9)
    expect_true(prarp$area_ratio > 0.9)
})

bruker_folder <- test_that("GLC works for bruker folder", {
    data_path <- metabodecon_file("bruker/sim_subset")
    x <- generate_lorentz_curves(
        data_path,
        sfr = c(3.55, 3.35),
        wshw = 0,
        ask = FALSE,
        verbose = FALSE
    )
    expect_identical(object = names(x), expected = c("sim_01", "sim_02"))
    expect_identical(object = class(x), expected = "decons1")
    expect_identical(object = class(x$sim_01), expected = "decon1")
    expect_identical(object = class(x$sim_02), expected = "decon1")
    prarps <- mapply(calc_prarp, x, list(sim[[1]]$meta$simpar, sim[[2]]$meta$simpar))
    expect_true(all(prarps > 0.8))
})

wrong_sfr <- test_that("GLC works when no peaks are filtered out", {
    x <- simulate_spectrum(ndp = 256, npk = 3)
    expect_error(generate_lorentz_curves(
        x, sfr = c(Inf, -Inf), wshw = 0, smopts = c(0, 3), ask = FALSE
    ))
    decon <- generate_lorentz_curves(
        x, sfr = c(Inf, -Inf), wshw = 0, smopts = c(0, 3), ask = FALSE, force = TRUE
    )
    expect_identical(length(decon), 32L)
})
