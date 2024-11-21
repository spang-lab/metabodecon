single_spectrum <- test_that("deconvolute works for 1 spectrum", {

    x1_old <- generate_lorentz_curves_sim(sim[[1]])
    x2_new <- deconvolute(sim[[1]])
    x1_new <- as_decon1(x2_new)

    prarp_x1_old <- calc_prarp(x1_old, truepar = sim[[1]]$meta$simpar)
    prarp_x2_new <- calc_prarp(x2_new, truepar = sim[[1]]$meta$simpar)
    prarp_x1_new <- calc_prarp(x1_new, truepar = sim[[1]]$meta$simpar)

    plot_spectrum(x1_old, main = "Old", foc_rgn = c(0.75, 0.25))
    plot_spectrum(x1_new, main = "Old", foc_rgn = c(0.75, 0.25))

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
