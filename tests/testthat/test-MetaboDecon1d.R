single_spectrum <- test_that("deconvolute works for 1 spectrum", {

    x1_old <- generate_lorentz_curves_sim(sim[[1]])
    x2_new <- deconvolute(sim[[1]])
    x1_new <- as_decon1(x2_new)

    prarp_x1_old <- calc_prarp(x1_old, truepar = sim[[1]]$meta$simpar)
    prarp_x2_new <- calc_prarp(x2_new, truepar = sim[[1]]$meta$simpar)
    prarp_x1_new <- calc_prarp(x1_new, truepar = sim[[1]]$meta$simpar)

    plot_spectrum(x1_old, main = "Old", foc_rgn = c(0.75, 0.25))
    plot_spectrum(x1_new, main = "Old", foc_rgn = c(0.75, 0.25))

    expect_identical(object = names(x2_new), expected = decon1_members)
    expect_identical(object = class(x2_new), expected = "decon1")
    expect_true(prarp$peak_ratio > 0.9)
    expect_true(prarp$area_ratio > 0.9)
})
