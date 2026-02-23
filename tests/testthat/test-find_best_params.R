testthat::test_that("optimize_settings aligns with find_best_params", {
    testthat::skip_if_not(metabodecon::check_mdrb())
    spec <- metabodecon::sim[[1]]
    sfr <- stats::quantile(spec$cs, c(0.9, 0.1))

    best <- find_best_params(
        spec,
        sfr = sfr,
        smopts = c(2, 5),
        delta = 6.4,
        nfit = 3,
        verbose = FALSE
    )
    mse_grid <- mdrb_mse_for_params(
        spec,
        sfr = sfr,
        smopts = best$smopts,
        delta = best$delta,
        nfit = best$nfit
    )

    dec <- mdrb::Deconvoluter$new()
    ref <- mdrb::Spectrum$new(spec$cs, spec$si, sfr)
    mse_opt <- dec$optimize_settings(ref)
    sm <- dec$smoothing_settings()
    sel <- dec$selection_settings()
    fit <- dec$fitting_settings()
    mse_opt2 <- mdrb_mse_for_params(
        spec,
        sfr = sfr,
        smopts = c(sm$iterations, sm$window_size),
        delta = sel$threshold,
        nfit = fit$iterations
    )

    testthat::expect_true(is.finite(mse_opt))
    testthat::expect_true(is.finite(mse_grid))
    testthat::expect_lte(mse_opt2, mse_grid + abs(mse_grid) * 0.05 + 1e-6)
})
