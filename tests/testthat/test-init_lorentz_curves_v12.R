# skip_if_slow_tests_disabled()

test_that("init_lorentz_curves_v10 works", {

    x <- MD1D(cache = TRUE)$rv

    spec <- list()
    spec$sdp <- x$x_values
    spec$y_smooth <- x$debuglist$smoothed$spectrum_y
    spec$peak$center <- x$debuglist$peak_scores_calc$filtered_peaks + 1
    spec$peak$right <- x$debuglist$peak_scores_calc$filtered_left_position + 1
    spec$peak$left <- x$debuglist$peak_scores_calc$filtered_right_position + 1
    spec$peak$high <- rep(TRUE, length(filtered_peaks))

    y2 <- init_lorentz_curves_v12(spec)$lc

    expect_equal(y2$A, x$debuglist$params_init$A)
    expect_equal(y2$lambda, x$debuglist$params_init$lambda)
    expect_equal(y2$w, x$debuglist$params_init$w)
})
