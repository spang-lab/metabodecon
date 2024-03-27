# skip_if_slow_tests_disabled()

test_that("init_lorentz_curves_v13 works", {

    obj <- MD1D(cache = TRUE)$rv

    x <- obj$x_values
    y <- obj$debuglist$smoothed$spectrum_y
    pc <- obj$debuglist$peak_scores_calc$filtered_peaks + 1
    pr <- obj$debuglist$peak_scores_calc$filtered_left_position + 1
    pl <- obj$debuglist$peak_scores_calc$filtered_right_position + 1

    y3 <- init_lorentz_curves_v13(x, y, pc, pr, pl)

    expect_equal(y3$A, obj$debuglist$params_init$A)
    expect_equal(y3$lambda, obj$debuglist$params_init$lambda)
    expect_equal(y3$w, obj$debuglist$params_init$w)
})
