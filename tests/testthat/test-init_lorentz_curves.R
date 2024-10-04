test_that("init_lorentz_curves_v14 works", {

    obj <- md1d("sim_01", cache = TRUE)$rv

    spec <- within(list(), {
        sdp <- obj$x_values
        y_smooth <- obj$debuglist$smoothed$spectrum_y
        peak <- within(list(), {
            center <- obj$debuglist$peak_scores_calc$filtered_peaks + 1
            right <- obj$debuglist$peak_scores_calc$filtered_left_position + 1
            left <- obj$debuglist$peak_scores_calc$filtered_right_position + 1
            high <- rep(TRUE, length(center))
        })
    })
    y14 <- init_lc(spec, verbose = FALSE)
    expect_equal(y14$A, obj$debuglist$params_init$A)
    expect_equal(y14$lambda, obj$debuglist$params_init$lambda)
    expect_equal(y14$w, obj$debuglist$params_init$w)
})
