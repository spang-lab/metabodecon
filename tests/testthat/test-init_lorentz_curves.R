test_that("init_lorentz_curves_v14 works", {
    sim <- metabodecon_file("sim_subset")
    obj <- MetaboDecon1D_silent(
        filepath = sim,
        filename ="sim_01",
        range_water_signal_ppm = 0,
        signal_free_region = c(3.52, 3.37),
        debug = TRUE
    )
    spec <- within(list(), {
        sdp <- obj$x_values
        y_smooth <- obj$debuglist$smooth$spectrum_y
        peak <- within(list(), {
            center <- obj$debuglist$peakscore$filtered_peaks + 1
            right <- obj$debuglist$peakscore$filtered_left_position + 1
            left <- obj$debuglist$peakscore$filtered_right_position + 1
            high <- rep(TRUE, length(center))
        })
    })
    lcpar <- init_lc(spec, verbose = FALSE)
    expect_equal(lcpar$A, obj$debuglist$parinit$A)
    expect_equal(lcpar$lambda, obj$debuglist$parinit$lambda)
    expect_equal(lcpar$w, obj$debuglist$parinit$w)
})
