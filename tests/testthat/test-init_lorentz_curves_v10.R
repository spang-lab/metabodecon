# skip_if_slow_tests_disabled()

test_that("init_lorentz_curves_v10 works", {

    x <- MD1D(cache = TRUE)$rv

    spectrum_x <- x$x_values
    spectrum_y <- x$debuglist$smoothed$spectrum_y
    filtered_peaks <- x$debuglist$peak_scores_calc$filtered_peaks
    filtered_left_position <- x$debuglist$peak_scores_calc$filtered_left_position
    filtered_right_position <- x$debuglist$peak_scores_calc$filtered_right_position
    save_scores <- x$debuglist$peak_scores_calc$save_scores

    y <- init_lorentz_curves_v10(spectrum_x, spectrum_y, filtered_peaks, filtered_left_position, filtered_right_position, save_scores)

    expect_equal(y$A, x$debuglist$params_init$A)
    expect_equal(y$lambda, x$debuglist$params_init$lambda)
    expect_equal(y$w, x$debuglist$params_init$w)
    # For a super strange reason, `MetaboDecon1D()` sometimes produces slightly different results (e.g. when called for the first time after package loading) on Windows.
    # The first difference occurs in entry `x1$debuglist$params_init$lambda`, which is a vector of length 1227.
    # Element 273 is sometimes `2.168404e-19` larger than in normal runs.
    # This difference then carries over to all following calculations.
    # Example: when following code is called directly after package loading, all elements are identical, except for the first entry (assuming `MetaboDecon1D()` has not been called before).
    #
    # >>> xs <- lapply(1:10, function(x) MetaboDecon1D_urine1_1010yy_ni1_dbg(overwrite = TRUE)$rv)
    # >>> pairwise_identical(xs)"
    #
    # The isolated calculations done in `init_lorentz_curves_v10` have the same problem. I.e., sometimes they produce `y$lambda[273] == -0.01896324208858969975755` an sometimes they produce `y$lambda[273] == -0.018963242088589703227`.
    # This makes testing very difficult, because we cannot use `identical` anymore.
    #
    # Workaround: use `all.equal` instead of `identical`.
    # We wanted to do that for later versions anyway.
    # But still it's very unsatisfying that we can't find the source of the differences and we cannot use identical anymore.
})
