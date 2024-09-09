library(testthat)

test_that("speaq_align works", {
    spectrum_data <- generate_lorentz_curves_sim("bruker/sim")
    feat <- gen_feat_mat(spectrum_data)
    maxShift <- 200
    si_size_real_spectrum <- length(spectrum_data[[1]]$y_values)

    # Check that the function returns a matrix with the correct dimensions
    M <- speaq_align(feat, maxShift, spectrum_data, si_size_real_spectrum, verbose = FALSE, show = TRUE)
    expect_equal(nrow(M), length(spectrum_data))
    expect_equal(ncol(M), length(spectrum_data$sim_01$y_values))

    # Check that the function returns the same result as the previous version (which was slower and more complex)
    output <- capture.output(
        M_v1 <- speaq_align_v1(feat, maxShift, spectrum_data, si_size_real_spectrum)
    )
    expect_equal(M_v1, M)
})
