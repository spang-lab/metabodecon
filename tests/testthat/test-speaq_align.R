library(testthat)

test_that("speaq_align works", {
    sim_subset <- metabodecon_file("bruker/sim_subset")
    decons <- generate_lorentz_curves_sim(sim_subset)
    feat_mat <- gen_feat_mat(decons)
    aligned_mat <<- speaq_align(
        feat = feat_mat,
        maxShift = 200,
        spectrum_data = decons,
        si_size_real_spectrum = length(spectrum_data[[1]]$y_values),
        verbose = FALSE,
        show = FALSE
    )
    expect_equal(nrow(aligned_mat), length(decons))
    expect_equal(ncol(aligned_mat), length(decons$sim_01$y_values))
})

test_that("speaq_align is backwards compatible to speaq_align_original", {
    withr::local_output_sink(nullfile())
    aligned_mat_original <- speaq_align_original(
        feat = feat_mat,
        maxShift = 200,
        spectrum_data = decons,
        si_size_real_spectrum = length(decons[[1]]$y_values)
    )
    expect_equal(aligned_mat_original, aligned_mat)
})
