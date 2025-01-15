library(testthat)

x <- generate_lorentz_curves(
    data_path = sim[1:2],
    sfr = c(3.55, 3.35),
    wshw = 0,
    ask = FALSE,
    verbose = FALSE
)
feat_mat <- gen_feat_mat(x)

test_that("speaq_align works", {
    aligned_mat <- speaq_align(
        feat = feat_mat,
        maxShift = 200,
        spectrum_data = x,
        si_size_real_spectrum = length(spectrum_data[[1]]$y_values),
        verbose = FALSE,
        show = FALSE
    )
    expect_equal(nrow(aligned_mat), length(x))
    expect_equal(ncol(aligned_mat), length(x$sim_01$y_values))
})

test_that("speaq_align is backwards compatible to speaq_align_original", {
    local_output_sink(nullfile())
    new <- speaq_align(
        feat = feat_mat,
        maxShift = 200,
        spectrum_data = x,
        si_size_real_spectrum = length(x[[1]]$y_values)
    )
    old <- speaq_align_original(
        feat = feat_mat,
        maxShift = 200,
        spectrum_data = x,
        si_size_real_spectrum = length(x[[1]]$y_values)
    )
    if (identical(environment(), globalenv())) deferred_run()
    expect_equal(colnames(old), colnames(new))
    expect_equal(rownames(old), rownames(new))
    expect_true(all.equal(old, new))
})
