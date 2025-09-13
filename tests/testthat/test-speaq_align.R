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

    skip_if_speaq_deps_missing()

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

    skip_if_speaq_deps_missing()

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
    old_vals <- old[!is.na(old)]
    new_vals <- new[!is.na(new)]
    expect_equal(old, new, tolerance = 1e-5) # (1)
    # (1) With version 1.5.0 of metabodecon, the integral calculation done in
    # speaq_align was refactored to use `A * pi` instead of `A * (atan((n - p) /
    # l) - atan((0 - p) / l))` (with limits from 0 to n). This means the values
    # between speaq_align and speaq_align_original are not identical anymore and
    # we need to set a higher tolerance (1e-5 in this case instead of the
    # appox. default of 1e-8).
})
