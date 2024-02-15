# Keep this test in a seperate file, because it's very slow. This is useful
# because testthat runs test files in parallel, i.e. by distributing slow tests
# across different files we ensure they can be also distributed across different
# cores.
library(testthat)

test_that("2. MetaboDecon1D with: 2 jcampdx, answers == y1yy", {

    ## Skip Conditions #####
    skip_on_cran()
    skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

    ## Call function #####
    x <- with(
        testdir = "MetaboDecon1D/2",
        output = "captured",
        message = "captured",
        plots = "plots.pdf",
        inputs = c(urine = "jcampdx/urine"),
        answers = c(
            "y", # Use the same parameters for all spectra?
            "1", # Choose number of file to adjust all parameters
            "y", # Signal free region borders correct selected?
            "y" # Water artefact fully inside red vertical lines?
        ),
        expr = {
            set.seed(123)
            MetaboDecon1D(
                filepath = "urine",
                file_format = "jcampdx",
                number_iterations = 1
            )
        }
    )

    ## Check return value #####
    expect_identical(capture.output(str(x$rv)), c(
        "List of 2",
        " $ urine.dx  :List of 19",
        "  ..$ number_of_files           : int 2",
        "  ..$ filename                  : chr \"urine.dx\"",
        "  ..$ x_values                  : num [1:131072] 131 131 131 131 131 ...",
        "  ..$ x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        "  ..$ y_values                  : num [1:131072] 0.000831 0.000783 0.000743 0.000717 0.00065 ...",
        "  ..$ spectrum_superposition    : num [1, 1:131072] 3.66e-05 3.66e-05 3.66e-05 3.66e-05 3.66e-05 ...",
        "  ..$ mse_normed                : num 4.46e-11",
        "  ..$ index_peak_triplets_middle: num [1:1227] 36159 37149 37419 37435 38943 ...",
        "  ..$ index_peak_triplets_left  : num [1:1227] 36161 37160 37423 37438 38949 ...",
        "  ..$ index_peak_triplets_right : num [1:1227] 36156 37140 37415 37432 38938 ...",
        "  ..$ peak_triplets_middle      : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ...",
        "  ..$ peak_triplets_left        : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ...",
        "  ..$ peak_triplets_right       : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ...",
        "  ..$ integrals                 : num [1, 1:1227] 0.000488 0.026444 0.000409 0.000404 0.007139 ...",
        "  ..$ signal_free_region        : num [1:2] 109.1 21.9",
        "  ..$ range_water_signal_ppm    : num 0.153",
        "  ..$ A                         : num [1:1227] -0.000155 -0.00842 -0.00013 -0.000129 -0.002273 ...",
        "  ..$ lambda                    : num [1:1227] -0.00768 -0.02186 -0.00703 -0.0059 -0.0129 ...",
        "  ..$ x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ...",
        " $ urine_2.dx:List of 19",
        "  ..$ number_of_files           : int 2",
        "  ..$ filename                  : chr \"urine_2.dx\"",
        "  ..$ x_values                  : num [1:131072] 131 131 131 131 131 ...",
        "  ..$ x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        "  ..$ y_values                  : num [1:131072] 0.00586 0.00578 0.00569 0.00557 0.00548 ...",
        "  ..$ spectrum_superposition    : num [1, 1:131072] 4.69e-05 4.69e-05 4.69e-05 4.69e-05 4.69e-05 ...",
        "  ..$ mse_normed                : num 3.35e-11",
        "  ..$ index_peak_triplets_middle: num [1:1393] 36290 37241 38346 38826 39025 ...",
        "  ..$ index_peak_triplets_left  : num [1:1393] 36297 37244 38349 38835 39028 ...",
        "  ..$ index_peak_triplets_right : num [1:1393] 36285 37234 38343 38823 39019 ...",
        "  ..$ peak_triplets_middle      : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ peak_triplets_left        : num [1:1393] 9.26 9.12 8.95 8.87 8.85 ...",
        "  ..$ peak_triplets_right       : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ integrals                 : num [1, 1:1393] 0.00679 0.00499 0.00317 0.01724 0.00267 ...",
        "  ..$ signal_free_region        : num [1:2] 109.1 21.9",
        "  ..$ range_water_signal_ppm    : num 0.153",
        "  ..$ A                         : num [1:1393] -0.002161 -0.001589 -0.001011 -0.005491 -0.000851 ...",
        "  ..$ lambda                    : num [1:1393] -0.0188 -0.0168 -0.0139 -0.029 -0.0145 ...",
        "  ..$ x_0                       : num [1:1393] 94.8 93.8 92.7 92.2 92 ..."
    ))

    ## Check created files #####
    expect_file_size(x$testdir, c(
        `plots.pdf` = 321364,
        `urine/urine.dx` = 1192696,
        `urine/urine.dx approximated_spectrum.txt` = 2581870,
        `urine/urine.dx parameters.txt` = 72101,
        `urine/urine_2.dx` = 1214431,
        `urine/urine_2.dx approximated_spectrum.txt` = 2571992,
        `urine/urine_2.dx parameters.txt` = 81481
    ))
})
