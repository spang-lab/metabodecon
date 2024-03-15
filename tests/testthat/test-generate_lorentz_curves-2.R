# Keep this test in a seperate file, because it's very slow. This is useful
# because testthat runs test files in parallel, i.e. by distributing slow tests
# across different files we ensure they can be also distributed across different
# cores.
library(testthat)

test_that("2. generate_lorentz_curves with: 2 jcampdx, answers == y1yy", {

    ## Skip conditions #####
    skip_on_cran()
    skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

    ## Call function #####
    x <- evalwith(
        testdir = "generate_lorentz_curves/2",
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
            generate_lorentz_curves(data_path = "urine", file_format = "jcampdx")
        }
    )

    ## Check return value #####
    expect_identical(capture.output(str(x$rv)), c(
        "List of 2",
        " $ urine_1.dx:List of 19",
        "  ..$ number_of_files           : int 2",
        "  ..$ filename                  : chr \"urine_1.dx\"",
        "  ..$ x_values                  : num [1:131072] 131 131 131 131 131 ...",
        "  ..$ x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        "  ..$ y_values                  : num [1:131072] 0.000831 0.000783 0.000743 0.000717 0.00065 ...",
        "  ..$ spectrum_superposition    : num [1, 1:131072] 3.51e-05 3.51e-05 3.51e-05 3.51e-05 3.52e-05 ...",
        "  ..$ mse_normed                : num 3.92e-11",
        "  ..$ index_peak_triplets_middle: num [1:1227] 36159 37149 37419 37435 38943 ...",
        "  ..$ index_peak_triplets_left  : num [1:1227] 36161 37160 37423 37438 38949 ...",
        "  ..$ index_peak_triplets_right : num [1:1227] 36156 37140 37415 37432 38938 ...",
        "  ..$ peak_triplets_middle      : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ...",
        "  ..$ peak_triplets_left        : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ...",
        "  ..$ peak_triplets_right       : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ...",
        "  ..$ integrals                 : num [1, 1:1227] 0.000501 0.026496 0.000402 0.000375 0.008274 ...",
        "  ..$ signal_free_region        : num [1:2] 109.1 21.9",
        "  ..$ range_water_signal_ppm    : num 0.153",
        "  ..$ A                         : num [1:1227] -0.00016 -0.008436 -0.000128 -0.000119 -0.002634 ...",
        "  ..$ lambda                    : num [1:1227] -0.00775 -0.02188 -0.00675 -0.00562 -0.01343 ...",
        "  ..$ x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ...",
        " $ urine_2.dx:List of 19",
        "  ..$ number_of_files           : int 2",
        "  ..$ filename                  : chr \"urine_2.dx\"",
        "  ..$ x_values                  : num [1:131072] 131 131 131 131 131 ...",
        "  ..$ x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        "  ..$ y_values                  : num [1:131072] 0.00586 0.00578 0.00569 0.00557 0.00548 ...",
        "  ..$ spectrum_superposition    : num [1, 1:131072] 4.21e-05 4.21e-05 4.21e-05 4.21e-05 4.21e-05 ...",
        "  ..$ mse_normed                : num 2.86e-11",
        "  ..$ index_peak_triplets_middle: num [1:1393] 36290 37241 38346 38826 39025 ...",
        "  ..$ index_peak_triplets_left  : num [1:1393] 36297 37244 38349 38835 39028 ...",
        "  ..$ index_peak_triplets_right : num [1:1393] 36285 37234 38343 38823 39019 ...",
        "  ..$ peak_triplets_middle      : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ peak_triplets_left        : num [1:1393] 9.26 9.12 8.95 8.87 8.85 ...",
        "  ..$ peak_triplets_right       : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ integrals                 : num [1, 1:1393] 0.00683 0.00504 0.00322 0.0174 0.00274 ...",
        "  ..$ signal_free_region        : num [1:2] 109.1 21.9",
        "  ..$ range_water_signal_ppm    : num 0.153",
        "  ..$ A                         : num [1:1393] -0.002176 -0.001604 -0.001025 -0.005541 -0.000872 ...",
        "  ..$ lambda                    : num [1:1393] -0.0189 -0.0168 -0.014 -0.0291 -0.0146 ...",
        "  ..$ x_0                       : num [1:1393] 94.8 93.8 92.7 92.2 92 ..."
    ))

    ## Check created files #####
    expect_file_size(x$testdir, c(
        `plots.pdf` = 321364,
        `urine/urine_1.dx` = 1192696,
        `urine/urine_1.dx approximated_spectrum.txt` = 2581870,
        `urine/urine_1.dx parameters.txt` = 72101,
        `urine/urine_2.dx` = 1214431,
        `urine/urine_2.dx approximated_spectrum.txt` = 2571992,
        `urine/urine_2.dx parameters.txt` = 81481
    ))
})
