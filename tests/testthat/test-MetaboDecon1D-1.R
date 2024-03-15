# Keep this test in a seperate file, because it's very slow. This is useful
# because testthat runs test files in parallel, i.e. by distributing slow tests
# across different files we ensure they can be also distributed across different
# cores.
library(testthat)

test_that("1. MetaboDecon1D with: 1 jcampdx, answers == yy", {

    ## Skip conditions #####
    skip_on_cran()
    skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

    ## Call function #####
    x <- evalwith(
        testdir = "MetaboDecon1D/1",
        output = "captured",
        message = "captured",
        plots = "plots.pdf",
        inputs = c(urine_1.dx = "jcampdx/urine/urine_1.dx"),
        answers = c(
            "y", # Signal free region borders correct selected?
            "y" # Water artefact fully inside red vertical lines?
        ),
        expr = {
            set.seed(123)
            MetaboDecon1D(
                filepath = ".",
                filename = "urine_1.dx",
                file_format = "jcampdx",
                number_iterations = 1
            )
        }
    )

    ## Check return value #####
    expect_identical(capture.output(str(x$rv)), c(
        "List of 17",
        " $ number_of_files           : num 1",
        " $ filename                  : chr \"urine_1.dx\"",
        " $ x_values                  : num [1:131072] 131 131 131 131 131 ...",
        " $ x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        " $ y_values                  : num [1:131072] 0.000831 0.000783 0.000743 0.000717 0.00065 ...",
        " $ spectrum_superposition    : num [1, 1:131072] 3.66e-05 3.66e-05 3.66e-05 3.66e-05 3.66e-05 ...",
        " $ mse_normed                : num 4.46e-11",
        " $ index_peak_triplets_middle: num [1:1227] 36159 37149 37419 37435 38943 ...",
        " $ index_peak_triplets_left  : num [1:1227] 36161 37160 37423 37438 38949 ...",
        " $ index_peak_triplets_right : num [1:1227] 36156 37140 37415 37432 38938 ...",
        " $ peak_triplets_middle      : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ...",
        " $ peak_triplets_left        : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ...",
        " $ peak_triplets_right       : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ...",
        " $ integrals                 : num [1, 1:1227] 0.000488 0.026444 0.000409 0.000404 0.007139 ...",
        " $ A                         : num [1:1227] -0.000155 -0.00842 -0.00013 -0.000129 -0.002273 ...",
        " $ lambda                    : num [1:1227] -0.00768 -0.02186 -0.00703 -0.0059 -0.0129 ...",
        " $ x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ..."
    ))

    ## Check created files #####
    expect_file_size(x$testdir, c(
        `plots.pdf` = 321364,
        `urine_1.dx` = 1192696,
        `urine_1.dx approximated_spectrum.txt` = 2581870,
        `urine_1.dx parameters.txt` = 72101
    ))
})
