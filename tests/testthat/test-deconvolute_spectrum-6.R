skip_on_cran()
skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

exp_str <- c(
    "List of 19",
    " $ number_of_files           : num 2",
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
    " $ signal_free_region        : num [1:2] 109.1 21.9",
    " $ range_water_signal_ppm    : num 0.153",
    " $ A                         : num [1:1227] -0.000155 -0.00842 -0.00013 -0.000129 -0.002273 ...",
    " $ lambda                    : num [1:1227] -0.00768 -0.02186 -0.00703 -0.0059 -0.0129 ...",
    " $ x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ..."
)

test_that("6. deconvolute_spectrum(format = jcampdx, sameparam = TRUE, nfile = 2/2)", {
    x <- with(
        testdir = "deconvolute_spectrum/6",
        inputs = c(urine.dx = "jcampdx/urine/urine_1.dx"),
        output = "captured", message = "captured", plots = "plots.pdf",
        answers = NULL,
        expr = {
            set.seed(1234)
            dspec <- deconvolute_spectrum(
                filepath = ".",
                name = "urine_1.dx",
                file_format = "jcampdx",
                same_parameter = TRUE,
                processing_value = 10,
                number_iterations = 1,
                range_water_signal_ppm = 0.1527692,
                signal_free_region = c(109.09458303373, 21.8529143006947), # When using the same parameters as for the last run, the signal_free_region is not interpreted as "signal free region in ppm" but as "signal free region in xxx" where "xxx" is the unit used by `deconvolute_spectrum()$signal_free_region`. From what I could see by skimming over the code, this unit is something like a "scaled rank". I.e. if signal_free_region_in_ppm[1] is 11.44 and 11.44 is the 109094th smallest value in the spectrum, then `deconvolute_spectrum()$signal_free_region[0]` will be `109094/scale_factor[1]`, which is 109.094 if scale_factor[1] is 1000.
                smoothing_param = c(2, 5),
                delta = 6.4,
                scale_factor = c(1000, 1000000),
                current_filenumber = 2,
                number_of_files = 2
            )
        }
    )
    expect_equal(capture.output(str(x$rv)), exp_str)
})
