library(testthat)

SKIP <- c()


test_that("1. MetaboDecon1D with: 1 jcampdx, answers == yy", {
    # Test 1 #####

    ## Skip conditions #####
    skip_on_cran()
    skip_if(1 %in% SKIP)
    skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

    ## Call function #####
    x <- with(
        testdir = "MetaboDecon1D/1",
        output = "captured",
        message = "captured",
        plots = "plots.pdf",
        inputs = c(urine.dx = "jcampdx/urine/urine.dx"),
        answers = c(
            "y", # Signal free region borders correct selected?
            "y" # Water artefact fully inside red vertical lines?
        ),
        expr = {
            set.seed(123)
            MetaboDecon1D(
                filepath = ".",
                filename = "urine.dx",
                file_format = "jcampdx",
                number_iterations = 1
            )
        }
    )

    ## Check return value #####
    expect_identical(capture.output(str(x$rv)), c(
        "List of 17",
        " $ number_of_files           : num 1",
        " $ filename                  : chr \"urine.dx\"",
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
        `urine.dx` = 1192696,
        `urine.dx approximated_spectrum.txt` = 2581870,
        `urine.dx parameters.txt` = 72101
    ))
})


test_that("2. MetaboDecon1D with: 2 jcampdx, answers == y1yy", {
    # Test 2 #####

    ## Skip Conditions #####
    skip_on_cran()
    skip_if(2 %in% SKIP)
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


test_that("3. MetaboDecon1D with: 2 jcampdx, answers == nyyn**y*n*y", {
    # Test 3 #####

    ## Skip conditions #####
    skip_on_cran()
    skip_if(3 %in% SKIP)
    skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

    ## Call function #####
    x <- with(
        testdir = "MetaboDecon1D/3",
        answers = c(
            "n", # Use the same parameters for all spectra?
            "y", # Signal free region borders correct selected?
            "y", # Water artefact fully inside red vertical lines?
            "n", # Signal free region borders correct selected?
            "11", # Choose another left border: [e.g. 12]
            "-1", # Choose another right border: [e.g. -2]
            "y", # Signal free region borders correct selected?
            "asdf", # Water artefact fully inside red vertical lines?
            "n", # Water artefact fully inside red vertical lines?
            "0.13", # Choose another half width range (in ppm) [e.g. 0.1222154]
            "y" # Water artefact fully inside red vertical lines?
        ),
        output = "captured",
        message = "captured",
        plots = "plots.pdf",
        inputs = c(urine = "jcampdx/urine"),
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
        " $ urine.dx  :List of 17",
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
        "  ..$ A                         : num [1:1227] -0.000155 -0.00842 -0.00013 -0.000129 -0.002273 ...",
        "  ..$ lambda                    : num [1:1227] -0.00768 -0.02186 -0.00703 -0.0059 -0.0129 ...",
        "  ..$ x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ...",
        " $ urine_2.dx:List of 17",
        "  ..$ number_of_files           : int 2",
        "  ..$ filename                  : chr \"urine_2.dx\"",
        "  ..$ x_values                  : num [1:131072] 131 131 131 131 131 ...",
        "  ..$ x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        "  ..$ y_values                  : num [1:131072] 0.00586 0.00578 0.00569 0.00557 0.00548 ...",
        "  ..$ spectrum_superposition    : num [1, 1:131072] 4.7e-05 4.7e-05 4.7e-05 4.7e-05 4.7e-05 ...",
        "  ..$ mse_normed                : num 3.33e-11",
        "  ..$ index_peak_triplets_middle: num [1:1402] 36290 37241 38346 38826 39025 ...",
        "  ..$ index_peak_triplets_left  : num [1:1402] 36297 37244 38349 38835 39028 ...",
        "  ..$ index_peak_triplets_right : num [1:1402] 36285 37234 38343 38823 39019 ...",
        "  ..$ peak_triplets_middle      : num [1:1402] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ peak_triplets_left        : num [1:1402] 9.26 9.12 8.95 8.87 8.85 ...",
        "  ..$ peak_triplets_right       : num [1:1402] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ integrals                 : num [1, 1:1402] 0.00679 0.00499 0.00317 0.01724 0.00267 ...",
        "  ..$ A                         : num [1:1402] -0.002161 -0.001589 -0.001011 -0.005491 -0.000851 ...",
        "  ..$ lambda                    : num [1:1402] -0.0188 -0.0168 -0.0139 -0.029 -0.0145 ...",
        "  ..$ x_0                       : num [1:1402] 94.8 93.8 92.7 92.2 92 ..."
    ))

    ## Check created files #####
    expect_file_size(x$testdir, c(
        `plots.pdf` = 961666,
        `urine/urine.dx approximated_spectrum.txt` = 2581870,
        `urine/urine.dx parameters.txt` = 72101,
        `urine/urine_2.dx approximated_spectrum.txt` = 2571332,
        `urine/urine_2.dx parameters.txt` = 82012
    ))

    ## Check output #####
    expect_identical(x$output$text, c("[1] 4.456557e-11", "[1] 3.332342e-11"))
    expect_identical(x$message$text, c(
        "",
        "    <MetaboDecon1D>  Copyright (C) <2021>  <Martina Haeckl>",
        "    This program comes with ABSOLUTELY NO WARRANTY.",
        "    This is free software, and you are welcome to redistribute it",
        "    under certain conditions; type `show_license()' for details.",
        "Do you want to use the same parameters (signal_free_region, range_water_signal_ppm) for all spectra? (y/n) n",
        "Start deconvolution of urine.dx:",
        "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): y",
        "Water artefact fully inside red vertical lines? (y/n): y",
        "",
        "Normed MSE value of iteration 1 is: ",
        "",
        "Saving parameters to txt documents...",
        "Start deconvolution of urine_2.dx:",
        "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): n",
        "Choose another left border: [e.g. 12] 11",
        "Choose another right border: [e.g. -2] -1",
        "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): y",
        "Water artefact fully inside red vertical lines? (y/n): asdf",
        "Error. Please type only y or n.",
        "Water artefact fully inside red vertical lines? (y/n): n",
        "Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154] 0.13",
        "Water artefact fully inside red vertical lines? (y/n): y",
        "",
        "Normed MSE value of iteration 1 is: ",
        "",
        "Saving parameters to txt documents..."
    ))
})

test_that("4. MetaboDecon1D with: 2 bruker, answers == nyyn**y*n*y", {
    # Test 4 #####

    ## Skip conditions #####
    skip_if(4 %in% SKIP)
    skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

    ## Call function #####
    x <- with(
        testdir = "MetaboDecon1D/4",
        answers = c(
            "10", # What is the name of the subfolder of your filepath? [e.g. 10]
            "10", # What is the name of the subsubsubfolder of your filepath? [e.g. 10]
            "n", # Use the same parameters for all spectra?
            "y", # Signal free region borders correct selected?
            "y", # Water artefact fully inside red vertical lines?
            "n", # Signal free region borders correct selected?
            "11", # Choose another left border: [e.g. 12]
            "-1", # Choose another right border: [e.g. -2]
            "y", # Signal free region borders correct selected?
            "asdf", # Water artefact fully inside red vertical lines?
            "n", # Water artefact fully inside red vertical lines?
            "0.13", # Choose another half width range (in ppm) [e.g. 0.1222154]
            "y" # Water artefact fully inside red vertical lines?
        ),
        output = "captured",
        message = "captured",
        plots = "plots.pdf",
        inputs = c(urine = "bruker/urine"),
        expr = {
            set.seed(123)
            MetaboDecon1D(
                filepath = "urine",
                file_format = "bruker",
                number_iterations = 1
            )
        }
    )

    ## Check return value #####
    expect_identical(capture.output(str(x$rv)), c(
        "List of 2",
        " $ urine  :List of 17",
        "  ..$ number_of_files           : int 2",
        "  ..$ filename                  : chr \"urine\"",
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
        "  ..$ A                         : num [1:1227] -0.000155 -0.00842 -0.00013 -0.000129 -0.002273 ...",
        "  ..$ lambda                    : num [1:1227] -0.00768 -0.02186 -0.00703 -0.0059 -0.0129 ...",
        "  ..$ x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ...",
        " $ urine_2:List of 17",
        "  ..$ number_of_files           : int 2",
        "  ..$ filename                  : chr \"urine_2\"",
        "  ..$ x_values                  : num [1:131072] 131 131 131 131 131 ...",
        "  ..$ x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        "  ..$ y_values                  : num [1:131072] 0.00586 0.00578 0.00569 0.00557 0.00548 ...",
        "  ..$ spectrum_superposition    : num [1, 1:131072] 4.7e-05 4.7e-05 4.7e-05 4.7e-05 4.7e-05 ...",
        "  ..$ mse_normed                : num 3.33e-11",
        "  ..$ index_peak_triplets_middle: num [1:1402] 36290 37241 38346 38826 39025 ...",
        "  ..$ index_peak_triplets_left  : num [1:1402] 36297 37244 38349 38835 39028 ...",
        "  ..$ index_peak_triplets_right : num [1:1402] 36285 37234 38343 38823 39019 ...",
        "  ..$ peak_triplets_middle      : num [1:1402] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ peak_triplets_left        : num [1:1402] 9.26 9.12 8.95 8.87 8.85 ...",
        "  ..$ peak_triplets_right       : num [1:1402] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ integrals                 : num [1, 1:1402] 0.00679 0.00499 0.00317 0.01724 0.00267 ...",
        "  ..$ A                         : num [1:1402] -0.002161 -0.001589 -0.001011 -0.005491 -0.000851 ...",
        "  ..$ lambda                    : num [1:1402] -0.0188 -0.0168 -0.0139 -0.029 -0.0145 ...",
        "  ..$ x_0                       : num [1:1402] 94.8 93.8 92.7 92.2 92 ..."
    ))

    ## Check created files #####
    expect_file_size(x$testdir, c(
        `plots.pdf` = 961664,
        `urine/urine approximated_spectrum.txt` = 2581865,
        `urine/urine parameters.txt` = 72104,
        `urine/urine_2 approximated_spectrum.txt` = 2571332,
        `urine/urine_2 parameters.txt` = 82012
    ))

    ## Check output #####
    expect_identical(x$output$text, c("[1] 4.456526e-11", "[1] 3.332342e-11"))
    expect_identical(x$message$text, c(
        "",
        "    <MetaboDecon1D>  Copyright (C) <2021>  <Martina Haeckl>",
        "    This program comes with ABSOLUTELY NO WARRANTY.",
        "    This is free software, and you are welcome to redistribute it",
        "    under certain conditions; type `show_license()' for details.",
        "What is the name of the subfolder of your filepath: ",
        "[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10] 10",
        "What is the name of the subsubsubfolder of your filepath: ",
        "[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10/pdata/10] 10",
        "Do you want to use the same parameters (signal_free_region, range_water_signal_ppm) for all spectra? (y/n) n",
        "Start deconvolution of urine:",
        "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): y",
        "Water artefact fully inside red vertical lines? (y/n): y",
        "",
        "Normed MSE value of iteration 1 is: ",
        "",
        "Saving parameters to txt documents...",
        "Start deconvolution of urine_2:",
        "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): n",
        "Choose another left border: [e.g. 12] 11",
        "Choose another right border: [e.g. -2] -1",
        "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): y",
        "Water artefact fully inside red vertical lines? (y/n): asdf",
        "Error. Please type only y or n.",
        "Water artefact fully inside red vertical lines? (y/n): n",
        "Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154] 0.13",
        "Water artefact fully inside red vertical lines? (y/n): y",
        "",
        "Normed MSE value of iteration 1 is: ",
        "",
        "Saving parameters to txt documents..."
    ))
})
