# Keep this test in a seperate file, because it's very slow. This is useful
# because testthat runs test files in parallel, i.e. by distributing slow tests
# across different files we ensure they can be also distributed across different
# cores.
library(testthat)

test_that("4. generate_lorentz_curves with: 2 bruker, answers == nyyn**y*n*y", {
    ## Skip conditions #####
    skip_on_cran()
    skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

    ## Call function #####
    x <- evalwith(
        testdir = "generate_lorentz_curves/4",
        output = "captured",
        message = "captured",
        plots = "plots.pdf",
        inputs = c(urine = "bruker/urine"),
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
        expr = {
            set.seed(123)
            generate_lorentz_curves(data_path = "urine", file_format = "bruker")
        }
    )

    ## Check return value #####
    expect_identical(capture.output(str(x$rv)), c(
        "List of 2",
        " $ urine_1:List of 17",
        "  ..$ number_of_files           : int 2",
        "  ..$ filename                  : chr \"urine_1\"",
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
        "  ..$ A                         : num [1:1227] -0.00016 -0.008436 -0.000128 -0.000119 -0.002634 ...",
        "  ..$ lambda                    : num [1:1227] -0.00775 -0.02188 -0.00675 -0.00562 -0.01343 ...",
        "  ..$ x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ...",
        " $ urine_2:List of 17",
        "  ..$ number_of_files           : int 2",
        "  ..$ filename                  : chr \"urine_2\"",
        "  ..$ x_values                  : num [1:131072] 131 131 131 131 131 ...",
        "  ..$ x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        "  ..$ y_values                  : num [1:131072] 0.00586 0.00578 0.00569 0.00557 0.00548 ...",
        "  ..$ spectrum_superposition    : num [1, 1:131072] 4.23e-05 4.23e-05 4.23e-05 4.23e-05 4.23e-05 ...",
        "  ..$ mse_normed                : num 2.87e-11",
        "  ..$ index_peak_triplets_middle: num [1:1402] 36290 37241 38346 38826 39025 ...",
        "  ..$ index_peak_triplets_left  : num [1:1402] 36297 37244 38349 38835 39028 ...",
        "  ..$ index_peak_triplets_right : num [1:1402] 36285 37234 38343 38823 39019 ...",
        "  ..$ peak_triplets_middle      : num [1:1402] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ peak_triplets_left        : num [1:1402] 9.26 9.12 8.95 8.87 8.85 ...",
        "  ..$ peak_triplets_right       : num [1:1402] 9.26 9.12 8.95 8.88 8.85 ...",
        "  ..$ integrals                 : num [1, 1:1402] 0.00683 0.00504 0.00322 0.0174 0.00274 ...",
        "  ..$ A                         : num [1:1402] -0.002176 -0.001604 -0.001025 -0.005541 -0.000872 ...",
        "  ..$ lambda                    : num [1:1402] -0.0189 -0.0168 -0.014 -0.0291 -0.0146 ...",
        "  ..$ x_0                       : num [1:1402] 94.8 93.8 92.7 92.2 92 ..."
    ))
    ## Check created files #####
    expect_file_size(x$testdir, c(
        `plots.pdf` = 961664,
        `urine/urine_1 approximated_spectrum.txt` = 2581865,
        `urine/urine_1 parameters.txt` = 72104,
        `urine/urine_2 approximated_spectrum.txt` = 2571332,
        `urine/urine_2 parameters.txt` = 82012
    ))

    ## Check output #####
    expect_identical(x$output$text, c(
        "[1] 4.456526e-11",
        "[1] 4.297345e-11",
        "[1] 4.104424e-11",
        "[1] 4.039879e-11",
        "[1] 3.936393e-11",
        "[1] 3.935647e-11",
        "[1] 3.893839e-11",
        "[1] 3.936722e-11",
        "[1] 3.903905e-11",
        "[1] 3.915808e-11",
        "[1] 3.332342e-11",
        "[1] 3.02236e-11",
        "[1] 2.972469e-11",
        "[1] 2.941173e-11",
        "[1] 2.909153e-11",
        "[1] 2.897956e-11",
        "[1] 2.866265e-11",
        "[1] 2.866066e-11",
        "[1] 2.859706e-11",
        "[1] 2.870111e-11"
    ))
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
        "Start deconvolution of urine_1:",
        "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): y",
        "Water artefact fully inside red vertical lines? (y/n): y",
        "",
        "Normed MSE value of iteration 1 is: ",
        "",
        "Normed MSE value of iteration 2 is: ",
        "",
        "Normed MSE value of iteration 3 is: ",
        "",
        "Normed MSE value of iteration 4 is: ",
        "",
        "Normed MSE value of iteration 5 is: ",
        "",
        "Normed MSE value of iteration 6 is: ",
        "",
        "Normed MSE value of iteration 7 is: ",
        "",
        "Normed MSE value of iteration 8 is: ",
        "",
        "Normed MSE value of iteration 9 is: ",
        "",
        "Normed MSE value of iteration 10 is: ",
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
        "Normed MSE value of iteration 2 is: ",
        "",
        "Normed MSE value of iteration 3 is: ",
        "",
        "Normed MSE value of iteration 4 is: ",
        "",
        "Normed MSE value of iteration 5 is: ",
        "",
        "Normed MSE value of iteration 6 is: ",
        "",
        "Normed MSE value of iteration 7 is: ",
        "",
        "Normed MSE value of iteration 8 is: ",
        "",
        "Normed MSE value of iteration 9 is: ",
        "",
        "Normed MSE value of iteration 10 is: ",
        "",
        "Saving parameters to txt documents..."
    ))
})
