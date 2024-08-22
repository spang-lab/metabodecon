library(testthat)

testdir_name <- "read_decon_params"
testdir_path <- testdir(testdir_name)
sim_subset_path <- file.path(testdir_path, "sim_subset")

#' @description Creates files `sim_01 approximated_spectrum.txt` and `sim_01 parameters.txt` inside folder `testdir_name` by calling the old `MetaboDecon1D()` function (which writes these files to the input folder). These files can be used as inputs to test the `read_decon_params()` function.
make_test_inputs <- function() {
    ewo <- evalwith(
        testdir = testdir_name,
        output = "captured", message = "captured",
        plot = "plots.pdf",
        inputs = "bruker/sim_subset",
        answers = c(
            10,   # Name of subfolder?
            10,   # Name of subsubsubfolder?
            "y",  # Use same parameters for all spectra?
            1,    # Number of spectra to process?
            "n",  # Signal free region borders correctly selected?
            3.58, # Left border
            3.42, # Right border
            "y",  # Signal free region borders correctly selected?
            "n",  # Water artefact inside red vertical lines?
            0,    # Choose half width range for the water artefact.
            "y",  # Water artefact inside red vertical lines?
            "y"   # Save results as text documents at `data_path`?
        ),
        expr = {
            MetaboDecon1D("sim_subset")
        }
    )
    stopifnot(all.equal(dir(sim_subset_path), c(
        "sim_01", "sim_01 approximated_spectrum.txt", "sim_01 parameters.txt",
        "sim_02", "sim_02 approximated_spectrum.txt", "sim_02 parameters.txt"
    )))
}


test_that("read_decon_params stops without input files", {
    clear(testdir_path)
    expect_error(x1 <- read_decon_params_v1(data_path = sim_subset_path))
    expect_error(x2 <- read_decon_params_v2(data_path = sim_subset_path))
})

test_that("read_decon_params fails for partly missing input files", {
    clear(testdir_path)
    make_test_inputs() # Creates the following files in `sim_subset_path`:
    # sim_01 approximated_spectrum.txt, sim_01 parameters.txt,
    # sim_02 approximated_spectrum.txt, sim_02 parameters.txt,
    unlink(file.path(sim_subset_path, "sim_01 approximated_spectrum.txt"))
    x1 <- read_decon_params_v1(data_path = sim_subset_path)
    expect_error(x2 <- read_decon_params_v2(data_path = sim_subset_path))
})

test_that("read_decon_params works for normal input", {
    clear(testdir_path)
    make_test_inputs() # Creates the following files in `sim_subset_path`:
    # sim_01 approximated_spectrum.txt, sim_01 parameters.txt,
    # sim_02 approximated_spectrum.txt, sim_02 parameters.txt,
    x1 <- read_decon_params_v1(data_path = sim_subset_path)
    x2 <- read_decon_params_v2(data_path = sim_subset_path, verbose = TRUE)
    # list(
    #     w =               list(sim_01 = c(+0.958, ...), sim_02 = c(+0.938, ...)),
    #     lambda =          list(sim_01 = c(-0.007, ...), sim_02 = c(-0.007, ...)),
    #     A =               list(sim_01 = c(-0.044, ...), sim_02 = c(-0.272, ...)),
    #     noise_threshold = list(sim_01 =   +0.000,       sim_02 =   +0.000      ),
    #     spectrum_approx = list(sim_01 = c(+0.013, ...), sim_02 = c(+0.024, ...))
    # )                             |____names only in v2____|
    x2_unnamed <- lapply(x2, function(x) { unname(x) })
    expect_true(all.equal(x1, x2_unnamed))
})
