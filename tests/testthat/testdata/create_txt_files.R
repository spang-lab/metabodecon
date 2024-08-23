
Usage <- quote({
    source("tests/testthat/testdata/create_txt_files.R")
})

devtools::load_all(quiet = TRUE)
testdata_dir <- pkg_file("tests/testthat/testdata")
testdir_name <- "test_data_create_txt_files"
testdir_path <- testdir(testdir_name)
sim_subset_path <- file.path(testdir_path, "sim_subset")
cat(sprintf("Clearing test directory: %s\n", testdir_path))
clear(testdir_path)
cat("Creating parameter files\n")
ewo <- evalwith(
    testdir = testdir_name,
    output = "captured", message = "captured", plot = "plots.pdf",
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
created_files <- dir(sim_subset_path)
stopifnot(
    "sim_01 approximated_spectrum.txt" %in% created_files &&
    "sim_02 approximated_spectrum.txt" %in% created_files &&
    "sim_01 parameters.txt" %in% created_files &&
    "sim_02 parameters.txt" %in% created_files
)
fetch <- function(x) {
    src <- file.path(sim_subset_path, x)
    dst <- file.path(testdata_dir)
    msg <- sprintf("Copying:\nFrom: %s\nTo:   %s\n", src, dst)
    cat(msg)
    file.copy(src, dst)
}
fetch("sim_01 approximated_spectrum.txt")
fetch("sim_02 approximated_spectrum.txt")
fetch("sim_01 parameters.txt")
fetch("sim_02 parameters.txt")
