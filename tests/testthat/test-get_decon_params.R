library(testthat)

testdir_name <- "get_decon_params"
testdir_path <- testdir(testdir_name)
testdata_path <- pkg_file("tests/testthat/testdata")

expected_structure <- c(
    "List of 4",
    " $ w                     :List of 2",
    "  ..$ sim_01: num [1:79] 1.3 1.3 1.28 1.28 1.25 ...",
    "  ..$ sim_02: num [1:79] 1.3 1.3 1.28 1.28 1.28 ...",
    " $ lambda                :List of 2",
    "  ..$ sim_01: num [1:79] -0.00154 -0.00289 -0.00432 -0.00359 -0.00459 ...",
    "  ..$ sim_02: num [1:79] -0.0028 -0.00422 -0.00297 -0.00575 -0.00563 ...",
    " $ A                     :List of 2",
    "  ..$ sim_01: num [1:79] -2.24e-07 -8.75e-07 -7.96e-07 -7.78e-07 -1.60e-06 ...",
    "  ..$ sim_02: num [1:79] -6.91e-07 -1.18e-06 -6.08e-07 -1.25e-06 -1.56e-06 ...",
    " $ spectrum_superposition:List of 2",
    "  ..$ sim_01: num [1:1309] 0.0136 0.0137 0.0138 0.0138 0.0139 ...",
    "  ..$ sim_02: num [1:1309] 0.0246 0.0247 0.0249 0.0251 0.0252 ..."
)

test_that("get_decon_params works", {
    clear(testdir_path)
    sim <- metabodecon_file("bruker/sim_subset")
    decons <- generate_lorentz_curves(sim, sfr = c(3.58, 3.42), wshw = 0, delta = 0.1, ask = FALSE, verbose = FALSE)
    obj <- get_decon_params(decons)
    expect_str(obj, expected_structure)
    expect_error(get_decon_params("asdf", warn = FALSE), "asdf does not exist.")
    expect_error(get_decon_params(testdir_path, warn = FALSE), "No parameter files found in the given directory.")
    write_parameters_txt(decons, testdir_path)
    expect_warning(obj2 <- get_decon_params(testdir_path, warn = TRUE), "You have provided a path.*")
    expect_str(obj2, expected_structure)
    expect_equal(obj[[1]], obj2[[1]], tolerance = 1e-6)
})
