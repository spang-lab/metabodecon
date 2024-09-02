library(testthat)

testdir_name <- "get_decon_params"
testdir_path <- testdir(testdir_name)
testdata_path <- pkg_file("tests/testthat/testdata")

expected_structure <- c(
    "List of 4",
    " $ w                     :List of 2",
    "  ..$ sim_01: num [1:79] 1.3 1.3 1.28 1.28 1.25 ...",
    "  ..$ sim_02: num [1:41] 1.3 1.25 1.21 1.16 1.15 ...",
    " $ lambda                :List of 2",
    "  ..$ sim_01: num [1:79] -0.00154 -0.00289 -0.00432 -0.00359 -0.00459 ...",
    "  ..$ sim_02: num [1:41] -0.00497 -0.00401 -0.00362 -0.0051 -0.00524 ...",
    " $ A                     :List of 2",
    "  ..$ sim_01: num [1:79] -2.24e-07 -8.75e-07 -7.96e-07 -7.78e-07 -1.60e-06 ...",
    "  ..$ sim_02: num [1:41] -5.00e-06 -6.51e-06 -2.76e-06 -2.51e-06 -4.61e-06 ...",
    " $ spectrum_superposition:List of 2",
    "  ..$ sim_01: num [1:1309] 0.0136 0.0137 0.0138 0.0138 0.0139 ...",
    "  ..$ sim_02: num [1:1309] 0.0114 0.0115 0.0116 0.0117 0.0119 ..."
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
