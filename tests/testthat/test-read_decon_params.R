library(testthat)

testdir_name <- "read_decon_params"
testdir_path <- testdir(testdir_name)
testdata_path <- pkg_file("tests/testthat/testdata")
path <- metabodecon_file("bruker/sim_subset")
decons <- generate_lorentz_curves_sim(path)

test_that("read_decon_params stops without input files", {
    clear(testdir_path)
    expect_error(x1 <- read_decon_params_original(data_path = testdir_path))
    expect_error(x2 <- read_decon_params(data_path = testdir_path))
})

test_that("read_decon_params fails for partly missing input files", {
    clear(testdir_path)
    write_parameters_txt(decons, testdir_path)
    unlink(file.path(testdir_path, "sim_01 approximated_spectrum.txt"))
    x1 <- read_decon_params_original(data_path = testdir_path)
    expect_error(x2 <- read_decon_params(data_path = testdir_path))
})

test_that("read_decon_params works for normal input", {
    clear(testdir_path)
    write_parameters_txt(decons, testdir_path)
    x1 <- read_decon_params_original(data_path = testdir_path)
    x2 <- read_decon_params(data_path = testdir_path)
    # x1 = list(
    #   list(sim_01 = c( 0.958, ...), sim_02 = c(...)),
    #   list(sim_01 = c(-0.007, ...), sim_02 = c(...)),
    #   list(sim_01 = c(-0.044, ...), sim_02 = c(...)),
    #   list(sim_01 = c( 0.013, ...), sim_02 = c(...))
    # )
    # names(x1) <- c("w", "lambda", "A", "spectrum_superposition")
    x2_unnamed <- lapply(x2, function(x) { unname(x) }) # (1)
    # (1) Names weren't present in the original function
    expect_true(all.equal(x1, x2_unnamed))
})
