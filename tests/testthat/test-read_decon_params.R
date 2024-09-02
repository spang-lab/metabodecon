library(testthat)

testdir_name <- "read_decon_params"
testdir_path <- testdir(testdir_name)
testdata_path <- pkg_file("tests/testthat/testdata")
decons <- glc_sim()

test_that("read_decon_params stops without input files", {
    clear(testdir_path)
    expect_error(x1 <- read_decon_params_v1(data_path = testdir_path))
    expect_error(x2 <- read_decon_params(data_path = testdir_path))
})

test_that("read_decon_params fails for partly missing input files", {
    clear(testdir_path)
    write_parameters_txt(decons, testdir_path)
    unlink(file.path(testdir_path, "sim_01 approximated_spectrum.txt"))
    x1 <- read_decon_params_v1(data_path = testdir_path)
    expect_error(x2 <- read_decon_params(data_path = testdir_path))
})

test_that("read_decon_params works for normal input", {
    clear(testdir_path)
    write_parameters_txt(decons, testdir_path)
    x1 <- read_decon_params_v1(data_path = testdir_path)
    x2 <- read_decon_params(data_path = testdir_path)
    # list(
    #     list(sim_01 = c(+0.958, ...), sim_02 = c(+0.938, ...)), # w
    #     list(sim_01 = c(-0.007, ...), sim_02 = c(-0.007, ...)), # lambda
    #     list(sim_01 = c(-0.044, ...), sim_02 = c(-0.272, ...)), # A
    #     list(sim_01 = c(+0.013, ...), sim_02 = c(+0.024, ...))  # spectrum_superposition
    # )           |____names only in v2____|
    x2_unnamed <- lapply(x2, function(x) { unname(x) })
    expect_true(all.equal(x1, x2_unnamed))
})
