library(testthat)

from_decon1 <- test_that("get_decon_params works for decon1 objects", {
    x <- generate_lorentz_curves(sim[1:2], sfr = c(3.55, 3.35), wshw = 0, ask = FALSE, verbose = FALSE)
    p <- get_decon_params(x)
    expect_equal(names(p), expected = c("w", "lambda", "A", "spectrum_superposition"))
    expect_equal(names(p$w), expected = c("sim_01", "sim_02"))
    expect_equal(names(p$A), expected = c("sim_01", "sim_02"))
    expect_true(is.numeric(p$lambda$sim_01))
    expect_true(is.numeric(p$spectrum_superposition$sim_02))
    expect_equal(p$w$sim_01, x$sim_01$x_0)
})

correct_error <- test_that("get_decon_params raises correct errors", {
    testdir_path <- testdir("get_decon_params")
    clear(testdir_path)
    expect_error(
        object = get_decon_params("asdf"),
        regexp = "asdf does not exist."
    )
    expect_error(
        object = get_decon_params(testdir_path, warn = FALSE),
        regexp = "No parameter files found in the given directory."
    )
})

from_files <- test_that("get_decon_params works for *.txt files", {
    testdir_path <- testdir("get_decon_params")
    clear(testdir_path)
    x <- generate_lorentz_curves(sim[1:2], sfr = c(3.55, 3.35), wshw = 0, ask = FALSE, verbose = FALSE)
    write_parameters_txt(x, testdir_path)
    expect_warning(p <- get_decon_params(testdir_path, warn = TRUE), "You have provided a path.*")
    expect_equal(names(p), expected = c("w", "lambda", "A", "spectrum_superposition"))
    expect_equal(names(p$w), expected = c("sim_01", "sim_02"))
    expect_equal(names(p$A), expected = c("sim_01", "sim_02"))
    expect_true(is.numeric(p$lambda$sim_01))
    expect_true(is.numeric(p$spectrum_superposition$sim_02))
    expect_equal(p$w$sim_01, x$sim_01$x_0)
})
