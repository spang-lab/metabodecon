skip_if_slow_tests_disabled()

test_that("GLC works with: dp=urine_1, an=yy, ni=3", {
    x <- glc_urine1_yy_ni3_dbg()$rv$urine_1
    y <- MetaboDecon1D_urine1_1010yy_ni3_dbg()$rv
    r <- compare_spectra(x, y, silent = TRUE)
    expect_true(sum(r == 0) >= 51 &&sum(r == 1) >= 9 &&sum(r %in% c(2, 3) == 0)) # >=51 identical, >=9 equal, no errors
})

skip()

test_that("generate_lorentz_curves works with 1 jcampdx, answers == y1yy", {
    x <- evalwith(
        testdir = "glc_urine1dx_1010yy_ni3_dbg", inputs = "jcampdx/urine_1.dx",
        answers = c("y", "1", "y", "y"),
        output = "captured", message = "captured", plot = "plots.pdf",
        expr = generate_lorentz_curves(data_path = ".", file_format = "jcampdx")
    )
    expect_str(x$rv, str_urine_1_deconvoluted())
    expect_file_size(x$testdir, c(`plots.pdf` = 321364, `urine_1.dx` = 1192696, `urine_1.dx approximated_spectrum.txt` = 2581870, `urine_1.dx parameters.txt` = 72101))
})

test_that("generate_lorentz_curves works with 2 jcampdx, answers == y1yy", {
    x <- evalwith(
        testdir = "generate_lorentz_curves/2", inputs = "jcampdx/urine", answers = c("y", "1", "y", "y"),
        output = "captured", message = "captured", plot = "plots.pdf",
        expr = generate_lorentz_curves(data_path = "urine", file_format = "jcampdx")
    )
    expect_str(x$rv, str_urine_1_deconvoluted())
    expect_file_size(x$testdir, c(`plots.pdf` = 321364, `urine_1.dx` = 1192696, `urine_1.dx approximated_spectrum.txt` = 2581870, `urine_1.dx parameters.txt` = 72101))
})

test_that("generate_lorentz_curves works with 2 jcampdx, answers == nyyn**y*n*y", {
    x <- evalwith(
        testdir = "generate_lorentz_curves/3", inputs = "jcampdx/urine", answers = "n",
        output = "captured", message = "captured", plot = "plots.pdf",
        expr = generate_lorentz_curves(data_path = "urine", file_format = "jcampdx")
    )
    expect_str(x$rv, str_urine_deconvoluted())
    expect_file_size(x$testdir, c(`plots.pdf` = 961666, `urine_2.dx approximated_spectrum.txt` = 2571332, `urine_2.dx parameters.txt` = 82012))
})

test_that("generate_lorentz_curves works with 2 bruker, answers == nyyn**y*n*y", {
    x4 <- evalwith(
        testdir = "generate_lorentz_curves/4", inputs = "bruker/urine",
        answers = c("10", "10", "n", "y", "y", "n", "11", "-1", "y", "asdf", "n", "0.13", "y"),
        output = "captured", message = "captured", plot = "plots.pdf",
        expr = generate_lorentz_curves(data_path = "urine", file_format = "bruker")
    )
    expect_str(x4$rv, str_urine_deconvoluted())
    expect_file_size(x$testdir, c(`plots.pdf` = 961664, `urine_1 approximated_spectrum.txt` = 2581865, `urine_1 parameters.txt` = 72104, `urine_2 approximated_spectrum.txt` = 2571332, `urine_2 parameters.txt` = 82012))
})
