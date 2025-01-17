library(testthat)

test_that("read_spectrum works for bruker spectra", {
    urine <- pkg_file("example_datasets/bruker/urine")
    urine_1 <- file.path(urine, "urine_1")
    urine_2 <- file.path(urine, "urine_2")
    X1 <- read_spectrum(urine_1)
    X2 <- read_spectrum(urine_2)
    XX1 <- read_spectra(urine_1)
    XX <- read_spectra(urine)
    expect_equal(X1, XX$urine_1)
    expect_equal(X2, XX$urine_2)
    expect_equal(XX1$urine_1, XX$urine_1)
    expect_equal(names(X1), c("si", "cs", "meta"))
})

skip_if_slow_tests_disabled()

test_that("read_spectrum works for jcampdx spectra", {
    urine_1 <- pkg_file("example_datasets/bruker/urine/urine_1")
    urine_1_dx <- pkg_file("example_datasets/jcampdx/urine/urine_1.dx")
    X1 <- read_spectrum(urine_1, raw = TRUE)
    system.time(X1_dx <- read_spectrum(urine_1_dx, raw = TRUE, file_format = "jcampdx"))
    expect_equal(X1$si, X1_dx$si)
    expect_equal(X1$cs, X1_dx$cs)
    expect_equal(X1$meta$fq, X1_dx$meta$fq)
    expect_equal(paste0(X1$meta$name, ".dx"), X1_dx$meta$name)
    expect_equal(paste0(X1$meta$path, ".dx"), gsub("jcampdx", "bruker", X1_dx$meta$path))
})
