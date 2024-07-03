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
    expect_equal(X1[1, ], structure(list(si = 316.25, cs = 14.80254, fq = 600243921.683815), row.names = 1L, class = "data.frame"))
    expect_equal(X1[131072, ], structure(list(si = 269.5, cs = -5.2210744338963, fq = 600255940.914584), row.names = 131072L, class = "data.frame"))
})

skip_if_slow_tests_disabled()

test_that("read_spectrum works for jcampdx spectra", {
    urine_1 <- pkg_file("example_datasets/bruker/urine/urine_1")
    urine_1_dx <- pkg_file("example_datasets/jcampdx/urine/urine_1.dx")
    X1 <- read_spectrum(urine_1)
    X1_dx <- read_spectrum(urine_1_dx, file_format = "jcampdx")
    expect_equal(X1, X1_dx)
})
