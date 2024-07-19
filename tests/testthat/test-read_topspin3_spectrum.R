library(testthat)

test_that("read_topspin3_spectrum works", {
    urine1_dir <- pkg_file("example_datasets/bruker/urine/urine_1")
    S <- read_topspin3_spectrum(urine1_dir)
    expect_equal(dim(S), c(131072, 3))
    expect_equal(colnames(S), c("si", "cs", "fq"))
})
