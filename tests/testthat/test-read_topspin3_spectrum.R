library(testthat)

test_that("read_topspin3_spectrum works", {
    blood1_dir <- pkg_file("example_datasets/bruker/blood/blood_01")
    S <- read_topspin3_spectrum(blood1_dir)
    expect_equal(dim(S), c(131072, 3))
    expect_equal(colnames(S), c("si", "cs", "fq"))
})
