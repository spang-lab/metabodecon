library(testthat)

test_that("read_bruker_spectrum works", {
    urine1_dir <- pkg_file("example_datasets/bruker/urine/urine_1")
    S <- read_bruker_spectrum(urine1_dir)
    expect_equal(dim(data.frame(S$si, S$cs, S$fq)), c(131072, 3))
    expect_equal(names(S), c("si", "cs", "fq", "name", "path", "type", "mfs"))
})
