library(testthat)

test_that("read_bruker works", {
    urine1_dir <- pkg_file("example_datasets/bruker/urine/urine_1")
    x <- read_bruker(urine1_dir)
    expect_equal(names(x), c("si", "cs", "meta"))
    expect_equal(length(x$si), 131072)
    expect_equal(length(x$cs), 131072)
})
