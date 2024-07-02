library(testthat)

test_that("read_acqus_file works", {
    blood1_dir <- pkg_file("example_datasets/bruker/blood/blood_01")
    acqus <- read_acqus_file(blood1_dir)
    expect_equal(as.numeric(acqus$SW), 20.0236139622347)
    expect_equal(as.numeric(acqus$SW_h), 12019.2307692308)
})
