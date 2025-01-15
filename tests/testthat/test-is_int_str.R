library(testthat)

test_that("is_int_str works", {
    fixed_point_notation <- c("2.0", "3.", ".4", "-5.", "-.6")
    scientific_notation <- c("0.45e+04", "66e-05", "0.2e-3", "-33.e-1")

    floats <- c(fixed_point_notation, scientific_notation)
    ints <- c("5", "-5", "1234")
    words <- c("Hello", "world", "!", "It was nice seeing you", ".")

    expect_true(all(is_int_str(floats) == FALSE))
    expect_true(all(is_int_str(ints) == TRUE))
    expect_true(all(is_int_str(words) == FALSE))
})
