library(testthat)

test_that("combine_scores works", {
    M <- rbind(
        c(2, 0, 2, 2, 0),
        c(2, 1, 0, 2, 0),
        c(0, 1, 0, 2, 0),
        c(0, 0, 3, 0, 1)
    )
    U <- M != 0
    uu <- colSums(U) # 2 2 2 3 1

    cc <- combine_scores(U, uu, j = 2, nn = c(1, 3))
    expect_true(length(cc) == 2)
    expect_true(cc[1] == 0) # M[, 1] and M[, 2] are not combinable
    expect_true(cc[2] == 2) # M[, 3] and M[, 2] are combinable and M[, 3] has two nonzero elements

    cc <- combine_scores(U, uu, j = 1, nn = 2:5)
    expect_true(length(cc) == 4)
    expect_true(cc[1] == 0) # M[, 2] and M[, 1] are not combinable
    expect_true(cc[2] == 0) # M[, 3] and M[, 1] are not combinable
    expect_true(cc[3] == 0) # M[, 4] and M[, 1] are not combinable
    expect_true(cc[4] == 1) # M[, 5] and M[, 1] are combinable and M[, 5] has one nonzero element
})
