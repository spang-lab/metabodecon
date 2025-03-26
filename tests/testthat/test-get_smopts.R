test_that("get_smopts works", {
    x <- get_smopts(sim[1:2], c(2,3))
    expect_equal(x, list(sim_01 = c(2,3), sim_02 = c(2,3)))

    x <- get_smopts(sim[1:2], list(c(2,3), c(3,5)))
    expect_equal(x, list(sim_01 = c(2,3), sim_02 = c(3,5)))
})
