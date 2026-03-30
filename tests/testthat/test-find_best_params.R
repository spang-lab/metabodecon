testthat::test_that("find_best_params returns best row for binary factors", {
    testthat::skip("find_best_params not yet implemented")
    testthat::skip_if_not_installed("e1071")

    set.seed(1)
    X <- matrix(rnorm(24 * 12), nrow = 24)
    y <- factor(rep(c("Control", "AKI"), each = 12),
      levels = c("Control", "AKI"))

    out <- find_best_params(
        X,
        y,
        costs = c(0.2, 0.5),
        gammas = 2^c(-11, -10),
        nfeats = 2:3,
        normalize = TRUE,
        k = 3
    )

    testthat::expect_true(is.list(out))
    testthat::expect_named(out, c("cost", "gamma", "nfeat", "grid"))
    testthat::expect_true(is.data.frame(out$grid))
    testthat::expect_true(all(c("cost", "gamma", "nfeat", "auc") %in%
      names(out$grid)))

    idx <- which.max(out$grid$auc)
    best <- out$grid[idx, ]
    testthat::expect_equal(out$cost, best$cost)
    testthat::expect_equal(out$gamma, best$gamma)
    testthat::expect_equal(out$nfeat, best$nfeat)
})
