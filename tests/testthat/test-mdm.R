# Minimal smoke tests for fit_mdm() and benchmark().
# Uses tiny simulated spectra to keep CI runtime short.

set.seed(1)
n <- 32
npk <- 3
cs <- seq(from = 3.6, length.out = 512, by = -0.0006)
x0 <- sort(runif(npk, 3.42, 3.56))
A <- runif(npk, 8, 14) * 1e3
lam <- runif(npk, 0.9, 1.3) / 1e3
y <- factor(rep(c("A", "B"), each = n / 2))
sp <- vector("list", n)
for (i in seq_len(n)) {
    xi <- x0 + rnorm(npk, sd = 0.0003)
    Ai <- A * runif(npk, 0.8, 1.2)
    li <- lam * runif(npk, 0.9, 1.1)
    Ai[1] <- Ai[1] * (if (y[i] == "A") 1.3 else 0.7)
    sp[[i]] <- simulate_spectrum(
        name = sprintf("s_%02d", i), cs = cs,
        x0 = sort(xi), A = Ai, lambda = li,
        noise = rnorm(length(cs), sd = 500)
    )
}
class(sp) <- "spectra"

testthat::test_that("fit_mdm returns mdm with scalar perf and resolved params", {
    m <- fit_mdm(
        sp, y,
        npmax=0L, maxShift=50L, maxCombine=20L,
        use_rust = 0.5, nworkers = 1, verbosity = 0
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_true(is.finite(m$acc))
    testthat::expect_true(is.finite(m$auc))
    testthat::expect_true(is.integer(m$params$maxShift) ||
                          is.numeric(m$params$maxShift))
})

testthat::test_that("benchmark returns predictions and performance", {
    res <- benchmark(
        sp, y,
        npmax=0L, maxShift=50L, maxCombine=20L, k = 4,
        use_rust = 0.5, nworkers = 1, verbosity = 0
    )
    testthat::expect_true(is.data.frame(res$predictions))
    testthat::expect_true("true" %in% names(res$predictions))
    testthat::expect_equal(nrow(res$predictions), length(y))
    testthat::expect_true(is.list(res$overall))
    testthat::expect_true(is.numeric(res$overall$acc))
})

testthat::test_that("fit_mdm with bin/identity2 returns mdm object", {
    m <- fit_mdm(
        sp, y,
        feat_fun = bin, decon_fun = identity2,
        align_fun = identity_align, snap_fun = identity_snap,
        npmax=0L, maxShift=0L, maxCombine=64L, igrs = list(),
        verbosity = 0
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_true(!is.null(m$model))
    testthat::expect_true("peakPos" %in% names(m$params))
})

testthat::test_that("fit_mdm with fit_ranger returns mdm with OOB scores", {
    testthat::skip_if_not_installed("ranger")
    m <- fit_mdm(
        sp, y,
        npmax=0L, maxShift=50L, maxCombine=20L,
        fit_fun = fit_ranger, predict_fun = predict_ranger,
        use_rust = 0.5, nworkers = 1, verbosity = 0
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_true(inherits(m$model, "ranger"))
    testthat::expect_true(is.finite(m$acc))
    testthat::expect_true(is.finite(m$auc))
    p <- stats::predict(m, sp[1:4], type = "prob", verbosity = 0)
    testthat::expect_length(p, 4L)
})

testthat::test_that("fit_mdm with snap_nw_blind predicts on held-out spectra", {
    m <- fit_mdm(
        sp, y,
        npmax=0L, maxShift=50L, maxCombine=20L,
        snap_fun = snap_nw_blind,
        fit_fun = fit_ranger, predict_fun = predict_ranger,
        use_rust = 0.5, nworkers = 1, verbosity = 0
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_identical(m$params$snap_fun, snap_nw_blind)
    testthat::expect_true(!is.null(m$ref$snap))
    p <- stats::predict(m, sp[1:4], type = "prob", verbosity = 0)
    testthat::expect_length(p, 4L)
    testthat::expect_true(all(is.finite(p)))
})
