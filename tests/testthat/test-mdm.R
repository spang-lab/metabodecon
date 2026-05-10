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

mog1 <- function(rows = 1) {
    data.frame(
        nfit = 1, smit = 1, smws = 3, delta = if (rows == 1) 8 else c(6, 8),
        npmax = 0, maxShift = 50, maxCombine = 20,
        stringsAsFactors = FALSE
    )
}

testthat::test_that("fit_mdm returns mdm with attached mog", {
    m <- fit_mdm(
        sp, y, mog = mog1(2),
        use_rust = 0.5, nworkers = 1, verbosity = 0,
        nfolds = 3
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_true(is.data.frame(m$mog))
    testthat::expect_equal(nrow(m$mog), 2)
    testthat::expect_true(all(c("acc", "auc") %in% names(m$mog)))
})

testthat::test_that("benchmark returns predictions and performance", {
    res <- benchmark(
        sp, y, mog = mog1(1), k = 4,
        use_rust = 0.5, nworkers = 1, verbosity = 0,
        nfolds = 3
    )
    testthat::expect_true(is.data.frame(res$predictions))
    testthat::expect_true("true" %in% names(res$predictions))
    testthat::expect_equal(nrow(res$predictions), length(y))
    testthat::expect_true(is.list(res$overall))
    testthat::expect_true(is.numeric(res$overall$acc))
})

testthat::test_that("get_mog produces required columns", {
    g <- get_mog("default")
    testthat::expect_true(is.data.frame(g))
    cols <- c("nfit", "smit", "smws", "delta", "npmax",
              "maxShift", "maxCombine")
    testthat::expect_true(all(cols %in% names(g)))
})

testthat::test_that("fit_mdm with bin/identity2 returns mdm object", {
    mog_bm <- data.frame(
        nfit = 0L, smit = 0L, smws = 0L, delta = 0,
        npmax = 0L, maxShift = 0L, maxCombine = 64L
    )
    m <- fit_mdm(
        sp, y,
        feat_mat = bin, decon_fun = identity2,
        align_fun = "identity_align",
        mog = mog_bm, igrs = list(),
        nfolds = 3, verbosity = 0
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_true(!is.null(m$model))
    testthat::expect_true("peakPos" %in% names(m$params))
})
