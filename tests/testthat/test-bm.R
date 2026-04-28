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

testthat::test_that("benchmark_bm returns predictions", {
    bm <- benchmark_bm(
        sp, y,
        regions = rbind(c(3.60, 3.30)),
        binwidth = 0.02,
        nfo = 4, nfl = 3, nwo = 1
    )
    testthat::expect_true(is.data.frame(bm$predictions))
    testthat::expect_equal(nrow(bm$predictions), length(y))
    testthat::expect_true(all(vapply(bm$models, inherits, logical(1), what = "bm")))
})