check_prarp <- test_that("deconvolute_ispec works", {
    x <- deconvolute_ispec(
        as_ispec(sim[[1]]),
        sfr = c(3.55, 3.35),
        wshw = 0,
        ask = FALSE,
        verbose = FALSE
    )
    obj <- calc_prarp(x, truepar = sim[[1]]$meta$simpar)
    expect_true(obj$prarp >= 0.836)
    expect_true(obj$prarpx >= 0.779)
})

check_prarp <- test_that("bwc2 > bwc0", {
    args0 <- args2 <- list(
        as_ispec(sim[[1]]), sfr = c(3.55, 3.35), wshw = 0, ask = FALSE,
        verbose = FALSE, bwc = 0
    )
    args2$bwc <- 2
    x <- do.call(deconvolute_ispec, args0)
    y <- do.call(deconvolute_ispec, args2)
    calc_prarp(x)$prarpx
    calc_prarp(y)$prarpx
})