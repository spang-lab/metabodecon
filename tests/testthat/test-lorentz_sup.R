library(testthat)

# Older (and slower) reference implementations #####

# v1: loop over evaluation points
lorentz_sup_v1 <- function(x, x0, Al, l2) {
    sapply(x, function(xi) {
        sum(Al / (l2 + (xi - x0)^2))
    })
}

# v2: loop over peaks
lorentz_sup_v2 <- function(x, x0, Al, l2) {
    y <- numeric(length(x))
    for (j in seq_along(x0)) {
        y <- y + (Al[j] / (l2[j] + (x - x0[j])^2))
    }
    y
}

# Tests #####

test_that("lorentz_sup == lorentz for a single peak", {
    x <- c(-2, -1, 0, 1, 2)
    y <- lorentz_sup(x, x0 = 0, Al = 1, l2 = 1)
    expect_equal(y, 1 / (1 + x^2))
})

test_that("lorentz_sup: A/lambda interface matches Al/l2 interface", {
    x  <- seq(0, 10, length.out = 50)
    x0 <- c(2, 5, 8); A <- c(1, -2, 0.5); lambda <- c(0.3, 0.5, 0.2)
    Al <- abs(A * lambda); l2 <- lambda^2
    expect_equal(
        lorentz_sup(x, x0, A = A, lambda = lambda),
        lorentz_sup(x, x0, Al = Al, l2 = l2)
    )
})

test_that("lorentz_sup: lcpar interface matches explicit args", {
    x     <- seq(0, 10, length.out = 50)
    x0    <- c(2, 5, 8); A <- c(1, -2, 0.5); lambda <- c(0.3, 0.5, 0.2)
    lcpar <- data.frame(x0 = x0, A = A, lambda = lambda)
    expect_equal(
        lorentz_sup(x, x0, A = A, lambda = lambda),
        lorentz_sup(x, x0, lcpar = lcpar)
    )
})

test_that("lorentz_sup_c == lorentz_sup_v1 == lorentz_sup_v2", {
    set.seed(42)
    x  <- seq(0, 10, length.out = 200)
    x0 <- runif(15, 1, 9); A <- runif(15, 0.5, 2); lambda <- runif(15, 0.1, 0.5)
    Al <- abs(A * lambda); l2 <- lambda^2
    yC <- lorentz_sup(x, x0, Al = Al, l2 = l2)
    expect_equal(yC, lorentz_sup_v1(x, x0, Al, l2))
    expect_equal(yC, lorentz_sup_v2(x, x0, Al, l2))
})
