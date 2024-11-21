
testthat("lorentz_int can produce backwards compatible results", {
    decon <- generate_lorentz_curves_sim(sim[1])
    A <- decon$A
    x0 <- decon$x_0
    lambda <- decon$lambda
    old <- decon$integrals[1, ]
    new <- lorentz_int(x0, A, lambda, limits = c(0, max(decon$x_values) + 0.001))
    # + 0.001 because we have an error in the orignal function where the
    # scaled datapoints are shifted by 0.001 during integral calculation
    expect_equal(old, new)
})

testthat("using `A * pi` instead of `integrals` produces similar results", {
    decon <- generate_lorentz_curves_sim(sim[1])
    A <- decon$A
    x0 <- decon$x_0
    lambda <- decon$lambda
    old <- decon$integrals[1, ]
    new <- (abs(A) * pi)
    err <- abs(old - new)
    ratio <- old / new
    ones <- rep(1, length(ratio))
    expect_equal(ratio, ones, tolerance = 0.01) # less than 1% difference
})
