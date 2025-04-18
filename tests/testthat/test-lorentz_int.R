# Inputs
sim_dir <- metabodecon_file("sim_subset")
decon <- MetaboDecon1D_silent_sim(filepath = sim_dir, filename = "sim_01")
A <- decon$A
x0 <- decon$x_0
lambda <- decon$lambda
old <- decon$integrals[1, ]

test_that("lorentz_int can produce backwards compatible results", {
    new <- lorentz_int(x0, A, lambda, limits = c(0, max(decon$x_values) + 0.001))
    # We add + 0.001 because there is an error in the orignal function where the
    # limits for calculating the integral are set one "scaled data point" (i.e.
    # 1/1000) too high. (In fact, setting any limits at all is inferior compared
    # to using A * pi for calculating the integrals).
    expect_equal(old, new)
})

test_that("using `A * pi` instead of `integrals` produces similar results", {
    old <- decon$integrals[1, ]
    new <- (abs(A) * pi)
    err <- abs(old - new)
    ratio <- old / new
    ones <- rep(1, length(ratio))
    expect_equal(ratio, ones, tolerance = 0.01) # less than 1% difference
})
