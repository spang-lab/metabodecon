test_that("smooth_signals_v12 works", {
    y <- c(1, 4, 2, 4, 3)
    spec <- list(y_pos = y)
    z1 <- evalwith(output = "captured", smooth_signals_v12(spec, k = 3, reps = 1))$rv$y_smooth
    z2 <- evalwith(output = "captured", smooth_signals_v12(spec, k = 3, reps = 2))$rv$y_smooth
    z3 <- evalwith(output = "captured", smooth_signals_v12(spec, k = 3, reps = 3))$rv$y_smooth
    expect_equal(z1, sapply(list(y[1:2], y[1:3], y[2:4], y[3:5], y[4:5]), mean))
    expect_equal(z2, sapply(list(z1[1:2], z1[1:3], z1[2:4], z1[3:5], z1[4:5]), mean))
    expect_equal(z3, sapply(list(z2[1:2], z2[1:3], z2[2:4], z2[3:5], z2[4:5]), mean))
})

test_that("smooth_signals_v20 works", {
    y <- c(1, 4, 2, 4, 3)
    spec <- list(y_pos = y)
    z1 <- evalwith(output = "captured", smooth_signals_v20(spec, k = 3, reps = 1))$rv$y_smooth
    z2 <- evalwith(output = "captured", smooth_signals_v20(spec, k = 3, reps = 2))$rv$y_smooth
    z3 <- evalwith(output = "captured", smooth_signals_v20(spec, k = 3, reps = 3))$rv$y_smooth
    expect_equal(z1, sapply(list(y[1:2], y[1:3], y[2:4], y[3:5], y[4:5]), mean))
    expect_equal(z2, sapply(list(z1[1:2], z1[1:3], z1[2:4], z1[3:5], z1[4:5]), mean))
    expect_equal(z3, sapply(list(z2[1:2], z2[1:3], z2[2:4], z2[3:5], z2[4:5]), mean))
})

test_that("runtime(smooth_signals_v20) < runtime(smooth_signals_v12)", {
    r1 <- n <- 0
    while (r1 == 0) {
        spec <- list(y_pos = rnorm(2^10 * 2^(n <- n + 1)))
        r1 <- system.time(smooth_signals_v12(spec, k = 5, reps = 10, verbose = FALSE))[["elapsed"]]
    }
    r2 <- system.time(smooth_signals_v20(spec, k = 5, reps = 10))[["elapsed"]]
    if (r2 >= r1) message("r1: ", r1, " r2: ", r2)
    expect_true(r2 < r1)
})
