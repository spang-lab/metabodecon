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
