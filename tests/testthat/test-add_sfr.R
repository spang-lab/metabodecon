test_that("add_sfr works", {
    # Given below PPM values as input borders (3 points left, 2 points right), the correct output would be DP1, but the expected output is DP2 to maintain backwards compatibility.
    # PPM: 4.7  3.4  2.1| 0.8 -0.5  -1.8| -3.1 -4.4
    # DP1: 7.0  6.0  5.0| 4.0  3.0   2.0|  1.0  0.0
    # DP2: 7.0 |6.0  5.0  4.0  3.0|  2.0   1.0  0.0
    spec <- within(list(), {
        ppm <- c(4.7, 3.4, 2.1, 0.8, -0.5, -1.8, -3.1, -4.4)
        dp <- c(7, 6, 5, 4, 3, 2, 1, 0)
        n <- length(ppm)
        ppm_min <- min(ppm)
        ppm_max <- max(ppm)
        ppm_step <- (max(ppm) - min(ppm)) / (n - 1)
        ppm_nstep <- (max(ppm) - min(ppm)) / (n)
        sfx <- 1000
        sfy <- 1000000
    })
    spec <- add_sfr(spec, sfr = c(2.0, -3.0))
    expect_equal(spec$sfr, list(
        right_sdp = 0.00223076923076923, right_dp = 2.23076923076923,
        left_sdp = 0.00662637362637363, left_dp = 6.62637362637363,
        right_ppm = -3, left_ppm = 2
    ))
})
