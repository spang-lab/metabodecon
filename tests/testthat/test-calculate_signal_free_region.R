test_that("1. determine_signal_free_region()", {

    # Inputs
    ppm <- c(4.7, 3.4, 2.1, 0.8, -0.5, -1.8, -3.1, -4.4)
    dp <-  c(7,   6,   5,   4,   3,    2,    1,    0)
    n <- length(ppm)
    ppm_step <- (max(ppm) - min(ppm)) / (n - 1)
    ppm_nstep <- (max(ppm) - min(ppm)) / (n)
    spectrum <- list(ppm = ppm, dp = dp, ppm_min = min(ppm), ppm_max = max(ppm), ppm_step = ppm_step, ppm_nstep = ppm_nstep)
    left_ppm <- 2.0 # i.e. 3 points left
    right_ppm <- -3.0 # i.e. 2 points right

    # Function call
    sfr <- determine_signal_free_region(spectrum, left_ppm, right_ppm, sf = 1000, bwc = TRUE)
    expect_equal(sfr, list(
        left_ppm = 2, #               c(4.7,   3.4, 2.1, | 0.8, -0.5, -1.8, -3.1, -4.4)
        left_dp = 6.62637362637363, # c(7,   | 6,   5,     4,   3,    2,    1,    0)
        left_sdp = 0.00662637362637363,
        right_ppm = -3, #              c(4.7, 3.4, 2.1, 0.8, -0.5,   -1.8, | -3.1, -4.4)
        right_dp = 2.23076923076923, # c(7,   6,   5,   4,   3,    | 2,      1,    0)
        right_sdp = 0.00223076923076923
    ))

    # Function call
    sfr <- determine_signal_free_region(spectrum, left_ppm, right_ppm, sf = 1000, bwc = FALSE)
    expect_equal(sfr, list(
        left_ppm = 2, #               c(4.7, 3.4, 2.1, | 0.8, -0.5, -1.8, -3.1, -4.4)
        left_dp = 4.92307692307692, # c(7,   6,   5,   | 4,   3,    2,    1,    0)
        left_sdp = 0.00492307692307692,
        right_ppm = -3, #              c(4.7, 3.4, 2.1, 0.8, -0.5, -1.8, | -3.1, -4.4)
        right_dp = 1.07692307692308, # c(7,   6,   5,   4,   3,    2,    | 1,    0)
        right_sdp = 0.00107692307692308
    ))
})
