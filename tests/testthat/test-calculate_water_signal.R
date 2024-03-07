test_that("1. determine_water_signal()", {

    # Inputs
    ppm <- c(4.7, 3.4, 2.1, 0.8, -0.5, -1.8, -3.1, -4.4)
    n <- length(ppm)
    ppm_step <- (max(ppm) - min(ppm)) / (n - 1)
    ppm_nstep <- (max(ppm) - min(ppm)) / (n)
    spectrum <- list(x_ppm = ppm, n = n, ppm_step = ppm_step, ppm_nstep = ppm_nstep)
    hwidth_ppm <- 0.2

    # Function call
    ws <- determine_water_signal(spectrum, hwidth_ppm, bwc = TRUE)
    expect_equal(ws, list(
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
