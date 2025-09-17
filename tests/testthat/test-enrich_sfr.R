test_that("enrich_sfr works", {

    # Given below PPM values as input borders (3 points left, 2 points right),
    # the correct output would be DP1, but the expected output is DP2 to
    # maintain backwards compatibility with MetaboDecon1D (v0.2.2).
    #
    # PPM: 4.7   3.4   2.1 | 0.8  -0.5   -1.8 | -3.1  -4.4
    # DP1: 7.0   6.0   5.0 | 4.0   3.0    2.0 |  1.0   0.0
    # DP2: 7.0 | 6.0   5.0   4.0   3.0 |  2.0    1.0   0.0

    ppm       <- c(4.7, 3.4, 2.1, 0.8, -0.5, -1.8, -3.1, -4.4)
    dp        <- c(7,   6,   5,   4,    3,    2,    1,    0)
    n         <- length(ppm)
    ppm_min   <- min(ppm)
    ppm_max   <- max(ppm)
    ppm_step  <- (max(ppm) - min(ppm)) / (n - 1)
    ppm_nstep <- (max(ppm) - min(ppm)) / (n)
    sf        <- c(1e3, 1e6)
    spec <- named(ppm, dp, n, ppm_min, ppm_max, ppm_step, ppm_nstep, sf)
    class(spec) <- "ispec"

    sfr <- enrich_sfr(sfr = c(2.0, -3.0), x = spec)
    expect_equal(sfr, list(
        left_ppm  = 2,
        right_ppm = -3,
        left_dp   = 6.62637362637363,
        right_dp  = 2.23076923076923,
        left_sdp  = 0.00662637362637363,
        right_sdp = 0.00223076923076923
    ))

    sfr2 <- enrich_sfr2(sfr = c(2.0, -3.0), cs = ppm)
    expect_equal(sfr, sfr2)

})
