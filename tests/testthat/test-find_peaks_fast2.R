library(testthat)

# Inputs (Second Derivative (d) and Peak Center Indices (pd))
x <- c(1,   2,   3,  4,  5,   6,   7,   8,   9,   10,  11,  12,  13)
y <- c(100, 100, 98, 92, 84,  77,  72,  68,  63,  56,  50,  42,  33)
# plot(x, y)
d <- calc_second_derivative(y)
expect_equal(d, c(NA, -2, -4, -2,   1,   2,  1,   -1,  -2,   1,  -2,  -1,  NA))
#
#        Nearest maximum or root right
#        right of PC1 (Peak Center 1).
#        ==> Right Border for PC1.
#                |
#                |     Nearest maximum or root right of
#                |     PC2 AND nearest maximum or root
#                |     left of PC3 ==> Right Border for
#                |     PC2 AND left border for PC3.
#                |                 |
#     ___________|_________________|____________
#  2 |           |    ###          |            |
#  1 |           | #########      R2L3          |
# -1 | __ #########         #L2###    ###### __ |
# -1 |    #########            PC2    ######    |
# -2 |    #######R1                   PC3       |
# -3 |       ###                                |
# -4 |       PC1                                |
#    |__________________________________________|
#       1  2  3  4  5  6  7  8  9  10  11 12 13
#

test_that("find_peaks works", {

    withr::local_output_sink(nullfile())

    # Individual steps
    pc <- get_peak_centers_fast(d)
    rb <- get_right_borders_fast(d, pc)
    lb <- get_left_borders_fast(d, pc)
    sc <- get_peak_scores_fast(d, pc, lb, rb)
    expect_equal(pc, c(3, 9, 11))
    expect_equal(rb, c(4, 10, NA))
    expect_equal(lb, c(NA, 8, 10))
    expect_equal(sc, c(0, 3, 0))

    # Final output
    pk <- data.frame(left = lb, center = pc, right = rb, score = sc)
    pk <- pk[!is.na(pk$left) & !is.na(pk$right), ]
    rownames(pk) <- NULL
    peaks <- find_peaks2(y)
    expect_equal(peaks, pk)

    # Compare with C
    peaks_c <- find_peaks2(y, use_c = TRUE)
    expect_equal(peaks_c, peaks)

    # Real data
    pk_r <- find_peaks2(sim[[1]]$si, use_c = FALSE)
    pk_c <- find_peaks2(sim[[1]]$si, use_c = TRUE)
    expect_equal(pk_c, pk_r)

    # Edge cases
    for (use_c in c(FALSE, TRUE)) {
        # Empty signal: no peaks
        expect_equal(nrow(find_peaks2(numeric(0), use_c=use_c)), 0L)

        # Signals too short to contain any peak
        expect_equal(nrow(find_peaks2(c(1), use_c=use_c)), 0L)
        expect_equal(nrow(find_peaks2(c(1, 2), use_c=use_c)), 0L)
        expect_equal(nrow(find_peaks2(c(1, 2, 3), use_c=use_c)), 0L)

        # Flat signal: no peaks
        expect_equal(nrow(find_peaks2(rep(1, 20), use_c=use_c)), 0L)

        # Single isolated peak (well-separated from edges)
        y_one <- c(rep(0, 5), 10, rep(0, 5))
        expect_equal(nrow(find_peaks2(y_one, use_c=use_c)), 1L)
    }

})
