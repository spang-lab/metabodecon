library(testthat)

# Reference signal:
#
# x: 1   2   3   4   5   6   7   8   9   10  11  12  13
# y: 100 100 98  92  84  77  72  68  63  56  50  42  33
#
# Second derivative d (computed as y[i-1] + y[i+1] - 2*y[i]) is:
# d: NA  -2  -4  -2  1   2   1   -1  -2  1   -2  -1  NA
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
# Peak centers (PC) are at indices 3, 9, 11. Only PC at index 9 has both a
# left and a right border (indices 8 and 10), so it is the only peak that
# survives the NA filter.

test_that("find_peaks2 detects the expected peak on a hand-crafted signal", {
    withr::local_output_sink(nullfile())
    y <- c(100, 100, 98, 92, 84, 77, 72, 68, 63, 56, 50, 42, 33)
    pk <- find_peaks2(y)
    expect_equal(pk$center, 9L)
    expect_equal(pk$left, 8L)
    expect_equal(pk$right, 10L)
    expect_equal(pk$score, 3)
})

test_that("find_peaks2 returns no peaks on degenerate inputs", {
    withr::local_output_sink(nullfile())

    # Empty signal
    expect_equal(nrow(find_peaks2(numeric(0))), 0L)

    # Signals too short to contain a peak
    expect_equal(nrow(find_peaks2(c(1))), 0L)
    expect_equal(nrow(find_peaks2(c(1, 2))), 0L)
    expect_equal(nrow(find_peaks2(c(1, 2, 3))), 0L)

    # Flat / monotonic signals
    expect_equal(nrow(find_peaks2(rep(1, 20))), 0L)
    expect_equal(nrow(find_peaks2(seq_len(20))), 0L)
    expect_equal(nrow(find_peaks2(rev(seq_len(20)))), 0L)
})

test_that("find_peaks2 detects a single isolated peak", {
    withr::local_output_sink(nullfile())
    y <- c(rep(0, 5), 10, rep(0, 5))
    pk <- find_peaks2(y)
    expect_equal(nrow(pk), 1L)
    expect_equal(pk$center, 6L)
})

test_that("find_peaks2 detects multiple peaks in a synthetic signal", {
    withr::local_output_sink(nullfile())
    # Two well-separated peaks
    y <- c(0, 0, 1, 5, 10, 5, 1, 0, 0, 1, 5, 10, 5, 1, 0, 0)
    pk <- find_peaks2(y)
    expect_equal(nrow(pk), 2L)
    expect_equal(pk$center, c(5L, 12L))
})

test_that("find_peaks2 works on real spectrum data", {
    withr::local_output_sink(nullfile())
    pk <- find_peaks2(sim[[1]]$si)
    expect_s3_class(pk, "data.frame")
    expect_setequal(names(pk), c("left", "center", "right", "score"))
    expect_gt(nrow(pk), 0L)
    # All borders must be valid indices (no NAs after filtering)
    expect_false(anyNA(pk$left))
    expect_false(anyNA(pk$right))
    expect_false(anyNA(pk$center))
    # Borders must bracket the center
    expect_true(all(pk$left < pk$center))
    expect_true(all(pk$center < pk$right))
    # Scores are non-negative
    expect_true(all(pk$score >= 0))
})
