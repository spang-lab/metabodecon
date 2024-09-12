f <- function(x, s) {
    s3 <- round(s/3, 4)
    metabodecon:::logf("Sleeping %s seconds", s3); Sys.sleep(s/3)
    metabodecon:::logf("Sleeping another %s seconds", s3); Sys.sleep(s/3)
    metabodecon:::logf("And once more %s seconds", s3); Sys.sleep(s/3)
    metabodecon:::logf("Done sleeping")
    x
}
rt <- system.time

test_that("mcmapply works", {
    s0.01 <- list(0.01, 0.01, 0.01, 0.01) # ==> 1s
    s2.00 <- list(2.00, 2.00) # ==> 1s
    rt(x1 <- mcmapply(1, f, 6:9, s0.01, log = TRUE))  # 0.06
    rt(x2 <- mcmapply(1, f, 6:9, s0.01, log = FALSE)) # 0.06
    rt(x3 <- mimapply(1, f, 6:9, s0.01))  # 0.17
    rt(x4 <- mcmapply(2, f, list(3, 4), s0.01))  # 0.61
    rt(x5 <- mcmapply(2, f, list(3, 4), s0.01, log = FALSE)) # 0.62
    rt(x6 <- mcmapply(2, f, list(3, 4), s0.01, loadpkg = FALSE)) # 0.32
    rt(x7 <- mcmapply(2, f, list("a", "b"), s2.00, log = TRUE)) # 2.63
    rt(x8 <- mimapply(2, f, list("a", "b"), s2.00)) # 2.67
    # Observations:
    # 1. Process creation takes approx 0.3 sec on my test machine.
    # 2. Loading metabodecon another 0.3 sec. I.e. in total 0.6s.
    # Conclusions:
    # For mapply calls that take less than 0.6 sec, it doesn't make sense
    # to use mcmapply. For mapply calls that take more than one 0.6 sec,
    # mcmapply CAN be faster.
    expect_equal(x1, as.list(6:9))
    expect_equal(x2, as.list(6:9))
    expect_equal(x3, as.list(6:9))
    expect_equal(x4, list(3, 4, 3, 4))
    expect_equal(x5, list(3, 4, 3, 4))
    expect_equal(x6, list(3, 4, 3, 4))
    expect_equal(x7, list("a", "b"))
    expect_equal(x8, list("a", "b"))
})
