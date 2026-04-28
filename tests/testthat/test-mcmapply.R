
# Observations:
# 1. Process creation on Windows takes approx 0.25 sec on my test machine. See x3.
# 2. Loading metabodecon takes another 0.2 sec. I.e. in total 0.45s. See x4.
# 3. On MacOS/Unix, process creation is much faster, because we can use forking (0.03 instead of 0.25).

# Conclusions:
# 1. For mapply calls that take less than 0.45 sec, it doesn't make sense to use mcmapply on Windows.
# 2. For mapply calls that take more than 0.45 sec, mcmapply CAN be faster.
# 3. On MacOS/Linux mcmapply can be used for even faster function.

test_that("mcmapply works", {
    f <- function(s, x) { message("Sleeping: ", s); Sys.sleep(s); x}
    out <- tmpfile()
    err <- tmpfile()
    local_message_sink(out)
    local_output_sink(err)
    system.time(x1 <- mcmapply(1, f, list(0, 0), 6:9,        log = FALSE, loadpkg = FALSE)) # Win: 0.020, Linux: 0.003
    system.time(x2 <- mcmapply(1, f, 0,          list(6, 7), log = FALSE, loadpkg = TRUE))  # Win: 0.000, Linux: 0.001
    system.time(x3 <- mcmapply(2, f, rep(0, 2),  3:4,        log = TRUE,  loadpkg = FALSE)) # Win: 0.250, Linux: 0.034
    system.time(x4 <- mcmapply(2, f, rep(0, 4),  3:4,        log = TRUE,  loadpkg = TRUE))  # Win: 0.440, Linux: 0.036
    system.time(x5 <- mcmapply(2, f, 0.5,        3:4,        log = TRUE,  loadpkg = FALSE)) # Win: 0.730, Linux: 0.537
    if (identical(environment(), globalenv())) deferred_run() # for interactive testing
    expect_equal(x1, as.list(6:9))
    expect_equal(x2, as.list(6:7))
    expect_equal(x3, list(3, 4))
    expect_equal(x4, list(3, 4, 3, 4))
    expect_equal(x5, list(3, 4))
})
