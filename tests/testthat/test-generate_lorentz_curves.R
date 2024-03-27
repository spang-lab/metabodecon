skip_if_slow_tests_disabled()

test_that("GLC works for 1 bruker", {
    new <- glc(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE)$rv$urine_1
    old <- MD1D(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE)$rv
    r <- compare_spectra(new, old, silent = TRUE)
    expect_true(sum(r %in% 0:1) >= 60 && sum(r %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
})

test_that("GLC work for n bruker", {
    x <- glc(dp = "urine_2", ff = "bruker", nfit = 3, simple = TRUE)$rv$urine_2
    y <- MD1D(dp = "urine_2", ff = "bruker", nfit = 3, simple = TRUE)$rv
    r <- compare_spectra(new = x, old = y, silent = TRUE)
    expect_true(sum(r %in% 0:1) >= 60 && sum(r %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
})

# We dont test for jcampdx files because after calling `read_spectrum()`, the data is the same as for bruker, which is tested in `test-read_spectrum.R`. Also, the calculations in the old `MetaboDecon1D()` function are slightly different for jcampdx and bruker files (in the jcampdx case it calculates with n-1 instead of n) and our compare function currently only accounts for bruker-case errors of MetaboDecon1D, but not jcampdx-case errors.

su1 <- MD1D(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE)$rv
su2 <- MD1D(dp = "urine_2", ff = "bruker", nfit = 3, simple = TRUE)$rv
mu <- MD1D(dp = "urine", ff = "bruker", nfit = 3, simple = TRUE)$rv
mu1 <- mu$urine_1
mu2 <- mu$urine_2
