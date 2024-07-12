# GLC v13 #####

# IMPORTANT: we dont test for jcampdx files because after calling `read_spectrum()`, the data is the same as for bruker, which is tested in `test-read_spectrum.R`. Also, the calculations in the old `MetaboDecon1D()` function are slightly different for jcampdx and bruker files (in the jcampdx case it calculates with n-1 instead of n) and our `compare_spectra_v13` function currently only accounts for bruker-type errors of MetaboDecon1D, but not jcampdx errors.

test_that("GLC works for 1 bruker", {
    new <- glc(dp = "sim_01", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE)$rv$sim_01
    old <- md1d(dp = "sim_01", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE)$rv
    r <- compare_spectra_v13(new, old, silent = TRUE)
    expect_true(sum(r %in% 0:1) >= 60 && sum(r %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
})

test_that("GLC works for bruker folder", {
    x <- glc(dp = "sim", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE, nworkers = 1)$rv
    y <- md1d(dp = "sim", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE)$rv
    r1 <- compare_spectra_v13(new = x$sim_01, old = y$sim_01, silent = TRUE)
    r2 <- compare_spectra_v13(new = x$sim_02, old = y$sim_02, silent = TRUE)
    expect_true(sum(r1 %in% 0:1) >= 60 && sum(r1 %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
    expect_true(sum(r2 %in% 0:1) >= 60 && sum(r2 %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
})
