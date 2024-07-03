# GLC v13 #####

# IMPORTANT: we dont test for jcampdx files because after calling `read_spectrum()`, the data is the same as for bruker, which is tested in `test-read_spectrum.R`. Also, the calculations in the old `MetaboDecon1D()` function are slightly different for jcampdx and bruker files (in the jcampdx case it calculates with n-1 instead of n) and our `compare_spectra_v13` function currently only accounts for bruker-type errors of MetaboDecon1D, but not jcampdx errors.

skip_if_slow_tests_disabled()

test_that("GLC v13 works for 1 bruker", {
    new <- glc_v13(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE)$rv$urine_1
    old <- MD1D(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE)$rv
    r <- compare_spectra_v13(new, old, silent = FALSE)
    expect_true(sum(r %in% 0:1) >= 60 && sum(r %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
})

test_that("GLC v13 works for n bruker", {
    x <- glc_v13(dp = "urine", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE)$rv
    y <- MD1D(dp = "urine", ff = "bruker", nfit = 3, simple = TRUE)$rv
    r1 <- compare_spectra_v13(new = x$urine_1, old = y$urine_1, silent = TRUE)
    r2 <- compare_spectra_v13(new = x$urine_2, old = y$urine_2, silent = TRUE)
    expect_true(sum(r1 %in% 0:1) >= 60 && sum(r1 %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
    expect_true(sum(r2 %in% 0:1) >= 60 && sum(r2 %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
})

test_that("GLC v13 works for 16 bruker blood samples", {
    x1 <- glc_v13(dp = "blood", nfit = 10, cache = FALSE, nworkers = 1) # 56s on R4
    x2 <- glc_v13(dp = "blood", nfit = 10, cache = FALSE, nworkers = 2) # 40s on R4
    x4 <- glc_v13(dp = "blood", nfit = 10, cache = FALSE, nworkers = 4) # 28s on R4
    x8 <- glc_v13(dp = "blood", nfit = 10, cache = FALSE, nworkers = 8) # 18s on R4
    x16 <- glc_v13(dp = "blood", nfit = 10, cache = FALSE, nworkers = 16) # 15s on R4
    xa <- glc_v13(dp = "blood", nfit = 10, cache = FALSE)
    y <- MD1D(dp = "blood", nfit = 10) # 490s on R4
    r <- lapply(1:16, function(i) compare_spectra_v13(new = xa$rv[[i]], old = y$rv[[i]], silent = TRUE))
    for (i in 1:16) {
        expect_true(sum(r[[i]] %in% 0:1) >= 60 && sum(r[[i]] %in% 2:3) == 0)
        # >=60 identical/equal && no diffs/errors
    }
})

# GLC v12 #####

skip()

test_that("GLC works for 1 bruker", {
    new <- glc_v12(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE)$rv$urine_1
    old <- MD1D(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE)$rv
    r <- compare_spectra(new, old, silent = TRUE)
    expect_true(sum(r %in% 0:1) >= 60 && sum(r %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
})

test_that("GLC work for n bruker", {
    x <- glc_v12(dp = "urine", ff = "bruker", nfit = 3, simple = TRUE, overwrite = FALSE)$rv
    y <- MD1D(dp = "urine", ff = "bruker", nfit = 3, simple = TRUE)$rv
    r1 <- compare_spectra(new = x$urine_1, old = y$urine_1, silent = TRUE)
    r2 <- compare_spectra(new = x$urine_2, old = y$urine_2, silent = TRUE)
    expect_true(sum(r1 %in% 0:1) >= 60 && sum(r1 %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
    expect_true(sum(r2 %in% 0:1) >= 60 && sum(r2 %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
})
