# Tests for the lightweight speaq replacements in R/speaq.R.
# Expected values were generated with the original speaq functions
# using the script in tmp/gen_expectations.R.

# find_ref -----------------------------------------------------------

test_that("find_ref picks the best reference spectrum", {
    # Generated with:
    #   pl <- list(c(10, 50, 100, 200), c(12, 48, 102, 198),
    #              c(11, 51, 99, 201))
    #   speaq::findRef(pl)
    #   => $refInd = 1, $orderSpec = c(1, 3, 2)
    pl <- list(
        c(10, 50, 100, 200),
        c(12, 48, 102, 198),
        c(11, 51, 99, 201)
    )
    res <- find_ref(pl)
    expect_equal(res$refInd, 1)
    expect_equal(res$orderSpec, c(1, 3, 2))
})

# fft_shift ----------------------------------------------------------

test_that("fft_shift detects a known shift", {
    # Generated with:
    #   set.seed(1)
    #   n <- 64
    #   ref <- sin(seq(0, 4*pi, length.out = n)) + rnorm(n, 0, 0.1)
    #   tar <- c(rep(0, 3), ref[1:(n-3)])
    #   speaq::findShiftStepFFT(ref, tar, maxShift = 10)
    #   => $stepAdj = -3, $corValue = 495781616782.328
    set.seed(1)
    n <- 64
    ref <- sin(seq(0, 4 * pi, length.out = n)) + rnorm(n, 0, 0.1)
    tar <- c(rep(0, 3), ref[1:(n - 3)])
    res <- fft_shift(ref, tar, maxShift = 10)
    expect_equal(res$stepAdj, -3L)
    expect_equal(res$corValue, 495781616782.328, tolerance = 1e-6)
})

test_that("fft_shift returns 0 when spectra are identical", {
    # Generated with:
    #   ref2 <- c(0, 0, 5, 0, 8, 0, 0)
    #   speaq::findShiftStepFFT(ref2, ref2, maxShift = 5)
    #   => $stepAdj = 0
    ref <- c(0, 0, 5, 0, 8, 0, 0)
    res <- fft_shift(ref, ref, maxShift = 5)
    expect_equal(res$stepAdj, 0L)
})

# do_shift -----------------------------------------------------------

test_that("do_shift shifts right (positive step)", {
    # Generated with:
    #   speaq::doShift(1:8, 2)
    #   => c(1, 1, 1, 2, 3, 4, 5, 6)
    expect_equal(do_shift(1:8, 2), c(1, 1, 1, 2, 3, 4, 5, 6))
})

test_that("do_shift shifts left (negative step)", {
    # Generated with:
    #   speaq::doShift(1:8, -3)
    #   => c(4, 5, 6, 7, 7, 7, 7, 7)
    expect_equal(do_shift(1:8, -3), c(4, 5, 6, 7, 7, 7, 7, 7))
})

test_that("do_shift with step 0 matches speaq quirk", {
    # Generated with:
    #   speaq::doShift(1:8, 0)
    #   => c(1, 2, 3, 4, 5, 6, 7, 7)  # last elem overwritten
    expect_equal(do_shift(1:8, 0), c(1, 2, 3, 4, 5, 6, 7, 7))
})

# hclust_align -------------------------------------------------------

test_that("hclust_align corrects a +2 shift across 3 peaks", {
    # Generated with:
    #   n2 <- 128
    #   refSpec <- double(n2)
    #   refSpec[30] <- 5; refSpec[60] <- 8; refSpec[100] <- 3
    #   tarSpec <- double(n2)
    #   tarSpec[32] <- 5; tarSpec[62] <- 8; tarSpec[102] <- 3
    #   peakList <- c(30, 60, 100, 32, 62, 102)
    #   peakLabel <- c(1, 1, 1, 0, 0, 0)
    #   speaq::hClustAlign(refSpec, tarSpec, peakList, peakLabel,
    #       startP = 1, endP = 128, maxShift = 10,
    #       acceptLostPeak = FALSE)
    #   => $peakList = c(30, 60, 100, 30, 60, 99)
    #   => nonzero tarSpec at indices 30, 60, 99:104
    #   => nonzero tarSpec values 5, 8, 3, 3, 3, 3, 3, 3
    n2 <- 128
    refSpec <- double(n2)
    refSpec[30] <- 5; refSpec[60] <- 8; refSpec[100] <- 3
    tarSpec <- double(n2)
    tarSpec[32] <- 5; tarSpec[62] <- 8; tarSpec[102] <- 3
    peakList <- c(30, 60, 100, 32, 62, 102)
    peakLabel <- c(1, 1, 1, 0, 0, 0)
    res <- hclust_align(refSpec, tarSpec, peakList, peakLabel,
        startP = 1, endP = n2, maxShift = 10)
    expect_equal(res$peakList, c(30, 60, 100, 30, 60, 99))
    nz <- which(res$tarSpec != 0)
    expect_equal(nz, c(30, 60, 99, 100, 101, 102, 103, 104))
    expect_equal(
        res$tarSpec[nz],
        c(5, 8, 3, 3, 3, 3, 3, 3)
    )
})

test_that("hclust_align does not shift already-aligned spectra", {
    # Generated with:
    #   n2 <- 128
    #   refSpec2 <- double(n2); refSpec2[30] <- 5; refSpec2[60] <- 8
    #   tarSpec2 <- double(n2); tarSpec2[30] <- 4; tarSpec2[60] <- 7
    #   peakList2 <- c(30, 60, 30, 60)
    #   peakLabel2 <- c(1, 1, 0, 0)
    #   speaq::hClustAlign(refSpec2, tarSpec2, peakList2, peakLabel2,
    #       startP = 1, endP = 128, maxShift = 10,
    #       acceptLostPeak = FALSE)
    #   => $peakList = c(30, 60, 30, 60)
    #   => nonzero tarSpec at 30, 60
    n2 <- 128
    refSpec <- double(n2); refSpec[30] <- 5; refSpec[60] <- 8
    tarSpec <- double(n2); tarSpec[30] <- 4; tarSpec[60] <- 7
    peakList <- c(30, 60, 30, 60)
    peakLabel <- c(1, 1, 0, 0)
    res <- hclust_align(refSpec, tarSpec, peakList, peakLabel,
        startP = 1, endP = n2, maxShift = 10)
    expect_equal(res$peakList, c(30, 60, 30, 60))
    expect_equal(which(res$tarSpec != 0), c(30, 60))
})
