library(testthat)

test_that("snap_to_peakPos snaps to reference positions", {
    M <- rbind(
        c(0, 2, 0, 0, 3, 0, 0),
        c(0, 0, 4, 0, 0, 5, 0)
    )
    obj <- snap_to_peakPos(M, peakPos=c(2, 5), maxCombine=1)
    exp <- rbind(
        c(0, 2, 0, 0, 3, 0, 0),
        c(0, 4, 0, 0, 5, 0, 0)
    )
    expect_equal(obj, exp)
})

test_that("snap_to_peakPos snaps midpoint ties to the left reference", {
    M <- matrix(c(0, 0, 0, 0, 7, 0, 0, 0, 0), nrow=1)
    obj <- snap_to_peakPos(M, peakPos=c(3, 7), maxCombine=2)
    exp <- matrix(c(0, 0, 7, 0, 0, 0, 0, 0, 0), nrow=1)
    expect_equal(obj, exp)
})

test_that("snap_to_peakPos drops peaks farther than maxCombine", {
    M <- matrix(c(0, 0, 0, 0, 7, 0, 0, 0, 0), nrow=1)
    obj <- snap_to_peakPos(M, peakPos=c(3, 7), maxCombine=1)
    exp <- matrix(0, nrow=1, ncol=ncol(M))
    expect_equal(obj, exp)
})

test_that("snap_to_peakPos sums multiple peaks mapping to same peakPos", {
    M <- matrix(c(0, 2, 0, 3, 0, 0, 0), nrow=1)
    obj <- snap_to_peakPos(M, peakPos=c(3), maxCombine=2)
    exp <- matrix(c(0, 0, 5, 0, 0, 0, 0), nrow=1)
    expect_equal(obj, exp)
})
