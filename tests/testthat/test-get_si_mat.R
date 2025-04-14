test_that("get_si_mat returns a matrix of the correct dimensions", {
    system.time({
        decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
        aligns <- align(decons)
        si_mat <- get_si_mat(aligns) # 2048 x 2 matrix (2048 datapoints, 2 spectra)
        expect_equal(dim(si_mat), c(2048, 2))
    })
})