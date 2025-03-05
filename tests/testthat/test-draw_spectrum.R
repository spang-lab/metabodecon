testthat::skip_on_cran()
testthat::skip_on_ci()

test_result <- test_that("draw_spectrum works", {
    tmp <- vdiffr::expect_doppelganger(
        title = "draw_spectrum",
        fig = test_draw_spectrum,
        writer = function(plot, file, title = "") {
            with_svg(file, plot(), width = 14, height = 14)
        }
    )
})
