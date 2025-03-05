testthat::skip_on_cran()
testthat::skip_on_ci()

test_result <- test_that("plot_spectrum works", {
    tmp <- vdiffr::expect_doppelganger(
        title = "plot_spectrum",
        fig = test_plot_spectrum,
        writer = function(plot, file, title = "") {
            with_svg(file, plot(), width = 12, height = 16)
        }
    )
})
