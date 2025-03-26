testthat::skip_on_cran()
testthat::skip_on_ci()

test_plot_sfr <- function() {
    plot_sfr(
        cs = sim[[1]]$cs,
        si = sim[[1]]$si,
        sfr = c(3.55, 3.35)
    )
}

test_result <- test_that("plot_sfr works", {
    tmp <- vdiffr::expect_doppelganger(
        title = "plot_sfr",
        fig = test_plot_sfr,
        writer = function(plot, file, title = "") {
            with_svg(file, plot())
        }
    )
})
