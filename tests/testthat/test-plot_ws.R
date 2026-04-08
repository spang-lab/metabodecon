testthat::skip_on_cran()
testthat::skip_on_ci()

test_plot_ws <- function() {
    plot_ws(
        cs = metabodecon::sim[[1]]$cs,
        si = metabodecon::sim[[1]]$si,
        wshw = 0.01
    )
}

test_result <- test_that("plot_ws works", {
    tmp <- vdiffr::expect_doppelganger(
        title = "plot_ws",
        fig = test_plot_ws,
        writer = metabodecon:::make_stable_svg_writer()
    )
})
