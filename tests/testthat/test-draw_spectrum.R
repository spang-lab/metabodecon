testthat::skip_on_cran()
testthat::skip_on_ci()

test_result <- test_that("draw_spectrum works", {
    withr::local_output_sink(nullfile())
    tmp <- vdiffr::expect_doppelganger(
        title = "draw_spectrum",
        fig = test_draw_spectrum,
        writer = metabodecon:::make_stable_svg_writer(
            width = 14,
            height = 14
        )
    )
})
