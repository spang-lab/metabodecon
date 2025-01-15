testthat::skip_on_cran()
testthat::skip_on_ci()

develop_draw_spectrum <- function() {

    # Achieve clean state
    while (dev.cur() %!=% c(`null device` = 1L)) dev.off()
    rm_all()
    deferred_run()
    untrace(draw_spectrum)

    # Define Input Object
    target <- c("sim1", "sap1")[1]
    obj <- switch(target,
        "sim1" = deconvolute(sim[[1]], sfr = c(3.55, 3.35)),
        "sap1" = deconvolute(sap[[1]], sfr = c(3.2, -3.2), smopts = c(1, 3), delta = 3)
    )

    # Get to state within [plot_spectrum()] where [draw_spectrum()] is called
    stub(plot_spectrum, obj = obj, ... = NULL)
    foc_frac <- foc_frac %||% get_foc_frac(obj, foc_rgn)
    foc_rgn <- foc_rgn %||% get_foc_rgn(obj, foc_frac)
    layout <- layout %||% get_ps_layout(obj, foc_rgn)
    local_par(mar = mar)
    plot_empty()
    args <- get_ds_arglists(obj, foc_rgn, foc_frac, layout, args1, args2, args3)

    # Stub [draw_spectrum()] with correspondings args as defined by [plot_spectrum()]
    invisible(do.call(stub, c(draw_spectrum, args[[1]])))

    # Not "Go to Definition" and start developing
    if (FALSE) draw_spectrum(obj)
}

test_draw_spectrum <- function() {
    decon <- deconvolute(sim[[1]], sfr = c(3.55, 3.35))
    local_par(mfrow = c(4, 2), mar = c(2, 2, 0.5, 0.5))
    plot_dummy()
    draw_spectrum(obj = decon)
    draw_spectrum(obj = decon, lgd = list(x = "top", bg = NA))
    draw_spectrum(obj = decon, foc_rgn = c(3.45, 3.37))
    plot_dummy()
    draw_spectrum(obj = decon, fig = c(0.1, 0.4, 0.30, 0.45), add = TRUE)
    plot_dummy()
    draw_spectrum(obj = decon, fig = c(0.1, 0.4, 0.05, 0.20), add = FALSE)
    draw_spectrum(obj = decon, lc_lines = NULL, lc_rects = NULL, foc_only = FALSE)
}

test_result <- test_that("draw_spectrum works", {
    tmp <- vdiffr::expect_doppelganger(
        title = "draw_spectrum",
        fig = test_draw_spectrum,
        writer = function(plot, file, title = "") {
            with_svg(file, plot(), width = 12, height = 16)
        }
    )
})
