develop_draw_spectrum <- function() {

    # Achieve clean state
    while (dev.cur() %!=% c(`null device` = 1L)) dev.off()
    rm_all()
    deferred_run()
    untrace(draw_spectrum)

    # Define Input Object
    target <- c("sim1", "sap2")[1]
    obj <- switch(target, "sim1" = get_sim1_decon1(), "sap2" = get_sap2_idecon())
    obj <- as_v2_obj(obj)

    # Get to state within [plot_spectrum()] where [draw_spectrum()] is called
    stub(plot_spectrum, obj = obj)
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
    path <- metabodecon_file("sim/sim_01")
    spec <- read_spectrum(path)
    decon <- generate_lorentz_curves(spec, sfr = c(3.55, 3.35), ws = 0, ask = FALSE, verbose = FALSE)
    plot_dummy <- function() {
        plot(0, 0, ylim = c(0, 1), xlim = c(0, 1), xaxs = "i", yaxs = "i")
        text(0.5, 0.5, "dummy")
    }
    leftmiddle <- c(0.1, 0.4, 0.30, 0.45)
    leftbottom <- c(0.1, 0.4, 0.05, 0.20)
    p <- local({
        local_par(mfrow = c(4, 2), mar = c(2, 2, 0.5, 0.5))
        p <- list()
        p[[1]] <- plot_dummy()
        p[[2]] <- draw_spectrum(obj = spec)
        p[[3]] <- draw_spectrum(obj = spec)
        p[[4]] <- draw_spectrum(obj = decon, foc_rgn = c(3.45, 3.37))
        p[[5]] <- plot_dummy()
        p[[5]] <- draw_spectrum(obj = decon, foc_rgn = c(3.45, 3.37), foc_only = TRUE, fig = leftmiddle, bg_fill = rgb(0, 0, 1, 0.1))
        p[[6]] <- plot_dummy()
        p[[7]] <- draw_spectrum(obj = decon, fig = leftbottom, add = FALSE)
        p[[8]] <- draw_spectrum(obj = decon, lc_show = FALSE, dp_show = FALSE)
        p
    })
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
