test_ps_internal <- function() {
    sim_01 <- metabodecon_file("sim/sim_01")
    decon <- generate_lorentz_curves(
        sim_01,
        sfr = c(3.42, 3.58), ws = 0,
        ask = FALSE, verbose = FALSE,
        delta = 0.1
    )
    plot_dummy <- function() {
        plot(0, 0, ylim = c(0, 1), xlim = c(0, 1), xaxs = "i", yaxs = "i")
        text(0.5, 0.5, "dummy")
    }
    leftmiddle <- c(0.1, 0.4, 0.4, 0.6)
    leftbottom <- c(0.1, 0.4, 0.1, 0.3)
    p <- local({
        p <- list()
        opar <- par(mfrow = c(3, 2), mar = c(2, 2, 0.5, 0.5))
        on.exit(par(opar))
        p[[1]] <- plot_dummy()
        p[[2]] <- ps_internal(decon, foc_rgn = c(3.55, 3.52))
        p[[3]] <- plot_dummy()
        p[[3]] <- ps_internal(decon,
            foc_rgn = c(3.55, 3.52),
            foc_only = TRUE, fig = leftmiddle,
            fill_col = rgb(0, 0, 1, 0.1)
        )
        p[[4]] <- plot_dummy()
        p[[5]] <- ps_internal(decon, fig = leftbottom, add = FALSE)
        p[[6]] <- ps_internal(decon, lc_show = FALSE, trp_show = FALSE)
        p
    })
}
