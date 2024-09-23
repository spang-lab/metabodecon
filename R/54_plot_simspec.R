
plot_si_mat <- function(Y = generate_lorentz_curves_sim("bruker/sim"),
                        lgdcex = "auto",
                        main = NULL,
                        mar = par("mar")) {
    n <- nrow(Y) # Number of Spectra
    p <- ncol(Y) # Number of Datapoints
    cols <- rainbow(n) # Colors
    dpis <- seq_len(p) # Datapoint Indices
    spis <- seq_len(n) # Spectrum Indices
    ltxt <- paste("Spectrum", spis)
    xlab <- "Datapoint Number"
    ylab <- "Signal Intensity"
    xlim <- c(1, p)
    ylim <- c(0, max(Y))
    args <- named(x = NA, type = "n", xlim, ylim, xlab, ylab)
    opar <- par(mar = mar)
    on.exit(par(opar))
    do.call(plot, args)
    for (i in spis) lines(x = dpis, y = Y[i, ], col = cols[i])
    if (lgdcex == "auto") lgdcex <- 1 / max(1, log(n, base = 8))
    legend(x = "topright", legend = ltxt, col = cols, lty = 1, cex = lgdcex)
    if (!is.null(main)) title(main = main)
}

#' @noRd
#' @title Plots a spectrum simulated with [get_sim_params()]
#' @param simspec A simulated spectrum as returned by [get_sim_params()].
#' @examples
#' simspec <- get_sim_params()
#' plot_sim_spec(simspec)
plot_sim_spec <- function(simspec = get_sim_params()) {
    X <- simspec$X
    top_ticks <- seq(from = min(X$cs), to = max(X$cs), length.out = 5)
    top_labels <- round(seq(from = min(X$fq), to = max(X$fq), length.out = 5))
    line1_text <- sprintf("Name: %s", gsub("blood", "sim", simspec$filename))
    line2_text <- sprintf("Base: %s ", simspec$filename)
    line3_text <- sprintf("Range: 3.6-3.4 ppm")
    plot(
        x = X$cs, y = X$si_raw * 1e-6, type = "l",
        xlab = "Chemical Shift [PPM]", ylab = "Signal Intensity [AU]",
        xlim = c(3.6, 3.4), col = "black"
    )
    lines(x = X$cs, y = X$si_smooth, col = "blue")
    lines(x = X$cs, y = X$si_sim, col = "red")
    legend(
        x = "topleft", lty = 1,
        col = c("black", "blue", "red"),
        legend = c("Original Raw / 1e6", "Original Smoothed", "Simulated")
    )
    axis(3, at = top_ticks, labels = top_labels, cex.axis = 0.75)
    mtext(sprintf("Frequency [Hz]"), side = 3, line = 2)
    mtext(line1_text, side = 3, line = -1.1, col = "red", cex = 1, adj = 0.99)
    mtext(line2_text, side = 3, line = -2.1, col = "red", cex = 1, adj = 0.99)
    mtext(line3_text, side = 3, line = -3.1, col = "red", cex = 1, adj = 0.99)
}
