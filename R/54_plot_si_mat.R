
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
