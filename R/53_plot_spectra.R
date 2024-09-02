plot_spectra <- function(ss = glc_sim(),
                         mar = c(4.1, 4.1, 1.1, 0.1),
                         peak_rng = get_ppm_range(ss, show = FALSE)) {
    if (is_decon_obj(ss)) ss <- list(ss)
    if (is_decon_list(ss)) {
        xrng <- range(c(sapply(ss, function(s) s$x_values_ppm)))
        ymax <- max(sapply(ss, function(s) max(s$y_values)))
    } else if (is.data.frame(ss)) {
        stop("ss must be a list of deconvoluted spectra.")
    }
    a <- peak_rng[1]
    b <- peak_rng[2]
    w <- (b - a) / 4
    y8 <- ymax * 0.8
    cols <- rainbow(length(ss))
    ltxt <- paste("Spectrum", 1:length(ss))
    opar <- par(mar = mar)
    on.exit(par(opar))
    plot(x = NA, type = "n", xlim = xrng[2:1], ylim = c(0, ymax), xlab = "Chemical Shift [ppm]", ylab = "Signal Intensity [au]")
    abline(v = peak_rng, lty = 2)
    for (i in 1:length(ss)) lines(x = ss[[i]]$x_values_ppm, y = ss[[i]]$y_values, col = cols[i])
    arrows(x0 = c(a + w, b - w), x1 = c(a, b), y0 = y8, y1 = y8, length = 0.1, lty = 2, col = "black")
    text(x = mean(c(a, b)), y = y8, labels = "ppm range")
    mtext(text = round(c(a, b), 4), side = 3, line = 0, at = c(a, b))
    legend(x = "topright", legend = ltxt, col = cols, lty = 1)
}
