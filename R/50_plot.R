#' @noRd
#' @title Plot peaks of a spectrum
#' @description  This function plots the peaks of a spectrum, including the smoothed and scaled signal intensity and the second derivative. It also allows for the specification of peak positions and the option to draw vertical lines at these positions.#'
#' @param spec A data frame containing the spectrum data. It should have columns 'ppm', 'Y', 'ip', 'ip_left', 'ip_right', and 'd'.
#' @param ppm A vector of length 2 specifying the range of ppm values to consider for the plot. Default is c(3.402, 3.437).
#' @param dp A vector specifying the positions of the peaks. If NULL (default), the function will determine the peak positions based on the 'ppm' range.
#' @param vlines A logical value indicating whether to draw vertical lines at the peak positions. Default is FALSE.
#' @return A data frame with columns 'x' (ppm values), 'y' (smoothed and scaled signal intensity), 'd' (second derivative), 'is_ip' (whether the position is a peak), and 'is_ip_left' (whether the position is to the left of a peak).
#' @examples
#' \dontrun{
#' plot_peaks(spec, ppm = c(3.402, 3.437), dp = NULL, vlines = FALSE) # region from 3.402 to 3.437 ppm
#' pdf("spec_3.400_3.500.pdf", width = 24, height=8); plot_peaks(spec, c(3.400, 3.500), vlines = FALSE);  dev.off()
#' plot_peaks(spec, dp = 1:200, vlines = FALSE) # first 200 data points
#' pdf("spec_n1_n500.pdf", width = 24, height=8); plot_peaks(spec, dp = 1:500, vlines = FALSE);  dev.off()
#' }
plot_peaks <- function(spec, ppm = c(3.402, 3.437), dp = NULL, vlines = FALSE) {
    if (is.null(dp)) dp <- which(spec$ppm > min(ppm) & spec$ppm < max(ppm))
    x <- spec$ppm[dp]
    y <- spec$y_smooth[dp]
    l <- which(dp %in% spec$left) %||% numeric()
    p <- which(dp %in% spec$peak)
    r <- which(dp %in% spec$right) %||% numeric()
    m <- which(!(dp %in% c(spec$left, spec$peak, spec$right)))
    d <- spec$d[dp]
    withr::with_par(list(mfrow = c(3, 1), mar = c(0, 6, 4, 2), las = 1), {
        plot_spectrum(spec, focus = c(min(x), max(x)))
        # Plot 2: x ~ y + peaks (focussed region)
        plot(x, y, type = "l", xlab = "ppm", ylab = "", xaxt = "n", xlim = c(max(x), min(x)))
        mtext("smoothed and scaled signal intensity", side = 2, line = 5, las = 0)
        points(x[m], y[m], type = "p", pch = 124) # vertical dash
        points(x[p], y[p], col = "red", pch = 17) # triangle
        points(x[l], y[l], col = "blue", pch = 0) # open square
        points(x[r], y[r], col = "blue", pch = 4) # x character
        if (vlines) {
            abline(v = x[p], col = "red")
            abline(v = x[l], col = "blue")
            abline(v = x[r], col = "blue")
        }
        for (i in seq_along(l)) {
            rect(
                x[l[i]], par("usr")[3], x[r[i]], par("usr")[4],
                col = rgb(0, 0, 0, alpha = 0.1),
                border = rgb(0, 0, 0, alpha = 0.2)
            )
        }
        axis(3, at = x, labels = dp)
        legend("topright", legend = c("peak", "left", "right", "other"), col = c("red", "blue", "blue", "black"), pch = c(2, 0, 4, 124))
        # Plot 3: x ~ d + peaks
        withr::with_par(list(mar = c(5, 6, 0, 2)), {
            plot(x, d, type = "l", xlab = "ppm", ylab = "", xlim = c(max(x), min(x)))
            mtext("second derivative", side = 2, line = 5, las = 0)
            points(x[m], d[m], type = "p", pch = "|")
            points(x[p], d[p], col = "red", pch = 17)
            points(x[l], d[l], col = "blue", pch = 0)
            points(x[r], d[r], col = "blue", pch = 4)
            if (vlines) {
                abline(v = x[p], col = "red")
                abline(v = x[l], col = "blue")
                abline(v = x[r], col = "blue")
            }
            abline(h = 0, col = "black")
        })
    })
    df <- data.frame(x = x, y = y, d = d, is_ip = dp %in% spec$peak, is_ip_left = dp %in% spec$ip_left)
    invisible(df)
}


#' @title Plot Spectrum
#' @description Plot a spectrum based on the provided deconvolution data and highlight a specific region of interest in the spectrum.
#' @param decon An object as returned by [generate_lorentz_curves()], containing the deconvolution data. Must include either `x_values_ppm` or `ppm` for the x-axis values, and either `y_values` or `y_smooth` for the y-axis values.
#' @param focus A numeric vector of length 2 specifying the region of interest to highlight on the plot. The region is defined by its start and end points on the x-axis (in ppm).
#' @return A plot is generated as a side effect, highlighting the specified focus region on the spectrum.
#' @examples
#' sim_01 <- system.file("example_datasets/bruker/sim/sim_01", package = "metabodecon")
#' decon <- generate_lorentz_curves(
#'     sim_01, sfr = c(3.42, 3.58), ws = 0, ask = FALSE,
#'     smopts = c(1, 5), delta = 0.1
#' )
#' plot_spectrum(decon)
plot_spectrum <- function(decon, focus = c(3.45, 3.55)) {
    x <- decon$x_values_ppm %||% decon$ppm
    y <- decon$y_values %||% decon$y_smooth
    plot(x, y, type = "l", xlab = "ppm", ylab = "", xlim = c(max(x), min(x)))
    rect(
        min(focus), par("usr")[3], max(focus), par("usr")[4],
        col = rgb(0, 0, 0, alpha = 0.1),
        border = rgb(0, 0, 0, alpha = 0.2)
    )
}

#' @noRd
#' @title Plot Signal Free Region
#' @description Draws the signal free region as green vertical lines into the given spectrum.
#' @param spec A list representing the spectrum as returned by [read_spectrum()] or [load_bruker_spectrum()].
#' @param left_ppm The left border of the signal free region in ppm.
#' @param right_ppm The right border of the signal free region in ppm.
#' @return NULL. Called for side effect of plotting the signal free region.
plot_sfr <- function(spec, left_ppm, right_ppm) {
    plot(
        x = spec$ppm,
        y = spec$y_scaled,
        type = "l",
        xlab = "[ppm]",
        ylab = "Intensity [a.u.]",
        xlim = c(spec$ppm_max, spec$ppm_min)
    )
    graphics::abline(v = c(left_ppm, right_ppm), col = "green")
}

#' @noRd
#' @title Plot Water Signal
#' @description Draws the water signal as red vertical lines into the given spectrum.
#' @param spec A list representing the spec as returned by [read_spectrum()].
#' @param hwidth_ppm The half width of the water signal in ppm.
#' @return NULL. Called for side effect of plotting the water signal.
plot_ws <- function(spec, hwidth_ppm) {
    center_ppm <- (spec$ppm_max + spec$ppm_min) / 2
    plot(
        spec$ppm,
        spec$y_scaled,
        type = "l",
        xlab = "[ppm]",
        ylab = "Intensity [a.u.]",
        xlim = c(center_ppm + 2 * hwidth_ppm, center_ppm - 2 * hwidth_ppm)
    )
    graphics::abline(v = center_ppm + hwidth_ppm, col = "red")
    graphics::abline(v = center_ppm - hwidth_ppm, col = "red")
}
