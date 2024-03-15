#' Plot peaks of a spectrum
#'
#' This function plots the peaks of a spectrum, including the smoothed and scaled signal intensity and the second derivative.
#' It also allows for the specification of peak positions and the option to draw vertical lines at these positions.
#'
#' @param spec A data frame containing the spectrum data. It should have columns 'ppm', 'Y', 'ip', 'ip_left', 'ip_right', and 'd'.
#' @param ppm A vector of length 2 specifying the range of ppm values to consider for the plot. Default is c(3.402, 3.437).
#' @param dp A vector specifying the positions of the peaks. If NULL (default), the function will determine the peak positions based on the 'ppm' range.
#' @param vlines A logical value indicating whether to draw vertical lines at the peak positions. Default is FALSE.
#' @return A data frame with columns 'x' (ppm values), 'y' (smoothed and scaled signal intensity), 'd' (second derivative), 'is_ip' (whether the position is a peak), and 'is_ip_left' (whether the position is to the left of a peak).
#' @examples \dontrun{
#' plot_peaks(spec, ppm = c(3.402, 3.437), dp = NULL, vlines = FALSE) # region from 3.402 to 3.437 ppm
#' pdf("spec_3.400_3.500.pdf", width = 24, height=8); plot_peaks(spec, c(3.400, 3.500), vlines = FALSE);  dev.off()
#' plot_peaks(spec, dp = 1:200, vlines = FALSE) # first 200 data points
#' pdf("spec_n1_n500.pdf", width = 24, height=8); plot_peaks(spec, dp = 1:500, vlines = FALSE);  dev.off()
#' }
#' @noRd
plot_peaks <- function(spec, ppm = c(3.402, 3.437), dp = NULL, vlines = FALSE) {
    if (is.null(dp)) dp <- which(spec$ppm > min(ppm) & spec$ppm < max(ppm))
    x <- spec$ppm[dp]
    y <- spec$Y$smooth[dp]
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

plot_spectrum <- function(spec, focus) {
    x <- spec$ppm
    y <- spec$Y$smooth
    plot(x, y, type = "l", xlab = "ppm", ylab = "", xlim = c(max(x), min(x)))
    rect(
        min(focus), par("usr")[3], max(focus), par("usr")[4],
        col = rgb(0, 0, 0, alpha = 0.1),
        border = rgb(0, 0, 0, alpha = 0.2)
    )
}
