# Internal #####

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
#' pdf("spec_3.400_3.500.pdf", width = 24, height = 8)
#' plot_peaks(spec, c(3.400, 3.500), vlines = FALSE)
#' dev.off()
#' plot_peaks(spec, dp = 1:200, vlines = FALSE) # first 200 data points
#' pdf("spec_n1_n500.pdf", width = 24, height = 8)
#' plot_peaks(spec, dp = 1:500, vlines = FALSE)
#' dev.off()
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
        xlab = "Chemical Shift [ppm]",
        ylab = "Signal Intensity [au]",
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
    hw <- hwidth_ppm
    x0 <- (spec$ppm_max + spec$ppm_min) / 2
    xlim <- if (hw > 0) c(x0 + 5 * hw, x0 - 5 * hw) else range(spec$ppm)[2:1]
    plot(
        spec$ppm,
        spec$y_scaled,
        type = "l",
        xlab = "Chemical Shift [ppm]",
        ylab = "Signal Intensity [au]",
        xlim = xlim
    )
    graphics::abline(v = x0, col = "red", lty = 1) # center
    mtext("center", side = 3, line = 0, at = x0, col = "red")
    if (hw > 0) {
        ya <- max(spec$y_scaled) * 0.8
        abline(v = c(x0 + hw, x0 - hw), col = "red") # left/right border
        arrows(x0, ya, x0 + c(hw, -hw), ya, col = "red", length = 0.1)
        text(x0 + c(hw, -hw) / 2, ya, labels = "wshw", pos = 3, col = "red")
    }
}


# Helpers #####

#' @description Returns TRUE if the current multi-figure gets filled by row, else FALSE.
#' @examples
#' callwith <- function(by_row = TRUE) {
#'     grid <- c(2, 2)
#'     opar <- if (by_row) par(mfrow = grid) else par(mfcol = grid)
#'     on.exit(opar)
#'     plot(1:10)
#'     by_row <- mf_filled_by_row()
#'     plot(1:5)
#'     by_row
#' }
#' x <- callwith(by_row = TRUE)
#' y <- callwith(by_row = FALSE)
#' stopifnot(isTRUE(x) && isFALSE(y))
mf_filled_by_row <- function() {
    mfg <- par("mfg")
    row <- mfg[1]
    col <- mfg[2]
    nrows <- mfg[3]
    ncols <- mfg[4]
    if (nrows == 1 || ncols == 1) {
        # In this case it doesn't matter, as the figure spots in the grid will be filled top-to-bottom / left-to-right in both cases.
        return(TRUE)
    } else {
        # We have at least two rows AND two cols. So what we can do is to set c(1, 1) as next figure, and then advance one frame. If we end up in c(1, 2) we are row-oriented, if we end up in c(2, 1) we are column-oriented. After doing this we can reset the current figure number to the original value.
        par(mfg = c(1, 1, nrows, ncols))
        on.exit({
            par(mfg = mfg)
            plot_empty() # When querying `mfg` we get the "current figure number". But when setting it, we set the "next figure number". I.e. we need to advance one frame, or we would set the "current figure number" as "next figure number".
        })
        plot_empty() # Draw into c(1, 1)
        plot_empty() # Draw into c(1, 2) or c(2, 1)
        mfg2 <- par("mfg") # Query current position
        return(if (mfg2[1] == 1) TRUE else FALSE)
    }
}

plot_empty <- function() {
    plot(
        x = 0.5, y = 0.5, type = "n",
        xlim = c(0, 1), xlab = "", xaxs = "i",
        ylim = c(0, 1), ylab = "", yaxs = "i",
        axes = FALSE, main = "", ann = FALSE
    )
}

#' @noRd
#' @title Set Figure Region
#' @description Calling `par(fig=xxyy)` resets the current multi-figure configuration (MFC) to one row and one column. `set_fig()` handles this scenario, by first storing the current MFC, then calling `par(fig=xxyy)` and finally returning a function that can be used to restore the MFC. See 'Details' for further information.
#' @param fig Region to draw into, given as normalized device coordinates.
#' @param add If TRUE, the new plot is added to the existing plot.
#' @return Function to reset the MFC.
#' @details
#' Note 1: Setting `par(fig=xxyy)` resets the current MFC to one row and one column. I.e., array layouts defined by setting `mfrow` or `mfcol`, must be saved before changing `fig` and restored after the plot has been drawn. The same applies to the current figure number `mfg` (see Note 2).
#'
#' Note 2: When restoring `mfg` it's important to additionally advance one frame, because when querying `mfg`, the "current figure number" is returned, but when setting `mfg`, the value is interpreted as "next figure number". I.e. without advancing one frame, the next figure would be drawn into the spot of the current figure.
#'
#' Note 3: If `fig=xxyy` and `add=FALSE`, it is still necessary to use `par(fig=xxyy, new=TRUE)` to prevent the device from being cleared (which would be bad in a multi-figure environment). I.e., the figure number must be advanced manually, as would have been done by a normal plot. This manual increment must be done before the MFC is reset by `par(fig=xxyy, new=TRUE)`.
#' @examples
#' plot_dummy <- function() {
#'     plot(0, 0, ylim = c(0, 1), xlim = c(0, 1), xaxs = "i", yaxs = "i")
#'     text(0.5, 0.5, "dummy")
#' }
#' p <- local({
#'     opar <- par(mfrow = c(2, 2), mar = c(2, 2, 0.5, 0.5))
#'     on.exit(par(opar))
#'
#'     topleft <- plot_dummy()
#'     reset_mfc <- set_fig(fig = c(0.25, 0.50, 0.50, 0.75), add = TRUE)
#'     topleft2 <- plot_dummy()
#'     reset_mfc()
#'
#'     topright <- plot_dummy()
#'
#'     reset_mfc <- set_fig(fig = c(0.1, 0.4, 0.1, 0.4), add = FALSE)
#'     bottom_left <- plot_dummy()
#'     reset_mfc()
#'
#'     bottom_right <- plot_dummy()
#' })
set_fig <- function(fig = NULL, add = TRUE) {
    if (is.null(fig)) {
        return(function() {})
    } # Nothing to do if region is NULL
    if (isFALSE(add)) plot_empty() # Advance one frame if `add=FALSE` (Note 3)
    op <- par(c("mar", "mfrow", "mfcol", "mfg", "fig")) # Store MF conf (Note 1)
    byrow <- mf_filled_by_row() # Store MF conf (Note 1)
    par(fig = fig, new = TRUE) # Set new figure region (Note 3)
    reset_mfc <- function() {
        if (byrow) par(mfrow = op$mfrow) else par(mfcol = op$mfcol) # Restore MF Layout (Note 1)
        par(mfg = op$mfg) # Restore current figure number (Note 1)
        plot_empty() # Advance one frame (Note 2)
    }
    reset_mfc
}

#' @noRd
#' @title Plot into specific figure region
#' @description For Details see [set_fig()].
#' @examples
#' plot_dummy <- function() {
#'     plot(0, 0, ylim = c(0, 1), xlim = c(0, 1), xaxs = "i", yaxs = "i")
#'     text(0.5, 0.5, "dummy")
#' }
#' p <- local({
#'     opar <- par(mfrow = c(2, 2), mar = c(2, 2, 0.5, 0.5))
#'     on.exit(par(opar))
#'     topleft <- plot_dummy()
#'     topleft2 <- with_fig(fig = c(0.25, 0.50, 0.50, 0.75), add = TRUE, plot_dummy())
#'     topright <- plot_dummy()
#'     bottom_left <- with_fig(fig = c(0.1, 0.4, 0.1, 0.4), add = FALSE, plot_dummy())
#'     bottom_right <- plot_dummy()
#' })
with_fig <- function(expr, fig = NULL, add = TRUE) {
    reset_mfc <- set_fig(fig = fig, add = add)
    on.exit(reset_mfc())
    expr
}

plot_dummy <- function() {
    plot(
        x = 0, y = 0, main = "dummy main",
        ylim = c(0, 1), xlim = c(0, 1),
        xaxs = "i", yaxs = "i",
        xlab = "dummy xlab", ylab = "dummy ylab"
    )
    text(0.5, 0.5, "dummy text")
}
