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


#' @export
#' @title Plot Spectrum
#' @description Plot a spectrum based on the provided deconvolution data and zoom in on a specific region of interest in the spectrum.
#'
#' `r lifecycle::badge("experimental")`
#' @param focus A numeric vector of length 2 specifying the region of interest to zoom in on the plot. The region is defined by its start and end points on the x-axis (in ppm).
#' @param decon An object as returned by [generate_lorentz_curves()], containing the deconvolution data. Must include either `x_values_ppm` or `ppm` for the x-axis values, and either `y_values` or `y_smooth` for the y-axis values.
#' @param foc_rgn Numeric vector specifying the start and end of focus region.
#' @param foc_unit Character string specifying the unit in which `foc_rgn` is given. Can be "fraction" or "ppm".
#' @param mar Number of lines below/left/above/right plot region (used for axis annotations).
#' @param line_col Color of the spectrum line.
#' @param axis_col Color of tickmarks and ticklabels.
#' @param rect_fill_col Backghround color of the rectangle around the focus region.
#' @param rect_border_col Border color of the rectangle around the focus region.
#' @param sub_fig Either NULL or a numeric vector of the form `c(x1, x2, y1, y2)` giving the left/right/bottom/top coordinates of the sub figure region in "normalized plot coordinates" (as described in [graphics::grconvertX()]). If provided, the focussed region is drawn as sub figure within the main plot region. Setting `sub_fig` to NULL will prevent the sub figure from being drawn. If used, active multi-figure configurations (as defined by setting `mfrow` or `mfcol`) will be reset. See 'Details' for more information.
#' @param sub_mar Margins of the sub figure.
#' @param sub_line_col Color of the lines in the sub figure.
#' @param sub_axis_col Color of the axis in the sub figure.
#' @param sub_fill_col Background color of the sub figure.
#' @param sub_box_col Border color of the box surrounding the sub figure.
#' @param connect_line_col Color of the lines connecting the main and sub figure.
#' @param ysquash Fraction of plot height to squash y-values into. Useful in combination with `sub_fig` to prevent the spectrum lines from overlapping with the sub figure showing the focussed region.
#' @details Setting `sub_fig` to a value unequal NULL resets a potential multi-figure configuration defined via `mfrow` or `mfcol`. In such cases you need to restore the array config (`mfrow`/`mfcol`) as well as the current figure number (`mfg`) manually yourself. E.g. if you have a 2 x 2 grid, defined via `par(mfrow = c(2, 2))` and draw into position (1, 2), then you need to call `par(mfrow = c(2, 2), mfg = c(1, 2))` after calling `.plot_spectrum(...)`, to draw into position (2, 1) in the next plot.
#' @return NULL. Called for side effect of plotting the spectrum, as sketched below.
#'
#' ```
#' 0__0.1_____________0.5_____________0.9__1
#' |    _______________________________    |0.9
#' |   |Focus Region  _/\_   _         |   |
#' |   |       _    _/    \_/ \__      |   |
#' |   |    __/ \__/             \_    |   |
#' |   |___/_______________________\___|   |0.5
#' |                    _______________    |
#' | Full /\           |Focus /\       |   |
#' | Spectrum          |Region  \/\    |   |
#' |    /    \     /\  |  /\/      \   |   |
#' |___/______\___/__\_|_/__________\__|___|0
#' ```
#'
#' @examples
#'
#' # Prepare a deconvoluted spectrum as input
#' sim_01 <- metabodecon_file("sim/sim_01")
#' decon <- generate_lorentz_curves(
#'     sim_01,
#'     sfr = c(3.42, 3.58), ws = 0,
#'     ask = FALSE, verbose = FALSE
#' )
#'
#' # Plot the full spectrum
#' plot_spectrum(decon, foc_rgn = NULL)
#'
#' # Focus on a specific region, specified in ppm or as fraction
#' plot_spectrum(decon, foc_rgn = c(3.49, 3.45), foc_unit = "ppm")
#' plot_spectrum(decon, foc_rgn = c(0.4, 0.3), foc_unit = "fraction")
#'
#' # Change color and position of the focus region
#' plot_spectrum(decon,
#'     rect_fill_col = rgb(1, 0, 0, alpha = 0.1),
#'     rect_border_col = "red",
#'     sub_fig = c(x1 = 0.1, x2 = 0.9, y1 = 0.5, y2 = 0.9)
#' )
#'
#' # Remove connecting lines and fill colors
#' plot_spectrum(decon,
#'     rect_fill_col = NULL,
#'     sub_fill_col = NULL
#'     connect_line_col = NULL
#' )
plot_spectrum <- function(decon,
                          foc_rgn = c(0.40, 0.35),
                          foc_unit = "fraction",
                          mar = c(4.1, 4.1, 0.1, 0.1),
                          line_col = "black",
                          axis_col = "black",
                          rect_fill_col = rgb(0, 0, 1, alpha = 0.05),
                          rect_border_col = "blue",
                          sub_fig = c(x1 = 0.05, x2 = 0.95, y1 = 0.2, y2 = 0.95),
                          sub_mar = c(2.1, 2.1, 0.1, 0.1),
                          sub_line_col = line_col,
                          sub_axis_col = axis_col,
                          sub_fill_col = rect_fill_col,
                          sub_box_col = rect_border_col,
                          connect_line_col = rect_border_col,
                          ysquash = if (!is.null(foc_rgn) && !is.null(sub_fig)) sub_fig[3] * 0.96 else 0.96
                          ) {

    # Check args.
    foc_unit <- match.arg(foc_unit, c("ppm", "fraction"))
    cs <- decon[["ppm"]] %||% decon[["x_values_ppm"]] %||% decon[["cs"]]
    si <- decon[["y_smooth"]] %||% decon[["y_values"]] %||% decon[["si"]]
    if (foc_unit == "fraction" && !is.null(foc_rgn)) foc_rgn <- quantile(cs, foc_rgn)
    if (identical(sort(names(sub_fig)), c("x1", "x2", "y1", "y2"))) {
        sub_fig <- sub_fig[order(names(sub_fig))]
    } else if (is.null(sub_fig)) {} else {
        stop("sub_fig must be NULL or a numeric vector with names 'x1', 'x2', 'y1', 'y2'")
    }

    # Draw full spectrum if sub_fig is NULL, else focussed region.
    main <- .plot_spectrum(
        cs, si, foc_rgn,
        foc_only = if (is.null(sub_fig)) TRUE else FALSE,
        mar = mar,
        line_col = line_col,
        box_col = "black",
        axis_col = axis_col,
        fill_col = NULL,
        rect_fill_col = rect_fill_col,
        rect_border_col = rect_border_col,
        ysquash = ysquash
    )

    # Draw focussed region if sub_fig is defined.
    sub <- NULL
    if (!is.null(sub_fig) && !is.null(foc_rgn)) {
        xfig <- convert_pos(sub_fig[1:2], c(0, 1), main$par$plt[1:2])
        yfig <- convert_pos(sub_fig[3:4], c(0, 1), main$par$plt[3:4])
        xndc <- convert_pos(xfig, c(0, 1), main$par$fig[1:2])
        yndc <- convert_pos(yfig, c(0, 1), main$par$fig[3:4])
        sub <- .plot_spectrum(
            cs, si, foc_rgn,
            foc_only = TRUE,
            mar = sub_mar,
            line_col = sub_line_col,
            box_col = sub_box_col,
            axis_col = sub_axis_col,
            fill_col = sub_fill_col,
            fig = c(xndc, yndc),
            add = TRUE
        )
    }

    # Draw connecting lines between main and sub figure.
    if (!is.null(sub_fig) && !is.null(foc_rgn)) {
        line1_x0 <- grconvertX(main$rct$xleft_ndc, "ndc", "user")
        line2_x0 <- grconvertX(main$rct$xright_ndc, "ndc", "user")
        line1_x1 <- grconvertX(sub$plt$xlim_ndc[1], "ndc", "user")
        line2_x1 <- grconvertX(sub$plt$xlim_ndc[2], "ndc", "user")
        y0 <- grconvertY(main$rct$ytop_ndc, "ndc", "user")
        y1 <- grconvertY(sub$plt$ylim_ndc[1], "ndc", "user")
        opar <- par(xpd = TRUE)
        on.exit(par(opar), add = TRUE)
        segments(x0 = line1_x0, y0 = y0, x1 = line1_x1, y1 = y1, col = connect_line_col)
        segments(x0 = line2_x0, y0 = y0, x1 = line2_x1, y1 = y1, col = connect_line_col)
    }

    # Return plot parameters
    invisible(named(main, sub))
}

#' @noRd
#' @title Plot Spectrum Internal
#' @description Plot a spectrum based on the provided chemical shift and signal intensity data. Should not be called directly by the user. Use [plot_spectrum()] instead.
#' @param cs Chemical shift of spectrum in ppm.
#' @param si Signal intensity of spectrum in arbitrary units.
#' @param foc_rgn Focus region in ppm.
#' @param foc_only If TRUE, only the focussed region is drawn. If FALSE, the full spectrum is drawn.
#' @param mar Number of lines below/left/above/right plot region.
#' @param line_col Color of spectrum line.
#' @param box_col Color of box around plot region.
#' @param axis_col Color of tickmarks and ticklabels.
#' @param fill_col Background color of plot region.
#' @param rect_fill_col Fill color of rectangle around foc_rgn region.
#' @param rect_border_col Border color of rectangle around foc_rgn region.
#' @param ysquash Fraction of plot height to squash y-values into.
#' @param fig Region to draw into, given as normalized device coordinates. Doesn't work for multi-figure configurations defined by setting `mfrow` or `mfcol`. See 'Details' for more information.
#' @param add If TRUE, the new plot is added to the existing plot.
#' @details Setting `fig` to a value unequal NULL resets a potential multi-figure configuration defined via `mfrow` or `mfcol`. In such cases you need to restore the array config (`mfrow`/`mfcol`) as well as the current figure number (`mfg`) manually yourself. E.g. if you have a 2 x 2 grid, defined via `par(mfrow = c(2, 2))` and draw into position (1, 2), then you need to call `par(mfrow = c(2, 2), mfg = c(1, 2))` after calling `.plot_spectrum(...)`, to draw into position (2, 1) in the next plot.
#' @examples
#' sim_01 <- metabodecon_file("sim/sim_01")
#' decon <- generate_lorentz_curves(
#'     sim_01,
#'     sfr = c(3.42, 3.58), ws = 0,
#'     ask = FALSE, verbose = FALSE
#' )
#' cs <- decon[["x_values_ppm"]]
#' si <- decon[["y_values"]]
#' po <- .plot_spectrum(cs, si)
.plot_spectrum <- function(cs,
                           si,
                           foc_rgn = NULL,
                           foc_only = FALSE,
                           mar = c(4.1, 4.1, 0.1, 0.1),
                           line_col = "black",
                           box_col = "black",
                           axis_col = "black",
                           fill_col = NULL,
                           rect_fill_col = rgb(1.0, 0.0, 0.0, alpha = 0.1),
                           rect_border_col = "red",
                           ysquash = 0.96,
                           fig = NULL,
                           add = FALSE
                           ) {

    if (add && !is.null(fig)) {
        opar <- par(mar = mar, new = add, fig = fig)
    } else {
        opar <- par(mar = mar, new = add)
    }
    on.exit(add = TRUE, par(opar))

    has_foc_rgn <- !is.null(foc_rgn)
    if (has_foc_rgn && length(foc_rgn) != 2) stop("foc_rgn must be a numeric vector of length 2")
    if (foc_only && !has_foc_rgn) stop("foc_only requires foc_rgn to be specified")
    in_focus <- if (has_foc_rgn) (cs >= min(foc_rgn) & cs <= max(foc_rgn)) else FALSE
    si_foc <- si[in_focus]
    cs_foc <- cs[in_focus]

    # Prepare plot region
    plt <- list(
        x = x <- if (foc_only) cs_foc else cs,
        y = y <- if (foc_only) si_foc else si,
        type = "n",
        axes = FALSE,
        xlab = "Chemical Shift [ppm]",
        xlim = c(max(x), min(x)),
        xaxs = "i",
        ylab = "Signal Intensity [au]",
        ylim = c(0, max(y) / ysquash),
        yaxs = "i"
    )
    do.call(plot, plt)

    # Color background of plotting region
    if (!is.null(fill_col)) rect(
        xleft = plt$xlim[1],
        ybottom = plt$ylim[1],
        xright = plt$xlim[2],
        ytop = plt$ylim[2],
        col = fill_col,
        border = NA
    )

    # Draw spectrum as line
    lines(plt$x, plt$y, type = "l", col = line_col)
    xticks <- seq(min(plt$x), max(plt$x), length.out = 5)
    yticks <- seq(0,          max(plt$y), length.out = 5)
    axis(1, at = xticks, labels = round(xticks, 2), col = NA, col.ticks = axis_col, col.axis = axis_col)
    axis(2, at = yticks, labels = round(yticks, 2), col = NA, col.ticks = axis_col, col.axis = axis_col, las = 1)
    box(col = box_col)

    # Draw rectangle around focus region
    rct <- NULL
    if (has_foc_rgn && !foc_only) {
        rct <- list(
            xleft = max(foc_rgn),
            xright = min(foc_rgn),
            ybottom = 0,
            ytop = max(si_foc) / 0.96,
            col = rect_fill_col,
            border = rect_border_col
        )
        do.call(rect, rct)
    }

    # Return plot parameters
    plt$xlim_ndc <- grconvertX(plt$xlim, from = "user", to = "ndc")
    plt$ylim_ndc <- grconvertY(plt$ylim, from = "user", to = "ndc")
    if (!is.null(rct)) {
        rct$xleft_ndc <- grconvertX(rct$xleft, from = "user", to = "ndc")
        rct$xright_ndc <- grconvertX(rct$xright, from = "user", to = "ndc")
        rct$ybottom_ndc <- grconvertY(rct$ybottom, from = "user", to = "ndc")
        rct$ytop_ndc <- grconvertY(rct$ytop, from = "user", to = "ndc")
    }
    invisible(list(plt = plt, rct = rct, par = par()))
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
