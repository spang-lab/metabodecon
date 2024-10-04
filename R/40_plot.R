# Public API #####

#' @export
#'
#' @title Plot Spectrum
#'
#' @description
#' Plot a spectrum and zoom in on a specific region.
#'
#' `r lifecycle::badge("experimental")`
#'
#' @param decon An object as returned by [generate_lorentz_curves()], containing
#' the deconvolution data. Must include either `x_values_ppm` or `ppm` for the
#' x-axis values, and either `y_values` or `y_smooth` for the y-axis values.
#'
#' @param foc_rgn Numeric vector specifying the start and end of focus region.
#'
#' @param foc_unit Character string specifying the unit in which `foc_rgn` is
#' given. Can be "fraction" or "ppm".
#'
#' @param foc_only Logical. If TRUE, only the focused region is drawn. If FALSE,
#' the full spectrum is drawn.
#'
#' @param sub_show Logical. If TRUE, a sub figure is drawn within the main plot
#' region. If FALSE, the sub figure is not drawn.
#'
#' @param sub_rgn Either NULL or a numeric vector of the form `c(x1, x2, y1,
#' y2)` giving the left/right/bottom/top coordinates of the sub figure region in
#' "normalized plot coordinates" (as described in [graphics::grconvertX()]). If
#' provided, the focused region is drawn as sub figure within the main plot
#' region. Setting `sub_rgn` to NULL will prevent the sub figure from being
#' drawn.
#'
#' @param verbose Logical. If TRUE, print additional information.
#'
#' @param rct_show Logical. If TRUE, the focus region is shown as a rectangle.
#'
#' @param rct_fill Background color of the rectangle around the focus region.
#'
#' @param rct_col Border color of the rectangle around the focus region.
#'
#' @param main Title of the plot.
#'
#' @param xlab Label for the x-axis.
#'
#' @param ylab Label for the y-axis.
#'
#' @param mar Number of lines below/left/above/right plot region (used for axis
#' annotations).
#'
#' @param line_col Color of raw signal intensities.
#'
#' @param axis_col Color of tickmarks and ticklabels.
#'
#' @param fill_col Background color of the plot region.
#'
#' @param box_col Border color of the box surrounding the plot region.
#'
#' @param ysquash Fraction of plot height to squash y-values into. Useful in
#' combination with `sub_rgn` to prevent the spectrum lines from overlapping
#' with the sub figure showing the focused region.
#'
#' @param trp_show Logical. If TRUE, the peak triplets are shown.
#'
#' @param lc_show Logical. If TRUE, the Lorentzian Curves are shown.
#'
#' @param sup_show Logical. If TRUE, the superposition of the Lorentzian Curves
#' is shown.
#'
#' @param sub_mar Margins of the sub figure.
#'
#' @param sub_line_col Color of the lines in the sub figure.
#'
#' @param sub_axis_col Color of the axis in the sub figure.
#'
#' @param sub_fill_col Background color of the sub figure.
#'
#' @param sub_box_col Border color of the box surrounding the sub figure.
#'
#' @param sub_trp_show Logical. If TRUE, the peak triplets are shown in the sub
#' figure.
#'
#' @param sub_lc_show Logical. If TRUE, the Lorentzian Curves are shown in the
#' sub figure.
#'
#' @param sub_sup_show Logical. If TRUE, the superposition of the Lorentzian
#' Curves is shown in the sub figure.
#'
#' @param sf_y_raw Numeric value to divide raw signal intensities by before
#' drawing.
#'
#' @param sm_col Color of smoothed signal intensities.
#'
#' @param lc_col Color of the Lorentzian Curves.
#'
#' @param lc_lty Line type of the Lorentzian Curves.
#'
#' @param lc_fill Color of the rectangles shown at center of each lorentzian
#' curve. The width of the rectangle equals the half width at half height.
#'
#' @param trp_col Vector of length 4 giving the colors for the peak triplets.
#' The first three colors specify the color used for each peak-center,
#' left-border and right-border and the fourth color used for any non-peak data
#' point.
#'
#' @param trp_pch Vector of length 4 giving the plotting characters for the peak
#' triplets.
#'
#' @param sup_col Color of the superposition of the Lorentzian Curves.
#'
#' @param sup_lty Line type of the superposition of the Lorentzian Curves.
#'
#' @param cnct_show Logical. If TRUE, connecting lines are drawn between the
#' main and sub figure.
#'
#' @param cnct_col Color of the lines connecting the main and sub figure.
#'
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
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' ## Prepare a deconvoluted spectrum as input
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' sim_01 <- metabodecon_file("sim/sim_01")
#' spec <- read_spectrum(sim_01)
#' decon <- generate_lorentz_curves_sim(spec)
#' opar <- par(mfrow = c(3, 2))
#' on.exit(par(opar))
#'
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' ## 1. Plot the full (non-deconvoluted) spectrum
#' ## 2. Focus on a specific region, specified in ppm
#' ## 3. Focus on a specific region, specified in as fraction of the full width
#' ## 4. Change margin, color and position of the focus region
#' ## 5. Remove connecting lines and fill colors
#' ## 6. Hide xlab and ylab
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' plot_spectrum(spec)
#' plot_spectrum(decon, foc_rgn = c(3.49, 3.45), foc_unit = "ppm")
#' plot_spectrum(decon, foc_rgn = c(0.40, 0.30))
#' plot_spectrum(decon,
#'     sub_mar = c(4, 4, 0, 0),
#'     rct_fill = rgb(0.9, 0.5, 0.9, alpha = 0.1),
#'     rct_col = "violet",
#'     sub_rgn = c(x1 = 0.1, x2 = 0.9, y1 = 0.4, y2 = 0.9)
#' )
#' plot_spectrum(decon, rct_fill = NULL, sub_fill_col = NULL, cnct_col = NULL)
#' plot_spectrum(decon, xlab = "", ylab = "", mar = c(2, 2, 0, 1))
plot_spectrum <- function(
    # Common Settings
    decon,
    main = "",
    foc_rgn = c(0.40, 0.35),
    foc_unit = "fraction",
    foc_only = FALSE,
    sub_show = if (foc_only) FALSE else TRUE,
    sub_rgn = c(x1 = 0.05, x2 = 0.95, y1 = 0.2, y2 = 0.95),
    verbose = FALSE,
    # Focus Rectangle
    rct_show = sub_show,
    rct_col = "black",
    rct_fill = transp("yellow"),
    # Settings for Main Figure
    xlab = "Chemical Shift [ppm]",
    ylab = "Signal Intensity [au]",
    mar = c(4, 4, if (main == "") 0 else 4, 1),
    line_col = "black",
    axis_col = "black",
    fill_col = NULL,
    box_col = "black",
    ysquash = if (sub_show) sub_rgn[3] * 0.96 else 0.96,
    trp_show = if (sub_show) FALSE else TRUE,
    lc_show = if (sub_show) FALSE else TRUE,
    sup_show = if (sub_show) FALSE else TRUE,
    # Settings for Sub Figure
    sub_mar = c(2, 2, 0, 0),
    sub_line_col = line_col,
    sub_axis_col = axis_col,
    sub_fill_col = rct_fill,
    sub_box_col = rct_col,
    sub_trp_show = TRUE,
    sub_lc_show = TRUE,
    sub_sup_show = TRUE,
    # Settings for both Figures
    sf_y_raw = 1e6,
    sm_col = "blue",
    lc_col = "darkgrey",
    lc_lty = 1,
    lc_fill = transp(lc_col, 0.25),
    trp_col = rep(sm_col, 4),
    trp_pch = c(17, 4, 4, NA),
    sup_col = "red",
    sup_lty = 1,
    # Connecting Lines
    cnct_show = if (sub_show) TRUE else FALSE,
    cnct_col = rct_col
) {
    main_args <- ps_get_main_args()
    main_plot <- do.call(ps_internal, main_args)
    sub_args <- if (sub_show) ps_get_sub_args(main_args, main_plot)
    sub_plot <- if (sub_show) do.call(ps_internal, sub_args)
    cnct_plot <- if (cnct_show) ps_connect_main_sub(main_plot, sub_plot, cnct_col)
    invisible(named(main_args, main_plot, sub_args, sub_plot, cnct_plot))
}

# Internal #####

#' @noRd
#' @title Plot Signal Free Region
#' @description Draws the SFR as green vertical lines into the given spectrum.
#' @param spec The spectrum object.
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

#' @noRd
#' @title Plot Aligned Spectra
#' @description Plots the aligned and unaligned spectra for comparison.
#' @param YA,YB Matrix.
#' `YA[i,j]` == Signal Intensity of spectrum i at index j AFTER alignmnent.
#' `YB[i,j]` == Signal Intensity of spectrum i at index j BEFORE alignment.
#' @param PA,PB List of vectors.
#' `PA[[i]][j]` == index of peak j of spectrum i AFTER alignment.
#' `PB[[i]][j]` == index of peak j of spectrum i BEFORE alignment.
#' @param mfcol Vector of two integers specifying the number of rows and columns
#' of the plot grid
#' @examples
#' x <- seq(1.5 * pi, 9.5 * pi, length.out = 90)
#' y <- 10 * sin(x) # y without noise
#' p <- sort(order(y, decreasing = TRUE)[1:4]) # peaks without noise
#' Y <- lapply(1:4, function(i) smooth(smooth(y + rnorm(90)))) # add noise
#' Y <- do.call(rbind, Y)
#' Y <- Y - min(Y)
#' YB <- rbind(
#'     c(rep(0, 1), Y[1, ], rep(0, 9)),
#'     c(rep(0, 4), Y[1, ], rep(0, 6)),
#'     c(rep(0, 8), Y[1, ], rep(0, 2)),
#'     c(rep(0, 3), Y[1, ], rep(0, 7))
#' )
#' YA <- rbind(
#'     c(rep(0, 5), Y[1, ], rep(0, 5)),
#'     c(rep(0, 5), Y[1, ], rep(0, 5)),
#'     c(rep(0, 5), Y[1, ], rep(0, 5)),
#'     c(rep(0, 5), Y[1, ], rep(0, 5))
#' )
#' PA <- list(p + 5, p + 5, p + 5, p + 5)
#' PB <- list(p + 1, p + 4, p + 8, p + 3)
#' plot_align(YA, YB, PA, PB)
plot_align <- function(YA, YB, PA, PB, mfcol = c(nrow(YA), 1)) {
    stop("Implementation of this function is not finished yet.")
    s <- nrow(YA)
    if (!is.null(mfcol)) {
        opar <- par(mfcol = mfcol, mar = c(0, 2, 0, 0), oma = c(4.1, 2.1, 0, 0))
        on.exit(par(opar), add = TRUE)
    }
    for (i in seq_len(s)) {
        plot(x = seq_len(ncol(YB)), y = YB[i, ],
            type = "l", lty = 1, col = "darkgrey",
            xlim = c(1, ncol(YB)), ylim = c(0, max(YB)),
            xaxt = if (i == s) "s" else "n",
            xlab = "Datapoint Number", ylab = "Signal Intensity"
        )
        lines(x = seq_len(ncol(YA)), y = YA[i, ], col = "blue", lty = 1)
        abline(v = PB[[i]], col = transp("darkgrey", 0.5), lty = 2)
        abline(v = PA[[i]], col = transp("blue", 0.5), lty = 2)
    }
    mtext("Datapoint Number", side = 1, outer = TRUE, line = 3)
    mtext("Signal Intensity", side = 2, outer = TRUE, line = 1)
}

plot_empty <- function() {
    plot(
        x = 0.5, y = 0.5, type = "n",
        xlim = c(0, 1), xlab = "", xaxs = "i",
        ylim = c(0, 1), ylab = "", yaxs = "i",
        axes = FALSE, main = "", ann = FALSE
    )
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

# Plot Spectrum Helpers #####

ps_get_main_args <- function(env = parent.frame()) {
    args <- sapply(names(formals(plot_spectrum)), get, envir = env)
    args$foc_unit <- match.arg(args$foc_unit, c("ppm", "fraction"))
    if (identical(sort(names(args$sub_rgn)), c("x1", "x2", "y1", "y2"))) {
        args$sub_rgn <- args$sub_rgn[order(names(args$sub_rgn))]
    } else if (is.null(args$sub_rgn)) {
        # do nothing
    } else {
        stop("sub_rgn must be NULL or a numeric vector with names x1, x2, y1, y2")
    }
    args$lgd <- if (args$sub_show) FALSE else TRUE
    args$sm_show <- if (args$sub_show) FALSE else TRUE
    invisible(args)
}

ps_get_sub_args <- function(main_args, main_plot) {
    sub_args <- within(main_args, {
        foc_only <- TRUE
        mar <- sub_mar
        box_col <- sub_box_col
        axis_col <- sub_axis_col
        fill_col <- sub_fill_col
        line_col <- sub_line_col
        trp_show <- sub_trp_show
        lc_show <- sub_lc_show
        sup_show <- sub_sup_show
        sm_show <- TRUE
        main <- ""
        xlab <- ""
        ylab <- ""
        ysquash <- 0.96
        fig <- c(
            convert_pos(sub_rgn[1:2], c(0, 1), main_plot$plt$ndc$xlim),
            convert_pos(sub_rgn[3:4], c(0, 1), main_plot$plt$ndc$ylim)
        )
        add <- TRUE
        lgd <- TRUE
    })
    invisible(sub_args)
}

ps_connect_main_sub <- function(main_plot, sub_plot, cnct_col) {
    rct <- main_plot$foc$ndc
    sfg <- sub_plot$plt$ndc
    xds <- list(
        c(x0 = rct$xleft, x1 = sfg$xlim[1]),
        c(x0 = rct$xright, x1 = sfg$xlim[2])
    )
    yd <- c(y0 = rct$ytop, y1 = sfg$ylim[1])
    y <- grconvertY(yd, "ndc")
    opar <- par(xpd = TRUE)
    on.exit(par(opar), add = TRUE)
    for (xd in xds) {
        x <- grconvertX(xd, "ndc")
        segments(x[1], y[1], x[2], y[2], col = cnct_col)
    }
}

# styler: off
#' @noRd
#' @title Helper for [plot_spectrum()].
#' @description
#' For arguments see [plot_spectrum()].
#' For examples see `test/testthat/test-ps_internal.R`
ps_internal <- function(
    # Mandatory
    decon, # Output of [generate_lorentz_curves()]
    # Figure Region
    fig = NULL, # Drawing region in NDC.
    add = !is.null(fig), # If TRUE, the plot is added to the current figure.
    # Plot Region
    main = "", # Title of the plot.
    lgd = TRUE, # If TRUE, a legend is drawn.
    xlab = "Chemical Shift [ppm]", # Label for the x-axis.
    ylab = "Signal Intensity [au]", # Label for the y-axis.
    mar = c(4.1, 4.1, 0.1, 0.1), # Lines below/left/above/right plot region.
    box_col = "black", # Color of box surrounding plot region.
    axis_col = "black", # Color of tickmarks and ticklabels.
    fill_col = NULL, # Background color of plot region.
    # Focus Region (FR)
    foc_rgn = NULL, # Start and end of FR.
    foc_unit = "ppm", # Unit of `foc_rgn`. Either "fraction" or "ppm".
    foc_only = FALSE, # If TRUE, draw focusregion, else full spectrum.
    # Focus Rectangle
    rct_show = !is.null(foc_rgn) && !foc_only, # Draw rectangle around FR?
    rct_fill = transp("yellow"), # Background color of rectangle around FR.
    rct_col = "black", # Border color of rectangle around FR.
    # Spectrum Lines
    line_col = "black", # Color of raw signal intensities.
    sm_show = TRUE, # If TRUE, smoothed signal intensities are shown.
    sm_col = "blue", # Color of smoothed signal intensities.
    ysquash = 0.96, # Fraction of plot height to squash y-values into.
    sf_y_raw = 1e6, # Divide raw SI by this factor before drawing.
    # Lorentzians
    lc_show = TRUE, # If TRUE, the Lorentzian Curves (LCs) are shown.
    lc_col = "black", # Color of the Lorentzian Curves.
    lc_lty = 1, # Line type of the Lorentzian Curves.
    lc_fill = transp("black"), # BG-color of rectangles shown at LC-center.
    sup_show = TRUE, # If TRUE, the LC-Superposition (LC-Sup) is shown.
    sup_col = "red", # Color of LC-Sup.
    sup_lty = 1, # Line type of LC-Sup.
    # Peak Triplets
    trp_show = TRUE, # If TRUE, the peak triplets are shown.
    trp_col = rep("black", 4), # Colors for center, left, right, non-peak DPs.
    trp_pch = c(17, 4, 4, NA), # Pchars for center, left, right, non-peak DPs.
    # Misc
    verbose = FALSE,
    # Unused
    ...
) {
    args <- psi_get_args()
    old <- par(mar = mar, new = add)
    on.exit(par(old), add = TRUE)
    rst <- set_fig(fig = fig, add = add)
    on.exit(rst(), add = TRUE)
    dat <- psi_get_dat(args)
    plt <- psi_init_plot_region(dat, verbose)
    bgr <- psi_draw_bg(dat, fill_col, verbose)
    lns <- psi_draw_lines(dat, line_col, foc_only, sm_show, sm_col, verbose)
    trp <- psi_draw_triplets(dat, trp_show, trp_pch, trp_col, verbose)
    lcs <- psi_draw_lorentz_curves(dat, args, verbose)
    foc <- psi_draw_focus_rectangle(dat, rct_show, rct_fill, rct_col, verbose)
    axs <- psi_draw_axis(dat, main, xlab, ylab, axis_col, box_col, verbose)
    lgd <- psi_draw_legend(args, verbose)
    named(dat, plt, bgr, axs, lns, trp, lcs, foc, par = par())
}
# styler: on

# Plot Spectrum Internal Helpers #####

psi_get_args <- function(env = parent.frame()) {
    args <- sapply(names(formals(ps_internal)), get, envir = env)
    if (!is.null(args$foc_rgn) && length(args$foc_rgn) != 2) {
        stop("foc_rgn must be a numeric vector of length 2")
    }
    if (args$foc_only && is.null(args$foc_rgn)) {
        stop("foc_only requires foc_rgn to be specified")
    }
    invisible(args)
}

psi_get_dat <- function(args) {
    decon <- args$decon
    foc_rgn <- args$foc_rgn
    foc_unit <- args$foc_unit
    foc_only <- args$foc_only
    sf_y_raw <- args$sf_y_raw
    ysquash <- args$ysquash

    # Chemical Shift and Signal Intensity
    cs <- decon[["ppm"]] %||% decon[["x_values_ppm"]] %||% decon[["cs"]]
    si <- decon[["y_raw"]] %||% decon[["y_values_raw"]] %||% decon[["si"]]
    sis <- decon[["y_smooth"]] %||% decon[["y_values"]]
    si <- si / sf_y_raw

    # Peak Indices
    ipc <- decon$index_peak_triplets_middle # Peak centers
    ipl <- decon$index_peak_triplets_left # Peak borders left
    ipr <- decon$index_peak_triplets_right # Peak borders right
    idp <- rev(seq_along(cs)) # Data points
    ipp <- c(ipc, ipl, ipr) # Peak points
    inp <- setdiff(idp, ipp) # Non-peak points

    # Focus Region Limits
    foc_lim <- {
        if (is.null(foc_rgn)) {
            NULL
        } else if (foc_unit == "fraction") {
            quantile(cs, foc_rgn)
        } else if (foc_unit == "ppm") {
            foc_rgn
        } else {
            stop("foc_unit must be 'fraction' or 'ppm'")
        }
    }
    ifp <- {
        if (is.null(foc_rgn)) {
            integer()
        } else {
            which(cs > min(foc_lim) & cs < max(foc_lim))
        }
    }

    # Lorentzians
    x_0 <- decon$x_0_ppm
    A <- decon$A_ppm
    lambda <- decon$lambda_ppm

    # Plot region
    x <- if (foc_only) cs[ifp] else cs
    y <- if (foc_only) si[ifp] else si
    xlim <- c(max(x), min(x))
    ylim <- c(0, max(y) / ysquash)

    locals(without = c("decon"))
}

psi_init_plot_region <- function(dat = psi_get_dat(),
                                 verbose = FALSE) {
    if (verbose) logf("Initializing plot region")
    usr <- list(
        x    = dat$x,    y    = dat$y,
        xlim = dat$xlim, ylim = dat$ylim,
        xaxs = "i",      yaxs = "i",
        xlab = "",       ylab = "",
        type = "n",      axes = FALSE
    )
    do.call(plot, usr)
    ndc <- list(
        # Must be done AFTER plotting, or wrong user coords will be used
        xlim = grconvertX(usr$xlim, from = "user", to = "ndc"),
        ylim = grconvertY(usr$ylim, from = "user", to = "ndc")
    )
    named(usr, ndc)
}

psi_draw_bg <- function(dat = psi_get_dat(),
                        fill_col = NULL,
                        verbose = FALSE) {
    if (verbose) logf("Drawing background")
    bgr <- within(list(), {
        xleft <- dat$xlim[1]
        xright <- dat$xlim[2]
        ybottom <- dat$ylim[1]
        ytop <- dat$ylim[2]
        col <- fill_col
        border <- NA
    })
    if (!is.null(fill_col)) do.call(rect, bgr)
    return(bgr)
}

psi_draw_axis <- function(dat,
                          main = "",
                          xlab = "Chemical Shift [ppm]",
                          ylab = "Signal Intensity [au]",
                          axis_col = "black",
                          box_col = "black",
                          verbose = FALSE) {
    if (verbose) logf("Drawing axis")
    xtks <- seq(dat$xlim[1], dat$xlim[2], length.out = 5)
    ytks <- seq(dat$ylim[1], max(dat$y), length.out = 5)
    for (i in 2:12) if (length(unique(xtklabs <- round(xtks, i))) >= 5) break
    for (i in 2:12) if (length(unique(ytklabs <- round(ytks, i))) >= 5) break
    axis(
        side = 1, at = xtks, labels = xtklabs,
        col = NA, col.ticks = axis_col, col.axis = axis_col
    )
    axis(
        side = 2, at = ytks, labels = ytklabs,
        col = NA, col.ticks = axis_col, col.axis = axis_col, las = 1
    )
    title(main = main, xlab = xlab, ylab = ylab)
    box(col = box_col)
    named(xtks, ytks, xlab, ylab)
}

psi_draw_lines <- function(dat,
                           line_col = "black",
                           foc_only = FALSE,
                           sm_show = TRUE,
                           sm_col = "blue",
                           verbose = FALSE) {
    x <- if (foc_only) dat$cs[dat$ifp] else dat$cs
    y <- if (foc_only) dat$si[dat$ifp] else dat$si
    ys <- if (foc_only) dat$sis[dat$ifp] else dat$sis
    if (verbose) logf("Drawing raw signal")
    lines(x, y, type = "l", col = line_col, lty = 1)
    if (sm_show && !is.null(ys)) {
        if (verbose) logf("Drawing smoothed signal")
        lines(x, ys, type = "l", col = sm_col, lty = 1)
    }
}

psi_draw_triplets <- function(dat,
                              show = TRUE,
                              pch = c(17, 4, 4, NA),
                              col = c("red", "blue", "blue", "black"),
                              verbose = FALSE) {
    if (isFALSE(show)) return(NULL) # styler: off
    if (verbose) logf("Drawing peak triplets")
    x <- dat$cs
    y <- dat$sis %||% dat$si
    p <- dat$ipc
    l <- dat$ipl
    r <- dat$ipr
    q <- dat$inp
    points(x[p], y[p], col = col[1], pch = pch[1]) # 017 = triangle
    points(x[l], y[l], col = col[2], pch = pch[2]) # 000 = open square
    points(x[r], y[r], col = col[3], pch = pch[3]) # 004 = x character
    points(x[q], y[q], col = col[4], pch = pch[4]) # 124 = vertical dash
}

psi_draw_focus_rectangle <- function(dat,
                                     show = TRUE,
                                     fill = rgb(0, 0, 1, alpha = 0.1),
                                     col = "blue",
                                     verbose = FALSE) {
    if (is.null(dat$foc_lim) || !show) return(NULL) # styler: off
    usr <- list(
        xleft = max(dat$foc_lim),
        xright = min(dat$foc_lim),
        ybottom = 0,
        ytop = max(dat$si[dat$ifp]) / 0.96,
        col = fill,
        border = col
    )
    ndc <- list(
        xleft = grconvertX(usr$xleft, to = "ndc"),
        xright = grconvertX(usr$xright, to = "ndc"),
        ybottom = grconvertY(usr$ybottom, to = "ndc"),
        ytop = grconvertY(usr$ytop, to = "ndc")
    )
    if (verbose) logf("Drawing focus rectangle")
    do.call(rect, usr)
    named(usr, ndc)
}

psi_draw_lorentz_curves <- function(dat, args = NULL, verbose = FALSE) {
    lc_show <- args$lc_show %||% TRUE
    lc_col <- args$lc_col %||% "black"
    lc_lty <- args$lc_lty %||% 1
    lc_fill <- args$lc_fill %||% transp(lc_col)
    sup_show <- args$sup_show %||% TRUE
    sup_col <- args$sup_col %||% "red"
    sup_lty <- args$sup_lty %||% 1
    foc_only <- args$foc_only %||% FALSE
    if ((is.null(lc_show) && is.null(sup_show)) || is.null(dat$A)) {
        return(invisible(NULL))
    }
    x <- if (foc_only) dat$cs[dat$ifp] else dat$cs
    Y <- matrix(nrow = length(x), ncol = length(dat$A))
    if (verbose) logf("Drawing individual Lorentzian Curves")
    for (i in seq_along(dat$A)) {
        Y[, i] <- psi_draw_lorentz_curve(
            x,
            x_0 = dat$x_0[i],
            A = dat$A[i],
            lambda = dat$lambda[i],
            lc_show = lc_show,
            lc_col = lc_col,
            lc_lty = lc_lty,
            lc_fill = lc_fill
        )
    }
    if (sup_show) {
        if (verbose) logf("Drawing superposition of Lorentzian Curves")
        lines(x = x, y = rowSums(Y), col = sup_col, lty = sup_lty, type = "l")
    }
    invisible(NULL)
}

psi_draw_lorentz_curve <- function(x,
                                   x_0,
                                   A,
                                   lambda,
                                   lc_show = TRUE,
                                   lc_col = "black",
                                   lc_lty = 1,
                                   lc_fill = NULL) {
    y <- lorentz(x, x_0, A, lambda)
    if (lc_show) {
        near_zero <- abs(y) < 0.01
        y_big <- y[!near_zero]
        x_big <- x[!near_zero]
        lines(
            x = x_big,
            y = y_big,
            col = lc_col,
            lty = lc_lty,
            type = "l"
        )
    }
    if (!is.null(lc_fill)) {
        rect(
            xleft = x_0 + lambda,
            xright = x_0 - lambda,
            ybottom = par("usr")[3],
            ytop = lorentz(x_0, x_0, A, lambda),
            col = lc_fill,
            border = NA
        )
    }
    return(y)
}

psi_draw_legend <- function(args, verbose = FALSE) {
    if (isFALSE(args$lgd)) return(invisible(NULL)) # styler: off
    if (verbose) logf("Drawing legend")
    legend(
        x = "topright",
        legend = c(
            "Raw Signal",
            "Smoothed Signal",
            "Peak Triplets",
            "Single Lorentzian",
            "Sum of Lorentzians"
        ),
        col = c(
            args$line_col, args$sm_col, args$trp_col[1], args$lc_col, args$sup_col
        ),
        lty = c(1, 1, NA, 1, 1),
        pch = c(NA, NA, args$trp_pch[1], NA, NA),
    )
}

psi_setup_dev_env <- function() {
    sim_01 <- metabodecon_file("sim/sim_01")
    decon <- generate_lorentz_curves(
        sim_01,
        sfr = c(3.42, 3.58), wshw = 0, delta = 0.1,
        ask = FALSE, verbose = FALSE
    )
    args <- stub("ps_internal", decon = decon, foc_rgn = c(3.55, 3.52))
    invisible(args)
}

# Helpers #####

#' @noRd
#'
#' @title Set Figure Region
#'
#' @description
#' Calling `par(fig=xxyy)` resets the current multi-figure configuration (MFC)
#' to one row and one column. `set_fig()` handles this scenario, by first
#' storing the current MFC, then calling `par(fig=xxyy)` and finally returning a
#' function that can be used to restore the MFC. See 'Details' for further
#' information.
#'
#' @param fig Region to draw into, given as normalized device coordinates.
#'
#' @param add If TRUE, the new plot is added to the existing plot.
#'
#' @return Function to reset the MFC.
#'
#' @details
#' Note 1: Setting `par(fig=xxyy)` resets the current MFC to one row and one
#' column. I.e., array layouts defined by setting `mfrow` or `mfcol`, must be
#' saved before changing `fig` and restored after the plot has been drawn. The
#' same applies to the current figure number `mfg` (see Note 2).
#'
#' Note 2: When restoring `mfg` it's important to additionally advance one
#' frame, because when querying `mfg`, the "current figure number" is returned,
#' but when setting `mfg`, the value is interpreted as "next figure number".
#' I.e. without advancing one frame, the next figure would be drawn into the
#' spot of the current figure.
#'
#' Note 3: If `fig=xxyy` and `add=FALSE`, it is still necessary to use
#' `par(fig=xxyy, new=TRUE)` to prevent the device from being cleared (which
#' would be bad in a multi-figure environment). I.e., the figure number must be
#' advanced manually, as would have been done by a normal plot. This manual
#' increment must be done before the MFC is reset by `par(fig=xxyy, new=TRUE)`.
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

#' @noRd
#' @description
#' Returns TRUE if the current multi-figure gets filled by row, else FALSE.
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
        # In this case it doesn't matter, as the figure spots in the grid will
        # be filled top-to-bottom / left-to-right in both cases.
        return(TRUE)
    } else {
        # We have at least two rows AND two cols. So what we can do is to set
        # c(1, 1) as next figure, and then advance one frame. If we end up in
        # c(1, 2) we are row-oriented, if we end up in c(2, 1) we are
        # column-oriented. After doing this we can reset the current figure
        # number to the original value.
        par(mfg = c(1, 1, nrows, ncols))
        on.exit({
            par(mfg = mfg) # Reset current figure number
            plot_empty() # (1)
            # (1) When querying `mfg` we get the "current figure number". But
            # when setting it, we set the "next figure number". I.e. we need to
            # advance one frame, or we would set the "current figure number" as
            # "next figure number".
        })
        plot_empty() # Draw into c(1, 1)
        plot_empty() # Draw into c(1, 2) or c(2, 1)
        mfg2 <- par("mfg") # Query current position
        return(if (mfg2[1] == 1) TRUE else FALSE)
    }
}

# Deprecated #####

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

plot_spectra <- function(ss = generate_lorentz_curves_sim(),
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

plot_noise_methods <- function(siRND, siSFR, n = 300, start = 5000) {
    ymin <- min(c(min(siRND), min(siSFR)))
    ymax <- max(c(max(siRND), max(siSFR)))
    ylim <- c(ymin, ymax)
    opar <- par(mfrow = c(5, 1), mar = c(3, 4, 0, 1))
    on.exit(par(opar), add = TRUE)
    for (i in 1:5) {
        redT <- rgb(1, 0, 0, 0.1)
        bluT <- rgb(0, 0, 1, 0.1)
        idx <- ((i - 1) * n + 1):(i * n) + start
        ysiRND <- siRND[idx]
        ysiSFR <- siSFR[idx]
        plot(1:n, ysiRND, type = "n", ylim = ylim, ylab = "", xlab = "", xaxt = "n")
        axis(1, at = seq(1, n, by = 50), labels = idx[seq(1, n, by = 50)])
        points(1:n, ysiRND, col = "red", pch = 20)
        points(1:n, ysiSFR, col = "blue", pch = 20)
        lines(1:n, ysiRND, col = "red")
        lines(1:n, ysiSFR, col = "blue")
        lines(1:n, rep(0, n), col = "black", lty = 2)
        polygon(c(1:n, n:1), c(ysiRND, rep(0, n)), col = redT, border = NA)
        polygon(c(1:n, n:1), c(ysiSFR, rep(0, n)), col = bluT, border = NA)
        legend("topleft", NULL, c("RND", "SFR"), col = c("red", "blue"), lty = 1)
    }
}
