# Exported #####

#' @export
#'
#' @title Plot Spectrum
#'
#' @description
#' Plot a spectrum and zoom in on a specific region.
#'
#' `r lifecycle::badge("experimental")`
#'
#' @param obj An object as returned by [generate_lorentz_curves()], containing
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
#' obj <- generate_lorentz_curves_sim(spec)
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
#' plot_spectrum(obj, foc_rgn = c(3.49, 3.45), foc_unit = "ppm")
#' plot_spectrum(obj, foc_rgn = c(0.40, 0.30))
#' plot_spectrum(obj,
#'     sub_mar = c(4, 4, 0, 0),
#'     rct_fill = rgb(0.9, 0.5, 0.9, alpha = 0.1),
#'     rct_col = "violet",
#'     sub_rgn = c(x1 = 0.1, x2 = 0.9, y1 = 0.4, y2 = 0.9)
#' )
#' plot_spectrum(obj, rct_fill = NULL, sub_fill_col = NULL, cnct_col = NULL)
#' plot_spectrum(obj, xlab = "", ylab = "", mar = c(2, 2, 0, 1))
plot_spectrum <- function(
    # Common Settings
    obj,
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
    main_plot <- do.call(draw_spectrum, main_args)
    sub_args <- if (sub_show) ps_get_sub_args(main_args, main_plot)
    sub_plot <- if (sub_show) do.call(draw_spectrum, sub_args)
    cnct_plot <- if (cnct_show) ps_connect_main_sub(main_plot, sub_plot, cnct_col)
    invisible(named(main_args, main_plot, sub_args, sub_plot, cnct_plot))
}


# Helpers #####

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

