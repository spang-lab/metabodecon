# Public #####

#' @export
#' @title Plot Spectrum
#' @description Plot a spectrum based on the provided deconvolution data and zoom in on a specific region of interest in the spectrum.
#'
#' `r lifecycle::badge("experimental")`
#' @param decon An object as returned by [generate_lorentz_curves()], containing the deconvolution data. Must include either `x_values_ppm` or `ppm` for the x-axis values, and either `y_values` or `y_smooth` for the y-axis values.
#' @param foc_rgn Numeric vector specifying the start and end of focus region.
#' @param foc_unit Character string specifying the unit in which `foc_rgn` is given. Can be "fraction" or "ppm".
#' @param sub_show Logical. If TRUE, a sub figure is drawn within the main plot region. If FALSE, the sub figure is not drawn.
#' @param sub_rgn Either NULL or a numeric vector of the form `c(x1, x2, y1, y2)` giving the left/right/bottom/top coordinates of the sub figure region in "normalized plot coordinates" (as described in [graphics::grconvertX()]). If provided, the focussed region is drawn as sub figure within the main plot region. Setting `sub_rgn` to NULL will prevent the sub figure from being drawn.
#' @param foc_only Logical. If TRUE, only the focussed region is drawn. If FALSE, the full spectrum is drawn.
#' @param foc_fill Background color of the rectangle around the focus region.
#' @param foc_col Border color of the rectangle around the focus region.
#' @param mar Number of lines below/left/above/right plot region (used for axis annotations).
#' @param line_col Color of the spectrum line.
#' @param axis_col Color of tickmarks and ticklabels.
#' @param fill_col Background color of the plot region.
#' @param box_col Border color of the box surrounding the plot region.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param ysquash Fraction of plot height to squash y-values into. Useful in combination with `sub_rgn` to prevent the spectrum lines from overlapping with the sub figure showing the focussed region.
#' @param trp_show Logical. If TRUE, the peak triplets are shown.
#' @param lc_show Logical. If TRUE, the Lorentzian Curves are shown.
#' @param sup_show Logical. If TRUE, the superposition of the Lorentzian Curves is shown.
#' @param sub_mar Margins of the sub figure.
#' @param sub_line_col Color of the lines in the sub figure.
#' @param sub_axis_col Color of the axis in the sub figure.
#' @param sub_fill_col Background color of the sub figure.
#' @param sub_box_col Border color of the box surrounding the sub figure.
#' @param sub_trp_show Logical. If TRUE, the peak triplets are shown in the sub figure.
#' @param sub_lc_show Logical. If TRUE, the Lorentzian Curves are shown in the sub figure.
#' @param sub_sup_show Logical. If TRUE, the superposition of the Lorentzian Curves is shown in the sub figure.
#' @param yscale10 Logical. If TRUE, scales the y-axis by a power of 10 so that all y-values are between 0 and 100.
#' @param lc_col Color of the Lorentzian Curves.
#' @param lc_lty Line type of the Lorentzian Curves.
#' @param trp_col Vector of length 4 giving the colors for the peak triplets. The first three colors specify the color used for each peak-center, left-border and right-border and the fourth color used for any non-peak data point.
#' @param trp_pch Vector of length 4 giving the plotting characters for the peak triplets.
#' @param sup_col Color of the superposition of the Lorentzian Curves.
#' @param sup_lty Line type of the superposition of the Lorentzian Curves.
#' @param cnct_show Logical. If TRUE, connecting lines are drawn between the main and sub figure.
#' @param cnct_col Color of the lines connecting the main and sub figure.
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
#' # Prepare a deconvoluted spectrum as input
#' sim_01 <- metabodecon_file("sim/sim_01")
#' spec <- read_spectrum(sim_01)
#' decon <- generate_lorentz_curves(
#'     sim_01,
#'     sfr = c(3.42, 3.58), ws = 0, delta = 0.1,
#'     ask = FALSE, verbose = FALSE
#' )
#' testfunc <- function(figs = 1:6) {
#'     n <- length(figs)
#'     nr <- ceiling(sqrt(n))
#'     nc <- if (nr^2 > n) nr - 1 else nr
#'     spec <- read_spectrum(metabodecon_file("sim/sim_01"))
#'     decon <- glc("sim_01", debug = FALSE)$rv
#'     opar <- par(mfrow = c(nr, nc))
#'     on.exit(par(opar))
#'
#'     # Plot the full (non-deconvoluted) spectrum
#'     if (1 %in% figs) plot_spectrum(spec)
#'
#'     # Focus on a specific region, specified in ppm or as fraction
#'     if (2 %in% figs) plot_spectrum(decon, foc_rgn = c(3.49, 3.45), foc_unit = "ppm")
#'     if (3 %in% figs) plot_spectrum(decon, foc_rgn = c(0.40, 0.30))
#'
#'     # Change margin, color and position of the focus region
#'     if (4 %in% figs) plot_spectrum(decon,
#'         sub_mar = c(4, 4, 0, 0),
#'         foc_fill = rgb(0.9, 0.5, 0.9, alpha = 0.1),
#'         foc_col = "violet",
#'         sub_rgn = c(x1 = 0.1, x2 = 0.9, y1 = 0.4, y2 = 0.9)
#'     )
#'
#'     # Remove connecting lines and fill colors
#'     if (5 %in% figs) plot_spectrum(decon,
#'         foc_fill = NULL,
#'         sub_fill_col = NULL,
#'         cnct_col = NULL
#'     )
#'
#'     # Hide xlab and ylab
#'     if (6 %in% figs) plot_spectrum(decon,
#'         xlab = "",
#'         ylab = "",
#'         mar = c(2, 2, 0, 1)
#'     )
#' }
#' testfunc(3:4)
#' testfunc(1:6)
plot_spectrum <- function(
    decon,
    # Common Settings
    foc_rgn = c(0.40, 0.35),
    foc_unit = "fraction",
    sub_show = TRUE,
    sub_rgn = c(x1 = 0.05, x2 = 0.95, y1 = 0.2, y2 = 0.95),
    # Focus Region
    foc_only = if (sub_show) FALSE else TRUE,
    foc_fill = rgb(1, 1, 0, alpha = 0.08),
    foc_col = "black",
    # Settings for Main Figure
    mar = c(4, 4, 0, 1),
    line_col = "black",
    axis_col = "black",
    fill_col = NULL,
    box_col = "black",
    xlab = "Chemical Shift [ppm]",
    ylab = "Signal Intensity [au]",
    ysquash = if (sub_show) sub_rgn[3] * 0.96 else 0.96,
    trp_show = if (sub_show) FALSE else TRUE,
    lc_show = if (sub_show) FALSE else TRUE,
    sup_show = if (sub_show) FALSE else TRUE,
    # Settings for Sub Figure
    sub_mar = c(2, 2, 0, 0),
    sub_line_col = line_col,
    sub_axis_col = axis_col,
    sub_fill_col = foc_fill,
    sub_box_col = foc_col,
    sub_trp_show = TRUE,
    sub_lc_show = TRUE,
    sub_sup_show = TRUE,
    # Settings for both Figures
    yscale10 = TRUE,
    trp_col = c("red", "red", "red", "black"),
    lc_col = "darkgrey",
    trp_pch = c(17, 4, 4, NA),
    lc_lty = 1,
    sup_col = "red",
    sup_lty = 1,
    # Connecting Lines
    cnct_show = if (sub_show) TRUE else FALSE,
    cnct_col = foc_col) {
    #
    # Function Body
    #
    # px_setup_dev_env()
    main_args <- px_get_main_args()
    main_plot <- do.call(plot_spec, main_args)
    sub_args <- if (sub_show) px_get_sub_args(main_args, main_plot)
    sub_plot <- if (sub_show) do.call(plot_spec, sub_args)
    cnct_plot <- if (cnct_show) px_connect_main_sub(main_plot, sub_plot, cnct_col)
    invisible(named(main_args, main_plot, sub_args, sub_plot, cnct_plot))
}

# Helpers #####

px_get_main_args <- function(env = parent.frame()) {
    args <- sapply(names(formals(plot_spectrum)), get, envir = env)
    args$foc_unit <- match.arg(args$foc_unit, c("ppm", "fraction"))
    if (identical(sort(names(args$sub_rgn)), c("x1", "x2", "y1", "y2"))) {
        args$sub_rgn <- args$sub_rgn[order(names(args$sub_rgn))]
    } else if (is.null(args$sub_rgn)) {
        # do nothing
    } else {
        stop("sub_rgn must be NULL or a numeric vector with names x1, x2, y1, y2")
    }
    invisible(args)
}

px_get_sub_args <- function(main_args, main_plot) {
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
        xlab <- ""
        ylab <- ""
        ysquash <- 0.96
        fig <- c(
            convert_pos(sub_rgn[1:2], c(0, 1), main_plot$plt$ndc$xlim),
            convert_pos(sub_rgn[3:4], c(0, 1), main_plot$plt$ndc$ylim)
        )
        add <- TRUE
    })
    invisible(sub_args)
}

px_connect_main_sub <- function(main_plot, sub_plot, cnct_col) {
    rct <- main_plot$foc$ndc
    sfg <- sub_plot$plt$ndc
    xds <- list(
        c(x0 = rct$xleft,  x1 = sfg$xlim[1]),
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


# Development #####

px_setup_dev_env <- function() {
    # Prepare a deconvoluted spectrum as input
    sim_01 <- metabodecon_file("sim/sim_01")
    spec <- read_spectrum(sim_01)
    decon <- generate_lorentz_curves(
        data_path = sim_01,
        sfr = c(3.42, 3.58),
        ws = 0,
        delta = 0.1,
        ask = FALSE,
        verbose = FALSE
    )
    args <- stub(
        func = plot_spectrum,
        decon = decon,
        foc_rgn = c(3.49, 3.45),
        foc_unit = "ppm",
        envir = .GlobalEnv
    )
    width <- options("width")
    line <- paste0(collapse(rep("-", width), ""), "\n")
    cat(line)
    cat("Assigned to .GlobalEnv:\n")
    cat(line)
    str(args, 1, give.attr = FALSE)
    cat(line)
    invisible(args)
}

#' @examples
#' px_test(1, 2) # first two plots
#' px_text(2:4)  # second to fourth plot
px_test <- function(figs = 1:5) {
    n <- length(figs)
    nr <- ceiling(sqrt(n))
    nc <- if (nr^2 > n) nr - 1 else nr
    spec <- read_spectrum(metabodecon_file("sim/sim_01"))
    decon <- glc("sim_01", debug = FALSE)$rv
    opar <- par(mfrow = c(nr, nc))
    on.exit(par(opar))

    # Plot the full (non-deconvoluted) spectrum
    if (1 %in% figs) plot_spectrum(spec)

    # Focus on a specific region, specified in ppm or as fraction
    if (2 %in% figs) plot_spectrum(decon, foc_rgn = c(3.49, 3.45), foc_unit = "ppm")
    if (3 %in% figs) plot_spectrum(decon, foc_rgn = c(0.40, 0.30))

    # Change margin, color and position of the focus region
    if (4 %in% figs) plot_spectrum(decon,
        sub_mar = c(4, 4, 0, 0),
        foc_fill = rgb(0.9, 0.5, 0.9, alpha = 0.1),
        foc_col = "violet",
        sub_rgn = c(x1 = 0.1, x2 = 0.9, y1 = 0.4, y2 = 0.9)
    )

    # Remove connecting lines and fill colors
    if (5 %in% figs) plot_spectrum(decon,
        foc_fill = NULL,
        sub_fill_col = NULL,
        cnct_col = NULL
    )

    # Hide xlab and ylab
    if (6 %in% figs) plot_spectrum(decon,
        xlab = "",
        ylab = "",
        mar = c(2, 2, 0, 1)
    )
}

# Throwaway #####

px_test_gr_convert <- function() {
    par(mfrow = c(1, 2), xpd = TRUE)
    plot_dummy()
    x <- c(0.25, 0.75)
    xds <- list()
    units <- c("in", "dev", "ndc", "nfc", "npc", "nic", "lines", "chars")
    for (i in seq_along(units)) {
        y <- 0.1 * i
        unit <- units[i]
        xu <- grconvertX(x, from = unit, to = "user")
        xds[[i]] <- grconvertX(x, from = unit, to = "ndc")
        lines(x = xu, y = rep(y, 2))
        xur <- collapse(round(xu, 2), " - ")
        label <- sprintf("0.25 - 0.75 %s == %s %s", unit, xur, "user")
        text(x = 0.5, y = y, labels = label, pos = 3)
    }
    with_fig(
        fig = c(0, 1, 0, 1),
        expr = {
            plot_empty() # important to setup user coordinates
            for (i in seq_along(units)) {
                y <- 0.1 * i + 0.05
                xu <- grconvertX(xds[[i]], from = "ndc", to = "user")
                lines(x = xu, y = rep(y, 2), col = "red")
            }
        }
    )
    par(mfrow = c(1, 1), xpd = FALSE)
}
