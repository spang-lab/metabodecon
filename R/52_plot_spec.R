# Plot_Spec #####

#' @noRd
#' @title Plot Spectrum Internal
#' @description Plot a spectrum based on the provided chemical shift and signal intensity data. Helper of public function [plot_spectrum()]. Should not be called directly by the user.
#' @param decon An object as returned by [generate_lorentz_curves()], containing the deconvolution data. Must include either `x_values_ppm` or `ppm` for the x-axis values, and either `y_values` or `y_smooth` for the y-axis values.
#' @param foc_rgn Focus region in ppm.
#' @param foc_only If TRUE, only the focussed region is drawn. If FALSE, the full spectrum is drawn.
#' @param mar Number of lines below/left/above/right plot region.
#' @param line_col Color of spectrum line.
#' @param box_col Color of box around plot region.
#' @param axis_col Color of tickmarks and ticklabels.
#' @param fill_col Background color of plot region.
#' @param foc_fill Fill color of rectangle around foc_rgn region.
#' @param foc_col Border color of rectangle around foc_rgn region.
#' @param ysquash Fraction of plot height to squash y-values into.
#' @param fig Region to draw into, given as normalized device coordinates
#' @param add If TRUE, the new plot is added to the existing plot.
#' @return A list with the following 3 elements:
#' - plt: List with 12 elements used for creating the main plot, in particular `xlim_ndc` and `ylim_ndc`, which are the normalized device coordinates of the x- and y-axis limits.
#' - rct: List of 10 elements used for creating the rectangle around the focus region, in particular `xleft_ndc`, `xright_ndc`, `ybottom_ndc`, and `ytop_ndc`, which are the normalized device coordinates of the rectangle.
#' - par: List of the 72 graphical parameters ([graphics::par()]) used for creating the plot.
#' @examples
#' sim_01 <- metabodecon_file("sim/sim_01")
#' decon <- generate_lorentz_curves(
#'     sim_01,
#'     sfr = c(3.42, 3.58), ws = 0,
#'     ask = FALSE, verbose = FALSE,
#'     delta = 0.1
#' )
#' plot_dummy <- function() {
#'     plot(0, 0, ylim = c(0, 1), xlim = c(0, 1), xaxs = "i", yaxs = "i")
#'     text(0.5, 0.5, "dummy")
#' }
#' leftmiddle <- c(0.1, 0.4, 0.4, 0.6)
#' leftbottom <- c(0.1, 0.4, 0.1, 0.3)
#' p <- local({
#'     p <- list()
#'     opar <- par(mfrow = c(3, 2), mar = c(2, 2, 0.5, 0.5))
#'     on.exit(par(opar))
#'     p[[1]] <- plot_dummy()
#'     p[[2]] <- plot_spec(decon, foc_rgn = c(3.55, 3.52))
#'     p[[3]] <- plot_dummy()
#'     p[[3]] <- plot_spec(decon, foc_rgn = c(3.55, 3.52),
#'                         foc_only = TRUE, fig = leftmiddle,
#'                         fill_col = rgb(0, 0, 1, 0.1))
#'     p[[4]] <- plot_dummy()
#'     p[[5]] <- plot_spec(decon, fig = leftbottom, add = FALSE)
#'     p[[6]] <- plot_spec(decon, lc_show = FALSE, trp_show = FALSE)
#'     p
#' })
plot_spec <- function(
    decon,
    # Figure Region
    fig = NULL,
    add = !is.null(fig),
    # Plot Region
    mar = c(4.1, 4.1, 0.1, 0.1),
    box_col = "black",
    axis_col = "black",
    fill_col = NULL,
    xlab = "Chemical Shift [ppm]",
    ylab = "Signal Intensity [au]",
    # Focus Region
    foc_rgn = NULL,
    foc_unit = "ppm",
    foc_only = FALSE,
    foc_fill = rgb(0, 0, 1, alpha = 0.1),
    foc_col = "blue",
    # Spectrum Lines
    line_col = "black",
    ysquash = 0.96,
    yscale10 = TRUE,
    d2_show = FALSE,
    # Peak Triplets
    trp_show = TRUE,
    trp_col = c("red", "blue", "blue", "black"),
    trp_pch = c(17, 4, 4, NA),
    # Lorentzians
    lc_show = TRUE,
    lc_col = "darkgrey",
    lc_lty = 1,
    sup_show = TRUE,
    sup_col = "red",
    sup_lty = 1,
    # Unused
    ...
    ) {

    ps_check_args(foc_rgn, foc_only)
    old <- par(mar = mar, new = add)
    on.exit(par(old), add = TRUE)
    rst <- set_fig(fig = fig, add = add)
    on.exit(rst(), add = TRUE)

    dat <- ps_get_dat(decon, foc_rgn, foc_unit, foc_only, yscale10, ysquash)
    plt <- ps_init_plot_region(dat)
    bgr <- ps_draw_bg(dat, fill_col)
    lns <- ps_draw_lines(dat, line_col, d2_show, foc_only)
    trp <- if (trp_show) ps_draw_triplets(dat, trp_pch, trp_col)
    lcs <- ps_draw_lorentz_curves(dat, lc_show, lc_col, lc_lty, sup_show,
    sup_col, sup_lty, foc_only)
    foc <- if (!foc_only) ps_draw_focus_region(dat,  foc_fill, foc_col)
    axs <- ps_draw_axis(dat, xlab, ylab, axis_col, box_col)
    named(dat, plt, bgr, axs, lns, trp, lcs, foc, par = par())
}

# Helpers #####

ps_check_args <- function(foc_rgn = NULL, foc_only = FALSE) {
    if (!is.null(foc_rgn) && length(foc_rgn) != 2) {
        stop("foc_rgn must be a numeric vector of length 2")
    }
    if (foc_only && is.null(foc_rgn)) {
        stop("foc_only requires foc_rgn to be specified")
    }
}

ps_get_dat <- function(decon = glc("sim_01", debug = FALSE)$rv,
                       foc_rgn = NULL,
                       foc_unit = "ppm",
                       foc_only = FALSE,
                       yscale10 = TRUE,
                       ysquash = 0.96) {

    # Chemical Shift and Signal Intensity
    cs <- decon[["ppm"]] %||% decon[["x_values_ppm"]] %||% decon[["cs"]]
    si <- decon[["y_raw"]] %||% decon[["y_values_raw"]] %||% decon[["si"]]
    sis <- decon[["y_smooth"]] %||% decon[["y_values"]]
    if (yscale10) {
        xp <- floor(log10(max(si))) - 1
        si <- si / 10^xp
        if (!is.null(sis)) {
            # Scale smoothed SI to same range as the raw SI
            xps <- floor(log10(max(si) / max(sis)))
            sis <- sis / 10^xps
        }
    }

    # Peak Indices
    ipc <- decon$index_peak_triplets_middle # Peak centers
    ipl <- decon$index_peak_triplets_left # Peak borders left
    ipr <- decon$index_peak_triplets_right # Peak borders right
    idp <- rev(seq_along(cs)) # Data points
    ipp <- c(ipc, ipl, ipr) # Peak points
    inp <- setdiff(idp, ipp) # Non-peak points

    # Focus Region Limits
    foc_lim <- {
        if (is.null(foc_rgn)) NULL
        else if (foc_unit == "fraction") quantile(cs, foc_rgn)
        else if (foc_unit == "ppm") foc_rgn
        else stop("foc_unit must be 'fraction' or 'ppm'")
    }
    ifp <- {
        if (is.null(foc_rgn)) integer()
        else which(cs > min(foc_lim) & cs < max(foc_lim))
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

ps_init_plot_region <- function(dat = ps_get_dat()) {
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

ps_draw_bg <- function(dat = ps_get_dat(), fill_col) {
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

ps_draw_axis <- function(dat, xlab, ylab, axis_col, box_col) {
    xtks <- seq(dat$xlim[1], dat$xlim[2], length.out = 5)
    ytks <- seq(dat$ylim[1], max(dat$y), length.out = 5)
    xtklabs <- round(xtks, 3)
    ytklabs <- round(ytks, 2)
    axis(
        side = 1, at = xtks, labels = xtklabs,
        col = NA, col.ticks = axis_col, col.axis = axis_col
    )
    axis(
        side = 2, at = ytks, labels = ytklabs,
        col = NA, col.ticks = axis_col, col.axis = axis_col, las = 1
    )
    title(xlab = xlab)
    title(ylab = ylab)
    box(col = box_col)
    named(xtks, ytks, xlab, ylab)
}

ps_draw_lines <- function(dat,
                          line_col = "black",
                          d2_show = FALSE,
                          foc_only = FALSE,
                          sm_col = "blue") {
    x <- if (foc_only) dat$cs[dat$ifp] else dat$cs # Chemical Shift
    y <- if (foc_only) dat$si[dat$ifp] else dat$si # Raw SI
    ys <- if (foc_only) dat$sis[dat$ifp] else dat$sis # Smoothed SI
    if (d2_show) message("TODO")
    lines(x, y, type = "l",  col = line_col, lty = 1)
    if (!is.null(ys)) lines(x, ys, type = "l", col = sm_col, lty = 1)
}

ps_draw_triplets <- function(dat,
                             pch = c(17, 4, 4, NA),
                             col = c("red", "blue", "blue", "black")) {
    x <- dat$cs; y <- dat$sis %||% dat$si
    p <- dat$ipc; l <- dat$ipl; r <- dat$ipr; q <- dat$inp
    points(x[p], y[p], col = col[1], pch = pch[1]) # 017 = triangle
    points(x[l], y[l], col = col[2], pch = pch[2]) # 000 = open square
    points(x[r], y[r], col = col[3], pch = pch[3]) # 004 = x character
    points(x[q], y[q], col = col[4], pch = pch[4]) # 124 = vertical dash
}

ps_draw_focus_region <- function(dat,
                                 foc_fill = rgb(0, 0, 1, alpha = 0.1),
                                 foc_col = "blue") {
    if (is.null(dat$foc_lim)) {
        return(NULL)
    }
    usr <- list(
        xleft = max(dat$foc_lim),
        xright = min(dat$foc_lim),
        ybottom = 0,
        ytop = max(dat$si[dat$ifp]) / 0.96,
        col = foc_fill,
        border = foc_col
    )
    ndc <- list(
        xleft = grconvertX(usr$xleft, to = "ndc"),
        xright = grconvertX(usr$xright, to = "ndc"),
        ybottom = grconvertY(usr$ybottom, to = "ndc"),
        ytop = grconvertY(usr$ytop, to = "ndc")
    )
    do.call(rect, usr)
    named(usr, ndc)
}

ps_draw_lorentz_curves <- function(dat = ps_get_dat(),
                                   lc_show = TRUE,
                                   lc_col = "darkgrey",
                                   lc_lty = 1,
                                   sup_show = TRUE,
                                   sup_col = "red",
                                   sup_lty = 1,
                                   foc_only = FALSE) {
    if ((is.null(lc_show) && is.null(sup_show)) || is.null(dat$A)) {
        return(invisible(NULL))
    }
    x <- if (foc_only) dat$cs[dat$ifp] else dat$cs
    Y <- matrix(nrow = length(x), ncol = length(dat$A))
    for (i in seq_along(dat$A)) {
        Y[, i] <- lc(x, dat$x_0[i], dat$A[i], dat$lambda[i])
        if (lc_show) {
            y <- Y[, i]
            near_zero <- abs(y) < 0.01
            y_big <- y[!near_zero]
            x_big <- x[!near_zero]
            lines(x = x_big, y = y_big, col = lc_col, lty = lc_lty, type = "l")
        }
    }
    if (sup_show) {
        lines(x = x, y = rowSums(Y), col = sup_col, lty = sup_lty, type = "l")
    }
    invisible(NULL)
}

ps_setup_dev_env <- function() {
    sim_01 <- metabodecon_file("sim/sim_01")
    decon <- generate_lorentz_curves(
        sim_01, sfr = c(3.42, 3.58), ws = 0, ask = FALSE,
        verbose = FALSE, delta = 0.1
    )
    args <- stub("plot_spec", decon = decon, foc_rgn = c(3.55, 3.52))
    invisible(args)
}
