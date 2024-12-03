# Plot Spectrum #####

#' @export
#'
#' @title Plot Spectrum
#'
#' @description
#' Plot a spectrum and zoom in on a specific region.
#' `r lifecycle::badge("experimental")`
#'
#' @param obj
#' An object of type spectrum, decon0, decon1 or decon2. For details see
#' [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Metabodecon-Classes.html).
#'
#' @param foc_rgn
#' A numeric vector specifying the start and end of the focus region in ppm. If
#' set to NULL, `foc_frac` is used to determine the focus region. If both
#' `foc_rgn` and are set to NULL, [get_foc_rgn()] is used to determine a
#' suitable focus region automatically. Takes precedence over `foc_frac`.
#'
#' @param foc_frac
#' A numeric vector specifying the start and end of the focus region as fraction
#' of the full spectrum width. Only used if `foc_rgn` is set to NULL.
#'
#' @param cnct_show
#' Logical. If TRUE, connecting lines are drawn between sub figure 1 and the
#' focus rectangle in sub figure 3.  See 'Details'.
#'
#' @param cnct_col
#' Color of the lines connecting the main and sub figure.
#'
#' @param layout
#' A list with three elements. Each element should be a numeric vector of the
#' form `c(x1, x2, y1, y2)`. The i'th element gives the borders of the i'th sub
#' figure in "normalized plot coordinates". Setting an element to NULL will
#' prevent the corresponding sub figure from being drawn. For a description of
#' the individual sub figures see 'Details'. For a description of normalized
#' plot coordinates see [graphics::grconvertX()]. If `layout` is set to `NULL`,
#' [get_ps_layout()] is used to determine a suitable layout automatically.
#'
#' @param args1,args2,args3
#' List of arguments passed to [draw_spectrum()] when drawing sub figure
#' 1-3. See 'Details'.
#'
#' @return
#' NULL. Called for side effect of plotting as sketched in 'Details'.
#'
#' @details
#' This function first calls [plot_empty()] to initalize a new plotting canvas.
#' After that it calls [draw_spectrum()] four times to draw the following sub
#' figures onto the plotting canvas:
#'
#' 1. The signal intensities in the focus region
#' 2. The second derivative in the focus region
#' 3. The signal intensities over all datapoints
#'
#' The argument lists for the individual calls to [draw_spectrum()] are
#' determined at runtime and depend on the arguments passed to [plot_spectrum()]
#' as well as the currently active graphics device. To customize the appearance
#' of the individual sub plots, you can overwrite each value passed to
#' [draw_spectrum()] by providing a corresponding named element in `args1`,
#' `args2` or `args3`.
#'
#' A sketch of the resulting figure is shown below.
#'
#' ```
#'  _________________________________________________________
#' | Plotting canvas                                         |
#' |        ________________________________________         |
#' |       | Sub1: Signal Intensity in Focus Region  |       |
#' |       |                   /\                    |       |
#' |       |                  /  \                   |       |
#' |       |                 /    \  /\              |       |
#' |       |                /      \/  \             |       |
#' |       |           /\  /            \            |       |
#' |       |          /  \/              \           |       |
#' |       |         /                    \          |       |
#' |       |________/______________________\_________|       |
#' |        _________________________________________        |
#' |       | Sub2: Second Derivative in Focus Region |       |
#' |       |_________________________________________|       |
#' |   ___________________________________________________   |
#' |  | Sub3: Signal Intensity over all Datapoints        |  |
#' |  |                         ________________          |  |
#' |  |                        |Focus Rectangle |         |  |
#' |  |          /\            |       /\       |         |  |
#' |  |         /  \           |      /  \/\    |         |  |
#' |  |        /    \     /\   |   /\/      \   |         |  |
#' |  |_______/______\___/__\__|__/__________\__|_________|  |
#' |_________________________________________________________|
#' ```
#'
#' @examples
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' ## Prepare a deconvoluted spectrum as input
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' sim_01 <- metabodecon_file("sim/sim_01")
#' spec <- read_spectrum(sim_01)
#' obj <- generate_lorentz_curves_sim(spec)
#' opar <- par(mfrow = c(2, 2))
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
#'     sub_pos = c(x1 = 0.1, x2 = 0.9, y1 = 0.4, y2 = 0.9)
#' )
#' plot_spectrum(obj, rct_fill = NULL, sub_fill_col = NULL, cnct_col = NULL)
#' plot_spectrum(obj, xlab = "", ylab = "", mar = c(2, 2, 0, 1))
plot_spectrum <- function(obj,
                          foc_rgn = NULL,
                          foc_frac = NULL,
                          cnct_show = TRUE,
                          cnct_col = "black",
                          mar = c(4.1, 4.1, 0.1, 0.1),
                          layout = NULL,
                          args1 = list(),
                          args2 = list(),
                          args3 = list())
{
    obj <- as_v2_obj(obj)
    foc_frac <- foc_frac %||% get_foc_frac(obj, foc_rgn)
    foc_rgn <- foc_rgn %||% get_foc_rgn(obj, foc_frac)
    layout <- layout %||% get_ps_layout(obj, foc_rgn)
    local_par(mar = mar)
    plot_empty()
    args <- get_ds_arglists(obj, foc_rgn, foc_frac, layout, args1, args2, args3)
    sub1 <- if (!is.null(layout[[1]])) do.call(draw_spectrum, args[[1]])
    sub2 <- if (!is.null(layout[[2]])) do.call(draw_spectrum, args[[2]])
    sub3 <- if (!is.null(layout[[3]])) do.call(draw_spectrum, args[[3]])
    cnct <- if (cnct_show) draw_connection_lines(sub3, sub_plot, cnct_col)
    invisible(named(args, sub1, sub2, sub3, cnct))
}

get_ps_layout <- function(obj, foc_rgn) {
    foc_only_layout <- list(
        sub1 = c(x1 = 0.00, x2 = 1.00, y1 = 0.25, y2 = 1.00),
        sub2 = c(x1 = 0.00, x2 = 1.00, y1 = 0.00, y2 = 0.25),
        sub3 = NULL
    )
    zoom_layout <- list(
        sub1 = c(x1 = 0.05, x2 = 0.95, y1 = 0.50, y2 = 0.95),
        sub2 = c(x1 = 0.05, x2 = 0.95, y1 = 0.35, y2 = 0.50),
        sub3 = c(x1 = 0.00, x2 = 1.00, y1 = 0.00, y2 = 0.50)
    )
    if (width(foc_rgn) < width(obj$cs)) {
        zoom_layout
    } else {
        foc_only_layout
    }
}

get_ds_arglists <- function(obj, foc_rgn, foc_frac, layout, args1, args2, args3) {
    defaults <- named(obj, foc_rgn, foc_frac)
    sub1 <- list(
        bg_fill = if (is.null(layout[[3]])) "white" else transp("yellow"),
        fig = npc_to_ndc(layout[[1]])
    )
    sub2 <- list(
        d2_show = TRUE,
        ylab = "Second Derivative [au]",
        fig = npc_to_ndc(layout[[2]])
    )
    sub3 <- list(
        fig = npc_to_ndc(layout[[3]])
    )
    arglists <- mapply(modifyList, list(defaults), list(sub1, sub2, sub3))
    arglists <- mapply(modifyList, arglists, list(sub1, sub2, sub3))
    names(arglists) <- c("sub1", "sub2", "sub3")
    arglists
}

get_foc_frac <- function(obj, foc_rgn = NULL) {
    stopifnot(is_num(obj$cs))
    stopifnot(is.null(foc_rgn) || (is_num(foc_rgn, 2)))
    if (is.null(foc_rgn)) {
        n <- length(obj$cs)
        width <- min(256 / n, 0.5)
        center <- if (n <= 2048) 0.5 else 0.75
        c(center - width, center + width)
    } else {
        convert_pos(foc_rgn, obj$cs, c(0, 1))
    }
}

get_foc_rgn <- function(obj, foc_frac = NULL) {
    stopifnot(is_num(obj$cs))
    stopifnot(is.null(foc_frac) || (is_num(foc_frac, 2)))
    if (is.null(foc_frac)) foc_frac <- get_foc_frac(obj)
    quantile(obj$cs, foc_frac)
}



# Draw Spectrum #####

#' @export
#' @title Draw Spectrum
#' @description
#' Draws a single spectrum.
#' Should not be used directly, instead use [plot_spectrum()].
#' For usage examples see [test/testthat/test-draw_spectrum.R](https://github.com/spang-lab/metabodecon/blob/main/test/testthat/test-draw_spectrum.R).
#' @param obj An object of type `spectrum` or `decon2`. For details see [Metabodecon Classes](https://spang-lab.github.io/metabodecon/articles/Metabodecon-Classes.html).
#' @param add If TRUE, draw into the currently open figure. If FALSE, start a new figure.
#' @param axis_col Color of tickmarks and ticklabels.
#' @param box_col Color of box surrounding plot region.
#' @param fig Drawing region in normalized device coordinates.
#' @param bg_fill Background color of plot region.
#' @param foc_only Logical. If TRUE, only the focused region is drawn. If FALSE, the full spectrum is drawn.
#' @param foc_rgn Numeric vector specifying the start and end of focus region.
#' @param foc_unit Character string specifying the unit in which `foc_rgn` is given. Can be "fraction" or "ppm".
#' @param lc_col Color of the individual Lorentzian Curves.
#' @param lc_fill Color of the rectangles shown at center of each lorentzian curve. The width of the rectangle equals the half width at half height.
#' @param lc_lty Line type of the individual Lorentzian Curves.
#' @param lc_show If TRUE, the individual Lorentzian Curves (LCs) are shown.
#' @param lgd If TRUE, a legend is drawn.
#' @param line_col Color of raw signal intensities.
#' @param main Title of the plot.
#' @param mar Number of lines below/left-of/above/right-of plot region.
#' @param rct_col Border color of the rectangle around the focus region.
#' @param rct_fill Background color of the rectangle around the focus region.
#' @param rct_show Logical. If TRUE, the focus region is shown as a rectangle.
#' @param sf_si_raw Numeric value to divide raw signal intensities by before drawing.
#' @param sm_col Color of smoothed signal intensities.
#' @param sm_show If TRUE, smoothed signal intensities are shown.
#' @param sup_col Color of the superposition of the Lorentzian Curves.
#' @param sup_lty Line type of the superposition of the Lorentzian Curves.
#' @param sup_show Logical. If TRUE, the superposition of the Lorentzian Curves is shown.
#' @param dp_col Vector of length 4 giving the color for each datapoint. The first three entries specify the color used for peak-centers, left-borders and right-borders. The fourth entry specifies the color used for non-peak data points.
#' @param dp_pch Vector of length 4 giving the plotting character for each datapoint.
#' @param dp_show Logical. If TRUE, individual point characters are drawn for each datapoint.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
draw_spectrum <- function(obj,
                          # Focus Region
                          foc_rgn   = NULL,
                          foc_frac  = NULL,
                          foc_only  = TRUE,
                          # Positioning
                          add = FALSE,
                          fig = NULL,
                          # Annotations
                          main = "",
                          xlab = "Chemical Shift [ppm]",
                          ylab = "Signal Intensity [au]",
                          lgd  = FALSE,
                          mar  = c(4.1, 4.1, 0.1, 0.1),
                          # Frame
                          axis_col = "black",
                          box_col  = "black",
                          bg_fill  = "white",
                          # Element config
                          si_line    = list(),
                          sm_line    = list(col = "blue"),
                          lc_lines   = list(col = "grey"),
                          lc_rect    = list(fill = transp("black")),
                          sup_line   = list(col = "red"),
                          d2_line    = NULL,
                          foc_rect   = list(fill = transp("yellow")),
                          peak_pts   = list(col = "blue", pch = 17),
                          border_pts = list(col = "blue", pch = 4),
                          normal_pts = NULL,
                          x0_vlines  = list())
{
    # Check inputs
    if (is.null(obj$cs)) stop("Chemical shifts missing")
    if (is.null(obj$si)) stop("Signal Intensities missing")
    obj <- as_v2_obj(obj)
    foc_frac <- foc_frac %||% get_foc_frac(obj, foc_rgn)
    foc_rgn <- foc_rgn %||% get_foc_rgn(obj, foc_frac)

    # Set graphical parameters
    local_par(mar = mar, new = add)
    local_fig(fig = fig, add = add)

    # Get lorentz parameters
    lcpar <- obj$lcpar
    x0 <- lcpar$x0
    A <- lcpar$A
    lambda <- lcpar$lambda

    # Get xy values over all data points
    cs_all  <- obj$cs
    si_all  <- obj$si
    sm_all  <- obj$sit$sm
    d2_all  <- d2_line %&&% calc_second_derivative(si_all)
    sup_all <- sup_line %&&% obj$sit$sup

    # Get some index vectors
    ifr <- which(cs_all >= min(foc_rgn) & cs_all <= max(foc_rgn)) # Focus points
    idp <- seq_along(cs_all) # All points
    inp <- setdiff(idp, ipp) # Non-peak points
    ipc <- obj$peak$center # Peak centers
    ipl <- obj$peak$left # Peak borders left
    ipr <- obj$peak$right # Peak borders right
    ipp <- c(ipc, ipl, ipr) # Peak points

    # Get xy values over the focus region
    cs_foc  <- cs_all[ifr]
    si_foc  <- si_all[ifr]
    sm_foc  <- sm_all[ifr]
    d2_foc  <- d2_all[ifr]
    sup_foc <- sup_all[ifr]

    # Get x and y limits
    cs <- if (foc_only) cs_foc else cs_all
    si <- if (foc_only) si_foc else si_all
    sm <- if (foc_only) sm_foc else sm_all
    d2 <- if (foc_only) d2_foc else d2_all
    sup <- if (foc_only) sup_foc else sup_all
    xlim <- c(max(cs), min(cs))
    ylim <- c(min(si, sm, d2), max(si, sm, d2))

    # Do the actual drawing
    draw_bg(xlim, ylim, bg_fill)
    draw_line(cs, si,  si_line)
    draw_line(cs, sm,  sm_line)
    draw_line(cs, sup, sup_line)
    draw_line(cs, d2,  d2_line)
    apply(lcpar, 1, draw_lc_line, cs, lc_lines)

    # CONTINUE HERE: calculate index vectors correctly
    draw_points(cs[ipc], y[ipc], peak_pts)
    draw_points(cs[ipb], y[ipb], border_pts)
    draw_points(cs[ipr], y[ipr], normal_pts)
    draw_foc_rect(env)
    draw_axis(env)
    draw_legend(env)
    ndc <- list(bg = c(), foc_rect = c())
    ndc
}

draw_line <- function(func, x, y, args) {
    if (is.null(args)) return(NULL)
    if (!is.list(args)) stop("args must a list or NULL")
    do.call(lines, c(list(x, y), args))
}

draw_points <- function(func, x, y, args) {
    if (is.null(args)) return(NULL)
    if (!is.list(args)) stop("args must a list or NULL")
    do.call(points, c(list(x, y), args))
}

draw_bg <- function(xlim, ylim, bgfill) {
    plot(xlim, ylim, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
         xlab = "", ylab = "", type = "n", axes = FALSE)
    if (is.null(bg_fill)) return(NULL)
    rect(xleft = xlim[1], xright = xlim[2],
         ybottom = ylim[1], ytop = ylim[2],
         col = bg_fill, border = NA)
}

draw_lc_line <- function(p, x, args) {
    if (is.null(args)) return(NULL)
    if (!is.list(args)) stop("args must a list or NULL")
    y <- lorentz(x, p[["x0"]], p[["A"]], p[["lambda"]])
    near_zero <- abs(y) < min(ylim) + 0.001 * diff(ylim)
    y <- y[!near_zero]
    x <- x[!near_zero]
    fnargs <- c(list(x, y), args)
    do.call(lines, fnargs)
}

draw_lc_rect <- function(x, x_0, A, lambda, lc_fill = NULL) {
    if (is.null(lc_fill)) return(y)
    rect(xleft = x_0 + lambda,
         xright = x_0 - lambda,
         ybottom = par("usr")[3],
         ytop = lorentz(x_0, x_0, A, lambda),
         col = lc_fill,
         border = NA
    )
    return(y)
}

draw_foc_rect <- function(dat,
                          show = TRUE,
                          fill = rgb(0, 0, 1, alpha = 0.1),
                          col = "blue") {
    if (is.null(dat$foc_rgn) || !show) return(NULL) # styler: off
    usr <- list(
        xleft = max(dat$foc_rgn),
        xright = min(dat$foc_rgn),
        ybottom = 0,
        ytop = max(dat$si[dat$ifr]) / 0.96,
        col = fill,
        border = col
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

draw_axis <- function(main, xlim, ylim, xlab, ylab, axis_col, box_col) {
    xtks <- seq(xlim[1], xlim[2], length.out = 5)
    ytks <- seq(ylim[1], ylim[2], length.out = 5)
    for (i in 2:12) if (length(unique(xtklabs <- round(xtks, i))) >= 5) break
    for (i in 2:12) if (length(unique(ytklabs <- round(ytks, i))) >= 5) break
    axis(side = 1, at = xtks, labels = xtklabs,
         col = NA, col.ticks = axis_col, col.axis = axis_col)
    axis(side = 2, at = ytks, labels = ytklabs,
         col = NA, col.ticks = axis_col, col.axis = axis_col, las = 1)
    title(main = main, xlab = xlab, ylab = ylab)
    box(col = box_col)
}

draw_legend <- function(dat, args, verbose = FALSE) {
    if (isFALSE(args$lgd)) return(invisible(NULL))
    dsc <- c("Raw Signal", "Smoothed Signal", "Single Lorentzian", "Sum of Lorentzians", "Peak Triplets")
    col <- c(args$line_col, args$sm_col, args$lc_col, args$sup_col, args$dp_col[1])
    lty <- c(1, 1, 1, 1, NA)
    pch <- c(NA, NA, NA, NA, args$dp_pch[1])
    keep <- c(
        !is.null(dat$si),
        args$sm_show && !is.null(dat$sm),
        args$lc_show && !is.null(dat$A),
        args$sup_show && !is.null(dat$A),
        args$dp_show && !is.null(dat$ipc)
    )
    df <- data.frame(dsc, col, lty, pch)[keep,]
    legend(x = "topright", legend = df$dsc, col = df$col, lty = df$lty, pch = df$pch)
}

draw_connection_lines <- function(main_plot, sub_plot, cnct_col) {
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

# Helpers #####

usr_to_ndc <- function(usr = par("usr")) {
    x_ndc <- grconvertX(usr[1:2], from = "user", to = "ndc")
    y_ndc <- grconvertY(usr[3:4], from = "user", to = "ndc")
    c(x_ndc, y_ndc)
}

npc_to_ndc <- function(npc = c(0, 1, 0, 1)) {
    x_ndc <- grconvertX(npc[1:2], from = "npc", to = "ndc")
    y_ndc <- grconvertY(npc[3:4], from = "npc", to = "ndc")
    c(x_ndc, y_ndc)
}

