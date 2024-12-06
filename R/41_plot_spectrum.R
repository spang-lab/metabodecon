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

# Draw #####

#' @export
#' @title Draw Spectrum
#' @description
#' Draws a single spectrum. Internally used by [plot_spectrum()], which is
#' usually the recommended way to plot spectra. For usage examples see
#' [test/testthat/test-draw_spectrum.R](https://github.com/spang-lab/metabodecon/blob/main/test/testthat/test-draw_spectrum.R).
#'
#' @param obj
#' An object of type `spectrum` or `decon2`. For details see [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Metabodecon-Classes.html).
#'
#' @param add
#' If TRUE, draw into the currently open figure. If FALSE, start a new figure.
#'
#' @param fig
#' Drawing region in normalized device coordinates.
#'
#' @param foc_only
#' Logical. If TRUE, only the focused region is drawn. If FALSE, the full
#' spectrum is drawn.
#'
#' @param foc_rgn
#' Numeric vector specifying the start and end of focus region.
#'
#' @param foc_frac
#' Numeric vector specifying the start and end of focus region as fraction of
#' the full spectrum width.
#'
#' @param mar
#' Number of lines below/left-of/above/right-of plot region.
#'
#' @param lgd
#' List of parameters passed to [legend()] when drawing the legend.
#'
#' @param main
#' Main title of the plot. Drawn via [title()].
#'
#' @param si_line,sm_line,lc_lines,sp_line,d2_line
#' List of parameters passed to [lines()] when drawing the raw signal
#' intensities, smoothed signal intensities, lorentzian curves, superposition of
#' lorentzian curves and second derivative.
#'
#' @param cent_pts,bord_pts,norm_pts
#' List of parameters passed to [points()] when drawing the peak center points,
#' peak border points and non-peak points.
#'
#' @param bg_rect,lc_rects,foc_rect
#' List of parameters passed to [rect()] when drawing the background, lorentzian
#' curves and focus rectangle.
#'
#' @param bt_axis,lt_axis,tp_axis,rt_axis
#' List of parameters used to overwrite the [draw_spectrum_defaults()]. The
#' provided parameters are passed to [axis()] when drawing the bottom, left, top
#' and right axis. In addition to the parameters of [axis()], the following
#' additional parameters are supported as well:
#' - `text`: Description for the axis. Drawn via [mtext()].
#' - `n`: Number of tickmarks.
#' - `digits`: Number of digits for rounding the labels. If a vector of numbers
#'   is provided, all numbers are tried, until `n` unique labels are found. See
#'   'Details'.
#' - `sf`: Scaling factor. Axis values are divided by this number before the
#'   labels are calculated. If you set this to anything unequal 1, you should
#'   also choose `text` in a way that reflects the scaling. E.g. if you set `sf
#'   = 1e6` you could change the text from "Signal Intensity [au]" to "Signal
#'   Intensity [Mau]" or "Signal Intensity [au] / 1e6", with "Mau" meaning
#'   "Mega-Arbitrary-Units".
#'
#' @details
#' Parameters `bt_axis`, `lt_axis`, `tp_axis` and `rt_axis` all support option
#' `n` and `digits`, where `n = 5` means "Draw 5 tickmarks over the full axis
#' range" and `digits = 3` means "round the label shown beside each tickmark to
#' 3 digits". If `n` is omitted, a suiteable value is chosen automatically using
#' [axTicks()]. If `digits` is omitted, a default of `2:12` is used. Providing a
#' vector of `digits` causes each digit to be tried as argument for [round()],
#' until a digit is encountered that results in `n` unique labels. Example:
#'
#' Assume we have `n = 4` and the corresponding calculated tickmark positions
#' are: 1.02421, 1.02542, 1.02663 and 1.02784. If we provide `digits = 1:5`, the
#' following roundings are tried:
#'
#' | digit  | label 1 | label 2 | label 3 | label 4 |
#' | ------ | ------- | ------- | ------- | ------- |
#' | 1      | 1       | 1       | 1       | 1       |
#' | 2      | 1.02    | 1.03    | 1.03    | 1.03    |
#' | 3      | 1.024   | 1.025   | 1.027   | 1.028   |
#' | 4      | 1.0242  | 1.0254  | 1.0266  | 1.0278  |
#' | 5      | 1.02421 | 1.02542 | 1.02663 | 1.02784 |
#'
#' In the above example the process would stop at `digit = 3`, because at this
#' point we have n = 4 unique labels (1.024, 1.025, 1.027 and 1.028).
draw_spectrum <- function(obj,
    foc_rgn  = NULL,   foc_frac = NULL,   foc_only = TRUE,
    add      = FALSE,  fig      = NULL,   mar      = c(4.1, 4.1, 0.1, 0.1),
    main     = NULL,   show_d2  = FALSE,  lgd      = list(),
    si_line  = list(), sm_line  = list(), lc_lines = list(),
    sp_line  = list(), d2_line  = list(),
    cent_pts = list(), bord_pts = list(), norm_pts = list(),
    bg_rect  = list(), foc_rect = list(), lc_rects = list(),
    bt_axis  = list(), lt_axis  = list(), tp_axis  = list(),
    rt_axis  = list())
{
    # Check inputs (278us)
    obj <- as_v2_obj(obj)
    eval(draw_spectrum_check_arguments)
    env <- environment()
    defaults <- get_draw_spectrum_defaults(show_d2, foc_only)
    for (x in names(defaults)) env[[x]] <- modifyList(defaults[[x]], env[[x]])
    foc_frac <- foc_frac %||% get_foc_frac(obj, foc_rgn)
    foc_rgn <- foc_rgn %||% get_foc_rgn(obj, foc_frac)

    # Set graphical parameters (7ms)
    local_par(mar = mar, new = add)
    local_fig(fig = fig, add = add)

    # Get xy values over all data points (13us)
    cs <- cs_all <- obj$cs
    si <- si_all <- obj$si
    sm <- sm_all <- obj$sit$sm
    d2 <- d2_all <- if (isFALSE(d2_line$show)) NULL else calc_second_derivative(si_all)
    sp <- sp_all <- obj$sit$sup

    # Get indices of important points relative all data points (611us)
    idp <- idp_all <- seq_along(cs_all) # Data points
    icp <- icp_all <- obj$peak$center # Center points
    ibp <- ibp_all <- unique(sort(c(obj$peak$left, obj$peak$right))) # Border points
    ipp <- ipp_all <- sort(c(icp_all, ibp_all)) # Peak points
    inp <- inp_all <- setdiff(idp_all, ipp_all) # Nonpeak points
    ifp <- ifp_all <- which(cs_all >= min(foc_rgn) & cs_all <= max(foc_rgn)) # Focus points

    # Get lorentz parameters across all data points (8us)
    lcpar  <- lcpar_all  <- obj$lcpar
    x0     <- x0_all     <- lcpar_all$x0
    A      <- A_all      <- lcpar_all$A
    lambda <- lambda_all <- lcpar_all$lambda

    # Get indices of important points relative to the vector of focus points (198us)
    if (foc_only) {
        offset <- min(ifp_all) - 1
        idp <- idp_foc <- ifp_all - offset # Data points
        icp <- icp_foc <- intersect(icp_all, ifp_all) - offset # Center points
        ibp <- ibp_foc <- intersect(ibp_all, ifp_all) - offset # Border points
        ipp <- ipp_foc <- intersect(ipp_all, ifp_all) - offset # Peak points
        inp <- inp_foc <- intersect(inp_all, ifp_all) - offset # Nonpeak points

        # Get xy values over the focus region (18us)
        cs <- cs_foc <- cs_all[ifp_all]
        si <- si_foc <- si_all[ifp_all]
        sm <- sm_foc <- sm_all[ifp_all]
        d2 <- d2_foc <- d2_all[ifp_all]
        sp <- sp_foc <- sp_all[ifp_all]

        # Get lorentzians affecting focus region (239us)
        y_foc_rgn_start   <- lorentz(x = min(foc_rgn), lcpar = lcpar_all)
        y_foc_rgn_end     <- lorentz(x = max(foc_rgn), lcpar = lcpar_all)
        y_tresh           <- 0.001 * diff(range(si))
        high_y_in_foc_rgn <- y_foc_rgn_start > y_tresh | y_foc_rgn_end > y_tresh
        x0_in_foc         <- x0_all > min(foc_rgn) & x0_all < max(foc_rgn)
        affects_foc       <- high_y_in_foc_rgn | x0_in_foc
        lcpar  <- lcpar_foc  <- lcpar_all[affects_foc, ]
        x0     <- x0_foc     <- x0_all[affects_foc]
        A      <- A_foc      <- A_all[affects_foc]
        lambda <- lambda_foc <- lambda_all[affects_foc]
    }

    # Get x and y limits (43us)
    xlim <- c(max(cs), min(cs))
    ymin <- if (show_d2) min(d2) else min(0, si)
    ymax <- if (show_d2) max(d2) else max(si)
    ylim <- c(ymin, ymax)
    ythresh <- 0.001 * diff(range(si))
    if (!is.null(foc_rgn)) {
        foc_rgn <- sort(foc_rgn)
        ymin_foc <- min(si[ifp])
        ymax_foc <- max(si[ifp])
    }
    ymin_foc <- if (foc_only) ymin else min(0, si_all[ifp_all])
    ymax_foc <- if (foc_only) ymax else min(0, si_all[ifp_all])
    ylim_foc <- c(ymin_foc, ymax_foc)

    # Do the actual drawing
    plot_empty(xlim, ylim, main = main)
    draw_rect(xlim, ylim, bg_rect)
    apply(lcpar_foc, 1, draw_lc_rect, cs, lc_rects, ymin) # 3ms
    apply(lcpar_foc, 1, draw_lc_line, cs, lc_lines, ythresh) # 13ms
    draw_line(cs, si, si_line)
    draw_line(cs, sm, sm_line)
    draw_points(cs[icp], sm[icp], cent_pts)
    draw_points(cs[ibp], sm[ibp], bord_pts)
    draw_points(cs[inp], sm[inp], norm_pts)
    draw_line(cs, sp, sp_line)
    draw_line(cs, d2, d2_line)
    draw_rect(foc_rgn, ylim_foc, foc_rect)
    draw_axis(xlim, side = 1, bt_axis)
    draw_axis(ylim, side = 2, lt_axis)
    draw_axis(xlim, side = 3, tp_axis)
    draw_axis(ylim, side = 4, rt_axis)
    draw_legend(lgd, si_line, sm_line, lc_lines, sp_line, d2_line, cent_pts, bord_pts, norm_pts)
    usr_to_ndc(c(foc_rgn, ylim))
}

draw_legend <- function(lgd, si_line, sm_line, lc_lines, sp_line, d2_line,
                        cent_pts, bord_pts, norm_pts) {
    if (isFALSE(pop(lgd, "show"))) return()
    dsc <- c("Raw Signal", "Smoothed Signal", "Single Lorentzian",
              "Sum of Lorentzians", "Second Derivative",
              "Peak Center", "Peak Border", "Non-Peak")
    col <- c(si_line$col %||% "black",  sm_line$col  %||% "black",
             lc_lines$col %||% "black", sp_line$col  %||% "black",
             d2_line$col %||% "black",  cent_pts$col %||% "black",
             bord_pts$col %||% "black", norm_pts$col %||% "black")
    lty <- c(1, 1, 1, 1, 1, NA, NA, NA)
    pch <- c(NA, NA, NA, NA, NA, cent_pts$pch %||% NA,
             bord_pts$pch %||% NA, norm_pts$pch %||% NA)
    keep <- c(!isFALSE(si_line$show), !isFALSE(sm_line$show),
              !isFALSE(lc_lines$show), !isFALSE(sp_line$show),
              !isFALSE(d2_line$show), !isFALSE(cent_pts$show),
              !isFALSE(bord_pts$show), !isFALSE(norm_pts$show))
    df <- data.frame(dsc, col, lty, pch)[keep, ]
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

draw_lc_line <- function(p, x, args = list(), threshold = 0) {
    if (isFALSE(pop(args, "show"))) return()
    y <- lorentz(x, p[["x0"]], p[["A"]], p[["lambda"]])
    if (threshold > 0) {
        large_enough <- abs(y) >= threshold
        y <- y[large_enough]
        x <- x[large_enough]
    }
    args$x <- x
    args$y <- y
    do.call(lines, args)
}

draw_lc_rect <- function(p, x, args = list(), ymin = 0) {
    if (isFALSE(pop(args, "show"))) return()
    A <- p[["A"]]
    x0 <- p[["x0"]]
    lambda <- p[["lambda"]]
    args$xleft <- x0 - lambda
    args$xright <- x0 + lambda
    args$ybottom <- ymin
    args$ytop <- A / lambda
    do.call(rect, args)
}

draw_axis <- function(lim, side = 1, args = list()) {
    if (isFALSE(pop(args, "show"))) return()
    n <- pop(args, "n", default = 5)
    digits <- pop(args, "digits", default = 2:12)
    sf <- pop(args, "sf", default = 1)
    text <- pop(args, "text")
    args$at <- if (is_num(n)) seq(min(lim), max(lim), length.out = n) else axTicks(side)
    labs_sc <- args$at / sf
    args$labels <- format(labs_sc, digits = 7) # default used by normal axis
    for (i in digits) {
        args$labels <- round(labs_sc, digits = i)
        if (length(unique(args$labels) >= 5)) break
    }
    if (!is.null(text)) mtext(text, side = side, line = 3)
    args$side <- side
    do.call(axis, args)
}

draw_rect <- function(xlim, ylim, args = list()) {
    if (isFALSE(pop(args, "show"))) return()
    args$xleft <- xlim[1]
    args$xright <- xlim[2]
    args$ybottom <- ylim[1]
    args$ytop <- ylim[2]
    do.call(rect, args)
}

draw_line <- function(x, y, args = list()) {
    if (isFALSE(pop(args, "show"))) return()
    args$x <- x
    args$y <- y
    do.call(lines, args)
}

draw_points <- function(x, y, args = list()) {
    if (isFALSE(pop(args, "show"))) return()
    args$x <- x
    args$y <- y
    do.call(points, args)
}

# Get #####

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
    defaults <- list(
        obj = obj,
        foc_rgn = foc_rgn,
        foc_frac = foc_frac,
        add = TRUE
    )
    sub1 <- list(
        fig = npc_to_ndc(layout[[1]])
    )
    sub2 <- list(
        fig = npc_to_ndc(layout[[2]]),
        bt_axis  = list(text = "Chemical Shift [ppm]"),
        lt_axis  = list(text = "Second Derivative [Mau]", las = 1, sf = 1e6),
        tp_axis  = list(skip = TRUE),
        rt_axis  = list(skip = TRUE),
        bg_rect  = list(skip = TRUE),
        si_line  = list(skip = TRUE),
        sm_line  = list(skip = TRUE),
        lc_lines = list(skip = TRUE),
        lc_rects = list(skip = TRUE),
        sp_line  = list(skip = TRUE),
        d2_line  = list(),
        foc_rect = list(skip = TRUE),
        cent_pts = list(skip = TRUE),
        bord_pts = list(skip = TRUE),
        norm_pts = list(skip = TRUE)
    )
    sub3 <- list(
        fig = npc_to_ndc(layout[[3]]),
        sm_line  = list(skip = TRUE),
        lc_lines = list(skip = TRUE),
        lc_rects = list(skip = TRUE),
        sp_line  = list(skip = TRUE),
        foc_rect = list(skip = TRUE),
        cent_pts = list(skip = TRUE),
        bord_pts = list(skip = TRUE),
        norm_pts = list(skip = TRUE)
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
        convert_pos(foc_rgn, range(obj$cs), c(0, 1))
    }
}

get_foc_rgn <- function(obj, foc_frac = NULL) {
    stopifnot(is_num(obj$cs))
    stopifnot(is.null(foc_frac) || (is_num(foc_frac, 2)))
    if (is.null(foc_frac)) foc_frac <- get_foc_frac(obj)
    quantile(obj$cs, foc_frac)
}

# Convert #####

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

get_draw_spectrum_defaults <- function(show_d2, foc_only = TRUE) {
    stopifnot(is_bool(show_d2), is_bool(foc_only))
    show_si  <- !show_d2
    ylab <- if (show_si) "Signal Intensity [Mau]" else"Second Derivative [au]"
    ysf <- if (show_si) 1e6 else 1
    list(
        lgd      = list(show = TRUE),
        d2_line  = list(show = show_d2),
        si_line  = list(show = show_si, col = "black"),
        sm_line  = list(show = show_si, col = "blue"),
        sp_line  = list(show = show_si, col = "red"),
        lc_lines = list(show = show_si && foc_only, col = "grey"),
        cent_pts = list(show = show_si && foc_only, col = "blue",  pch = 124),
        bord_pts = list(show = show_si && foc_only, col = "blue",  pch = 124),
        norm_pts = list(show = FALSE,               col = "black", pch = 124),
        bg_rect  = list(show = TRUE, col = "white"),
        foc_rect = list(show = TRUE, col = transp("yellow")),
        lc_rects = list(show = TRUE, col = transp("black"), border = NA),
        bt_axis  = list(show = TRUE, text = "Chemical Shift [ppm]"),
        lt_axis  = list(show = TRUE, text = ylab, sf = ysf, las  = 1),
        tp_axis  = list(show = FALSE),
        rt_axis  = list(show = FALSE)
    )
}
