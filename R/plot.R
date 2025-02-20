# Public #####

#' @export
#'
#' @title Plot Spectra
#'
#' @description
#' Plot a set of deconvoluted spectra.
#'
#' @param obj
#' An object of type `decons0`, `decons1` or `decons2`.
#' For details see [metabodecon_classes].
#'
#' @param ...
#' Additional arguments passed to the conversion function.
#'
#' @param sfy
#' Scaling factor for the y-axis.
#'
#' @param xlab
#' Label for the x-axis.
#'
#' @param ylab
#' Label for the y-axis.
#'
#' @param mar
#' A numeric vector of length 4, which specifies the margins of the plot.
#'
#' @return
#' A plot of the deconvoluted spectra.
#'
#' @seealso
#' [plot_spectrum()] for a much more sophisticated plotting routine
#' suitable for plotting a single spectrum.
#'
#' @examples
#' obj <- deconvolute(sim[1:4], sfr = c(3.55, 3.35))
#' plot_spectra(obj)
plot_spectra <- function(obj,
                         ...,
                         sfy = 1e6,
                         xlab = "Chemical Shift [ppm]",
                         ylab = paste("Signal Intensity [au] /", sfy),
                         mar = c(4.1, 4.1, 1.1, 0.1)) {
    decons <- as_v2_objs(obj, ...)
    sis <- lapply(decons, function(x) x$sit$al %||% x$sit$sup %||% x$si)
    x0s <- lapply(decons, function(x) x$lcpar$x0)
    css <- lapply(decons, function(x) x$cs)
    si_min <- 0
    si_max <- max(sapply(sis, max))
    cs_min <- min(sapply(css, min))
    cs_max <- max(sapply(css, max))
    x0_min <- min(sapply(x0s, min))
    x0_max <- max(sapply(x0s, max))
    x0_width <- x0_max - x0_min
    x0_quart <- x0_width / 4
    line_colors <- rainbow(length(decons))
    legend_text <- paste("Spectrum", 1:length(decons))
    local_par(mar = mar)
    plot(
        x = NA,
        type = "n",
        xlab = xlab,
        ylab = ylab,
        xlim = c(cs_max, cs_min),
        ylim = c(si_min, si_max)
    )
    abline(
        v = c(x0_min, x0_max),
        lty = 2
    )
    for (i in seq_along(decons)) {
        lines(x = css[[i]], y = sis[[i]], col = line_colors[[i]])
    }
    arrows(
        x0 = c(x0_min + x0_quart, x0_max - x0_quart),
        x1 = c(x0_min, x0_max),
        y0 = si_max * 0.8,
        y1 = si_max * 0.8,
        length = 0.2,
        lty = 2,
        col = "black"
    )
    text(
        x = mean(c(x0_min, x0_max)),
        y = 0.8 * si_max,
        labels = "ppm range"
    )
    mtext(
        text = round(c(x0_min, x0_max), 4),
        side = 3,
        line = 0,
        at = c(x0_min, x0_max)
    )
    legend(
        x = "topright",
        legend = legend_text,
        col = line_colors,
        lty = 1
    )
}

#' @export
#'
#' @title Plot Spectrum
#'
#' @description
#' Plot a spectrum and zoom in on a specific region.
#' `r lifecycle::badge("experimental")`
#'
#' @param x
#' An object of type `spectrum`, `decon0`, `decon1` or `decon2`. For details see
#' [metabodecon_classes].
#'
#' @param ...
#' Additional arguments passed to [draw_spectrum()] for **every** sub figure.
#' See 'Details'.
#'
#' @param obj
#' An object of type `spectrum` or `decon2`. Usually auto  generated  from  `x`,
#' but can be set manually in case the default conversion is not sufficient.
#'
#' @param foc_rgn
#' A numeric vector specifying the start and end of the focus region in ppm.  If
#' set to NULL, `foc_frac` is used  to  determine  the  focus  region.  If  both
#' `foc_rgn`  and  are  set  to  NULL,  a  suitable  focus  region  is   chosen
#' automatically. Takes precedence over `foc_frac`.
#'
#' @param foc_frac
#' A numeric vector specifying the start and end of the focus region as fraction
#' of the full spectrum width. Only used if `foc_rgn` is set to NULL.
#'
#' @param sub1,sub2,sub3
#' List of arguments passed to [draw_spectrum()] when drawing sub figure
#' 1-3. See 'Details'.
#'
#' @param mar
#' A numeric vector of length 4 passed, which specifies the margins of the plot.
#' Passed to [par()]. If set to `NULL`, a suitable value is chosen automatically.
#'
#' @param frame
#' A list of values passed to [box()] when drawing the frame around plot region.
#' If set to `NULL`, no frame is drawn.
#'
#' @param con_lines
#' A list of values passed to [lines()] when drawing the connecting lines between
#' sub figure 1 and the focus rectangle in sub figure 3. See 'Details'.
#' If set to `NULL`, the connecting lines are not drawn.
#'
#' @return
#' NULL. Called for side effect of plotting as sketched in 'Details'.
#'
#' @details
#' This function first initializes a new plotting canvas.  After  that  it  calls
#' [draw_spectrum()] multiple times to draw the following sub figures  onto  the
#' plotting canvas:
#'
#' 1. The signal intensities in the focus region
#' 2. The second derivative in the focus region
#' 3. The signal intensities over all datapoints
#'
#' The  argument  lists  for  the  individual  calls  to  [draw_spectrum()]  are
#' determined at runtime and depend on the arguments passed to [plot_spectrum()]
#' as well as the currently active graphics device. To customize the  appearance
#' of the  individual  sub  plots,  you  can  overwrite  each  value  passed  to
#' [draw_spectrum()] by providing  a  corresponding  named  element  in  `sub1`,
#' `sub2` or `sub3`.
#'
#' A sketch of the resulting figure is shown below.
#'
#' ```
#'  __________________________________________
#' |        ______________1_____________      |
#' |       | Sub1: Signal Intensity in  |     |
#' |       | Focus Region               |     |
#' |       |             /\             |     |
#' |       |            /  \            |     |
#' |       |           /    \  /\       |     |
#' |     11|          /      \/  \      |7    |
#' |       |     /\  /            \     |     |
#' |       |    /  \/              \    |     |
#' |       |   /                    \   |     |
#' |       |__/___________0__________\__|     |
#' |       | Sub2: Second Derivative    |     |
#' |     11| in Focus Region            |7    |
#' |       |____________________________|     |
#' |                      3                   |
#' |    __________________3_________________  |
#' |   |  Sub3: Signal Intensity over all   | |
#' |   |  Datapoints     ________________   | |
#' | 5 |                | Focus Rectangle|  |1|
#' |   |     /\         |       /\       |  | |
#' |   |    /  \        |      /  \/\    |  | |
#' |   |   /    \   /\  |   /\/      \   |  | |
#' |   |__/______\_/__\_|__/__________\__|__| |
#' |______________________5___________________|
#' ```
#'
#' Note  that  the  figure  created  by  `plot_spectrum()`  can  be  part  of  a
#' multi-figure configuration as created when setting  `mfrow`  or  `mfcol`  via
#' [par()]. Example:
#'
#' ```
#' _______________________________________
#'| Plot Spectrum with   | Other Figure  |
#'| sub3 = TRUE          | Other Figure  |
#'|      ___________     |  ___________  |
#'|     | Sub Fig 1 |    | | x      x  | |
#'|     |___________|    | |      x    | |
#'|     |_Sub_Fig_2_|    | |      x    | |
#'|   _________________  | |   x     x | |
#'|  |    Sub Fig 3    | | |      x    | |
#'|  |_________________| | |___________| |
#'|______________________|_______________|
#'| Some other Figure    | Plot Spectrum |
#'|                      | sub3 = FALSE  |
#'|  _________________   |  ___________  |
#'| |     ___         |  | | Sub Fig 1 | |
#'| | ___/   \___     |  | |           | |
#'| |/           \____|  | |___________| |
#'| |                 |  | | Sub Fig 2 | |
#'| |_________________|  | |___________| |
#'|______________________|_______________|
#' ```
#'
#' @examples
#' ## 1. Prepare a deconvoluted spectrum as input
#'
#' spec <- read_spectrum(metabodecon_file("sim/sim_01"))
#' decon <- generate_lorentz_curves_sim(spec)
#'
#' ## 2.1. Plot the full (non-deconvoluted) spectrum
#' ## 2.2. Remove connecting lines, and focus on a specific region specified in ppm
#' ## 2.3. Show second derivative and focus on a specific region specified as fraction
#' ## 2.4. Change color of focus rectangle and margins of sub figure 1
#' ## 2.5. Hide xlab and show second derivative
#' ## 2.6. Change the figure region for sub figure 1
#'
#' plot_spectrum(spec, sub1 = FALSE)
#' plot_spectrum(decon, foc_rgn = c(3.49, 3.45), con_lines = FALSE)
#' plot_spectrum(decon, sub2 = TRUE, foc_frac = c(0.40, 0.30))
#' plot_spectrum(decon,
#'     sub1 = list(mar = c(3, 6, 3, 6), lt_axis = list(col = "violet")),
#'     foc_rect = list(border = "violet", col = transp("violet")),
#'     con_lines = list(col = "violet")
#' )
#' plot_spectrum(decon,
#'     sub2 = TRUE,
#'     sub3 = list(bt_axis = list(text = "")),
#'     frame = TRUE,
#'     con_lines = FALSE
#' )
#' plot_spectrum(decon, sub1 = list(fig_rgn_npc = c(0,1,.3,1), mar = c(0,5,0,0)))
#'
plot_spectrum <- function(x,
                          ...,
                          obj = as_v2_obj(x),
                          foc_frac = get_foc_frac(obj),
                          foc_rgn = get_foc_rgn(obj, foc_frac),
                          sub1 = TRUE,
                          sub2 = FALSE,
                          sub3 = width(foc_rgn) < width(obj$cs),
                          mar  = NULL,
                          frame = FALSE,
                          con_lines = TRUE)
{
    # Parse Inputs
    dot_args <- if (identical(environment(), globalenv())) list() else list(...)
    sub1 <- combine(list(show = TRUE), sub1)
    sub2 <- combine(list(show = FALSE), sub2)
    sub3 <- combine(list(show = width(foc_rgn) < width(obj$cs)), sub3)
    frame <- combine(list(show = FALSE), frame)
    con_lines <- combine(list(show = TRUE), con_lines)

    # Setup Plotting Canvas
    default_mar_left <- if (sub3$show && any(sub1$show, sub2$show)) 3 else 6
    default_mar <- c(5, default_mar_left, 2, 2)
    local_par(mar = mar %||% default_mar)
    plot_empty()

    # Prepare Sub-Figure Arguments
    args <- get_sub_fig_args(obj, foc_frac, foc_rgn, sub1, sub2, sub3, dot_args)

    # Draw Sub-Figures
    fig1 <- do.call(draw_spectrum, args$sub1)
    fig2 <- do.call(draw_spectrum, args$sub2)
    fig3 <- do.call(draw_spectrum, args$sub3)
    figC <- draw_con_lines(fig1 %||% fig2, fig3, con_lines)
    figF <- draw_box(frame)

    invisible(NULL)
}

#' @export
#' @title Draw Spectrum
#' @description
#' Draws a single spectrum.  Internally  used  by  [plot_spectrum()],  which  is
#' usually  the  recommended  way  to  plot  spectra.  For  usage  examples  see
#' [test/testthat/test-draw_spectrum.R](https://github.com/spang-lab/metabodecon/blob/main/tests/testthat/test-draw_spectrum.R).
#'
#' @param obj
#' An object of type `spectrum` or `decon2`. For details see
#' [metabodecon_classes].
#'
#' @param add
#' If TRUE, draw into the currently open figure. If FALSE, start a new figure.
#'
#' @param fig_rgn
#' Drawing region in normalized device coordinates as vector of the form `c(x1,
#' x2, y1, y2)`.
#'
#' @param foc_only
#' Logical. If TRUE, only the focused region is drawn. If FALSE, the full
#' spectrum is drawn.
#'
#' @param foc_rgn
#' Numeric vector specifying the start and end of focus region in ppm.
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
#' @param show
#' Logical. If FALSE, the function returns without doing anything.
#'
#' @param show_d2
#' Logical. If TRUE, the second derivative of the spectrum is drawn. Setting
#' this to TRUE changes most of the defaults for the drawing, e.g. by disabling
#' the drawing of anything related to signal intensities and by changing the
#' y-axis label to "Second Derivative".
#'
#' @param truepar
#' Data frame with columns x0, A and lambda containing the true lorentzian that
#' were used to simulate the spectrum. Required if any `tp_*` argument is set.
#'
#' @param si_line,sm_line,lc_lines,sp_line,d2_line,tp_lines
#' List  of  parameters  passed  to  [lines()]  when  drawing  the  raw signal
#' intensities, smoothed signal intensities, lorentzian curves found by
#' deconvolution, superposition of lorentzian curves, second derivative and/or
#' true lorentzian curves.
#'
#' @param cent_pts,bord_pts,norm_pts
#' List of parameters passed to [points()] when drawing the peak center  points,
#' peak border points and non-peak points.
#'
#' @param bg_rect,lc_rects,foc_rect,tp_rects
#' List of parameters passed to [rect()] when drawing the background, lorentzian
#' curve substitutes, focus rectangle and/or true lorentzian curve substitutes.
#'
#' @param bt_axis,lt_axis,tp_axis,rt_axis
#' List of parameters used to overwrite the default values  passed  to  [axis()]
#' when drawing the bottom, left,  top  and  right  axis.  In  addition  to  the
#' parameters of [axis()], the following additional parameters are supported  as
#' well:
#'
#' - `text`:   Description for the axis. Drawn via [mtext()].
#' - `n`:      Number of tickmarks.
#' - `digits`: Number of digits for rounding the labels. If a vector of numbers
#'             is provided, all numbers are tried, until `n` unique labels are
#'             found. See 'Details'.
#' - `sf`:     Scaling factor. Axis values are divided by this number before the
#'             labels are calculated. If you set this to anything unequal 1, you
#'             should also choose `text` in a way that reflects the scaling.
#'             E.g. if you set `sf = 1e6` you could change the text from
#'             `"Signal Intensity [au]"` to `"Signal Intensity [Mau]"` or
#'             `"Signal Intensity [au] / 1e6"`, with `"Mau"` meaning
#'             "Mega-Arbitrary-Units".
#'
#' @param lc_verts,tp_verts
#' List of parameters passed to [segments()] when drawing vertical lines at the
#' centers of estimated, true or aligned lorentzian curves. Setting
#' `tp_verts$show` to TRUE requires `truepar` to be set. In addition to the
#' parameters of [segments()], the following additional parameters are supported
#' as well:
#'
#' - `height`: Height of the vertical line. The following options are supported:
#'   - `"full"`: The line is drawn from the bottom to the top of the plot region.
#'   - `"peak"`: The line is drawn from the bottom to the peak maximum.
#'   - `"int"`: The height of the line equals the peak integral. If this option
#'              is used, it is advisable to hide all other lines, as the
#'              integral values are usually on a different scale than the signal
#'              intensities. Two show both, signal intensities and integral
#'              values in the same plot, you can use the `sint` option.
#'   - `"sint"`: The height of the line equals the "scaled peak integral". This
#'               scaling is done in a way that the highest scaled integral equals
#'               the height of the plot region.
#'
#' @return
#' NULL. Called for side effect of plotting.
#'
#' @details
#' Parameters `bt_axis`, `lt_axis`, `tp_axis` and `rt_axis` all  support  option
#' `n` and `digits`, where `n = 5` means "Draw 5 tickmarks over  the  full  axis
#' range" and `digits = 3` means "round the label shown beside each tickmark  to
#' 3 digits". If `n` is omitted, a suitable value is chosen automatically  using
#' [axTicks()]. If `digits` is omitted, a default of `2:12` is used. Providing a
#' vector of `digits` causes each digit to be tried as argument  for  [round()],
#' until a digit is encountered that results in `n` unique labels. Example:
#'
#' Assume we have `n = 4` and the corresponding  calculated  tickmark  positions
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
#' In the above example the process would stop at `digit = 3`, because  at  this
#' point we have n = 4 unique labels (1.024, 1.025, 1.027 and 1.028).
#'
#' @examples
#' decon <- deconvolute(sim[[1]], sfr = c(3.55, 3.35))
#' draw_spectrum(obj = decon)
#' draw_spectrum(obj = decon, lgd = list(x = "top", bg = NA))
#' draw_spectrum(obj = decon, foc_rgn = c(3.45, 3.37))
#' draw_spectrum(obj = decon, fig = c(0.1, 0.4, 0.30, 0.45), add = TRUE)
#' draw_spectrum(obj = decon, fig = c(0.1, 0.4, 0.05, 0.20), add = FALSE)
#' draw_spectrum(obj = decon, lc_lines = NULL, lc_rects = NULL, foc_only = FALSE)
draw_spectrum <- function(
    obj,
    foc_rgn  = NULL,   foc_frac = NULL,   foc_only = TRUE,
    add      = FALSE,  fig_rgn  = NULL,   main     = NULL,
    show     = TRUE,   show_d2  = FALSE,  truepar  = NULL,
    mar      = c(4.1, 5.1, .1, .1),
    si_line  = list(), sm_line  = list(), sp_line  = list(),
    d2_line  = list(), lc_lines = list(), tp_lines = list(),
    cent_pts = list(), bord_pts = list(), norm_pts = list(),
    bg_rect  = list(), foc_rect = list(), lc_rects = list(), tp_rects = list(),
    bt_axis  = list(), lt_axis  = list(), tp_axis  = list(), rt_axis  = list(),
    tp_verts = list(), lc_verts = list(),
    lgd      = list()
    )
{
    # [ ]  Add `lc_int` argument.
    # [ ]  Add `fill` option to lc_lines and tp_lines
    # [ ]  Add `al_verts` argument.
    # [ ]  Add `al_arows` argument.
    # [ ]  Add `height` option to draw_verts.
    #      If "full", the line is drawn from bottom to top
    #      If "peak", the line is drawn from bottom to peak maximum
    #      If "sint", the line is drawn from bottom to

    # Check and enrich inputs (278us)
    if (isFALSE(show)) return()
    obj <- as_v2_obj(obj)
    stopifnot(
        is_num(foc_rgn, 2) || is.null(foc_rgn),
        is_num(foc_frac, 2) || is.null(foc_frac),
        is_bool(foc_only),
        is_bool(add),
        is_num(fig_rgn, 4) || is.null(fig_rgn),
        is_num(mar, 4) || is.null(mar),
        is_str(main) || is.null(main)
    )
    foc_frac <- foc_frac %||% get_foc_frac(obj, foc_rgn)
    foc_rgn <- foc_rgn %||% get_foc_rgn(obj, foc_frac)
    env <- environment()
    defaults <- get_ds_def_args(show_d2, foc_only)
    for (k in names(defaults)) {
        d <- defaults[[k]]
        x <- env[[k]]
        env[[k]] <- {
            if (is.null(x)) d
            else if (isFALSE(x)) modifyList(d, list(show = FALSE))
            else if (isTRUE(x)) modifyList(d, list(show = TRUE))
            else if (is.list(x)) modifyList(d, x)
            else stopf("Argument '%s' must be a list")
        }
    }
    decon_only <- c("sm_line", "lc_lines", "sp_line", "lc_rects", "cent_pts", "bord_pts")
    if (is_spectrum(obj)) for (k in decon_only) env[[k]]$show <- FALSE

    # Set graphical parameters (7ms)
    local_par(mar = mar, new = add)
    local_fig(fig = fig_rgn, add = add)

    # Get xy values over all data points (13us)
    cs <- cs_all <- obj$cs
    si <- si_all <- obj$si
    sm <- sm_all <- obj$sit$sm # NULL for spectrum objects (NFS)
    d2 <- d2_all <- NULL
    sp <- sp_all <- obj$sit$sup # NFS
    if (!isFALSE(d2_line$show)) {
        if (is.null(sm_all)) {
            warning("Smoothed SI is missing. Calculating second derivative from raw SI.")
            d2 <- d2_all <- calc_second_derivative(si_all)
        } else {
            d2 <- d2_all <- calc_second_derivative(sm_all)
        }
    }

    # Get indices of important points relative to all data points (611us)
    idp <- idp_all <- seq_along(cs_all) # Data points
    icp <- icp_all <- obj$peak$center # Center points (NFS)
    ibp <- ibp_all <- unique(sort(c(obj$peak$left, obj$peak$right))) # Border points (NFS)
    ipp <- ipp_all <- sort(c(icp_all, ibp_all)) # Peak points (NFS)
    inp <- inp_all <- setdiff(idp_all, ipp_all) # Nonpeak points (NFS)
    ifp <- ifp_all <- which(cs_all >= min(foc_rgn) & cs_all <= max(foc_rgn)) # Focus points

    # Get estimated lorentz parameters across all data points (8us)
    n0vec   <- numeric(0)
    n0df    <- data.frame(x0 = n0vec, A = n0vec, lambda = n0vec)
    lcpar   <- lcpar_all  <- obj$lcpar %||% n0df # (NFS)
    x0      <- x0_all     <- lcpar_all$x0 %||% n0df # (NFS)
    A       <- A_all      <- lcpar_all$A %||% n0df # (NFS)
    lambda  <- lambda_all <- lcpar_all$lambda %||% n0df # (NFS)

    # Get true lorentz parameters across all data points
    truepar <- truepar %||% obj$meta$simpar
    if (is.null(truepar)) {
        warnmsg <- "True params missing. Provide 'truepar' or unset 'tp_*' arguments."
        if (any(tp_rects$show, tp_lines$show, tp_verts$show)) warning(warnmsg, call. = FALSE)
        tp_rects$show <- FALSE
        tp_lines$show <- FALSE
        tp_verts$show <- FALSE
        trpar <- trpar_all <- data.frame(x0 = numeric(0), A = numeric(0), lambda = numeric(0))
    } else {
        trpar <- trpar_all <- as.data.frame(truepar[c("x0", "A", "lambda")])
    }

    # Get indices of important points relative to the vector of focus points (198us)
    if (foc_only) {
        offset <- min(ifp_all) - 1
        idp <- idp_foc <- ifp_all - offset # Data points
        icp <- icp_foc <- intersect(icp_all, ifp_all) - offset # Center points (NFS)
        ibp <- ibp_foc <- intersect(ibp_all, ifp_all) - offset # Border points (NFS)
        ipp <- ipp_foc <- intersect(ipp_all, ifp_all) - offset # Peak points (NFS)
        inp <- inp_foc <- intersect(inp_all, ifp_all) - offset # Nonpeak points (NFS)
    }

    # Get xy values over the focus region (18us)
    if (foc_only) {
        cs <- cs_foc <- cs_all[ifp_all]
        si <- si_foc <- si_all[ifp_all]
        sm <- sm_foc <- sm_all[ifp_all]
        d2 <- d2_foc <- d2_all[ifp_all]
        sp <- sp_foc <- sp_all[ifp_all]
    }

    # Get estimated lorentzians affecting focus region (239us)
    if (foc_only) {
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

    # Get true lorentzians affecting focus region (8us)
    if (foc_only) {
        y_foc_rgn_start    <- lorentz(x = min(foc_rgn), lcpar = trpar_all)
        y_foc_rgn_end      <- lorentz(x = max(foc_rgn), lcpar = trpar_all)
        y_tresh            <- 0.001 * diff(range(si))
        high_y_in_foc_rgn  <- y_foc_rgn_start > y_tresh | y_foc_rgn_end > y_tresh
        x0_in_foc          <- trpar$x0 > min(foc_rgn) & trpar$x0 < max(foc_rgn)
        affects_foc        <- high_y_in_foc_rgn | x0_in_foc
        trpar <- trpar_foc <- trpar_all[affects_foc, ]
    }

    # Get x and y limits (43us)
    xlim <- c(max(cs), min(cs))
    ymin <- if (show_d2) min(d2, na.rm = TRUE) else min(0, si, na.rm = TRUE)
    ymax <- if (show_d2) max(d2, na.rm = TRUE) else max(si, na.rm = TRUE)
    ylim <- c(ymin, ymax)
    ythresh <- 0.001 * diff(range(si))
    if (!is.null(foc_rgn)) {
        foc_rgn <- sort(foc_rgn)
        ymin_foc <- min(si[ifp])
        ymax_foc <- max(si[ifp])
    }
    ymin_foc <- if (foc_only) ymin else min(0, si_all[ifp_all])
    ymax_foc <- if (foc_only) ymax else max(si_all[ifp_all])
    ylim_foc <- c(ymin_foc, ymax_foc)

    # Do the actual drawing
    plot_empty(xlim, ylim, main = main)
    draw_rect(xlim, ylim, bg_rect)
    apply(trpar, 1, draw_lc_rect, cs, tp_rects, ymin) # 3ms
    apply(lcpar, 1, draw_lc_rect, cs, lc_rects, ymin) # 3ms
    apply(trpar, 1, draw_lc_line, cs, tp_lines, ythresh) # 13ms
    apply(lcpar, 1, draw_lc_line, cs, lc_lines, ythresh) # 13ms
    draw_line(cs, si, si_line)
    draw_line(cs, sm, sm_line)
    draw_verts(trpar, tp_verts)
    draw_verts(lcpar, lc_verts)
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
    list(
        plt_rgn_ndc = usr_to_ndc(c(xlim, ylim)),
        foc_rgn_ndc = usr_to_ndc(c(foc_rgn, ylim_foc))
    )
}

# Private #####

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
    in_xlim <- which(spec$ppm <= max(xlim) & spec$ppm >= min(xlim))
    plot(
        spec$ppm[in_xlim],
        spec$y_scaled[in_xlim],
        type = "l",
        xlab = "Chemical Shift [ppm]",
        ylab = "Signal Intensity [au]",
        xlim = xlim
    )
    mtext("Water Signal Region", side = 3, line = 0, at = x0, col = "blue")
    rect(
        xleft = x0 - hw,
        xright = x0 + hw,
        ybottom = par("usr")[3],
        ytop = par("usr")[4],
        col = rgb(0, 0, 1, alpha = 0.2),
        border = NA
    )
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

plot_empty <- function(xlim = c(0, 1),
                       ylim = c(0, 1),
                       xlab = "",
                       ylab = "",
                       main = "",
                       ann = TRUE,
                       axes = FALSE,
                       xaxs = "i",
                       yaxs = "i",
                       type = "n") {
    plot(x = xlim, y = ylim, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = main, ann = ann,
         axes = axes, xaxs = xaxs, yaxs = yaxs, type = type)
}

plot_dummy <- function() {
    plot(
        x = 0, y = 0, main = "",
        ylim = c(0, 1), xlim = c(0, 1),
        xaxs = "i", yaxs = "i",
        xlab = "dummy xlab", ylab = "dummy ylab"
    )
    text(0.5, 0.5, "dummy text")
}

# Draw Helpers (Private) #####

draw_legend <- function(args, si_line, sm_line, lc_lines, sp_line, d2_line,
                        cent_pts, bord_pts, norm_pts) {
    if (isFALSE(pop(args, "show"))) return()
    dscs <- c("Raw Signal", "Smoothed Signal", "Single Lorentzian",
              "Sum of Lorentzians", "Second Derivative",
              "Peak Center", "Peak Border", "Non-Peak")
    objs <- list(si_line, sm_line, lc_lines, sp_line, d2_line, cent_pts, bord_pts, norm_pts)
    lins <- list(si_line, sm_line, lc_lines, sp_line, d2_line)
    pnts <- list(cent_pts, bord_pts, norm_pts)
    ltys <- sapply(lins, function(obj) obj$lty %||% 1)
    pchs <- sapply(pnts, function(obj) obj$pch %||% 1)
    cols <- sapply(objs, function(obj) obj$col %||% "black")
    keep <- sapply(objs, function(obj) !isFALSE(obj$show))
    args$legend <- dscs[keep]
    args$col <- cols[keep]
    args$lty <- c(ltys, NA, NA, NA)[keep]
    args$pch <- c(NA, NA, NA, NA, NA, pchs)[keep]
    args$x <- args$x %||% "topright"
    do.call(legend, args)
}

draw_con_lines <- function(top_fig, bot_fig, con_lines) {
    if (isFALSE(pop(con_lines, "show"))) return()
    if (is.null(top_fig$plt_rgn_ndc)) return()
    if (is.null(bot_fig$foc_rgn_ndc)) return()
    plt_rgn_usr <- ndc_to_usr(top_fig$plt_rgn_ndc)
    foc_rgn_usr <- ndc_to_usr(bot_fig$foc_rgn_ndc)
    local_par(xpd = NA)
    con_lines$y0 <- max(foc_rgn_usr[3:4])
    con_lines$y1 <- min(plt_rgn_usr[3:4])
    con_lines$x0 <- min(foc_rgn_usr[1:2])
    con_lines$x1 <- min(plt_rgn_usr[1:2])
    do.call(segments, con_lines)
    con_lines$x0 <- max(foc_rgn_usr[1:2])
    con_lines$x1 <- max(plt_rgn_usr[1:2])
    do.call(segments, con_lines)
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
    n          <- pop(args, "n", default = 5)
    digits     <- pop(args, "digits", default = 2:12)
    sf         <- pop(args, "sf", default = 1)
    text       <- pop(args, "text")
    skip_first <- pop(args, "skip_first", FALSE)
    skip_last  <- pop(args, "skip_last", FALSE)
    args$at <- if (is_num(n)) seq(min(lim), max(lim), length.out = n) else axTicks(side)
    labs_sc <- args$at / sf
    if (!isFALSE(args$labels)) {
        args$labels <- format(labs_sc, digits = 7) # default used by normal axis
        for (i in digits) {
            args$labels <- format(labs_sc, digits = i)
            if (length(unique(args$labels)) >= n) break
        }
        if (skip_first) args$labels[1] <- ""
        if (skip_last) args$labels[length(args$labels)] <- ""
    }
    if (!is.null(text)) mtext(text, side = side, line = c(3, 4, 3, 4)[side])
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

draw_box <- function(args = list()) {
    if (isFALSE(pop(args, "show"))) return()
    do.call(box, args)
}

draw_line <- function(x, y, args = list()) {
    if (isFALSE(pop(args, "show"))) return()
    args$x <- x
    args$y <- y
    do.call(lines, args)
}

draw_verts <- function(p, args = list()) {
    if (isFALSE(pop(args, "show"))) return()
    if (nrow(p) == 0) return()
    args$x0 <- p$x0
    args$x1 <- p$x0
    args$y0 <- rep(par("usr")[3], length(p$x0))
    args$y1 <- p$A / p$lambda
    do.call(segments, args)
}

draw_points <- function(x, y, args = list()) {
    if (isFALSE(pop(args, "show"))) return()
    args$x <- x
    args$y <- y
    do.call(points, args)
}

# Get Helpers (Private) #####

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

get_sub_fig_args <- function(obj, foc_frac, foc_rgn, sub1, sub2, sub3, dot_args) {
    if ((n <- length(dot_args)) != (k <- length(names(dot_args)))) {
        if (n == 1) dot_args <- dot_args[[1]]
        else stop("... must contain tag=value pairs or a single list of such pairs")
    }
    usr_args <- lapply(list(sub1, sub2, sub3), function(x)
        modifyList(dot_args, x)
    )
    layout <- paste0(
        if (isFALSE(sub1$show)) "_" else "1",
        if (isFALSE(sub2$show)) "_" else "2",
        if (isFALSE(sub3$show)) "_" else "3"
    )
    fig_rgns <- switch(layout,           # x1 x2 y1    y2
        "___" = list(c(0.00, 0.00, 0.00, 0.00), c(0.00, 0.00, 0.00, 0.00), c(0.00, 0.00, 0.00, 0.00)),
        "1__" = list(c(0.00, 1.00, 0.00, 1.00), c(0.00, 0.00, 0.00, 0.00), c(0.00, 0.00, 0.00, 0.00)),
        "_2_" = list(c(0.00, 0.00, 0.00, 0.00), c(0.00, 1.00, 0.00, 1.00), c(0.00, 0.00, 0.00, 0.00)),
        "__3" = list(c(0.00, 1.00, 0.00, 1.00), c(0.00, 0.00, 0.00, 0.00), c(0.00, 1.00, 0.00, 1.00)),
        "12_" = list(c(0.00, 1.00, 0.20, 1.00), c(0.00, 1.00, 0.00, 0.20), c(0.00, 0.00, 0.00, 0.00)),
        "1_3" = list(c(0.00, 1.00, 0.40, 1.00), c(0.00, 0.00, 0.00, 0.00), c(0.00, 1.00, 0.00, 0.20)),
        "_23" = list(c(0.00, 0.00, 0.00, 0.00), c(0.00, 1.00, 0.50, 1.00), c(0.00, 1.00, 0.00, 0.30)),
        "123" = list(c(0.00, 1.00, 0.50, 1.00), c(0.00, 1.00, 0.30, 0.50), c(0.00, 1.00, 0.00, 0.20)),
        stop("Invalid layout")
    )
    mars <- switch(layout,          # b  l  t  r
        "1_3" = list(c(1, 6, 2, 2), c(0, 0, 0, 0), c(0, 0, 1, 0)),
        "_23" = list(c(0, 0, 0, 0), c(1, 2, 2, 2), c(0, 0, 1, 0)),
        "123" = list(c(0, 6, 2, 2), c(1, 6, 0, 2), c(0, 0, 1, 0)),
        list(c(0, 0, 0, 0), c(0, 0, 0, 0), c(0, 0, 0, 0))
    )
    for (i in 1:3) {
        fig_rgns[[i]] <- usr_args[[i]]$fig_rgn_npc %||% fig_rgns[[i]]
        usr_args[[i]]$fig_rgn_npc <- NULL
    }
    def_args <- lapply(1:3, function(i) list(
        obj = obj,
        foc_rgn = foc_rgn,
        foc_frac = foc_frac,
        add = TRUE,
        fig_rgn = npc_to_ndc(fig_rgns[[i]]),
        mar = mars[[i]]
    ))
    def_args[[1]]$bt_axis <- if (sub2$show) FALSE else if (sub3$show) list(text = "")
    def_args[[2]]$lt_axis <- if (layout == "_2_") list() else FALSE
    def_args[[2]]$bt_axis <-if (sub3$show) list(text = "")
    def_args[[2]]$show_d2 <- TRUE
    def_args[[3]]$bg_rect <- if (layout != "__3") list(border = NA)
    def_args[[3]]$foc_only <- FALSE
    def_args[[3]]$foc_rect <- if (layout == "__3") FALSE
    def_args[[3]]$lgd <- FALSE
    def_args[[3]]$lt_axis <- if (layout != "__3") list(text = ""  )
    structure(
        mapply(modifyList, def_args, usr_args, SIMPLIFY = FALSE),
        names = c("sub1", "sub2", "sub3")
    )
}

get_ds_def_args <- function(show_d2 = FALSE, foc_only = TRUE) {
    stopifnot(is_bool(show_d2), is_bool(foc_only))
    show_si  <- !show_d2
    ylab <- if (show_si) "Signal Intensity [au] / 1e6" else "Second Derivative / 1e6"
    lgdx <- if (show_si) "topright" else "bottomright"
    list(
        lgd       = list(show = TRUE, x = lgdx, bg = "white"),
        d2_line   = list(show = show_d2),
        si_line   = list(show = show_si, col = "black"),
        sm_line   = list(show = show_si, col = "blue"),
        sp_line   = list(show = show_si, col = "red"),
        lc_lines  = list(show = show_si && foc_only, col = "darkgrey"),
        tp_lines  = list(show = FALSE, col = "darkgreen"),
        tp_verts  = list(show = FALSE, col = "darkgreen"),
        lc_verts  = list(show = FALSE, col = "darkgrey"),
        cent_pts  = list(show = show_si && foc_only, col = "blue", bg = "blue", pch = 24),
        bord_pts  = list(show = FALSE, col = "blue", pch = 124),
        norm_pts  = list(show = FALSE, col = "black", pch = 124),
        bg_rect   = list(show = TRUE, col = "white"),
        foc_rect  = list(show = TRUE, col = transp("yellow")),
        lc_rects  = list(show = FALSE, col = transp("black"), border = NA),
        tp_rects  = list(show = FALSE, col = transp("darkgreen", 0.12), border = NA),
        bt_axis   = list(show = TRUE, text = "Chemical Shift [ppm]"),
        lt_axis   = list(show = TRUE, text = ylab, sf = 1e6, las  = 1),
        tp_axis   = list(show = FALSE),
        rt_axis   = list(show = FALSE)
    )
}

combine <- function(defaults, x) {
    name <- deparse(substitute(x))
    x <- if (is.null(x)) list()
        else if (isFALSE(x)) list(show = FALSE)
        else if (isTRUE(x)) list(show = TRUE)
        else if (is.list(x)) x
        else stopf("%s must be a bool or list", name)
    modifyList(defaults, x)
}

# Convert Helpers (Private) #####

usr_to_ndc <- function(usr = par("usr")) {
    x_ndc <- grconvertX(usr[1:2], from = "user", to = "ndc")
    y_ndc <- grconvertY(usr[3:4], from = "user", to = "ndc")
    c(x_ndc, y_ndc)
}

ndc_to_usr <- function(ndc) {
    x_usr <- grconvertX(ndc[1:2], from = "ndc", to = "user")
    y_usr <- grconvertY(ndc[3:4], from = "ndc", to = "user")
    c(x_usr, y_usr)
}

npc_to_ndc <- function(npc = c(0, 1, 0, 1)) {
    if (is.null(npc)) return(NULL)
    x_ndc <- grconvertX(npc[1:2], from = "npc", to = "ndc")
    y_ndc <- grconvertY(npc[3:4], from = "npc", to = "ndc")
    c(x_ndc, y_ndc)
}

# Multifigure Helpers (Private) #####

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
    if (is.null(fig)) return(function() {}) # Nothing to do if figure region is NULL
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
with_fig <- function(expr, fig = NULL, pos = NULL, add = TRUE) {
    reset_mfc <- set_fig(fig = fig, add = add)
    on.exit(reset_mfc())
    expr
}

local_fig <- function(fig = NULL, add = TRUE, envir = parent.frame()) {
  reset_mfc <- set_fig(fig = fig, add = add)
  defer(reset_mfc(), envir = envir)
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
