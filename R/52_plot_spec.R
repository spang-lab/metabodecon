# Internal #####

#' @noRd
#' @title Plot Spectrum Internal
#' @description Helper for [plot_spectrum()]. For arguments see [plot_spectrum()].
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
#'     p[[3]] <- plot_spec(decon,
#'         foc_rgn = c(3.55, 3.52),
#'         foc_only = TRUE, fig = leftmiddle,
#'         fill_col = rgb(0, 0, 1, 0.1)
#'     )
#'     p[[4]] <- plot_dummy()
#'     p[[5]] <- plot_spec(decon, fig = leftbottom, add = FALSE)
#'     p[[6]] <- plot_spec(decon, lc_show = FALSE, trp_show = FALSE)
#'     p
#' })
plot_spec <- function(
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
    rct_fill = trans("yellow"), # Background color of rectangle around FR.
    rct_col = "black", # Border color of rectangle around FR.
    # Spectrum Lines
    line_col = "black", # Color of raw signal intensities.
    sm_show = TRUE, # If TRUE, smoothed signal intensities are shown.
    sm_col = "blue", # Color of smoothed signal intensities.
    ysquash = 0.96, # Fraction of plot height to squash y-values into.
    sf_y_raw = 1e6, # Divide raw SI by this factor before drawing.
    d2_show = FALSE, # If TRUE, shows the second derivative of SIs.
    # Lorentzians
    lc_show = TRUE, # If TRUE, the Lorentzian Curves (LCs) are shown.
    lc_col = "black", # Color of the Lorentzian Curves.
    lc_lty = 1, # Line type of the Lorentzian Curves.
    lc_fill = trans("black"), # BG-color of rectangles shown at LC-center.
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
    ...) {
    #
    # Start Function Body
    #
    args <- psi_get_args()

    old <- par(mar = mar, new = add)
    on.exit(par(old), add = TRUE)
    rst <- set_fig(fig = fig, add = add)
    on.exit(rst(), add = TRUE)

    dat <- psi_get_dat(args)
    plt <- psi_init_plot_region(dat, verbose)
    bgr <- psi_draw_bg(dat, fill_col, verbose)
    lns <- psi_draw_lines(dat, line_col, d2_show, foc_only, sm_show, sm_col, verbose)
    trp <- psi_draw_triplets(dat, trp_show, trp_pch, trp_col, verbose)
    lcs <- psi_draw_lorentz_curves(dat, args, verbose)
    foc <- psi_draw_focus_rectangle(dat, rct_show, rct_fill, rct_col, verbose)
    axs <- psi_draw_axis(dat, main, xlab, ylab, axis_col, box_col, verbose)
    lgd <- psi_draw_legend(args, verbose)
    named(dat, plt, bgr, axs, lns, trp, lcs, foc, par = par())
}

# Helpers #####

psi_get_args <- function(env = parent.frame()) {
    args <- sapply(names(formals(plot_spec)), get, envir = env)
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
                           d2_show = FALSE,
                           foc_only = FALSE,
                           sm_show = TRUE,
                           sm_col = "blue",
                           verbose = FALSE) {
    x <- if (foc_only) dat$cs[dat$ifp] else dat$cs
    y <- if (foc_only) dat$si[dat$ifp] else dat$si
    ys <- if (foc_only) dat$sis[dat$ifp] else dat$sis
    if (verbose) logf("Drawing raw signal")
    lines(x, y, type = "l", col = line_col, lty = 1)
    if (d2_show) {
        if (verbose) logf("Drawing second derivative (TODO)")
    }
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
    y <- lc(x, x_0, A, lambda)
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
            ytop = lc(x_0, x_0, A, lambda),
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
    args <- stub("plot_spec", decon = decon, foc_rgn = c(3.55, 3.52))
    invisible(args)
}
