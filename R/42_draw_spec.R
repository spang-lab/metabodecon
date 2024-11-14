# Internal #####

# styler: off
#' @noRd
#' @title Plot Spectrum Internal
#' @description
#' Helper for [plot_spectrum()].
#' For argument descriptions see [plot_spectrum()].
#' For examples see `test/testthat/test-draw_spectrum.R`
draw_spectrum <- function(
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
    args <- ds_get_args()
    old <- par(mar = mar, new = add)
    on.exit(par(old), add = TRUE)
    rst <- set_fig(fig = fig, add = add)
    on.exit(rst(), add = TRUE)
    dat <- ds_get_dat(args)
    plt <- ds_init_plot_region(dat, verbose)
    bgr <- ds_draw_bg(dat, fill_col, verbose)
    lns <- ds_draw_lines(dat, line_col, foc_only, sm_show, sm_col, verbose)
    trp <- ds_draw_triplets(dat, trp_show, trp_pch, trp_col, verbose)
    lcs <- ds_draw_lorentz_curves(dat, args, verbose)
    foc <- ds_draw_focus_rectangle(dat, rct_show, rct_fill, rct_col, verbose)
    axs <- ds_draw_axis(dat, main, xlab, ylab, axis_col, box_col, verbose)
    lgd <- ds_draw_legend(args, verbose)
    named(dat, plt, bgr, axs, lns, trp, lcs, foc, par = par())
}
# styler: on

# Helpers #####

ds_get_args <- function(env = parent.frame()) {
    keys <- names(formals(draw_spectrum))
    keys <- keys[keys != "..."]
    args <- sapply(keys, get, envir = env)
    if (!is.null(args$foc_rgn) && length(args$foc_rgn) != 2) {
        stop("foc_rgn must be a numeric vector of length 2")
    }
    if (args$foc_only && is.null(args$foc_rgn)) {
        stop("foc_only requires foc_rgn to be specified")
    }
    invisible(args)
}

ds_get_dat <- function(args) {
    decon <- args$decon
    foc_rgn <- args$foc_rgn
    foc_unit <- args$foc_unit
    foc_only <- args$foc_only
    sf_y_raw <- args$sf_y_raw
    ysquash <- args$ysquash

    # Chemical Shift
    cs <- {                          # Available in:
        decon[["ppm"]] %||%          # - idecon, ispec
        decon[["x_values_ppm"]] %||% # - decon0, decon1
        decon[["cs"]] %||%           # - decon2, spectrum
        stop("chemical Shifts missing", call. = FALSE)
        # Use `[[` instead of `$` to prevent partial matching
    }

    # Raw Signal Intensity
    si <- {                          # Available in:
        decon[["y_raw"]] %||%        # - idecon, ispec
        decon[["y_values_raw"]] %||% # - decon1
        decon[["si"]] %||%           # - decon2, spectrum
        stop("raw signal intensities missing", call. = FALSE)
    }

    # Smoothed Signal Intensity
    sis <- {
        decon[["y_smooth"]] %||%
        decon[["y_values"]] %||%
        stop("smoothed signal intensities missing", call. = FALSE)
    }
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

    # Indices of Focus Region
    ifr <- {
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
    x <- if (foc_only) cs[ifr] else cs
    y <- if (foc_only) si[ifr] else si
    xlim <- c(max(x), min(x))
    ylim <- c(0, max(y) / ysquash)

    locals(without = c("decon"))
}

ds_init_plot_region <- function(dat = ds_get_dat(),
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

ds_draw_bg <- function(dat = ds_get_dat(),
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

ds_draw_axis <- function(dat,
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

ds_draw_lines <- function(dat,
                           line_col = "black",
                           foc_only = FALSE,
                           sm_show = TRUE,
                           sm_col = "blue",
                           verbose = FALSE) {
    x <- if (foc_only) dat$cs[dat$ifr] else dat$cs
    y <- if (foc_only) dat$si[dat$ifr] else dat$si
    ys <- if (foc_only) dat$sis[dat$ifr] else dat$sis
    if (verbose) logf("Drawing raw signal")
    lines(x, y, type = "l", col = line_col, lty = 1)
    if (sm_show && !is.null(ys)) {
        if (verbose) logf("Drawing smoothed signal")
        lines(x, ys, type = "l", col = sm_col, lty = 1)
    }
}

ds_draw_triplets <- function(dat,
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

ds_draw_focus_rectangle <- function(dat,
                                     show = TRUE,
                                     fill = rgb(0, 0, 1, alpha = 0.1),
                                     col = "blue",
                                     verbose = FALSE) {
    if (is.null(dat$foc_lim) || !show) return(NULL) # styler: off
    usr <- list(
        xleft = max(dat$foc_lim),
        xright = min(dat$foc_lim),
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
    if (verbose) logf("Drawing focus rectangle")
    do.call(rect, usr)
    named(usr, ndc)
}

ds_draw_lorentz_curves <- function(dat, args = NULL, verbose = FALSE) {
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
    x <- if (foc_only) dat$cs[dat$ifr] else dat$cs
    Y <- matrix(nrow = length(x), ncol = length(dat$A))
    if (verbose) logf("Drawing individual Lorentzian Curves")
    for (i in seq_along(dat$A)) {
        Y[, i] <- ds_draw_lorentz_curve(
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

ds_draw_lorentz_curve <- function(x,
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

ds_draw_legend <- function(args, verbose = FALSE) {
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

ds_setup_dev_env <- function() {
    sim_01 <- metabodecon_file("sim/sim_01")
    decon <- generate_lorentz_curves(
        sim_01,
        sfr = c(3.35, 3.55), wshw = 0, delta = 0.1,
        ask = FALSE, verbose = FALSE
    )
    args <- stub("draw_spectrum", decon = decon, foc_rgn = c(3.55, 3.52))
    invisible(args)
}
