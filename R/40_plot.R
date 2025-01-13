# Public #####

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
#' [get_sub_fig_rgns()] is used to determine a suitable layout automatically.
#'
#' @param sub1,sub2,sub3
#' List of arguments passed to [draw_spectrum()] when drawing sub figure
#' 1-3. See 'Details'.
#'
#' @return
#' NULL. Called for side effect of plotting as sketched in 'Details'.
#'
#' @details
#' This function first calls [plot_empty()] to initalize a new plotting canvas.
#' After that it calls [draw_spectrum()] multiple times to draw the following sub
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
#' [draw_spectrum()] by providing a corresponding named element in `sub1`,
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
#' Note that the figure created by `plot_spectrum()` can be part of a
#' multi-figure configuration as created when setting `mfrow` or `mfcol` via
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
#' @title Plot Spectra
#'
#' @description
#' Plot a set of deconvoluted spectra.
#'
#' @param x
#' An object of type decons0, decons1 or decons2. For details see
#' [Metabodecon Classes](https://spang-lab.github.io/metabodecon/articles/Metabodecon-Classes.html).
#'
#' @param mar A numeric vector of length 4, which specifies the margins of the plot.
#'
#' @param peak_rng A numeric vector of length 2, specifying the ppm range across
#' all peaks of the provided deconvoluted spectra.
#'
#' @return
#' A plot of the deconvoluted spectra.
#'
#' @seealso [plot_spectrum()] for a much more sofisticated plotting routine
#' suitable to plot a single spectrum.
#'
#' @examples
#' decons <- generate_lorentz_curves(
#'      data_path = sim,
#'      sfr = c(3.55, 3.35),
#'      wshw = 0,
#'      ask = FALSE,
#'      verbose = FALSE
#' )
#' plot_spectra(decons)
#'
plot_spectra <- function(x,
                         mar = c(4.1, 4.1, 1.1, 0.1),
                         peak_rng = get_ppm_range(x, show = FALSE)) {
    # CONTINUE HERE: Refactor to:
    # 1. Use decons_v2 object
    # 2. Use x$sit$al if available
    # 3. Else use x$sit$sup if available
    # 4. Else use x$si
    x <- as_decons1(x)
    if (is_decon_list(x)) {
        xrng <- range(c(sapply(x, function(s) s$x_values_ppm)))
        ymax <- max(sapply(x, function(s) max(s$y_values)))
    } else if (is.data.frame(x)) {
        stop("x must be a list of deconvoluted spectra.")
    }
    a <- peak_rng[1]
    b <- peak_rng[2]
    w <- (b - a) / 4
    y8 <- ymax * 0.8
    cols <- rainbow(length(x))
    ltxt <- paste("Spectrum", 1:length(x))
    opar <- par(mar = mar)
    on.exit(par(opar))
    plot(x = NA, type = "n", xlim = xrng[2:1], ylim = c(0, ymax), xlab = "Chemical Shift [ppm]", ylab = "Signal Intensity [au]")
    abline(v = peak_rng, lty = 2)
    for (i in 1:length(x)) lines(x = x[[i]]$x_values_ppm, y = x[[i]]$y_values, col = cols[i])
    arrows(x0 = c(a + w, b - w), x1 = c(a, b), y0 = y8, y1 = y8, length = 0.1, lty = 2, col = "black")
    text(x = mean(c(a, b)), y = y8, labels = "ppm range")
    mtext(text = round(c(a, b), 4), side = 3, line = 0, at = c(a, b))
    legend(x = "topright", legend = ltxt, col = cols, lty = 1)
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
