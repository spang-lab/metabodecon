develop_draw_spectrum <- function() {
    # Mandatory
    obj <- generate_lorentz_curves(sim[[1]], sfr = c(3.35, 3.55), wshw = 0, ask = FALSE, verbose = FALSE)
    # Figure Region
    fig <- NULL # Drawing region in NDC.
    add <- !is.null(fig) # If TRUE, the plot is added to the current figure.
    # Plot Region
    main <- "" # Title of the plot.
    lgd <- TRUE # If TRUE, a legend is drawn.
    xlab <- "Chemical Shift [ppm]" # Label for the x-axis.
    ylab <- "Signal Intensity [au]" # Label for the y-axis.
    mar <- c(4.1, 4.1, 0.1, 0.1) # Lines below/left/above/right plot region.
    box_col <- "black" # Color of box surrounding plot region.
    axis_col <- "black" # Color of tickmarks and ticklabels.
    fill_col <- NULL # Background color of plot region.
    # Focus Region (FR)
    foc_rgn <- c(3.50, 3.40) # Start and end of FR.
    foc_unit <- "ppm" # Unit of `foc_rgn`. Either "fraction" or "ppm".
    foc_only <- FALSE # If TRUE, draw focusregion, else full spectrum.
    # Focus Rectangle
    rct_show <- !is.null(foc_rgn) && !foc_only # Draw rectangle around FR?
    rct_fill <- transp("yellow") # Background color of rectangle around FR.
    rct_col <- "black" # Border color of rectangle around FR.
    # Spectrum Lines
    line_col <- "black" # Color of raw signal intensities.
    sm_show <- TRUE # If TRUE, smoothed signal intensities are shown.
    sm_col <- "blue" # Color of smoothed signal intensities.
    ysquash <- 0.96 # Fraction of plot height to squash y-values into.
    sf_y_raw <- 1e6 # Divide raw SI by this factor before drawing.
    # Lorentzians
    lc_show <- TRUE # If TRUE, the Lorentzian Curves (LCs) are shown.
    lc_col <- "black" # Color of the Lorentzian Curves.
    lc_lty <- 1 # Line type of the Lorentzian Curves.
    lc_fill <- transp("black") # BG-color of rectangles shown at LC-center.
    sup_show <- TRUE # If TRUE, the LC-Superposition (LC-Sup) is shown.
    sup_col <- "red" # Color of LC-Sup.
    sup_lty <- 1 # Line type of LC-Sup.
    # Peak Triplets
    trp_show <- TRUE # If TRUE, the peak triplets are shown.
    trp_col <- rep("black", 4) # Colors for center, left, right, non-peak DPs.
    trp_pch <- c(17, 4, 4, NA) # Pchars for center, left, right, non-peak DPs.
    # Misc
    verbose <- FALSE
}

test_draw_spectrum <- function() {
    path <- metabodecon_file("sim/sim_01")
    spec <- read_spectrum(path)
    decon <- generate_lorentz_curves(spec, sfr = c(3.55, 3.35), ws = 0, ask = FALSE, verbose = FALSE)
    plot_dummy <- function() {
        plot(0, 0, ylim = c(0, 1), xlim = c(0, 1), xaxs = "i", yaxs = "i")
        text(0.5, 0.5, "dummy")
    }
    leftmiddle <- c(0.1, 0.4, 0.30, 0.45)
    leftbottom <- c(0.1, 0.4, 0.05, 0.20)
    p <- local({
        local_par(mfrow = c(4, 2), mar = c(2, 2, 0.5, 0.5))
        p <- list()
        p[[1]] <- plot_dummy()
        p[[2]] <- draw_spectrum(obj = spec)
        p[[3]] <- draw_spectrum(obj = spec)
        p[[4]] <- draw_spectrum(obj = decon, foc_rgn = c(3.45, 3.37))
        p[[5]] <- plot_dummy()
        p[[5]] <- draw_spectrum(obj = decon, foc_rgn = c(3.45, 3.37), foc_only = TRUE, fig = leftmiddle, fill_col = rgb(0, 0, 1, 0.1))
        p[[6]] <- plot_dummy()
        p[[7]] <- draw_spectrum(obj = decon, fig = leftbottom, add = FALSE)
        p[[8]] <- draw_spectrum(obj = decon, lc_show = FALSE, trp_show = FALSE)
        p
    })
}

test_result <- test_that("draw_spectrum works", {
    tmp <- vdiffr::expect_doppelganger(
        title = "draw_spectrum",
        fig = test_draw_spectrum,
        writer = function(plot, file, title = "") {
            with_svg(file, plot(), width = 12, height = 16)
        }
    )
})
