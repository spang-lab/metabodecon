# Figures #####

#' @noRd
#' @examples
#' @param s # Scaling Factor for the created figure. By making the figure `s` times bigger than the textwidth you can make the font, lines, etc. `s` times smaller, thinner, etc. because it will get scaled down when included in the LaTeX document.
#' @param w # Figure width. Should be the textwidth of the Latex document in inches.
#' @param h # Figure height. Should be the textheight of the LaTeX document in inches minus some space for the fiure caption. Example: if the textheight is 9.5 inches, h = 8 could be a good choice.
#' @examples
#' mkfig_nmr_challenges()
mkfig_nmr_challenges <- function(show = TRUE,
                                 path = "challenges.pdf",
                                 s = 2,
                                 w = 5.45,
                                 h = 8,
                                 ...) {
    if (show) {
        if (dev.cur() == 1) dev.new(width = w * s, height = h * s, rescale = "fixed")
        plot_nmr_challenges()
    }
    if (is_str(path)) {
        withr::local_pdf(path, width = w * s, height = h * s)
        plot_nmr_challenges()
    }
}

# Plot #####
plot_nmr_challenges <- function(lwd = 0.1,
                                s = 1,
                                w = 5.45,
                                h = 8) {

    init_dev(s, w, h)
    local_par(mfrow = c(5, 2), mar = c(0, 0, 2, 0))
    clear_dev()
    plot_nmr_experiment()
    plot_yt_spectra()
    plot_preprocessing()
    plot_yf_spectra()
    plot_deconvolution()
    plot_deconvoluted_spectra()
    plot_alignment()
    plot_aligned_spectra()
    plot_annotation()
    plot_anntated_spectra()
}

plot_nmr_experiment <- function() {
    plot_empty(main = "NMR Experiment")
    marbox()
    text(.2, .75, "Samples")
    draw_vial(x1= .1, y1 = .6, h = .2)
    draw_vial(x1= .2, y1 = .3, h = .2)
    draw_vial(x1= .1, y1 = .1, h = .2)
    text(.6, .75, "NMR Spectrometer")
    draw_nmr_spectrometer(x1 = .5, y1 = .1, h = .6)
}

plot_yt_spectra <- function() {
    local_par(mar = c(1, 0, 2, 0))
    plot_empty(main = "FID Signals")
    marbox()
    ndc <- npc_to_ndc()
    fig_height <- diff(ndc[3:4])
    sub_height <- fig_height / 3
    local_par(mar = c(0, 0, 2, 0))
    with_fig(fig = c(0.50, 1.00, ndc[3] + 2 * sub_height, ndc[3] + 3 * sub_height), plot_yt_spectrum(1))
    with_fig(fig = c(0.50, 1.00, ndc[3] + 1 * sub_height, ndc[3] + 2 * sub_height), plot_yt_spectrum(2))
    with_fig(fig = c(0.50, 1.00, ndc[3] + 0 * sub_height, ndc[3] + 1 * sub_height), plot_yt_spectrum(3))
}

plot_yt_spectrum <- function(seed = 1) {
    set.seed(seed)
    time <- seq(0, 1, by = 0.001)
    decay_constant <- 5 + rnorm(1, 0, 0.5)
    frequency <- 50 + rnorm(1, 0, 2)
    fid_signal <- exp(-decay_constant * time) * sin(2 * pi * frequency * time)
    plot(time, fid_signal, type = "l", xlab = "", ylab = "", axes = FALSE)
    # main = "Simulated FID Signal"
    # xlab = "Time (s)"
    # ylab = "Signal Intensity"
}

# Plots the yf spectra on the left side, with an arrow going into a
# pipeline. In the pipeline, the following steps are written as text:,
# Apodization, Zero-Filling, Fourier Transform, Phase Correction, Baseline
# Correction, Referencing
plot_preprocessing <- function() {
    local_par(mar = c(1, 0, 2, 0))
    plot_empty(main = "Preprocessing")
    marbox()
    ndc <- npc_to_ndc()
    x1 <- ndc[1]; x2 <- ndc[2]; y1 <- ndc[3]; y2 <- ndc[4]
    H <- diff(ndc[3:4]); h <- H / 100
    W <- diff(ndc[1:2]); w <- W / 100
    # y~f spectra
    local_par(mar = c(0, 0, 0, 0))
    with_fig(fig = c(x1 + 5 * w, x1 + 40 * w, y1 + 66 * h, y1 + 99 * h), plot_yt_spectrum(1))
    with_fig(fig = c(x1 + 5 * w, x1 + 40 * w, y1 + 33 * h, y1 + 66 * h), plot_yt_spectrum(2))
    with_fig(fig = c(x1 + 5 * w, x1 + 40 * w, y1 +  0 * h, y1 + 33 * h), plot_yt_spectrum(3))

    # Write the preprocessing steps to the right
    local_par(mar = c(0, 0, 0, 0))
    with_fig(fig = c(x1 + 50 * w, x1 + 90 * w, y1, y1 + 98 * h), {
        plot_empty()
        lh <- 1.5 # Line height
        cex <- 0.8 # Chracter Size
        box()
        rect(0, 0, 1, 1, col = "lightyellow")
        mtext("Apodization", side = 3, line = -1 * lh, cex = cex)
        mtext("Zero-Filling", side = 3, line = -2 * lh, cex = cex)
        mtext("Fourier Transform", side = 3, line = -3 * lh, cex = cex)
        mtext("Phase Correction", side = 3, line = -4 * lh, cex = cex)
        mtext("Baseline Correction", side = 3, line = -5 * lh, cex = cex)
        mtext("Referencing", side = 3, line = -6 * lh, cex = cex)
    })
}

plot_deconvolution <- function() {
    plot_empty(main = "Deconvolution")
    marbox()
    metabodecon::draw_spectrum(
        metabodecon::deconvolute(metabodecon::sim[[1]]),
        foc_frac = c(0.48, 0.38),
        mar = c(0, 0, 0, 0),
        lgd = FALSE,
        bt_axis = FALSE,
        bg_rect = FALSE,
        foc_rect = FALSE,
        lc_lines = TRUE,
        sp_line = FALSE,
        cent_pts = FALSE,
        sm_line = FALSE,
        add = TRUE
    )
}

plot_yf_spectra <- function() {
    plot_empty(main = "Raw Data in Frequency Domain")
    marbox()
}

plot_deconvoluted_spectra <- function() {
    plot_empty(main = "Deconvoluted Signal Intensities")
    marbox()
}

plot_alignment <- function() {
    plot_empty(main = "Alignment")
    marbox()
}

plot_aligned_spectra <- function() {
    plot_empty(main = "Aligned Signal Intensities")
    marbox()
}

plot_annotation <- function() {
    plot_empty(main = "Identification")
    marbox()
}

plot_anntated_spectra <- function() {
    plot_empty(main = "Metabolite Concentrations")
    marbox()
}


# Draw #####

#' @noRd
#' @details See `misc/sketches/blood_vial.drawio` for a sketch of the vial
draw_vial_2 <- function(x, y, width = 0.2, border = "black") {
    xp <- grconvertWidth(width / 100, "inches", "user")   # 1 xp == 1/100 of full x-width given in user coords
    yp <- grconvertHeight(height / 100, "inches", "user") # 1 yp == 1/100 of full y-height given in user coords
    cap <- list(
        x1 = x - 50 * xu,
        x2 = x + 50 * xu,
        y1 = y + 40 * yu,
        y2 = y + 50 * yu,
    )
    vial <- sticker <- circle <- cross <- bottom <- list()
}

draw_vial <- function(x1 = 10, y1 = 10, h = 20, w = NULL) {
    # See `misc/sketches/blood_vial.drawio` for a sketch of the vial
    if (is.null(w)) w <- (h/3) * ipu_ratio()
    x2 <- x1 + w
    y2 <- y1 + h
    x <- seq(x1, x2, length.out = 13)
    y <- seq(y1, y2, length.out = 37)
    cross_width <- x[9] - x[5]
    cross_width_inch <- grconvertW(cross_width, "user", to = "inches")
    cross_height_inch <- cross_width_inch
    cross_height <- grconvertH(cross_width_inch, "inches", to = "user")
    cross_thickness <- cross_width / 6
    vial_height <- y[29] - y[5]
    if (2 * cross_height > vial_height) stop("Sticker height > vial height. Increase height/width ratio of vial.")
    sticker_bottom <- y[17] - cross_height
    sticker_top <- y[17] + cross_height
    symbols(x[7], y[5], circles = x[7] - x[3], inches = FALSE, add = TRUE, bg = "red") # bottom
    rect(x[3], y[5], x[11], y[29], col = "red", border = NA) # vial
    segments(x[3], y[5], x[3], y[29], col = "black") # vial
    segments(x[3], y[29], x[11], y[29], col = "black") # vial
    segments(x[11], y[5], x[11], y[29], col = "black") # vial
    rect(x[3], sticker_bottom, x[11], sticker_top, col = "white", border = "black") # sticker
    rect(x[1], y[29], x[13], y[37], col = "white", border = "black") # cap
    symbols(x[7], y[17], circles = x[7] - x[4], inches = FALSE, add = TRUE, bg = "red") # circle
    rect(x[7] - cross_width / 2, y[17] - cross_height / 6,
         x[7] + cross_width / 2, y[17] + cross_height / 6,
         col = "white", border = NA) # vial
    rect(x[7] - cross_width / 6, y[17] - cross_height / 2,
         x[7] + cross_width / 6, y[17] + cross_height / 2,
        col = "white", border = NA) # vial
}

draw_nmr_spectrometer <- function(x1 = 10, y1 = 20, h = 60, w = NULL) {
    if (is.null(w)) w <- ((2 * h) / 3) * ipu_ratio()
    x2 <- x1 + w
    y2 <- y1 + h
    xc <- (x1 + x2) / 2
    yc <- (y1 + y2) / 2

    # Magnet (large rectangle)
    rect(x1, y1, x2, y2, border = "black", col = "gray90")

    # Sample holder (small circle inside the magnet)
    circles(xc, yc, radii = w * 0.125, fg = "black", bg = "white")

    # Console/electronics (rectangle to the right)
    console_w <- w * 0.375
    console_h <- h * 0.5
    rect(x2 + .05, y1,
         x2 + console_w + .05,
         y1 + console_h,
         border = "black",
         col = "gray80"
    )

    # Draw the control panel (small rectangle near console)
    rect(
        x2 + .05,
        y1 + console_h + .05,
        x2 + console_w + .05,
        y1 + console_h + h * 0.166,
        border = "black",
        col = "gray70"
    )

    # Add labels
    text((x1 + x2) / 2, (y1 + y2) / 2, "Sample")
    text(x2 + console_w / 2 + .05, y1 + console_h / 2, "Console")
    text((x1 + x2) / 2, y2 - h * 0.05, "Magnet")
}

# Plot Helpers #####

init_dev <- function(s = 2, w = 5.45, h = 8) {
    if (dev.cur() == 1) dev.new(width = w * s, height = h * s, rescale = "fixed")
}

clear_dev <- function() {
    mfg <- par("mfg")
    while (mfg[1] != mfg[3] || mfg[2] != mfg[[4]]) {
        plot_empty()
        mfg <- par("mfg")
    }
    plot_empty()
    mfg[1:2] <- c(1, 1)
    par(mfg = mfg)
}

get_ndc <- function() {
    usr <- par("usr")  # returns c(x1, x2, y1, y2)
    xndc <- grconvertX(usr[1:2], from = "user", to = "ndc")
    yndc <- grconvertY(usr[3:4], from = "user", to = "ndc")
    ndc <- c(x1 = xndc[1], x2 = xndc[2], y1 = yndc[1], y2 = yndc[2])
}

marbox <- function() {
    withr::with_par(list(mar = c(0,0,0,0)), box())
}

circles <- function(x, y, radii, inches = FALSE, add = TRUE, ...) {
    symbols(x, y, circles = radii, inches = inches, add = add, ...)
}

figsize <- function() {
    inches <- par("pin")
    cm <- inches * 2.54
    dpi <- 72
    dp <- inches * dpi
    cat("Size in incheses:", inches, "\n")
    cat("Size in centimeters:", cm, "\n")
    cat("Size in pixels:", dp, "\n")
}

# Inch per Unit Ratio
ipu_ratio <- function() {
    d_inch <- par("pin") # plot dimension in inches
    w_inch <- d_inch[1]  # plot width in inches
    h_inch <- d_inch[2]  # plot height in inches
    x_user <- par("usr") # plot extremes (x1,x2,y1,y2) in user coordinates
    w_user <- diff(x_user[1:2]) # plot width in user coordinates
    h_user <- diff(x_user[3:4]) # plot height in user coordinates
    inch_per_xstep <- w_inch / w_user
    inch_per_ystep <- h_inch / h_user
    ipu_ratio <- inch_per_ystep / inch_per_xstep
    ipu_ratio
}

grconvertW <- function(x, from = "user", to = "inches") {
    grconvertX(x, from, to) - grconvertX(0, from, to)
}

grconvertH <- function(y, from = "user", to = "inches") {
    grconvertY(y, from, to) - grconvertY(0, from, to)
}

# Datasets #####

mkdat_sim2 <- function() {
    set.seed(815)
    nA <- 4
    nB <- 4
    n <- nA + nB
    sim2_data <- data.frame(
        classes    = c(rep("A", nA), rep("B", nB)),
        water      = rnorm(n, 8, 0.2),
        ethanol    = rnorm(n, 2, 0.1),
        aceticacid = c(rnorm(nA, 2, 0.1), rnorm(nB, 3, 0.1))
    )
}
