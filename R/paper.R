# Figures #####

cacheenv <- environment()

#' @noRd
#' @title Plot typical NMR Challenges
#' @param s Scaling Factor for the created figure. By making the figure `s` times bigger than the textwidth you can make the font, lines, etc. `s` times smaller, thinner, etc. because it will get scaled down when included in the LaTeX document.
#' @param w Figure width. Should be the textwidth of the Latex document in inches.
#' @param h Figure height. Should be the textheight of the LaTeX document in inches minus some space for the fiure caption. Example: if the textheight is 9.5 inches, h = 8 could be a good choice.
#' @author Tobias Schmidt, 2024-2025: initial version.
#' @examples
#' spectra <- sim[1:3]
#' decons <- deconvolute(sim[1:3])
#' aligns <- align(decons, maxShift = 100, maxCombine = 10)
#' mkfig_nmr_challenges(spectra, decons, aligns)
mkfig_nmr_challenges <- function(spectra = metabodecon::sim[1:3],
                                 decons = deconvolute(spectra),
                                 aligns = align(decons, maxShift = 100, maxCombine = 10),
                                 show = FALSE,
                                 store = TRUE,
                                 path = tmpdir(),
                                 name = c("challenges"),
                                 exts = c("pdf", "svg"),
                                 s = 1,
                                 w = 5.45,
                                 h = 8) {
    if (show) {
        init_dev(s, w, h)
        fill_dev()
        plot_nmr_challenges(spectra, decons, aligns)
    }
    if (store) {
        R.devices::devEval(
            path = path,
            name = name,
            type = exts,
            width = w * s,
            height = h * s,
            expr = plot_nmr_challenges(spectra, decons, aligns)
        )
    }
}

# Plot #####

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
test_plot_nmr_challenges <- function(interactive = TRUE) {
    spectra <- metabodecon::sim[1:3]
    decons <- deconvolute(spectra)
    aligns <- align(decons, maxShift = 100, maxCombine = 10)
    if (interactive) {
        .GlobalEnv$spectra <- spectra
        .GlobalEnv$decons <- decons
        .GlobalEnv$aligns <- aligns
        catf("%sInitialized:%s spectra, decons, aligns\n", esc$green, esc$reset)
        catf("%sTestcall:%s plot_nmr_challenges(spectra, decons, aligns)", esc$green, esc$reset)
    } else {
        plot_nmr_challenges(spectra, decons, aligns)
    }
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_nmr_challenges <- function(spectra = metabodecon::sim[1:3],
                                decons = deconvolute(spectra),
                                aligns = align(decons)) {
    local_par(mfrow = c(5, 2), mar = c(0, 0, 2, 0))
    plot_1_nmr_experiment()
    plot_2_raw_fids()
    plot_3_blackbox_preprocessing()
    plot_4_raw_spectra(spectra)
    plot_5_deconvolutions(decons)
    plot_6_deconvoluted_spectra(decons)
    plot_7_alignments(aligns)
    plot_8_aligned_spectra(aligns)
    plot_9_annotations(aligns)
    plot_10_annotated_spectra(aligns)
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_1_nmr_experiment <- function() {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "NMR Experiment")
    marbox()
    text(0, 1, "Samples", adj = c(0, 1))
    draw_vial(x1= .00, y1 = .48, h = .3)
    draw_vial(x1= .16, y1 = .25, h = .3)
    draw_vial(x1= .08, y1 = .02, h = .3)
    text(0.5, 1, "NMR Spectrometer", adj = c(0, 1))
    draw_nmr_spectrometer(x1 = .5, y1 = .0, h = .8)
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_2_raw_fids <- function() {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "FID Signals")
    marbox()
    box()
    # ndc <- npc_to_ndc()
    # fig_height <- diff(ndc[3:4])
    # sub_height <- fig_height / 3
    # local_par(mar = c(0, 0, 2, 0))
    # with_fig(fig = c(0.50, 1.00, ndc[3] + 2 * sub_height, ndc[3] + 3 * sub_height), plot_2_raw_fid(1))
    # with_fig(fig = c(0.50, 1.00, ndc[3] + 1 * sub_height, ndc[3] + 2 * sub_height), plot_2_raw_fid(2))
    # with_fig(fig = c(0.50, 1.00, ndc[3] + 0 * sub_height, ndc[3] + 1 * sub_height), plot_2_raw_fid(3))
    ndc <- npc_to_ndc()
    fig_width  <- diff(ndc[1:2])
    fig_height <- diff(ndc[3:4])
    sub_height <- fig_height / 3
    with_fig(fig = npc_to_ndc(c(0, 1, 0/3, 1/3)), plot_2_raw_fid(1))
    with_fig(fig = npc_to_ndc(c(0, 1, 1/3, 2/3)), plot_2_raw_fid(2))
    with_fig(fig = npc_to_ndc(c(0, 1, 2/3, 3/3)), plot_2_raw_fid(3))
    mtext2(1, "Time")
    mtext2(2, "Signal Intensity")
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_2_raw_fid <- function(seed = 1) {
    set.seed(seed)
    ymar <- abs(grconvertH(0.1, "nfc", "lines"))
    local_par(mar = c(ymar, 0, ymar, 0))
    time <- seq(0, 1, by = 0.001)
    decay_constant <- 5 + rnorm(1, 0, 0.5)
    frequency <- 50 + rnorm(1, 0, 2)
    fid_signal <- exp(-decay_constant * time) * sin(2 * pi * frequency * time)
    plot(time, fid_signal, type = "l", xlab = "", ylab = "", axes = FALSE)
    marbox()
}

#' @noRd
#' @description
#' Plots the yf spectra on the left side, with an arrow going into a
#' pipeline. In the pipeline, the following steps are written as text:,
#' Apodization, Zero-Filling, Fourier Transform, Phase Correction, Baseline
#' Correction, Referencing
#' @param cex Character Expansion.
#' @param lex Line Height Expansion. Must be greater than cex.
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_3_preprocessing <- function(cex = par("cex"),
                                 lex = cex * 1.5) {
    local_par(mar = c(1, 0, 2, 0))
    plot_empty(main = "Preprocessing")
    marbox()
    ndc <- npc_to_ndc()
    x1 <- ndc[1]; x2 <- ndc[2]
    y1 <- ndc[3]; y2 <- ndc[4]
    W <- diff(ndc[1:2]); w <- W / 100
    H <- diff(ndc[3:4]); h <- H / 100
    local_par(mar = c(0, 0, 0, 0))
    with_fig(fig = c(x1 + 10*w, x1 + 40*w, y1 + 66*h, y1 + 99*h), plot_2_raw_fid(1))
    with_fig(fig = c(x1 + 10*w, x1 + 40*w, y1 + 33*h, y1 + 66*h), plot_2_raw_fid(2))
    with_fig(fig = c(x1 + 10*w, x1 + 40*w, y1 + 00*h, y1 + 33*h), plot_2_raw_fid(3))
    local_par(mar = c(0, 0, 0, 0))
    local_fig(fig = c(x1 + 50 * w, x1 + 90 * w, y1 + 10 * h, y1 + 90 * h))
    plot_empty()
    lh <- grconvertH(1, "chars", "user") * lex # Line Height
    nl <- 1 / lh # Number of Lines per Subfigure
    cl <- (nl / 2) + 0.5 # Center Line
    if (lh > 1/7) warning("Text too large. Reduce `cex` or `lex`.")
    rect(0.00, 1 - 2 * lh, 1.00, 1 - 0 * lh, col = "lightyellow")
    rect(0.08, 0 + 2 * lh, 0.92, 1 - 2 * lh, col = "lightyellow")
    rect(0.00, 0 + 0 * lh, 1.00, 0 + 2 * lh, col = "lightyellow")
    mtext("Apodization",        3, (cl + 3) * lex, cex = cex)
    mtext("Zero-Filling",       3, (cl + 2) * lex, cex = cex)
    mtext("Fourier Transform",  3, (cl + 1) * lex, cex = cex)
    mtext("Phase Correction",   3, (cl + 0) * lex, cex = cex)
    mtext("Baseline Correction",3, (cl - 1) * lex, cex = cex)
    mtext("Referencing",        3, (cl - 2) * lex, cex = cex)
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_3_blackbox_preprocessing <- function(cex = par("cex"),
                                          lex = cex * 1.5) {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "Preprocessing")
    marbox()
    box()
    lh <- grconvertH(1, "chars", "user") * lex # Line Height in R does NOT depend on `cex`
    nl <- 1 / lh # Number of Lines per Subfigure
    cl <- (nl / 2) + 0.5 # Center Line
    if (7 * lh > 1) warning("Text does not fit in plot area. Reduce `cex` or `lex`.")
    rect(0, 0, 1, 1, col = "darkgrey")
    m3text <- function(line, text) mtext(text, line = (-cl + line) * lex, cex = cex)
    m3text(+3, "Apodization")
    m3text(+2, "Zero-Filling")
    m3text(+1, "Fourier Transform")
    m3text(+0, "Phase Correction")
    m3text(-1, "Baseline Correction")
    m3text(-2, "Referencing")
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_4_raw_spectra <- function(spectra) {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "Raw Data in Frequency Domain")
    marbox()
    ndc <- npc_to_ndc()
    fig_width  <- diff(ndc[1:2])
    fig_height <- diff(ndc[3:4])
    sub_height <- fig_height / 3
    box()
    with_fig(fig = npc_to_ndc(c(0, 1, 0/3, 1/3)), plot_4_raw_spectrum(spectra[[1]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 1/3, 2/3)), plot_4_raw_spectrum(spectra[[2]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 2/3, 3/3)), plot_4_raw_spectrum(spectra[[3]]))
    mtext2(1, "Chemical Shift")
    mtext2(2, "Signal Intensity")
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_4_raw_spectrum <- function(spectrum) {
    args <- draw_spectrum_plain_args
    args$obj <- spectrum
    do.call(draw_spectrum, args)
    with_par(list(mar = c(0, 0, 0, 0)), box())
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_5_deconvolutions <- function(decons) {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "Deconvolution")
    marbox()
    ndc <- npc_to_ndc()
    fig_width  <- diff(ndc[1:2])
    fig_height <- diff(ndc[3:4])
    sub_height <- fig_height / 3
    box()
    with_fig(fig = npc_to_ndc(c(0, 1, 0/3, 1/3)), plot_5_deconvolution(decons[[1]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 1/3, 2/3)), plot_5_deconvolution(decons[[2]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 2/3, 3/3)), plot_5_deconvolution(decons[[3]]))
    mtext2(1, "Chemical Shift")
    mtext2(2, "Signal Intensity")
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_5_deconvolution <- function(decon) {
    args <- draw_spectrum_plain_args
    args$obj <- decon
    do.call(draw_spectrum, args)
    with_par(list(mar = c(0, 0, 0, 0)), box())
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_6_deconvoluted_spectra <- function(decons) {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "Deconvoluted Signal Integrals")
    marbox()
    ndc <- npc_to_ndc()
    fig_width  <- diff(ndc[1:2])
    fig_height <- diff(ndc[3:4])
    sub_height <- fig_height / 3
    box()
    with_fig(fig = npc_to_ndc(c(0, 1, 0/3, 1/3)), plot_6_deconvoluted_spectrum(decons[[1]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 1/3, 2/3)), plot_6_deconvoluted_spectrum(decons[[2]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 2/3, 3/3)), plot_6_deconvoluted_spectrum(decons[[3]]))
    mtext2(1, "Chemical Shift")
    mtext2(2, "Signal Integral")
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_6_deconvoluted_spectrum <- function(decon) {
    args <- draw_spectrum_plain_args
    args$obj <- decon
    args$si_line <- FALSE
    args$lc_lines$col <- NA
    do.call(draw_spectrum, args)
    with_par(list(mar = c(0, 0, 0, 0)), box())
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_7_alignments <- function(aligns) {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "Alignment")
    marbox()
    ndc <- npc_to_ndc()
    fig_width  <- diff(ndc[1:2])
    fig_height <- diff(ndc[3:4])
    sub_height <- fig_height / 3
    box()
    with_fig(fig = npc_to_ndc(c(0, 1, 0/3, 1/3)), plot_7_alignment(aligns[[1]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 1/3, 2/3)), plot_7_alignment(aligns[[2]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 2/3, 3/3)), plot_7_alignment(aligns[[3]]))
    mtext2(1, "Chemical Shift")
    mtext2(2, "Signal Intensity")
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_7_alignment <- function(align) {
    args <- draw_spectrum_plain_args
    args$obj <- align
    args$si_line <- FALSE
    args$lc_lines$col <- NA
    do.call(draw_spectrum, args)
    with_par(list(mar = c(0, 0, 0, 0)), box())
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_8_aligned_spectra <- function(aligns) {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "Aligned Signal Integrals")
    marbox()
    ndc <- npc_to_ndc()
    fig_width  <- diff(ndc[1:2])
    fig_height <- diff(ndc[3:4])
    sub_height <- fig_height / 3
    box()
    with_fig(fig = npc_to_ndc(c(0, 1, 0/3, 1/3)), plot_8_aligned_spectrum(aligns[[1]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 1/3, 2/3)), plot_8_aligned_spectrum(aligns[[2]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 2/3, 3/3)), plot_8_aligned_spectrum(aligns[[3]]))
    mtext2(1, "Chemical Shift")
    mtext2(2, "Signal Intensity")

}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_8_aligned_spectrum <- function(align) {
    args <- draw_spectrum_plain_args
    args$obj <- align
    args$si_line <- FALSE
    args$lc_lines <- FALSE
    args$lc_verts <- FALSE
    args$al_arrows <- FALSE
    args$al_lines$col <- NA
    do.call(draw_spectrum, args)
    with_par(list(mar = c(0, 0, 0, 0)), box())
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_9_annotations <- function(aligns) {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "Identification")
    marbox()
    ndc <- npc_to_ndc()
    fig_width  <- diff(ndc[1:2])
    fig_height <- diff(ndc[3:4])
    sub_height <- fig_height / 3
    box()
    with_fig(fig = npc_to_ndc(c(0, 1, 0/3, 1/3)), plot_9_annotation(aligns[[1]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 1/3, 2/3)), plot_9_annotation(aligns[[2]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 2/3, 3/3)), plot_9_annotation(aligns[[3]]))
    mtext2(1, "Chemical Shift")
    mtext2(2, "Signal Intensity")
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_9_annotation <- function(align) {
    # Name        Formula  Peaks  Description
    # Methane     CH4      1      Simple alkane with one equivalent proton environment.
    # Ethane      C2H6     1      Alkane with all protons equivalent, producing one peak.
    # Chloroform  CHCl3    1      Common NMR solvent with a single proton.
    # Glycine     C2H5NO2  2      Simplest amino acid with alpha proton and NH3+ group.
    # Alanine     C3H7NO2  3      Contains CH3, alpha proton, and NH3+ group.
    # Serine      C3H7NO3  4      Has CH2, alpha proton, OH, and NH3+ groups.
    # Glucose     C6H12O6  5-6    Key sugar in metabolism, altered in diabetes.
    # Lactate     C3H6O3   3      Diagnostic metabolite for lactic acidosis.
    args <- draw_spectrum_plain_args
    args$obj <- align
    args$si_line <- FALSE
    args$lc_lines <- FALSE
    args$lc_verts <- FALSE
    args$al_arrows <- FALSE
    args$al_lines$col <- NA
    do.call(draw_spectrum, args)
    draw_spectrum(align)
    coords_text <- align$lcpar$x0_al
    segments

    with_par(list(mar = c(0, 0, 0, 0)), box())
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
plot_10_annotated_spectra <- function(aligns) {
    local_par(mar = c(2, 2, 2, 1))
    plot_empty(main = "Metabolite Concentrations")
    marbox()
    ndc <- npc_to_ndc()
    fig_width  <- diff(ndc[1:2])
    fig_height <- diff(ndc[3:4])
    sub_height <- fig_height / 3
    box()
    with_fig(fig = npc_to_ndc(c(0, 1, 0/3, 1/3)), plot_8_aligned_spectrum(aligns[[1]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 1/3, 2/3)), plot_8_aligned_spectrum(aligns[[2]]))
    with_fig(fig = npc_to_ndc(c(0, 1, 2/3, 3/3)), plot_8_aligned_spectrum(aligns[[3]]))
    mtext2(1, "Chemical Shift")
    mtext2(2, "Signal Intensity")

}


# Draw #####

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
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

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
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
    # text((x1 + x2) / 2, (y1 + y2) / 2, "Sample")
    # text(x2 + console_w / 2 + .05, y1 + console_h / 2, "Console")
    # text((x1 + x2) / 2, y2 - h * 0.1, "Magnet")
}

# Plot Helpers #####

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
init_dev <- function(s = 1, w = 5.45, h = 8) {
    if (dev.cur() == 1) dev.new(
        width = w * s,
        height = h * s,
        rescale = "fixed"
    )
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
fill_dev <- function() {
    mfg <- par("mfg")
    while (mfg[1] != mfg[3] || mfg[2] != mfg[[4]]) {
        plot_empty()
        mfg <- par("mfg")
    }
    # plot_empty()
    # mfg[1:2] <- c(1, 1)
    # par(mfg = mfg)
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
get_ndc <- function() {
    usr <- par("usr")  # returns c(x1, x2, y1, y2)
    xndc <- grconvertX(usr[1:2], from = "user", to = "ndc")
    yndc <- grconvertY(usr[3:4], from = "user", to = "ndc")
    ndc <- c(x1 = xndc[1], x2 = xndc[2], y1 = yndc[1], y2 = yndc[2])
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
marbox <- function() {
    withr::with_par(list(mar = c(0,0,0,0)), box())
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
circles <- function(x, y, radii, inches = FALSE, add = TRUE, ...) {
    symbols(x, y, circles = radii, inches = inches, add = add, ...)
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
figsize <- function() {
    inches <- par("pin")
    cm <- inches * 2.54
    dpi <- 72
    dp <- inches * dpi
    cat("Size in incheses:", inches, "\n")
    cat("Size in centimeters:", cm, "\n")
    cat("Size in pixels:", dp, "\n")
}

#' @noRd
#' @title Inch per Unit Ratio
#' @author Tobias Schmidt, 2024-2025: initial version.
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

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
grconvertW <- function(x, from = "user", to = "inches") {
    grconvertX(x, from, to) - grconvertX(0, from, to)
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
grconvertH <- function(y, from = "user", to = "inches") {
    grconvertY(y, from, to) - grconvertY(0, from, to)
}

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
mtext2 <- function(side, text, line = 0.5, cex = par("cex"), ...) {
    mtext(text, side, line, cex = cex, ...)
}

# Datasets #####

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
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

# Constants

#' @noRd
#' @author Tobias Schmidt, 2024-2025: initial version.
draw_spectrum_plain_args <- local({
    decon_dark  <- "#6C8EBFFF" # transp("#6c8ebf", 1.0)
    decon_light <- "#DAE8FC80" # transp("#dae8fc", 0.5)
    align_dark  <- "#D6B656FF" # transp("#d6b656", 1.0)
    align_light <- "#FFF2CC80" # transp("#fff2cc", 0.5)
    arrow_dark  <- "#BEBEBEFF" # transp("grey", 1.0)
    list(
        foc_rgn  = NULL,
        foc_frac = c(0.48, 0.34),
        foc_only = TRUE,
        add      = TRUE,
        fig_rgn  = NULL,
        main     = NULL,
        show     = TRUE,
        show_d2  = FALSE,
        truepar  = NULL,
        mar      = c(0, 0, 0, 0),
        # Spectrum Lines
        si_line  = list(show = TRUE),
        sm_line  = list(show = FALSE),
        d2_line  = list(show = FALSE),
        # Deconvolution Lines
        sp_line  = list(show = FALSE),
        lc_lines = list(show = TRUE, col = decon_dark, fill = decon_light),
        lc_verts = list(show = TRUE, col = decon_dark),
        tp_lines = list(show = FALSE),
        tp_verts = list(show = FALSE),
        # Alignment Lines
        al_line  = list(show = FALSE),
        al_lines = list(show = TRUE, col = align_dark, fill = align_light),
        al_verts = list(show = TRUE, col = align_dark),
        al_arrows = list(show = TRUE, hsf = 0.5),
        # Points
        cent_pts = list(show = FALSE),
        bord_pts = list(show = FALSE),
        norm_pts = list(show = FALSE),
        # Rectangles
        bg_rect  = list(show = FALSE),
        foc_rect = list(show = FALSE),
        lc_rects = list(show = FALSE),
        tp_rects = list(show = FALSE),
        # Axes
        bt_axis  = list(show = FALSE),
        lt_axis  = list(show = FALSE),
        tp_axis  = list(show = FALSE),
        rt_axis  = list(show = FALSE),
        bt_text  = list(show = FALSE),
        lt_text  = list(show = FALSE),
        tp_text  = list(show = FALSE),
        rt_text  = list(show = FALSE),
        # Legend
        lgd      = list(show = FALSE)
    )
})
