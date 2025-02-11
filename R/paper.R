# Figures #####

mkfig_nmr_challenges <- function() {
    withr::local_pdf("challenges.pdf", width = 5.45, height = 8)
    plot_nmr_challenges()
}

# Plot #####

plot_nmr_challenges <- function() {
    withr::local_par(mar = c(0, 0, 0, 0), mfrow = c(5, 2))

    newplot("NMR Experiment")
    draw_vial(x1= 10, x2 = 15, y1 = 60, y2 = 90)
    draw_vial(x1= 20, x2 = 25, y1 = 30, y2 = 60)
    draw_vial(x1= 10, x2 = 15, y1 = 10, y2 = 40)

    # # Label the vials group
    # text(1.5, 4.5, "Blood Plasma Vials", cex = 1)

    # # Draw an arrow to the right
    # arrows(2.5, 2, 5.5, 2, length = 0.1, lwd = 2)

    # # Draw NMR spectrometer box
    # rect(6, 1.5, 9, 2.5, col = 'gray90', border = 'black')
    # text(7.5, 2, "NMR\nSpectrometer", cex = 1)
    newplot("Raw Data in Time Domain")
    newplot("Preprocessing")
    newplot("Raw Data in Frequency Domain")

    # 3.1 Deconvolution
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
        sm_line = FALSE
    )
    box(bty = "o")
    mtext("Deconvolution", side = 3, line = -1.2)

    # 3.2 Deconvoluted Signal Intensities
    newplot("Deconvoluted Signal Intensities")
    newplot("Alignment")
    newplot("Aligned Signal Intensities")
    newplot("Identification")
    newplot("Metabolite Concentrations")
}

# Draw #####

draw_vial <- function(x1 = 10, x2 = 15, y1 = 10, y2 = 40) {
    # See `misc/sketches/blood_vial.drawio` for a sketch of the vial
    x <- seq(x1, x2, length.out = 13)
    y <- seq(y1, y2, length.out = 37)
    w <- x[13] - x[1]
    h <- y[37] - y[1]
    cross_width_user    <- x[9] - x[5]
    cross_width_inch    <- grconvertW(cross_width_user, "user", to = "inches")
    cross_height_inch   <- cross_width_inch
    cross_height_user   <- grconvertH(cross_width_inch, "inches", to = "user")
    vial_height_user    <- y[29] - y[5]
    if (2 * cross_height_user > vial_height_user) {
        stop("Sticker height > vial height. Increase height/width ratio of vial.")
    }
    sticker_bottom <- y[17] - cross_height_user
    sticker_top    <- y[17] + cross_height_user
    symbols(x[7], y[5], circles = x[7] - x[3], inches = FALSE, add = TRUE, bg = "red") # bottom
    rect(x[3], y[5],  x[11], y[29], col = 'red', border = NA) # vial
    segments(x[3], y[5], x[3], y[29], col = 'black') # vial
    segments(x[3], y[29], x[11], y[29], col = 'black') # vial
    segments(x[11], y[5], x[11], y[29], col = 'black') # vial
    rect(x[3], sticker_bottom,  x[11], sticker_top,  col = 'white', border = 'black') # sticker
    rect(x[1], y[29], x[13], y[37], col = 'white', border = 'black') # cap
    symbols(x[7], y[17], circles = x[7] - x[4], inches = FALSE, add = TRUE, bg = "red") # circle
    rect(x[9], y[17] - cross_height_user / 4,  x[5], y[17] + cross_height_user / 4, col = 'white', border = NA) # vial
    rect(x[6], y[17] - cross_height_user / 2,  x[8], y[17] + cross_height_user / 2, col = 'white', border = NA) # vial
}


# Plot Helpers #####

newplot <- function(text,
                    xlim = c(0, 100),
                    ylim = c(0, 100),
                    axes = FALSE,
                    ann = FALSE,
                    line = -1.2,
                    side = 3,
                    ...) {
    plot_empty(xlim = xlim, ylim = ylim, axes = axes, ann = ann)
    box(bty = "o")
    mtext(text, line = line, side = side)
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
