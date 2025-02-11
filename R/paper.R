# Figures #####

#' @noRd
#' @examples
#' @param s # Scaling Factor for the created figure. By making the figure `s` times bigger than the textwidth you can make the font, lines, etc. `s` times smaller, thinner, etc. because it will get scaled down when included in the LaTeX document.
#' dev.new(width = 5.45, height = 8)
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
plot_nmr_challenges <- function(lwd = 0.001) {
    withr::local_par(
        mar = c(0, 0, 0, 0),
        mfrow = c(5, 2),
        lwd = lwd
        # Assuming a h=15, w=10 inch device, we get h=3, w=2 inch per sub-figure.
        # Assuming further 100 DPI, we have h=300, w=200 pixels per sub-figure.
    )
    newplot("NMR Experiment")
    text(20, 75, "Samples", col = "red")
    draw_vial(x1= 10, y1 = 60)
    draw_vial(x1= 20, y1 = 30)
    draw_vial(x1= 10, y1 = 10)
    text(60, 75, "NMR Spectrometer")
    draw_nmr_spectrometer(x1 = 50, y1 = 10, h = 60)

    # # Draw an arrow to the right
    # arrows(2.5, 2, 5.5, 2, length = 0.1)

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

draw_vial <- function(x1 = 10, y1 = 10, h = 20, w = NULL) {
    # See `misc/sketches/blood_vial.drawio` for a sketch of the vial
    if (is.null(w)) w <- (h/3) * asp()
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
    if (2 * cross_height > vial_height) {
        stop("Sticker height > vial height. Increase height/width ratio of vial.")
    }
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
         col = "white", border = NA
    ) # vial
    rect(x[7] - cross_width / 6, y[17] - cross_height / 2,
         x[7] + cross_width / 6, y[17] + cross_height / 2,
        col = "white", border = NA
    ) # vial
}

draw_nmr_spectrometer <- function(x1 = 10, y1 = 20, h = 60, w = NULL) {
    if (is.null(w)) w <- ((2 * h) / 3) * asp()
    message(w)
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
    rect(x2 + 5, y1,
         x2 + console_w + 5,
         y1 + console_h,
         border = "black",
         col = "gray80"
    )

    # Draw the control panel (small rectangle near console)
    rect(
        x2 + 5,
        y1 + console_h + 5,
        x2 + console_w + 5,
        y1 + console_h + h * 0.166,
        border = "black",
        col = "gray70"
    )

    # Add labels
    text((x1 + x2) / 2, (y1 + y2) / 2, "Sample")
    text(x2 + console_w / 2 + 5, y1 + console_h / 2, "Console")
    text((x1 + x2) / 2, y2 - h * 0.05, "Magnet")
}

# Plot Helpers #####

circles <- function(x, y, radii, inches = FALSE, add = TRUE, ...) {
    symbols(x, y, circles = radii, inches = inches, add = add, ...)
}

newplot <- function(text,
                    xlim = c(0, 100),
                    ylim = c(0, 100),
                    axes = FALSE,
                    ann = FALSE,
                    line = -1.5,
                    side = 3,
                    ...) {
    plot_empty(xlim = xlim, ylim = ylim, axes = axes, ann = ann)
    box(bty = "o")
    mtext(text, line = line, side = side)
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

asp <- function() {
    d_inch <- par("pin") # plot dimension in inches
    w_inch <- d_inch[1]  # plot width in inches
    h_inch <- d_inch[2]  # plot height in inches
    x_user <- par("usr") # plot extremes (x1,x2,y1,y2) in user coordinates
    w_user <- diff(x_user[1:2]) # plot width in user coordinates
    h_user <- diff(x_user[3:4]) # plot height in user coordinates
    inch_per_xstep <- w_inch / w_user
    inch_per_ystep <- h_inch / h_user
    asp <- inch_per_ystep / inch_per_xstep
    asp
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
