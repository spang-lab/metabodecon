# Plot #####

plot_si_mat <- function(Y = generate_lorentz_curves_sim("bruker/sim"),
                        lgdcex = "auto",
                        main = NULL,
                        mar = par("mar")) {
    n <- nrow(Y) # Number of Spectra
    p <- ncol(Y) # Number of Datapoints
    cols <- rainbow(n) # Colors
    dpis <- seq_len(p) # Datapoint Indices
    spis <- seq_len(n) # Spectrum Indices
    ltxt <- paste("Spectrum", spis)
    xlab <- "Datapoint Number"
    ylab <- "Signal Intensity"
    xlim <- c(1, p)
    ylim <- c(0, max(Y))
    args <- named(x = NA, type = "n", xlim, ylim, xlab, ylab)
    opar <- par(mar = mar)
    on.exit(par(opar))
    do.call(plot, args)
    for (i in spis) lines(x = dpis, y = Y[i, ], col = cols[i])
    if (lgdcex == "auto") lgdcex <- 1 / max(1, log(n, base = 8))
    legend(x = "topright", legend = ltxt, col = cols, lty = 1, cex = lgdcex)
    if (!is.null(main)) title(main = main)
}

#' @noRd
#' @title Plots a spectrum simulated with [get_sim_params()]
#' @param simspec A simulated spectrum as returned by [get_sim_params()].
#' @examples
#' simspec <- get_sim_params()
#' plot_sim_spec(simspec)
plot_sim_spec <- function(simspec = get_sim_params()) {
    X <- simspec$X
    top_ticks <- seq(from = min(X$cs), to = max(X$cs), length.out = 5)
    top_labels <- round(seq(from = min(X$fq), to = max(X$fq), length.out = 5))
    line1_text <- sprintf("Name: %s", gsub("blood", "sim", simspec$filename))
    line2_text <- sprintf("Base: %s ", simspec$filename)
    line3_text <- sprintf("Range: 3.6-3.4 ppm")
    plot(
        x = X$cs, y = X$si_raw * 1e-6, type = "l",
        xlab = "Chemical Shift [PPM]", ylab = "Signal Intensity [AU]",
        xlim = c(3.6, 3.4), col = "black"
    )
    lines(x = X$cs, y = X$si_smooth, col = "blue")
    lines(x = X$cs, y = X$si_sim, col = "red")
    legend(
        x = "topleft", lty = 1,
        col = c("black", "blue", "red"),
        legend = c("Original Raw / 1e6", "Original Smoothed", "Simulated")
    )
    axis(3, at = top_ticks, labels = top_labels, cex.axis = 0.75)
    mtext(sprintf("Frequency [Hz]"), side = 3, line = 2)
    mtext(line1_text, side = 3, line = -1.1, col = "red", cex = 1, adj = 0.99)
    mtext(line2_text, side = 3, line = -2.1, col = "red", cex = 1, adj = 0.99)
    mtext(line3_text, side = 3, line = -3.1, col = "red", cex = 1, adj = 0.99)
}

plot_noise_methods <- function(siRND, siSFR, n = 300, start = 5000) {
    ymin <- min(c(min(siRND), min(siSFR)))
    ymax <- max(c(max(siRND), max(siSFR)))
    ylim <- c(ymin, ymax)
    opar <- par(mfrow = c(5, 1), mar = c(3, 4, 0, 1))
    on.exit(par(opar), add = TRUE)
    for (i in 1:5) {
        redT <- rgb(1, 0, 0, 0.1)
        bluT <- rgb(0, 0, 1, 0.1)
        idx <- ((i - 1) * n + 1):(i * n) + start
        ysiRND <- siRND[idx]
        ysiSFR <- siSFR[idx]
        plot(1:n, ysiRND, type = "n", ylim = ylim, ylab = "", xlab = "", xaxt = "n")
        axis(1, at = seq(1, n, by = 50), labels = idx[seq(1, n, by = 50)])
        points(1:n, ysiRND, col = "red", pch = 20)
        points(1:n, ysiSFR, col = "blue", pch = 20)
        lines(1:n, ysiRND, col = "red")
        lines(1:n, ysiSFR, col = "blue")
        lines(1:n, rep(0, n), col = "black", lty = 2)
        polygon(c(1:n, n:1), c(ysiRND, rep(0, n)), col = redT, border = NA)
        polygon(c(1:n, n:1), c(ysiSFR, rep(0, n)), col = bluT, border = NA)
        legend("topleft", NULL, c("RND", "SFR"), col = c("red", "blue"), lty = 1)
    }
}

# Align #####

#' @noRd
#' @author Initial version written by from Wolfram Gronwald as part of
#' `gen_feat_mat`. Moved into seperate function in 2024 by Tobias Schmidt.
read_decon_params_original <- function(data_path) {
    files <- list.files(data_path, ".txt", full.names = TRUE)
    num_spectra <- length(files) / 2
    spectrum_superposition <- vector(mode = "list", length = num_spectra)
    data_info <- vector(mode = "list", length = num_spectra)
    numerator_info <- 1
    numerator_spectrum <- 1
    for (i in 1:length(files)) {
        if (grepl(" parameters", files[i])) {
            data_info[[numerator_info]] <- as.matrix(data.table::fread(files[i], header = FALSE))
            numerator_info <- numerator_info + 1
        }
        if (grepl(" approximated_spectrum", files[i])) {
            spectrum_superposition[[numerator_spectrum]] <- as.vector(unlist(data.table::fread(files[i], header = FALSE)))[-1]
            numerator_spectrum <- numerator_spectrum + 1
        }
    }
    w <- vector(mode = "list", length = num_spectra)
    lambda <- vector(mode = "list", length = num_spectra)
    A <- vector(mode = "list", length = num_spectra)
    for (i in 1:num_spectra) {
        w[[i]] <- as.numeric(data_info[[i]][1, ][-1])
        lambda[[i]] <- as.numeric(data_info[[i]][2, ][-1])
        A[[i]] <- as.numeric(data_info[[i]][3, ][-1])
   }
    return(named(w, lambda, A, spectrum_superposition))
}

#' @noRd
#'
#' @title Align Signals using 'speaq'
#'
#' @description
#' Performs signal alignment across the individual spectra using the 'speaq'
#' package (Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S, Pieters L, et al.
#' (2018) speaq 2.0: A complete workflow for high-throughput 1D NMRspectra
#' processing and quantification. PLoS Comput Biol 14(3): e1006018.
#' https://www.doi.org/10.1371/journal.pcbi.1006018). The spectra deconvolution
#' process yields the signals of all spectra. Due to slight changes in
#' measurement conditions, e.g. pH variations, signal positions may vary
#' slightly across spectra. As a consequence, prior to further analysis signals
#' belonging to the same compound have to be aligned across spectra. This is the
#' purpose of the 'speaq' package.
#'
#' @param feat
#' Output of `gen_feat_mat()`.
#'
#' @param maxShift
#' Maximum number of points along the "ppm-axis" which a value can be moved by
#' speaq package e.g. 50. 50 is a suitable starting value for plasma spectra
#' with a digital resolution of 128K. Note that this parameter has to be
#' individually optimized depending on the type of analyzed spectra and the
#' digital resolution. For urine which is more prone to chemical shift
#' variations this value most probably has to be increased.
#'
#' @param spectrum_data
#' Output of `generate_lorentz_curves()`.
#'
#' @param si_size_real_spectrum
#' Number of real data points in your original spectra.
#'
#' @return
#' A matrix containing the aligned integral values of all spectra. Each row
#' contains the data of each spectrum and each column corresponds to one data
#' point. Each entry corresponds to the integral of a deconvoluted signal with
#' the signal center at this specific position after alignment by speaq.
#'
#' @author Wolfram Gronwald
#'
#' @examples
#' spectrum_data <- generate_lorentz_curves_sim("bruker/sim")
#' feat <- gen_feat_mat(spectrum_data)
#' maxShift <- 100
#' si_size_real_spectrum <- length(spectrum_data[[1]]$y_values)
#' after_speaq_mat <- speaq_align(feat, maxShift, spectrum_data, si_size_real_spectrum)
speaq_align_original <- function(feat = gen_feat_mat(spectrum_data),
                           maxShift = 50,
                           spectrum_data = generate_lorentz_curves_sim(),
                           si_size_real_spectrum = length(spectrum_data[[1]]$y_values)) {

    # Find reference spectrum, i.e the spectrum to which all others will be aligned to
    resFindRef <- speaq::findRef(feat$peakList)
    refInd <- resFindRef$refInd

    # Use overwritten CluPA function for multiple spectra
    data_result_aligned <- dohCluster(
        X = feat$data_matrix,
        peakList = feat$peakList,
        refInd = refInd,
        maxShift = maxShift,
        acceptLostPeak = FALSE,
        verbose = TRUE
    )
    # > str(data_result_aligned, 2)
    # List of 2
    # $ Y: num [1:16, 1:1309] 0.013 0.024 ... (matrix of aligned spectra)
    # $ new_peakList: List of 16              (aligned peaks of each spectrum)
    # ..$ sim_01: num [1:13] 350 416 457 ...
    # ..$ sim_02: num [1:11] 416 458 542 ...
    # ...

    # Aligned spectra values
    data_matrix_aligned <- data_result_aligned$Y
    # new peak list values
    new_peakList <- data_result_aligned$new_peakList
    # Recalculate peakList for the original format
    num_spectra <- dim(data_matrix_aligned)[1]
    peakList_result <- vector(mode = "list", length = num_spectra)
    for (i in 1:length(new_peakList)) {
        peakList_result[[i]] <- (dim(feat$data_matrix)[2]) - new_peakList[[i]]
    }

    # Generate feature matrix and integral matrix by using distance matrix
    # message("Warning: Generation of feature matrix and integral matrix have a high time duration (about 4h!)")
    # Calculate integral results for each positions
    integral_results <- vector(mode = "list", length = num_spectra)
    # Calculate integral_results
    for (spec in 1:length(peakList_result)) {
        for (ent in 1:length(peakList_result[[spec]])) {
            # Calculate integral value for each entry which is unequal 0
            if (peakList_result[[spec]][ent] != 0) {
                integral_results[[spec]][ent] <- feat$A[[spec]][ent] * (atan((-peakList_result[[spec]][ent] + si_size_real_spectrum) / feat$lambda[[spec]][ent]) - atan((-peakList_result[[spec]][ent]) / feat$lambda[[spec]][ent]))
            }
        }
    }

    # Generate a matrix containing the integrals of all spectra at their positions after alignment by speaq
    after_speaq_mat <- matrix(nrow = num_spectra, ncol = si_size_real_spectrum)
    # set ppm values for data points of matrix based on original data
    # Note signals will be shifted across data points, data points itself will keep their positions.
    colnames(after_speaq_mat) <- spectrum_data[[1]]$x_values_ppm
    mid_spec <- round(si_size_real_spectrum / 2)
    for (i in 1:num_spectra)
    {
        for (j in 1:length(peakList_result[[i]]))
        {
            k <- round(peakList_result[[i]][j])
            k1 <- mid_spec - (k - mid_spec) # mirror imaging of data at mid data point
            after_speaq_mat[i, k1] <- integral_results[[i]][j]
        }
    }
    return(after_speaq_mat)
}
