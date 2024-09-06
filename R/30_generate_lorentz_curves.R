# Exported #####

#' @export
#' @title Deconvolute the Sim Dataset
#' @description The simulated 'Sim' dataset is much smaller than real NMR spectra (1309 datapoints instead of 131072) and lacks a water signal. This makes it ideal for use in examples or as a default value for functions. However, the default values for `sfr`, `wshw`, and `delta` in the `generate_lorentz_curves()` function are not optimal for this dataset. To avoid repeatedly defining the optimal parameters in examples, this function is provided to deconvolute the "Sim" dataset with suitable parameters. The actual parameters used for the deconvolution can be found in the 'Examples' section.
#' @param name Name of the spectra to deconvolute. Use `'bruker/sim'` to deconvolute all spectra of the Sim dataset. Use `'bruker/sim_subset'` to deconvolute only the first two spectra of the Sim dataset. Use `'bruker/sim/sim_xx'`, with `xx` being a number between `01` and `16`, to deconvolute a single spectrum of the Sim dataset. The provided name is passed to `metabodecon_file()` to get the path to the spectra.
#' @param spectra List of spectra or path to spectra directory. If `NULL`, the spectra are read from the path obtained by passing `name` to `metabodecon_file()`.
#' @param ask Whether to ask for user input during the deconvolution process.
#' @param verbose Whether to print log messages during the deconvolution process.
#' @return A list representing the deconvoluted spectra. For details see the return value of `generate_lorentz_curves()`.
#' @examples
#' # Detailed version:
#' name = "bruker/sim_subset"
#' path <- metabodecon_file(name)
#' decons1 <- generate_lorentz_curves(
#'     data_path = path,    # Path to directory containing spectra
#'     sfr = c(3.42, 3.58), # Borders of signal free region (SFR) in ppm
#'     wshw = 0,            # Half width of water signal (WS) in ppm
#'     delta = 0.1,         # Threshold for peak filtering
#'     ask = FALSE,         # Don't ask for user input
#'     verbose = FALSE      # Suppress status messages
#' )
#'
#' # Short version with `glc_sim()`:
#' decons2 <- glc_sim("bruker/sim_subset")
#'
#' # Comparison of results:
#' all.equal(decons1, decons2)
glc_sim <- function(name = "bruker/sim_subset", spectra = NULL, ask = FALSE, verbose = FALSE) {
    if (is.null(spectra)) {
        pattern <- "^bruker/(sim|sim_subset|sim/sim_0[1-9]|sim/sim_1[0-6])$"
        errmsg <- "Invalid name. Use 'bruker/sim', 'bruker/sim_subset' or 'bruker/sim/sim_xx' with xx being a number between 01 and 16."
        if (!(is.character(name) && length(name) == 1 && grepl(pattern, name))) stop(errmsg)
        spectra <- metabodecon_file(name)
    }
    generate_lorentz_curves(
        data_path = spectra, # List of spectra or path to spectra directory
        sfr = c(3.42, 3.58), # Borders of signal free region (SFR) in ppm
        wshw = 0,            # Half width of water signal (WS) in ppm
        delta = 0.1,         # Configure threshold for peak filtering
        ask = ask,           # Don't ask for user input
        verbose = verbose    # Suppress status messages
    )
}

# Private Helpers #####

deconvolute_spectrum <- function(spec, smopts, delta, nfit, nfiles, force) {
    logf("Starting deconvolution of %s", spec$name)
    spec <- rm_water_signal(spec)
    spec <- rm_negative_signals(spec)
    spec <- smooth_signals(spec, reps = smopts[1], k = smopts[2])
    spec <- find_peaks(spec)
    spec <- filter_peaks(spec, delta, force)
    spec <- fit_lorentz_curves(spec, nfit)
    spec <- add_return_list(spec, nfiles, spec$name)
    logf("Finished deconvolution of %s", spec$name)
    spec
}

rm_water_signal <- function(spec) {
    logf("Removing water signal")
    y <- spec$y_scaled
    left <- spec$wsr$left_dp
    right <- spec$wsr$right_dp
    y[right:left] <- 0.01 / spec$sf[2] # Order, i.e. `right:left` instead of `left:right`, is important here, because `right` and `left` are floats. Example: `right <- 3.3; left <- 1.4` ==> `right:left == c(3.3, 2.3)` and `left:right == c(1.4, 2.4)`.
    spec$y_nows <- y
    spec
}

rm_negative_signals <- function(spec) {
    logf("Removing negative signals")
    if (is.null(spec$y_nows)) stop("Water signal not removed yet. Please call `rm_water_signal()` first.")
    spec$y_pos <- abs(spec$y_nows)
    spec
}

#' @inherit smooth_signals
#' @details New and fast version for smoothing of signals. Implements the same algorithm as [smooth_signal_v12()] using different R functions (e.g. [stats::filter()]), causing a massive speedup but also numeric differences compared to the old version.
#' @noRd
smooth_signals_v20 <- function(spec, reps = 2, k = 5) {
    if (k %% 2 == 0) stop("k must be odd")
    Z <- vector("list", length = reps)
    y <- spec$y_pos
    n <- length(y)
    for (i in seq_len(reps)) {
        filter <- rep(1 / k, k)
        z <- stats::filter(y, filter, sides = 2) # (1)
        q <- (k - 1) / 2 # (2)
        for (j in seq_len(q)) {
            z[j] <- mean(y[1:(q + j)]) # (3)
            z[n - j + 1] <- mean(y[(n - q - j + 1):n]) # (4)
        }
        y <- Z[[i]] <- as.numeric(z)
        # Calling (1) gives NAs at both sides of vector, as there are not enough values for the moving average. The number of NAs at each side is given by (2). Example: if n==100 and k==5, then q==2, so z[1]==NA, z[2]==NA, z[99]==NA and z[100]==NA. To stay backwards compatible, these values must be filled with the mean of the values that are available. To do so, we iterate from 1:q, i.e. j==1 and j==2 and set
        # z[1]   <- mean(y[1:3])    # 3 == 2+1 == q+j            # (3)
        # z[2]   <- mean(y[1:4])    # 4 == 2+2 == q+j            # (3)
        # z[99]  <- mean(y[97:100]) # 97 == 100-2-2+1 == n-q-j+1 # (4)
        # z[100] <- mean(y[98:100]) # 98 == 100-2-1+1 == n-q-j+1 # (4)
        # Note: we could also think of leaving the NAs as they are, which would be more correct I think and even faster, but would break compatibility with the old version completely. So not even `all.equal(v1, v2)` would be TRUE anymore.
    }
    spec$Z <- Z
    spec$y_smooth <- Z[[reps]]
    spec
}

#' @title Smooth signal intensities using a moving average
#' @description Smoothes signal intensities by applying a [moving average](https://en.wikipedia.org/wiki/Moving_average) filter with a window size of k.
#' @param spec A list representing the spectrum, which should include the scaled signal intensities, after removal of the water artefact and negative values (`spec$y_pos`).
#' @param reps The number of times to apply the moving average.
#' @param k The number of points within the moving average window. Must be odd, so the smoothed point is in the middle of the window.
#' @return A numeric vector of the smoothed values.
#' @details Old and slow version producing the same results as the implementation within `deconvolution` from `MetaboDecon1D_deconvolution.R`.
#' @noRd
smooth_signals <- function(spec, reps = 2, k = 5, verbose = TRUE) {
    if (verbose) logf("Smoothing signals")
    if (k %% 2 == 0) stop("k must be odd")
    n <- length(spec$y_pos) # number of data points in total
    ws <- floor(k / 2) # window size left/right of center
    ct <- seq_len(n) # center positions
    lb <- pmax(ct - ws, 1) # left borders
    rb <- pmin(ct + ws, n) # right borders
    nw <- rb - lb + 1 # number of data points in window
    Z <- data.frame(spec$y_pos)
    for (j in seq_len(reps)) {
        zj <- Z[[j]]
        zk <- sapply(seq_len(n), function(i) sum(zj[lb[i]:rb[i]]))
        zk <- (1 / nw) * zk
        Z[[j + 1]] <- zk
    }
    spec[c("Z", "y_smooth")] <- list(Z[, -1], Z[, reps + 1])
    spec
}

#' @noRd
#' @title Filter Peaks with Low Scores Outside Signal-Free Region
#' @description Calculates a score for each peak in the spectrum. Peaks with scores that are lower than `mu + delta * sigma` are filtered out, with `mu` being the mean and `sigma` the standard deviation of the scores of peaks within the signal-free region (SFR).
#' @param spec A list containing the spectrum data, including peaks and their scores, and the signal-free region definition.
#' @param delta A numeric value specifying how many standard deviations `s` a score needs to be above `mu` to not get filtered out. Here `s` denotes the standard deviation of peaks scores in the signal free region (SFR) and `mu` denotes the average peak score within the SFR.
#' @param force If no peaks are found in the SFR, the function stops with an error message by default. If `force` is TRUE, the function instead proceeds without filtering any peaks, potentially increasing runtime.
#' @return Returns the modified `spec` list with the `peak` component updated to indicate which peaks are considered significant based on their score relative to the SFR and the `delta` parameter. Peaks within the SFR are marked with a specific region code ("sfrl" for peaks in the left SFR and "sfrr" for peaks in the right SFR). Peaks in the normal region have the region-code "norm".
#' @details The function first identifies peaks within the SFR by comparing their center positions against the SFR boundaries. If peaks are found within the SFR, it calculates the mean and standard deviation of their scores to establish a filtering threshold. Peaks with scores below this threshold are considered low and filtered out. If no peaks are found within the SFR and `force` is FALSE, the function stops and issues an error message. If `force` is TRUE, the function proceeds without filtering, potentially increasing runtime.
#' @examples
#' spec <- list(
#'     peak = list(score = c(1, 5, 2), center = c(1, 2, 3)),
#'     sdp = c(3, 2, 1),
#'     sfr = list(left_sdp = 2.8, right_sdp = 1.2)
#' )
#' rm3 <- filtered_spec <- filter_peaks(spec)
#' rm2 <- filtered_spec <- filter_peaks(spec, delta = 1)
filter_peaks <- function(spec, delta = 6.4, force = FALSE) {
    logf("Removing peaks with low scores")
    peak_score <- spec$peak$score
    peak_defined <- !is.na(spec$peak$left) & !is.na(spec$peak$center) & !is.na(spec$peak$right)
    in_left_sfr <- spec$sdp[spec$peak$center] >= spec$sfr$left_sdp
    in_right_sfr <- spec$sdp[spec$peak$center] <= spec$sfr$right_sdp
    in_sfr <- in_left_sfr | in_right_sfr
    if (any(in_sfr)) {
        mu <- mean(peak_score[in_sfr])
        sigma <- sd(peak_score[in_sfr])
    } else {
        if (force) {
            logf("No signals found in signal free region. This is a clear indiciation that the deconvolution parameters are not set correctly. Continuing anyways without dynamic peak filtering, because `force` is TRUE. Note that this might increase runtime drastically.")
            mu <- 0
            sigma <- 0
        } else {
            stop("No signals found in signal free region. Please double check deconvolution parameters.")
        }
    }
    spec$peak$high <- peak_defined & (peak_score > mu + delta * sigma)
    spec$peak$region <- "norm"
    spec$peak$region[in_left_sfr] <- "sfrl"
    spec$peak$region[in_right_sfr] <- "sfrr"
    logf("Removed %d peaks", sum(!spec$peak$high))
    spec
}

filter_peaks_v13 <- function(ppm, # x values in ppm
                             pc, # peak center indices
                             ps, # peak scores
                             sfrl, # signal free region left in ppm
                             sfrr, # signal free region right in ppm
                             delta = 6.4 # peak filter threshold parameter
                             ) {
    if (any(is.na(ps))) stop("Peak scores must never be NA")
    logf("Removing peaks with low scores")
    in_sfr <- which(ppm[pc] >= ppm || ppm[pc] <= ppm)
    if (!any(in_sfr)) stop("No signals found in signal free region. Please double check deconvolution parameters.")
    mu <- mean(ps[in_sfr])
    sd <- sd(ps[in_sfr])
    threshold <- mu + delta * sd
    above_threshold <- ps >= threshold
    logf("Removed %d peaks", sum(!above_threshold))
    list(in_sfr, above_threshold)
}

#' @noRd
#' @title create backwards compatible return list
#' @param spec Deconvoluted spectrum as returned by [refine_lorentz_curves_v12()].
#' @param nfiles Number of deconvoluted spectrum.
#' @param nam Name of current spectrum.
#' @return The input spectrum with an additional list element `ret` containing the deconvolution results in a backwards compatible format.
add_return_list <- function(spec = glc(), nfiles = 1, nam = "urine_1") {
    # Prepare some shortcuts for later calculations
    sdp <- spec$sdp
    ppm <- spec$ppm
    hz <- spec$hz
    dp <- spec$dp
    y_raw <- spec$y_raw
    y_smooth <- spec$y_smooth
    A <- spec$lcr$A
    lambda <- spec$lcr$lambda
    x_0 <- spec$lcr$w

    # Calculate spectrum superposition
    s <- sapply(sdp, function(x_i) sum(abs(A * (lambda / (lambda^2 + (x_i - x_0)^2))))) # takes approx. 2.2 seconds for urine_1
    s_normed <- s / sum(s)

    # Calculate MSE_normed and MSE_normed_raw
    y_normed <- y_smooth / sum(y_smooth)
    y_raw_normed <- y_raw / sum(y_raw)
    mse_normed <- mean((y_normed - s_normed)^2)
    mse_normed_raw <- mean((y_raw_normed - s_normed)^2)

    # Create and return list
    spec$ret <- list(
        number_of_files = nfiles,
        filename = nam,
        x_values = spec$sdp,
        x_values_ppm = spec$ppm,
        y_values = spec$y_smooth,
        spectrum_superposition = s,
        mse_normed = mse_normed,
        index_peak_triplets_middle = spec$peak$center[spec$peak$high],
        index_peak_triplets_left = spec$peak$right[spec$peak$high],
        index_peak_triplets_right = spec$peak$left[spec$peak$high],
        peak_triplets_middle = spec$ppm[spec$peak$center[spec$peak$high]],
        peak_triplets_left = spec$ppm[spec$peak$right[spec$peak$high]],
        peak_triplets_right = spec$ppm[spec$peak$left[spec$peak$high]],
        integrals = spec$lcr$integrals,
        signal_free_region = c(spec$sfr$left_sdp, spec$sfr$right_sdp),
        range_water_signal_ppm = spec$wsr$hwidth_ppm,
        A = spec$lcr$A,
        lambda = spec$lcr$lambda,
        x_0 = spec$lcr$w,
        # New fields (available since v1.2.0)
        y_values_raw = spec$y_raw,
        x_values_hz = spec$hz,
        mse_normed_raw = mse_normed_raw,
        x_0_hz = convert_pos(x_0, sdp, hz),
        x_0_dp = convert_pos(x_0, sdp, dp),
        x_0_ppm = convert_pos(x_0, sdp, ppm),
        A_hz = convert_width(A, sdp, hz),
        A_dp = convert_width(A, sdp, dp),
        A_ppm = convert_width(A, sdp, ppm),
        lambda_hz = convert_width(lambda, sdp, hz),
        lambda_dp = convert_width(lambda, sdp, dp),
        lambda_ppm = convert_width(lambda, sdp, ppm)
    )
    class(spec$ret) <- "decon"
    spec
}

#' @noRd
#' @title Calculate Lorentz Curve values
#' @description Calculates the values of a Lorentz Curve for a vector of input values `x`. The Lorentz Curve is defined as \mjeqn{A \cdot \frac{\lambda}{\lambda^2 + (x_i - x_0)^2}}.
#' @param x Numeric vector of x values.
#' @param x0 Center of the Lorentz curve.
#' @param A Amplitude parameter of the Lorentz curve.
#' @param lambda Half width at half height of the Lorentz curve.
#' @return Numeric vector of y values.
#' @examples
#' x <- 1:10
#' x0 <- 5
#' A <- 10
#' lambda <- 2
#' y1 <- lorentz(x, x0, A, lambda)
#' y2 <- A * pi * dcauchy(x, location = x0, scale = lambda)
#' stopifnot(all.equal(y1, y2))
lorentz <- function(x, x0, A, lambda) {
    # For details see [Wikipedia > Cauchy_distribution > #Properties_of_PDF](https://en.wikipedia.org/wiki/Cauchy_distribution#Properties_of_PDF), in particular the formula below sentence "In physics, a three-parameter Lorentzian function is often used".
    A * (lambda / (lambda^2 + (x - x0)^2))
}

#' @noRd
#' @description Before version 1.2 of 'metabodecon', the deconvolution functions `generate_lorentz_curves` and `MetaboDecon1D` wrote their output partially as txt files to their input folder. The txt files were named "SPEC_NAME parameter.txt" and "SPEC_NAME approximated_spectrum.txt". Since version 1.2 these txt files are no longer created by default, to prevent accidental modifications of the input folders. However, to stay backwards compatible, functions that used to read "SPEC_NAME parameter.txt" and "SPEC_NAME approximated_spectrum.txt" still accept them as input (e.g. `gen_feat_mat()`). I.e., in order to test this functionality, we still need a way to create the corresponding txt files (which is no longer done by `generate_lorentz_curves()`). That's the purpose of this function: it takes the output of `generate_lorentz_curves()` as input and creates the (now deprecated) "SPEC_NAME parameter.txt" and "SPEC_NAME approximated_spectrum.txt" in folder `outdir`.
write_parameters_txt <- function(decon, outdir, verbose = FALSE) {
    if (is_decon_list(decon)) {
        for (obj in decon) write_parameters_txt(obj, outdir)
        return(invisible(NULL))
    }
    name <- decon$filename
    w_new <- decon$x_0
    lambda_new <- decon$lambda
    A_new <- decon$A
    noise_threshold <- rep(NA, length(A_new))
    pardf <- data.frame(rbind(w_new, lambda_new, A_new, noise_threshold))
    supdf <- data.frame(t(decon$spectrum_superposition))
    parfile <- file.path(outdir, paste(name, "parameters.txt"))
    supfile <- file.path(outdir, paste(name, "approximated_spectrum.txt"))
    if (verbose) cat(sprintf("Creating: %s\n", parfile))
    utils::write.table(pardf, parfile, sep = ",", col.names = FALSE, append = FALSE)
    if (verbose) cat(sprintf("Creating: %s\n", parfile))
    utils::write.table(supdf, supfile, sep = ",", col.names = FALSE, append = FALSE)
}

