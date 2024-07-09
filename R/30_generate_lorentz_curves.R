# Exported #####

#' @export
#' @title Generate Lorentz Curves from NMR Spectra
#' @description Deconvolutes NMR spetra and generates a Lorentz curve for each detected signal within a spectra.
#' @param data_path Either the path to an existing directory containing measured NMR spectra or a dataframe with columns `ppm` (parts per million) and `si` (signal intensity) or a list of such dataframes.
#' @param file_format Format of the spectra files. Either `"bruker"` or `"jcampdx"`. Only relevant if `data_path` is a directory.
#' @param make_rds Store results as rds file on disk? Should be set to TRUE if many spectra are evaluated to decrease computation time.
#' @param expno The experiment number for the spectra files. E.g. `"10"`. Only relevant if `data_path` is a directory and `file_format` is `"bruker"`.
#' @param procno The processing number for the spectra. E.g. `"10"`. Only relevant if `data_path` is a directory and `file_format` is `"bruker"`.
#' @param nfit Number of iterations for the approximation of the parameters for the Lorentz curves.
#' @param wshw Half width of the water artefact in ppm.
#' @param sfr Row vector with two entries consisting of the ppm positions for the left and right border of the signal free region of the spectrum.
#' @param smopts Vector with two entries consisting of the number of smoothing iterations and the number of data points to use for smoothing (must be uneven).
#' @param delta Threshold value to distinguish between signal and noise.
#' @param sf Vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis.
#' @param ask  Whether to ask for user input during the deconvolution process. If set to FALSE, the provided default values will be used.
#' @param debug Whether to return additional intermediate results for debugging purposes.
#' @param nworkers Number of workers to use for parallel processing. If set to `"auto"`, the number of workers will be determined automatically. If set to a number greater than 1, the number of workers will be limited to the number of spectra or 1 if the operating system is Windows.
#' @return A list of deconvoluted spectra. Each list element contains a list with the following elements:
#' * `number_of_files`: Number of deconvoluted spectra.
#' * `filename`: Name of the analyzed spectrum.
#' * `x_values`: Scaled datapoint numbers (SDP). Datapoints are numbered in descending order going from N to 0, where N equals the . Scaled data point numbers are obtained by dividing these numbers by the x-axis scale factor `sf[1]`. I.e., for a spectrum with 131072 datapoints and a scale factor of 1000, the first scale datapoint has value 131.071 and the last one has value 0.
#' * `x_values_ppm`: The chemical shift of each datapoint in ppm (parts per million).
#' * `y_values`: The scaled signal intensity (SSI) of each datapoint. Obtained by reading the raw intensity values from the provided `data_path` as integers and dividing them by the y-axis scale factor `sf[2]`.
#' * `spectrum_superposition`: Scaled signal intensity obtained by calculating the sum of all estimated Lorentz curves for each data point.
#' * `mse_normed`: Normalized mean squared error. Calculated as \eqn{\frac{1}{n} \sum_{i=1}^{n} (z_i - \hat{z}_i)^2} where \eqn{z_i} is the normalized, smoothed signal intensity of data point i and \eqn{\hat{z}_i} is the normalized superposition of Lorentz curves at data point i. Normalized in this context means that the vectors were scaled so the sum over all data points equals 1.
#' * `index_peak_triplets_middle`: Datapoint numbers of peak centers.
#' * `index_peak_triplets_left`: Datapoint numbers of left peak borders.
#' * `index_peak_triplets_right`: Datapoint numbers of right peak borders.
#' * `peak_triplets_middle`: Chemical shift of peak centers in ppm .
#' * `peak_triplets_left`: Chemical shift of left peak borders in ppm .
#' * `peak_triplets_right`: Chemical shift of right peak borders in ppm .
#' * `integrals`: Integrals of the Lorentz curves.
#' * `signal_free_region`: Borders of the signal free region of the spectrum in scaled datapoint numbers. Left of the first element and right of the second element no signals are expected.
#' * `range_water_signal_ppm`: Half width of the water signal in ppm. Potential signals in this region are ignored.
#' * `A`: Amplitude parameter of the Lorentz curves. Provided as negative number to maintain backwards compatiblity with MetaboDecon1D. The area under the Lorentz curve is calculated as \eqn{A \cdot \pi}.
#' * `lambda`: Half width of the Lorentz curves in scaled data points. Provided as negative value to maintain backwards compatiblity with MetaboDecon1D. Example: a value of -0.00525 corresponds to 5.25 data points. With a spectral width of 12019 Hz and 131072 data points this corresponds to a halfwidth at half height of approx. 0.48 Hz. The corresponding calculation is: (12019 Hz / 131071 dp) * 5.25 dp.
#' * `x_0`: Center of the Lorentz curves in scaled data points.
#' @details First, an automated curvature based signal selection is performed. Each signal is represented by 3 data points to allow the determination of initial Lorentz curves. These Lorentz curves are then iteratively adjusted to optimally approximate the measured spectrum.
#' @examples
#' \donttest{
#' example_datasets <- download_example_datasets()
#' urine <- file.path(example_datasets, "bruker/urine")
#' urine_1 <- file.path(urine, "urine_1")
#' deconv_urine <- generate_lorentz_curves(urine, ask = FALSE, nfit = 1)
#' deconv_urine_1 <- generate_lorentz_curves(urine_1, ask = FALSE)[[1]]
#' }
generate_lorentz_curves <- function(data_path = pkg_file("example_datasets/bruker/urine"),
                                    file_format = "bruker",
                                    make_rds = FALSE,
                                    expno = 10,
                                    procno = 10,
                                    nfit = 10,
                                    wshw = 0.1527692,
                                    sfr = c(11.44494, -1.8828),
                                    smopts = c(2, 5),
                                    delta = 6.4,
                                    sf = c(1000, 1000000),
                                    ask = TRUE,
                                    debug = FALSE,
                                    nworkers = "auto") {

    # Read spectra and ask user for parameters
    spectra_ds <- read_spectra(data_path, file_format, expno, procno, raw = TRUE)
    spectra <- lapply(spectra_ds, convert_spectrum, sfx = sf[1], sfy = sf[2])
    adjno <- get_adjno(spectra, sfr, wshw, ask)
    spectra <- get_sfrs(spectra, sfr, ask, adjno)
    spectra <- get_wsrs(spectra, wshw, ask, adjno)

    # Deconvolute spectra
    n <- length(spectra)
    nams <- names(spectra)
    if (nworkers == "auto") nworkers <- min(ceiling(parallel::detectCores() / 2), length(spectra))
    starttime <- Sys.time()

    if (nworkers == 1) {
        logf("Starting deconvolution of %d spectra using 1 worker", n)
        spectra <- lapply(seq_len(n), function(i) {
            deconvolute_spectrum(nams[[i]], spectra[[i]], smopts, delta, nfit, n, debug)
        })
    } else {
        logf("Starting %d worker processes", nworkers)
        cl <- parallel::makeCluster(nworkers, outfile = "")
        on.exit(parallel::stopCluster(cl), add = TRUE)
        logf("Exporting required functions and data to workers")
        exportvars <- c("logf", "fg", "deconvolute_spectrum", "nams", "spectra", "smopts", "delta", "nfit", "n", "debug")
        parallel::clusterExport(cl, exportvars, envir = environment())
        logf("Starting deconvolution of %d spectra using %d workers", n, nworkers)
        spectra <- parallel::parLapply(cl, seq_len(n), function(i) {
            opts <- options(toscutil.logf.sep1 = sprintf(" CORE %d ", Sys.getpid()))
            on.exit(options(opts), add = TRUE)
            deconvolute_spectrum(nams[[i]], spectra[[i]], smopts, delta, nfit, n, debug)
        })
    }

    endtime <- Sys.time()
    duration <- endtime - starttime
    logf("Finished deconvolution of %d spectra in %s", n, format(duration))
    names(spectra) <- nams

    # Prepare, store and return results
    ret <- if (debug) spectra else lapply(spectra, function(s) s$ret)
    if (isTRUE(make_rds)) {
        rdspath <- file.path(data_path, "spectrum_data.rds")
        if (interactive()) {
            yes <- get_yn_input(sprintf("Save results as '%s'?", rdspath))
            if (yes) saveRDS(ret, rdspath)
        } else {
            logf("Skipping RDS save: confirmation required but not in interactive mode. For details see `help('generate_lorentz_curves')`.")
        }
    } else if (is.character(make_rds)) {
        cat("Saving results as", make_rds, "\n")
        saveRDS(ret, make_rds)
    }
    ret
}

# Private Helpers #####

deconvolute_spectrum <- function(nam, spec, smopts, delta, nfit, n, debug) {
    logf("Starting deconvolution of %s", nam)
    spec <- rm_water_signal_v12(spec)
    spec <- rm_negative_signals_v12(spec)
    spec <- smooth_signals_v12(spec, reps = smopts[1], k = smopts[2])
    spec <- find_peaks_v12(spec)
    spec <- filter_peaks_v12(spec, delta)
    spec <- fit_lorentz_curves(spec, nfit)
    spec <- add_return_list_v13(spec, n, nam, debug)
    logf("Finished deconvolution of %s", nam)
    spec
}

rm_water_signal_v12 <- function(spec) {
    logf("Removing water signal")
    y <- spec$y_scaled
    left <- spec$wsr$left_dp
    right <- spec$wsr$right_dp
    y[right:left] <- 0.01 / spec$sfy # Order, i.e. `right:left` instead of `left:right`, is important here, because `right` and `left` are floats. Example: `right <- 3.3; left <- 1.4` ==> `right:left == c(3.3, 2.3)` and `left:right == c(1.4, 2.4)`.
    spec$y_nows <- y
    spec
}

rm_negative_signals_v12 <- function(spec) {
    logf("Removing negative signals")
    if (is.null(spec$y_nows)) stop("Water signal not removed yet. Please call `rm_water_signal_v12()` first.")
    spec$y_pos <- abs(spec$y_nows)
    spec
}

#' @inherit smooth_signals_v12
#' @details New and fast version for smoothing of signals. Implements the same algorithm as [smooth_signal_v12()] using different R functions (e.g. [stats::filter()]), causing a massive speedup but also numeric differences compared to the old version.
#' @noRd
smooth_signals_v20 <- function(spec, reps = 2, k = 5) {
    if (k %% 2 == 0) stop("k must be odd")

    Z <- vector("list", length = reps)
    y <- spec$y_pos
    n <- length(y)

    for (i in 1:reps) {
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
        # z[100] <- mean(y[98:100]) # 98 == 100-2-1+1 == n-q-j+1 # (4)
        # z[99]  <- mean(y[97:100]) # 97 == 100-2-2+1 == n-q-j+1 # (4)
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
smooth_signals_v12 <- function(spec, reps = 2, k = 5) {
    logf("Smoothing signals")
    if (k %% 2 == 0) stop("k must be odd")
    Z <- vector("list", length = reps)
    y <- spec$y_pos
    n <- length(y)
    for (i in 1:reps) {
        z <- y
        for (j in 1:(n)) {
            left_border <- j - floor(k / 2)
            right_border <- j + floor(k / 2)
            if (left_border <= 0) {
                left_border <- 1
                z[j] <- (1 / right_border) * sum(y[left_border:right_border])
            } else if (right_border >= n) {
                right_border <- n
                z[j] <- (1 / (right_border - left_border + 1)) * sum(y[left_border:right_border])
            } else {
                z[j] <- (1 / k) * sum(y[left_border:right_border])
            }
        }
        y <- Z[[i]] <- as.numeric(z)
    }
    spec$Z <- Z
    spec$y_smooth <- Z[[reps]]
    spec
}

filter_peaks_v13 <- function(ppm, # x values in ppm
                             pc, # peak center indices
                             ps, # peak scores
                             sfrl, # signal free region left in ppm
                             sfrr, # signal free region right in ppm
                             delta = 6.4) { # threshold parameter to distinguish between "real" peaks from noise
    if (any(is.na(ps))) stop("Peak scores must never be NA")
    logf("Removing peaks with low scores")
    in_sfr <- which(ppm[pc] >= ppm || ppm[pc] <= ppm)
    mu <- mean(ps[in_sfr]) # mean (greek letter mu)
    sigma <- sd(ps[in_sfr]) # standard deviation (greek letter sigma)
    gt_tau <- (ps >= mu + delta * sigma) # greater than threshold value (greek letter tau)
    logf("Removed %d peaks", sum(!gt_tau))
    list(in_sfr, gt_tau) # gttho
}

#' @title create backwards compatible return list
#' @param spec Deconvoluted spectrum as returned by [refine_lorentz_curves_v12()].
#' @param n Number of deconvoluted spectrum.
#' @param nam Name of current spectrum.
#' @param debug Add debug info to the return list
#' @noRd
add_return_list_v13 <- function(spec = glc_v13(), n = 1, nam = "urine_1", debug = TRUE) {

    # Prepare some shortcuts for later calculations
    x <- spec$sdp
    p <- spec$ppm
    r <- spec$y_raw
    y <- spec$y_smooth
    n <- length(x)
    A <- spec$lcr$A
    lambda <- spec$lcr$lambda
    x_0 <- spec$lcr$w

    # Calculate spectrum superposition
    s <- sapply(x, function(x_i) sum(abs(A * (lambda / (lambda^2 + (x_i - x_0)^2))))) # takes approx. 2.2 seconds for urine_1
    s_normed <- s / sum(s)

    # Calculate MSE_normed and MSE_normed_raw
    y_normed <- y / sum(y)
    y_raw_normed <- r / sum(r)
    mse_normed <- mean((y_normed - s_normed)^2)
    mse_normed_raw <- mean((y_raw_normed - s_normed)^2)

    # Prepare variables for converting between units in the next step
    width_ppm <- p[1] - p[n]
    width_sdp <- x[1] - x[n]
    width_hz  <- spec$hz[1] - spec$hz[n]
    width_dp  <- spec$dp[1] - spec$dp[n]
    ref_ppm   <- p[1]
    ref_sdp   <- x[1]
    ref_hz    <- spec$hz[1]
    ref_dp    <- spec$dp[1]

    # Create and return list
    spec$ret <- list(
        number_of_files = n,
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
        x_0_hz  = convert_pos(x_0, width_sdp, width_hz, ref_sdp, ref_hz),
        x_0_dp  = convert_pos(x_0, width_sdp, width_dp, ref_sdp, ref_dp),
        x_0_ppm = convert_pos(x_0, width_sdp, width_ppm, ref_sdp, ref_ppm),
        A_hz  = convert_width(A, width_sdp, width_hz),
        A_dp  = convert_width(A, width_sdp, width_dp),
        A_ppm = convert_width(A, width_sdp, width_ppm),
        lambda_hz  = convert_width(lambda, width_sdp, width_hz),
        lambda_dp  = convert_width(lambda, width_sdp, width_dp),
        lambda_ppm = convert_width(lambda, width_sdp, width_ppm)
    )
    spec
}


# Private Deprecated #####

#' @title Generate Lorentz Curves from NMR Spectra
#' @description Deconvolutes NMR spetra and generates a Lorentz curve for each detected signal within a spectra.
#' @param data_path Either the path to an existing directory containing measured NMR spectra or a dataframe with columns `ppm` (parts per million) and `si` (signal intensity) or a list of such dataframes.
#' @param file_format Format of the spectra files. Either `"bruker"` or `"jcampdx"`. Only relevant if `data_path` is a directory.
#' @param make_rds Store results as rds file on disk? Should be set to TRUE if many spectra are evaluated to decrease computation time.
#' @param expno The experiment number for the spectra files. E.g. `"10"`. Only relevant if `data_path` is a directory and `file_format` is `"bruker"`.
#' @param procno The processing number for the spectra. E.g. `"10"`. Only relevant if `data_path` is a directory and `file_format` is `"bruker"`.
#' @param nfit Number of iterations for the approximation of the parameters for the Lorentz curves.
#' @param wshw Half width of the water artefact in ppm.
#' @param sfr Row vector with two entries consisting of the ppm positions for the left and right border of the signal free region of the spectrum.
#' @param smopts Vector with two entries consisting of the number of smoothing iterations and the number of data points to use for smoothing (must be uneven).
#' @param delta Threshold value to distinguish between signal and noise.
#' @param sf Vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis.
#' @param ask  Whether to ask for user input during the deconvolution process. If set to FALSE, the provided default values will be used.
#' @details First, an automated curvature based signal selection is performed. Each signal is represented by 3 data points to allow the determination of initial Lorentz curves. These Lorentz curves are then iteratively adjusted to optimally approximate the measured spectrum.
#' @examples
#' \dontrun{
#' xds_path <- download_example_datasets()
#' data_path <- file.path(xds_path, "bruker/urine")
#' file_format <- "bruker"
#' ask <- FALSE
#' nfit <- 2
#' generate_lorentz_curves_v12(data_path, file_format, ask = ask, nfit = nfit)
#' }
#' @noRd
generate_lorentz_curves_v12 <- function(data_path = file.path(download_example_datasets(), "bruker/urine"),
                                        file_format = "bruker",
                                        make_rds = FALSE,
                                        expno = 10,
                                        procno = 10,
                                        nfit = 10,
                                        wshw = 0.1527692,
                                        sfr = c(11.44494, -1.8828),
                                        smopts = c(2, 5),
                                        delta = 6.4,
                                        sf = c(1000, 1000000),
                                        ask = TRUE,
                                        debug = FALSE,
                                        nworkers = 2) {

    # Read spectra and ask user for parameters
    spectra_ds <- read_spectra(data_path, file_format, expno, procno, ask, sf, raw = TRUE)
    spectra <- lapply(spectra_ds, convert_spectrum, sfx = sf[1], sfy = sf[2])
    adjno <- get_adjno(spectra, sfr, wshw, ask)
    spectra <- get_sfrs(spectra, sfr, ask, adjno)
    spectra <- get_wsrs(spectra, wshw, ask, adjno)

    # Deconvolute spectra
    n <- length(spectra)
    nams <- names(spectra)
    spectra <- lapply(seq_len(n), function(i) {
        nam <- nams[i]
        logf("Starting deconvolution of %s", nam)
        spec <- spectra[[i]]
        spec <- rm_water_signal_v12(spec)
        spec <- rm_negative_signals_v12(spec)
        spec <- smooth_signals_v12(spec, reps = smopts[1], k = smopts[2])
        spec <- find_peaks_v12(spec)
        spec <- filter_peaks_v12(spec, delta)
        spec <- init_lorentz_curves_v12(spec)
        spec <- refine_lorentz_curves_v12(spec, nfit)
        spec <- add_return_list_v12(spec, n, nam, debug)
    })
    names(spectra) <- nams

    # Create return lists
    if (debug) {
        return(spectra)
    } else {
        rets <- lapply(spectra, function(spec) spec$ret)
        names(rets) <- nams
        return(rets)
    }
}

filter_peaks_v12 <- function(spec, delta = 6.4) {
    logf("Removing peaks with low scores")
    score <- spec$peak$score
    l <- which(spec$sdp[spec$peak$center] >= spec$sfr$left_sdp)
    r <- which(spec$sdp[spec$peak$center] <= spec$sfr$right_sdp)
    mu <- mean(score[c(l, r)])
    sigma <- sd(c(score[l], score[r]))
    spec$peak$high <- score >= mu + delta * sigma
    spec$peak$region <- "norm"
    spec$peak$region[l] <- "sfrl"
    spec$peak$region[r] <- "sfrr"
    logf("Removed %d peaks", sum(!spec$peak$high))
    spec
}

#' @noRd
#' @title create backwards compatible return list
#' @param spec Deconvoluted spectrum as returned by [refine_lorentz_curves_v12()].
#' @param n Number of deconvoluted spectrum.
#' @param nam Name of current spectrum.
#' @param debug Add debug info to the return list
add_return_list_v12 <- function(spec = glc_v13()$rv, n = 1, nam = "urine_1", debug = TRUE) {
    spec$ret <- list(
        number_of_files = n,
        filename = nam,
        x_values = spec$sdp,
        x_values_ppm = spec$ppm,
        y_values = spec$y_smooth,
        spectrum_superposition = spec$lcr$spectrum_approx[1, ],
        mse_normed = spec$lcr$mse_normed,
        index_peak_triplets_middle = spec$peak$center[spec$peak$high],
        index_peak_triplets_left = spec$peak$right[spec$peak$high],
        index_peak_triplets_right = spec$peak$left[spec$peak$high],
        peak_triplets_middle = spec$ppm[spec$peak$center[spec$peak$high]],
        peak_triplets_left = spec$ppm[spec$peak$right[spec$peak$high]],
        peak_triplets_right = spec$ppm[spec$peak$left[spec$peak$high]],
        integrals = spec$lcr$integrals[1, ],
        signal_free_region = c(spec$sfr$left_sdp, spec$sfr$right_sdp),
        range_water_signal_ppm = spec$wsr$hwidth_ppm,
        A = spec$lcr$A_new,
        lambda = spec$lcr$lambda_new,
        x_0 = spec$lcr$w_new
    )
    spec
}
