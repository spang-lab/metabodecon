# Internal #####

#' @name deconvolute
#'
#' @title
#'  Generate Lorentz Curves from NMR Spectra
#'
#' @description
#'  Deconvolutes NMR spectra by modeling each detected signal within a spectrum as Lorentz Curve.
#'  \loadmathjax
#'
#' @inheritParams read_spectrum
#' @param ... Further arguments to be passed to internal subfunctions.
#' @param ask Logical. Whether to ask for user input during the deconvolution process. If FALSE, the provided default values will be used.
#' @param bwc Whether to produce results backwards compatible with [MetaboDecon1D()] and [generate_lorentz_curves()]. If `bwc < 2`, a `decon1` object is returned instead of a `decon2` object. If `bwc < 1`, fixes/improvements introduced after version 'metabodecon v1.0' are not used. Support for `bwc < 1` will be removed in 'metabodecon v2.0' and will result in a warning for 'metabodecon v1.x' with `x >= 1.2`.
#' @param data_path Either the path to a directory containing measured NMR spectra, a dataframe as returned by [read_spectrum()], or a list of such dataframes.
#' @param debug Logical. Whether to return additional intermediate results for debugging purposes.
#' @param delta Threshold for peak filtering. Higher values result in more peaks being filtered out. A peak is filtered if its score is below \mjeqn{\mu + \sigma \cdot \delta}{mu + s * delta}, where \mjeqn{\mu}{mu} is the average peak score in the signal-free region (SFR), and \mjeqn{\sigma}{s} is the standard deviation of peak scores in the SFR.
#' @param force If FALSE, the function stops with an error message if no peaks are found in the signal free region (SFR), as these peaks are required as a reference for peak filtering. If TRUE, the function instead proceeds without peak filtering, potentially increasing runtime and memory usage significantly.
#' @param make_rds Logical or character. If TRUE, stores results as an RDS file on disk. If a character string, saves the RDS file with the specified name. Should be set to TRUE if many spectra are evaluated to decrease computation time.
#' @param nfit Number of iterations for approximating the parameters for the Lorentz curves.
#' @param nworkers Number of workers to use for parallel processing. If `"auto"`, the number of workers will be determined automatically. If a number greater than 1, it will be limited to the number of spectra.
#' @param sf Numeric vector with two entries: the factors to scale the x-axis and y-axis.
#' @param sfr Numeric vector with two entries: the ppm positions for the left and right border of the signal-free region of the spectrum.
#' @param smopts Numeric vector with two entries: the number of smoothing iterations and the number of data points to use for smoothing (must be odd).
#' @param verbose Logical. Whether to print log messages during the deconvolution process.
#' @param wshw Half-width of the water artifact in ppm.
#'  Further arguments to be passed to [generate_lorentz_curves()].
#'
#' @return
#'  A 'GLCDecon' as described in [Metabodecon Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
#' @details
#'  First, an automated curvature based signal selection is performed. Each signal is represented by 3 data points to allow the determination of initial Lorentz curves. These Lorentz curves are then iteratively adjusted to optimally approximate the measured spectrum.
#'
#'  [generate_lorentz_curves_sim()] is identical to [generate_lorentz_curves()] except for the defaults, which are optimized for deconvoluting the 'Sim' dataset, shipped with 'metabodecon'. The 'Sim' dataset is a simulated dataset, which is much smaller than real NMR spectra (1309 datapoints instead of 131072) and lacks a water signal. This makes it ideal for use in examples or as a default value for functions. However, the default values for `sfr`, `wshw`, and `delta` in the "normal" [generate_lorentz_curves()] function are not optimal for this dataset. To avoid having to define the optimal parameters repeatedly in examples, this function is provided to deconvolute the "Sim" dataset with suitable parameters.
#'
#' @examples
#'
#' ## Define the paths to the example datasets we want to deconvolute:
#' ## `sim_dir`: directory containing 16 simulated spectra
#' ## `sim_01`: path to the first spectrum in the `sim` directory
#' ## `sim_01_spec`: the first spectrum in the `sim` directory as a dataframe
#' sim_dir <- metabodecon_file("sim_subset")
#' sim_1_dir <- file.path(sim_dir, "sim_01")
#' sim_2_dir <- file.path(sim_dir, "sim_02")
#' sim_1_spectrum <- read_spectrum(sim_1_dir)
#' sim_2_spectrum <- read_spectrum(sim_2_dir)
#' sim_spectra <- list(sim_1_spectrum, sim_2_spectrum)
#'
#' ## Show that `generate_lorentz_curves()` and `generate_lorentz_curves_sim()`
#' ## produce the same results:
#' sim_1_decon0 <- generate_lorentz_curves(
#'     data_path = path, # Path to directory containing spectra
#'     sfr = c(3.58, 3.42), # Borders of signal free region (SFR) in ppm
#'     wshw = 0, # Half width of water signal (WS) in ppm
#'     delta = 0.1, # Threshold for peak filtering
#'     ask = FALSE, # Don't ask for user input
#'     verbose = FALSE # Suppress status messages
#' )
#' sim_1_decon1 <- generate_lorentz_curves_sim(path)
#' stopifnot(all.equal(sim_1_decon0, sim_1_decon1))
#'
#' ## Show that passing a spectrum produces the same results as passing the
#' ## the corresponding directory:
#' decon_from_spectrum_dir <- generate_lorentz_curves_sim(sim_1_dir)
#' decon_from_spectrum_obj <- generate_lorentz_curves_sim(sim_1_spectrum)
#' decons_from_spectra_obj <- generate_lorentz_curves_sim(sim_spectra)
#' decons_from_spectra_dir <- generate_lorentz_curves_sim(sim_dir)
#' compare <- function(d1, d2) {
#'     ignore <- which(names(d1) %in% c("number_of_files", "filename"))
#'     equal <- all.equal(d1[-ignore], d2[-ignore])
#'     stopifnot(isTRUE(equal))
#' }
#' compare(decon_from_spectrum_dir, decon_from_spectrum_obj)
#' compare(decon_from_spectrum_dir, decons_from_spectra_obj[[1]])
#' compare(decon_from_spectrum_dir, decons_from_spectra_dir[[1]])
#'
#' ## Below example uses data from a real NMR experiment, instead of (small)
#' ## simulated datasets and therefor requires multiple seconds to run. Because
#' ## `ask` is TRUE in this example (the default value), the user will be asked
#' ## for input during the deconvolution. To avoid this, set `ask = FALSE`.
#' \dontrun{
#' example_datasets <- download_example_datasets()
#' urine_1 <- file.path(example_datasets, "bruker/urine/urine_1")
#' decon_urine_1 <- generate_lorentz_curves(urine_1)
#' }
NULL

deconvolute_gspec <- function(gspec, nfit, smopts, delta, sfr, wshw, force, bwc) {
    check_args_deconvolute_gspec()
    logf("Starting deconvolution of %s", gspec$name)
    gspec <- rm_water_signal(gspec, wshw, bwc)
    gspec <- rm_negative_signals(gspec, sfr, bwc)
    gspec <- smooth_signals(gspec, reps = smopts[1], k = smopts[2])
    gspec <- find_peaks(gspec)
    gspec <- filter_peaks(gspec, delta, force)
    gspec <- fit_lorentz_curves(gspec, nfit)
    logf("Finished deconvolution of %s", gspec$name)
    structure(gspec, class = "decon2")
}

deconvolute_gspecs <- function(gspecs,
                               nfit, smopts, delta, sfr, wshw,
                               ask, force, verbose, bwc) {
    check_args_deconvolute_gspecs()
    opts <- if (!verbose) options(toscutil.logf.file = nullfile())
    on.exit(options(opts), add = TRUE)

    adjno <- get_adjno(gspecs, sfr, wshw, ask)
    sfr <- get_sfr(gspecs, sfr, ask, adjno)
    wshw <- get_wshw(gspecs, wshw, ask, adjno)

    ns <- length(gspecs)
    nw <- if (nworkers == "auto") ceiling(parallel::detectCores() / 2) else nworkers
    nw <- min(nw, length(gspecs))
    nfstr <- if (ns == 1) "1 spectrum" else sprintf("%d spectra", ns)
    nwstr <- if (nw == 1) "1 worker" else sprintf("%d workers", nw)

    logf("Starting deconvolution of %s using %s", nfstr, nwstr)
    starttime <- Sys.time()

    decons <- if (nw == 1) {
        lapply(seq_along(gspecs), function(i) {
            deconvolute_gspec(
                gspecs[[i]], nfit, smopts, delta, sfr[[i]],
                wshw[[i]], force, bwc
            )
        })
    } else {
        cl <- parallel::makeCluster(nw, outfile = nullfile())
        on.exit(parallel::stopCluster(cl), add = TRUE)
        export <- c(
            "logf", "fg", "deconvolute_gspec", "gspecs", "nfit",
            "smopts", "delta", "sfr", "wshw", "force", "bwc"
        )
        parallel::clusterExport(cl, export, envir = environment())
        parallel::parLapply(cl, seq_along(gspecs), function(i) {
            opts <- options(toscutil.logf.sep1 = sprintf(" PID %d ", Sys.getpid()))
            on.exit(options(opts), add = TRUE)
            deconvolute_gspec(
                gspecs[[i]], nfit, smopts, delta, sfr[[i]],
                wshw[[i]], force, bwc
            )
        })
    }

    logf("Converting deconvolution results to class '%s'", rtyp)
    convert_func <- get("as_%s", rtyp)
    do.call(convert_func, c(decons, cnvpar))

    duration <- format(round(Sys.time() - starttime, 3))
    logf("Finished deconvolution of %s in %s", nfstr, duration)

    obj
}

# Public API #####

#' @export
#' @rdname deconvolute
generate_lorentz_curves <- function(# Parameters for getting spectra objects
                                    data_path,
                                    file_format = "bruker",
                                    make_rds = FALSE,
                                    expno = 10,
                                    procno = 10,
                                    # Parameters for deconvolution
                                    sfr = c(11.44494, -1.8828),
                                    wshw = 0.1527692,
                                    nfit = 10,
                                    smopts = c(2, 5),
                                    delta = 6.4,
                                    sf = c(1e3, 1e6),
                                    ask = TRUE,
                                    nworkers = 1,
                                    force = FALSE,
                                    verbose = TRUE) {
    args <- check_args()
    args <- rm(data_path, file_format, expno, procno, envir = args)
    args$gspecs <- as_gspecs(data_path, file_format, expno, procno)
    gdecons <- do.call(deconvolute_gspecs, args)
    decons2 <- as_decons2(gdecons)
    store_as_rds(decons2, make_rds, data_path)
    if (length(decons2) == 1) decons2[[1]] else decons2
}

#' @export
#' @rdname deconvolute
generate_lorentz_curves_sim <- function(# Parameters for getting spectra objects
                                        data_path,
                                        file_format = "bruker",
                                        make_rds = FALSE,
                                        expno = 10,
                                        procno = 10,
                                        # Parameters for deconvolution
                                        sfr = c(3.58, 3.42),
                                        wshw = 0,
                                        nfit = 10,
                                        smopts = c(2, 5),
                                        delta = 0.1,
                                        sf = c(1e3, 1e6),
                                        ask = FALSE,
                                        nworkers = 1,
                                        force = FALSE,
                                        verbose = TRUE) {
    generate_lorentz_curves(
        data_path, file_format, make_rds, expno, procno, sfr, wshw,
        nfit, smopts, delta, sf, ask, nworkers, force, verbose
    )
}


# Private Helpers #####

rm_water_signal <- function(x, wshw, bwc) {
    if (!is_gspec(x)) stop("Input must be a gspec object, not ", class(x))
    if (!version %in% 1:2) stop("version be 1 or 2")
    logf("Removing water signal")
    if (bwc) {
        wsr <- enrich_wshw(x, wshw)
        left <- wsr$left_dp
        right <- wsr$right_dp
        idx_wsr <- right:left # Order, i.e. `right:left` instead of `left:right`, is important here, because `right` and `left` are floats. Example: `right <- 3.3; left <- 1.4` ==> `right:left == c(3.3, 2.3)` and `left:right == c(1.4, 2.4)`.
    } else {
        ppm_center <- (x$ppm[1] + x$ppm[length(x$ppm)]) / 2
        idx_wsr <- which(min(wshw) < x$ppm & x$ppm < max(wshw))
    }
    x$y_nows <- x$y_scaled
    x$y_nows[idx_wsr] <- 0.01 / x$sf[2]
    x
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
        # Calling (1) `z <- stats::filter(y, filter, sides = 2)` gives NAs at both sides of vector, as there are not enough values for the moving average. The number of NAs at each side is given by (2) `q <- (k - 1) / 2` . Example: if n==100 and k==5, then q==2, so z[1]==NA, z[2]==NA, z[99]==NA and z[100]==NA. To stay backwards compatible, these values must be filled with the mean of the values that are available. To do so, we iterate from 1:q, i.e. j==1 and j==2 and set
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
filter_peaks <- function(gspec, delta = 6.4, force = FALSE, version = 1) {
    check_args_filter_peaks()
    logf("Removing peaks with low scores")
    peak_score <- gspec$peak$score
    peak_defined <- !is.na(gspec$peak$left) & !is.na(gspec$peak$center) & !is.na(gspec$peak$right)
    if (version == 1) {
        sfr <- enrich_sfr(gspec, gspec$sfr)
    }

    in_left_sfr <- gspec$sdp[gspec$peak$center] >= gspec$sfr$left_sdp
    in_right_sfr <- gspec$sdp[gspec$peak$center] <= gspec$sfr$right_sdp
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
    gspec$peak$high <- peak_defined & (peak_score > mu + delta * sigma)
    gspec$peak$region <- "norm"
    gspec$peak$region[in_left_sfr] <- "sfrl"
    gspec$peak$region[in_right_sfr] <- "sfrr"
    logf("Removed %d peaks", sum(!gspec$peak$high))
    gspec
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

store_as_rds <- function(decons, make_rds, data_path) {
    if (is.character(make_rds)) {
        cat("Saving results as", make_rds, "\n")
        saveRDS(decons, make_rds)
    } else if (isTRUE(make_rds)) {
        rdspath <- file.path(data_path, "spectrum_data.rds")
        if (interactive()) {
            yes <- get_yn_input(sprintf("Save results as '%s'?", rdspath))
            if (yes) saveRDS(decons, rdspath)
        } else {
            logf("Skipping RDS save: confirmation required but not in interactive mode. For details see `help('generate_lorentz_curves')`.")
        }
    }
}
