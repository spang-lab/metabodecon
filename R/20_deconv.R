# Public API #####

#' @name deconvolute
#'
#' @title Deconvolute one or more NMR spectra
#'
#' @description Deconvolutes NMR spectra by modeling each detected signal within
#' a spectrum as Lorentz Curve.
#'
#' @inheritParams read_spectrum
#'
#' @param ask Logical. Whether to ask for user input during the deconvolution
#' process. If FALSE, the provided default values will be used.
#'
#' @param x A `spectrum` or `spectra` object as described in
#' [metabodecon_classes].
#'
#' @param data_path Either the path to a directory containing measured NMR
#' spectra, a dataframe as returned by [read_spectrum()], or a list of such
#' dataframes.
#'
#' @param delta Threshold for peak filtering. Higher values result in more peaks
#' being filtered out. A peak is filtered if its score is below \eqn{\mu +
#' \sigma \cdot \delta}{mu + s * delta}, where \eqn{\mu}{mu} is the average
#' peak score in the signal-free region (SFR), and \eqn{\sigma}{s} is the
#' standard deviation of peak scores in the SFR. See 'Details'.
#'
#' @param force If FALSE, the function stops with an error message if no peaks
#' are found in the signal free region (SFR), as these peaks are required as a
#' reference for peak filtering. If TRUE, the function instead proceeds without
#' peak filtering, potentially increasing runtime and memory usage
#' significantly.
#'
#' @param make_rds Logical or character. If TRUE, stores results as an RDS file
#' on disk. If a character string, saves the RDS file with the specified name.
#' Should be set to TRUE if many spectra are evaluated to decrease computation
#' time.
#'
#' @param nfit Integer. Number of iterations for approximating the parameters
#' for the Lorentz curves. See 'Details'.
#'
#' @param nworkers Number of workers to use for parallel processing. If
#' `"auto"`, the number of workers will be determined automatically. If a number
#' greater than 1, it will be limited to the number of spectra.
#'
#' @param sfr Numeric vector with two entries: the ppm positions for the left
#' and right border of the signal-free region of the spectrum. See 'Details'.
#'
#' @param smopts Numeric vector with two entries: the number of smoothing
#' iterations and the number of data points to use for smoothing (must be odd).
#' See 'Details'.
#'
#' @param verbose Logical. Whether to print log messages during the
#' deconvolution process.
#'
#' @param wshw Half-width of the water artifact in ppm.  See 'Details'.
#'
#' @return A 'decon2' object as described in [metabodecon_classes].
#'
#' @details
#'
#' First, an automated curvature based signal selection is performed. Each
#' signal is represented by 3 data points to allow the determination of initial
#' Lorentz curves. These Lorentz curves are then iteratively adjusted to
#' optimally approximate the measured spectrum.
#'
#' [generate_lorentz_curves_sim()] is identical to [generate_lorentz_curves()]
#' except for the defaults, which are optimized for deconvoluting the 'Sim'
#' dataset, shipped with 'metabodecon'. The 'Sim' dataset is a simulated
#' dataset, which is much smaller than a real NMR spectra and lacks a water
#' signal. This makes it ideal for use in examples or as a default value for
#' functions. However, the default values for `sfr`, `wshw`, and `delta` in the
#' "normal" [generate_lorentz_curves()] function are not optimal for this
#' dataset. To avoid having to define the optimal parameters repeatedly in
#' examples, this function is provided to deconvolute the "Sim" dataset with
#' suitable parameters.
#'
#' In [generate_lorentz_curves()] the parameters `nfit`, `smopts`, `delta`,
#' `sfr` and `wshw` must be fully specified. In [deconvolute()], these
#' parameters can be set to `NULL` (the default value). In this case, the
#' function will try to determine the optimal values for these parameters
#' automatically. The values chosen are stored in field `args` of the returned
#' `decon2` object.
#'
#' @examples
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' ## Define the paths to the example datasets we want to deconvolute:
#' ## `sim_dir`: directory containing 16 simulated spectra
#' ## `sim_01`: path to the first spectrum in the `sim` directory
#' ## `sim_01_spec`: the first spectrum in the `sim` directory as a dataframe
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' sim_dir <- metabodecon_file("sim_subset")
#' sim_1_dir <- file.path(sim_dir, "sim_01")
#' sim_2_dir <- file.path(sim_dir, "sim_02")
#' sim_1_spectrum <- read_spectrum(sim_1_dir)
#' sim_2_spectrum <- read_spectrum(sim_2_dir)
#' sim_spectra <- read_spectra(sim_dir)
#'
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' ## Show that `generate_lorentz_curves()` and `generate_lorentz_curves_sim()`
#' ## produce the same results:
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' sim_1_decon0 <- generate_lorentz_curves(
#'     data_path = sim_1_dir, # Path to directory containing spectra
#'     sfr = c(3.55, 3.35), # Borders of signal free region (SFR) in ppm
#'     wshw = 0, # Half width of water signal (WS) in ppm
#'     ask = FALSE, # Don't ask for user input
#'     verbose = FALSE # Suppress status messages
#' )
#' sim_1_decon1 <- generate_lorentz_curves_sim(sim_1_dir)
#' stopifnot(all.equal(sim_1_decon0, sim_1_decon1))
#'
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' ## Show that passing a spectrum produces the same results as passing the
#' ## the corresponding directory:
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' decon_from_spectrum_dir <- generate_lorentz_curves_sim(sim_1_dir)
#' decon_from_spectrum_obj <- generate_lorentz_curves_sim(sim_1_spectrum)
#' decons_from_spectra_obj <- generate_lorentz_curves_sim(sim_spectra)
#' decons_from_spectra_dir <- generate_lorentz_curves_sim(sim_dir)
#'
#' most.equal <- function(x1, x2) {
#'     ignore <- which(names(x1) %in% c("number_of_files", "filename"))
#'     equal <- all.equal(x1[-ignore], x2[-ignore])
#'     invisible(stopifnot(isTRUE(equal)))
#' }
#'
#' all.equal(  decon_from_spectrum_dir, decon_from_spectrum_obj     )
#' all.equal(  decons_from_spectra_dir, decons_from_spectra_obj     )
#' most.equal( decon_from_spectrum_dir, decons_from_spectra_obj[[1]])
#' most.equal( decon_from_spectrum_dir, decons_from_spectra_dir[[1]])
#'
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' ## Below example uses data from a real NMR experiment, instead of (small)
#' ## simulated datasets and therefor requires multiple seconds to run. Because
#' ## `ask` is TRUE in this example (the default value), the user will be asked
#' ## for input during the deconvolution. To avoid this, set `ask = FALSE`.
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' \dontrun{
#' example_datasets <- download_example_datasets()
#' urine_1 <- file.path(example_datasets, "bruker/urine/urine_1")
#' decon_urine_1 <- generate_lorentz_curves(urine_1)
#' }
NULL

#' @export
#' @rdname deconvolute
deconvolute <- function(x,
                        nfit = 3,
                        smopts = c(2, 5),
                        delta = 6.4,
                        sfr = NULL,
                        wshw = 0,
                        ask = FALSE,
                        force = FALSE,
                        verbose = FALSE,
                        nworkers = 1) {
    check_args_deconvolute()
    ispecs <- as_ispecs(x)
    if (is.null(sfr)) sfr <- quantile(ispecs[[1]]$ppm, c(0.9, 0.1))
    idecons <- deconvolute_ispecs(
        ispecs,
        nfit, smopts, delta, sfr, wshw,
        ask, force, verbose,
        bwc = 2, nworkers = nworkers
    )
    decons2 <- as_decons2(idecons)
    if (length(decons2) == 1) decons2[[1]] else decons2
}

#' @export
#' @rdname deconvolute
generate_lorentz_curves <- function(data_path,
                                    file_format = "bruker",
                                    make_rds = FALSE,
                                    expno = 10,
                                    procno = 10,
                                    raw = TRUE,
                                    nfit = 10,
                                    smopts = c(2, 5),
                                    delta = 6.4,
                                    sfr = c(11.44494, -1.8828),
                                    wshw = 0.1527692,
                                    ask = TRUE,
                                    force = FALSE,
                                    verbose = TRUE,
                                    nworkers = 1) {
    # Check and convert input
    eval(generate_lorentz_curves_type_checks) # 99 us
    spectra <- as_spectra(
        data_path, file_format, expno, procno, raw,
        silent = !verbose, force = force
    )

    # Deconvolute
    idecons <- deconvolute_ispecs(
        ispecs = as_ispecs(x = spectra),
        nfit, smopts, delta, sfr, wshw,
        ask, force, verbose,
        bwc = 1, nworkers = nworkers
    )

    # Convert, store and return
    if (verbose) logf("Formatting deconvolution results as 'decon1' objects")
    decons1 <- as_decons1(idecons)
    store_as_rds(decons1, make_rds, data_path)
    if (length(decons1) == 1) decons1[[1]] else decons1
}

#' @export
#' @rdname deconvolute
generate_lorentz_curves_sim <- function(data_path,
                                        file_format = "bruker",
                                        make_rds = FALSE,
                                        expno = 10,
                                        procno = 10,
                                        raw = TRUE,
                                        nfit = 10,
                                        smopts = c(2, 5),
                                        delta = 6.4,
                                        sfr = c(3.55, 3.35),
                                        wshw = 0,
                                        ask = FALSE,
                                        force = FALSE,
                                        verbose = FALSE,
                                        nworkers = 1) {
    generate_lorentz_curves(
        data_path, file_format, make_rds, expno, procno, raw,
        nfit, smopts, delta, sfr, wshw,
        ask, force, verbose, nworkers
    )
}

# Internal main functions #####

#' @noRd
#' @inheritParams deconvolute
#' @param bwc Whether to produce results backwards compatible with
#' [MetaboDecon1D()]. If `bwc == 0`, bug fixes introduced after version 0.2.2 of
#' Metabodecon are not used. If `bwc == 1`, new features introduced after
#' version 0.2.2 of Metabodecon (e.g. faster algorithms) are not used. If `bwc
#' == 2`, all bug fixes and features introduced after version 0.2.2 are used.
#'
#' Support for `bwc == 0` will be removed in 'metabodecon v2.0'.
deconvolute_ispecs <- function(ispecs,
                               nfit = 3,
                               smopts = c(2, 5),
                               delta = 6.4,
                               sfr = c(3.55, 3.35),
                               wshw = 0,
                               ask = FALSE,
                               force = FALSE,
                               verbose = TRUE,
                               bwc = 2,
                               nworkers = 1) {

    # Check args & configure logging
    args <- check_args_deconvolute_ispecs()
    ispecs <- as_ispecs(ispecs)
    opts <- if (!verbose) options(toscutil.logf.file = nullfile())
    on.exit(options(opts), add = TRUE)

    # Init local vars
    ns <- length(ispecs)
    nw <- if (nworkers == "auto") ceiling(detectCores() / 2) else nworkers
    nw <- min(nw, length(ispecs))
    nfstr <- if (ns == 1) "1 spectrum" else sprintf("%d spectra", ns)
    nwstr <- if (nw == 1) "1 worker" else sprintf("%d workers", nw)
    adjno <- get_adjno(ispecs, sfr, wshw, ask)
    sfrlist <- get_sfr(ispecs, sfr, ask, adjno)
    wshwlist <- get_wshw(ispecs, wshw, ask, adjno)
    smoptslist <- get_smopts(ispecs, smopts)

    # Deconvolute spectra
    logf("Starting deconvolution of %s using %s", nfstr, nwstr)
    starttime <- Sys.time()
    idecon_list <- mcmapply(nw, deconvolute_ispec, ispecs,
        nfit, smoptslist, delta, sfrlist, wshwlist,
        ask, force, verbose, bwc, nworkers)
    idecons <- as_idecons(idecon_list)
    duration <- format(round(Sys.time() - starttime, 3))
    logf("Finished deconvolution of %s in %s", nfstr, duration)
    idecons
}

deconvolute_ispec <- function(ispec,
                              nfit = 3,
                              smopts = c(2, 5),
                              delta = 6.4,
                              sfr = c(3.55, 3.35),
                              wshw = 0,
                              ask = FALSE,
                              force = FALSE,
                              verbose = TRUE,
                              bwc = 2,
                              nworkers = 1) {

    # Check args & configure logging
    ispec <- as_ispec(ispec)
    args <- check_args_deconvolute_ispec()
    ispec$args <- args[names(args) != "ispec"]
    reps <- smopts[1]
    k <- smopts[2]
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())

    # Deconvolute ispec
    logf("Starting deconvolution of %s", ispec$name)
    ispec <- rm_water_signal(ispec, wshw, bwc)
    ispec <- rm_negative_signals(ispec)
    ispec <- smooth_signals(ispec, reps, k, bwc)
    ispec <- find_peaks(ispec)
    ispec <- filter_peaks(ispec, sfr, delta, force, bwc)
    ispec <- fit_lorentz_curves(ispec, nfit, bwc)
    idecon <- as_idecon(ispec)
    logf("Finished deconvolution of %s", ispec$name)
    idecon
}

# Helpers for deconvolute_ispecs #####

#' @description Return number of spectrum to adjust all others or 0 if every spectrum should be adjusted individually.
#' @noRd
get_adjno <- function(ispecs, sfr, wshw, ask) {
    if (isFALSE(ask) || length(ispecs) == 1) {
        return(0)
    }
    same_param <- get_yn_input("Use same parameters for all spectra?")
    if (!same_param) {
        return(0)
    }
    namestr <- paste(seq_along(ispecs), names(ispecs), sep = ": ", collapse = ", ")
    prompt <- sprintf("Number of spectrum for adjusting parameters? (%s)", namestr)
    get_num_input(prompt, min = 1, max = length(ispecs), int = TRUE)
}

#' @description Convert to list of correct length (one per spectrum) and let user confirm each entry.
#' @noRd
get_sfr <- function(ispecs, sfr, ask, adjno) {
    n <- length(ispecs)
    if (is_num(sfr, 2)) sfr <- rep(list(sfr), n)
    if (!is_list_of_nums(sfr, n, 2)) stop("sfr should be a [list of] num(2)")
    if (ask && adjno == 0) { # Different SFR for each spectrum.
        sfr <- lapply(seq_len(n), function(i) confirm_sfr(ispecs[[i]], sfr[[i]]))
    }
    if (ask && adjno >= 1) { # Same SFR for each spectrum.
        sfr_adjno <- confirm_sfr(ispecs[[adjno]], sfr[[adjno]])
        sfr <- rep(list(sfr_adjno), n)
    }
    names(sfr) <- names(ispecs)
    sfr
}

#' @description Convert to list of correct length (one per spectrum) and let user confirm each entry.
#' @noRd
get_wshw <- function(ispecs, wshw, ask, adjno) {
    n <- length(ispecs)
    if (is_num(wshw, 1)) wshw <- rep(list(wshw), n)
    if (is_num(wshw, n)) wshw <- as.list(wshw)
    if (!is_list_of_nums(wshw, n, 1)) stop("wshw should be a [list of] num(1)")
    if (ask && adjno == 0) wshw <- mapply(confirm_wshw, ispecs, wshw)
    if (ask && adjno >= 1) wshw <- rep(list(confirm_wshw(ispecs[[adjno]], wshw[[adjno]])), n)
    names(wshw) <- names(ispecs)
    wshw
}

#' @description Convert to list of correct length (one per spectrum).
#' @noRd
get_smopts <- function(ispecs, smopts) {
    n <- length(ispecs)
    if (is_int(smopts, 2)) smopts <- rep(list(smopts), n)
    if (!is_list_of_nums(smopts, n, 2)) stop("smopts should be a [list of] int(2)")
    names(smopts) <- names(ispecs)
    smopts
}

# Helpers for deconvolute_ispec #####

rm_water_signal <- function(x, wshw, bwc) {
    check_args_rm_water_signal()
    logf("Removing water signal")
    if (bwc >= 1) {
        ppm_center <- (x$ppm[1] + x$ppm[length(x$ppm)]) / 2
        idx_wsr <- which(x$ppm > ppm_center - wshw & x$ppm < ppm_center + wshw)
        x$y_nows <- x$y_scaled
        x$y_nows[idx_wsr] <- min(x$y_nows)
    } else {
        wsr <- enrich_wshw(wshw, x)
        left <- wsr$left_dp
        right <- wsr$right_dp
        idx_wsr <- right:left # (1)
        x$y_nows <- x$y_scaled
        x$y_nows[idx_wsr] <- 0.01 / x$sf[2]
        # (1) Order is important here, because right and left are floats. Example:
        # right <- 3.3; left <- 1.4
        # right:left == c(3.3, 2.3)
        # left:right == c(1.4, 2.4)
    }
    x
}

rm_negative_signals <- function(spec) {
    logf("Removing negative signals")
    errmsg <- "Water signal not removed yet. Call `rm_water_signal()` first."
    if (is.null(spec$y_nows)) stop(errmsg)
    spec$y_pos <- abs(spec$y_nows)
    spec
}

#' @noRd
#'
#' @title Smooth Signal Intensities using a Moving Average
#'
#' @description
#' Smoothes signal intensities by applying a [moving average](
#' https://en.wikipedia.org/wiki/Moving_average) filter with a window size of k.
#'
#' @param spec A list representing the spectrum, which should include the scaled
#' signal intensities, after removal of the water artefact and negative values
#' (`spec$y_pos`).
#'
#' @param reps The number of times to apply the moving average.
#'
#' @param k The number of points within the moving average window. Must be odd,
#' so the smoothed point is in the middle of the window.
#'
#' @return
#' A numeric vector of the smoothed values.
#'
#' @details
#' Old and slow version producing the same results as the
#' implementation within `deconvolution` from `MetaboDecon1D_deconvolution.R`.
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
#' @inherit smooth_signals
#' @details
#' New and fast version for smoothing of signals. Implements the same algorithm
#' as [smooth_signal_v12()] using different R functions (e.g.
#' [stats::filter()]), causing a massive speedup but also numeric differences
#' compared to the old version.
#'
#' WORK IN PROGRESS
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
        # Calling (1) gives NAs at both sides of vector, as there are not enough
        # values for the moving average. The number of NAs at each side is given
        # by (2). Example: if n==100 and k==5, then q==2, so z[1]==NA, z[2]==NA,
        # z[99]==NA and z[100]==NA. To stay backwards compatible, these values
        # must be filled with the mean of the values that are available. To do
        # so, we iterate from 1:q, i.e. j==1 and j==2 and set
        #
        # >>> z[1]   <- mean(y[1:3])    # 3 == 2+1 == q+j            # (3)
        # >>> z[2]   <- mean(y[1:4])    # 4 == 2+2 == q+j            # (3)
        # >>> z[99]  <- mean(y[97:100]) # 97 == 100-2-2+1 == n-q-j+1 # (4)
        # >>> z[100] <- mean(y[98:100]) # 98 == 100-2-1+1 == n-q-j+1 # (4)
        #
        # Note: we could also think of leaving the NAs as they are, which would
        # be more correct I think and even faster, but would break compatibility
        # with the old version completely. So not even `all.equal(v1, v2)` would
        # be TRUE anymore.
    }
    spec$Z <- Z
    spec$y_smooth <- Z[[reps]]
    spec
}

find_peaks <- function(spec) {
    logf("Starting peak selection")
    d <- spec$d <- calc_second_derivative(y = spec$y_smooth)
    a <- abs(d)
    m <- length(d)
    dl <- c(NA, d[-m]) # dl[i] == d[i-1]
    dr <- c(d[-1], NA) # dr[i] == d[i+1]
    center <- which(d < 0 & d <= dl & d < dr)
    spec$peak <- data.frame(left = NA, center = center, right = NA, score = NA)
    for (i in seq_along(center)) {
        j <- center[i]
        l <- spec$peak$left[i] <- get_left_border(j, d)
        r <- spec$peak$right[i] <- get_right_border(j, d, m)
        spec$peak$score[i] <- get_peak_score(j, l, r, a)
    }
    logf("Detected %d peaks", length(center))
    return(spec)
}

# WORK IN PROGRESS
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
#' @title Filter Peaks with Low Scores Outside Signal-Free Region
#'
#' @description
#' Calculates a score for each peak in the spectrum. Peaks with scores that are
#' lower than `mu + delta * sigma` are filtered out, with `mu` being the mean
#' and `sigma` the standard deviation of the scores of peaks within the
#' signal-free region (SFR).
#'
#' @param spec A list containing the spectrum data, including peaks and their
#' scores, and the signal-free region definition.
#'
#' @param delta A numeric value specifying how many standard deviations `s` a
#' score needs to be above `mu` to not get filtered out. Here `s` denotes the
#' standard deviation of peak scores in the signal free region (SFR) and `mu`
#' denotes the average peak score within the SFR.
#'
#' @param force If no peaks are found in the SFR, the function stops with an
#' error message by default. If `force` is TRUE, the function instead proceeds
#' without filtering any peaks, potentially increasing runtime.
#'
#' @return
#' Returns the modified `spec` list with the `peak` component updated to
#' indicate which peaks are considered significant based on their score relative
#' to the SFR and the `delta` parameter. Peaks within the SFR are marked with a
#' specific region code ("sfrl" for peaks in the left SFR and "sfrr" for peaks
#' in the right SFR). Peaks in the normal region have the region-code "norm".
#'
#' @details
#' The function first identifies peaks within the SFR by comparing their center
#' positions against the SFR boundaries. If peaks are found within the SFR, it
#' calculates the mean and standard deviation of their scores to establish a
#' filtering threshold. Peaks with scores below this threshold are considered
#' low and filtered out. If no peaks are found within the SFR and `force` is
#' FALSE, the function stops and issues an error message. If `force` is TRUE,
#' the function proceeds without filtering, potentially increasing runtime.
#'
#' @examples
#'
#' peak <- list(score = c(1, 5, 2), center = c(1, 2, 3))
#' sdp <- c(3, 2, 1)
#' ispec <- named(peak, sdp)
#' sfr <- list(left_sdp = 2.8, right_sdp = 1.2)
#' rm3 <- filtered_ispec <- filter_peaks(ispec, sfr)
#' rm2 <- filtered_ispec <- filter_peaks(ispec, sfr, delta = 1)
filter_peaks <- function(ispec, sfr, delta = 6.4, force = FALSE, bwc = 1) {
    stopifnot(is_ispec(ispec))
    logf("Removing peaks with low pscores")
    sdp <- ispec$sdp
    ppm <- ispec$ppm
    plb <- ispec$peak$left
    prb <- ispec$peak$right
    pct <- ispec$peak$center
    psc <- ispec$peak$score
    pok <- !is.na(plb) & !is.na(pct) & !is.na(prb)
    if (bwc < 1) {
        sfr <- enrich_sfr(sfr, ispec)
        in_left_sfr <- sdp[pct] >= sfr$left_sdp
        in_right_sfr <- sdp[pct] <= sfr$right_sdp
    } else {
        in_left_sfr <- ppm[pct] >= max(sfr)
        in_right_sfr <- ppm[pct] <= min(sfr)
    }
    in_sfr <- in_left_sfr | in_right_sfr
    if (sum(in_sfr) > 1) {
        mu <- mean(psc[in_sfr])
        sigma <- sd(psc[in_sfr])
    } else {
        if (!force) stop("Not enough signals found in signal free region. Please double check deconvolution parameters.")
        logf("No enough signals found in signal free region. This is a clear indication that the deconvolution parameters are not set correctly. Continuing anyways without dynamic peak filtering, because `force` is TRUE. Note that this might increase runtime drastically.")
        mu <- 0
        sigma <- 0
    }
    ispec$peak$high <- pok & (psc > mu + delta * sigma)
    ispec$peak$region <- "norm"
    ispec$peak$region[in_left_sfr] <- "sfrl"
    ispec$peak$region[in_right_sfr] <- "sfrr"
    logf("Removed %d peaks", sum(!ispec$peak$high))
    ispec
}

fit_lorentz_curves <- function(spec, nfit = 3, bwc = 1) {
    logf("Initializing Lorentz curves")
    spec$lci <- lc <- init_lc(spec) # Lorentz Curves Initialized
    spec$lca <- vector("list", length = nfit) # Lorentz Curves Approximated
    logf("Refining Lorentz Curves")
    for (i in 1:nfit) spec$lca[[i]] <- lc <- refine_lc_v14(spec, lc$Z)
    A <- lc$A
    lambda <- lc$lambda
    w <- lc$w
    limits <- range(spec$sdp)
    if (bwc < 2) limits[2] <- limits[2] + (1 / spec$sf[1])
    integrals <- lorentz_int(w, A, lambda, limits = limits)
    spec$lcr <- list(A = A, lambda = lambda, w = w, integrals = integrals)
    spec
}

# Helpers for get_sfr and get_wshw #####

#' @description Repeatedly ask the user to confirm/refine SFR borders.
#' @noRd
confirm_sfr <- function(x, sfr = c(11.44494, -1.8828)) {
    plot_sfr(x, sfr[1], sfr[2])
    sfr_ok <- get_yn_input("Signal free region correctly selected?")
    while (!sfr_ok) {
        get_border <- function(msg) get_num_input(msg, x$ppm_min, x$ppm_max)
        sfr[1] <- get_border("Choose another left border: [e.g. 12]")
        sfr[2] <- get_border("Choose another right border: [e.g. -2]")
        plot_sfr(x, sfr[1], sfr[2])
        sfr_ok <- get_yn_input("Signal free region correctly selected?")
    }
    sfr
}

#' @description Repeatedly ask the user to confirm/refine the WSHW.
#' @noRd
confirm_wshw <- function(x, wshw) {
    plot_ws(x, wshw)
    ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    while (!ws_ok) {
        wshw <- get_num_input("Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154]")
        plot_ws(x, wshw)
        ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    }
    wshw
}

# Helpers for filter_peaks #####

#' @noRd
#' @description
#' Takes the SFR in PPM and returns the SFR in PPM, DP and SDP.
#' @note
#' Because the conversion from PPM to DP/SDP is slightly off (by 1-2 data
#' points), the SFR borders in DP/SDP returned by this function are also
#' incorrect. However, to maintain backwards compatibility with the old
#' MetaboDecon1D function, the behaviour is not changed in this function.
#' Instead, to only work with the correct ppm values, set `bwc = 2` in
#' [filter_peaks()]. For details see `CHECK-2: signal free region calculation`
#' in `TODOS.md`.
enrich_sfr <- function(sfr, x) {
    stopifnot(is_ispec(x) || is_idecon(x))
    left_ppm <- sfr[1]
    right_ppm <- sfr[2]
    left_dp <- (x$n + 1) - (x$ppm_max - left_ppm) / x$ppm_nstep
    left_sdp <- left_dp / x$sf[1]
    right_dp <- (x$n + 1) - (x$ppm_max - right_ppm) / x$ppm_nstep
    right_sdp <- right_dp / x$sf[1]
    named(left_ppm, right_ppm, left_dp, right_dp, left_sdp, right_sdp)
}

#' @noRd
#' @description
#' Calculates the WSR in dp and ppm from the WSHW in ppm.
#' @note
#' Because the conversion from PPM to DP/SDP is slightly off (by 1-2 data
#' points), the SFR borders in DP/SDP returned by this function are also
#' incorrect. However, to maintain backwards compatibility with the old
#' MetaboDecon1D function, the behaviour is not changed in this function.
#' Instead, to only work with the correct ppm values, set `bwc = 2` in
#' [rm_water_signal()]. For details see `CHECK-3: water signal calculation` in
#' `TODOS.md`.
enrich_wshw <- function(wshw, x) {
    stopifnot(is_ispec(x) || is_idecon(x))
    x <- as_ispec(x)
    hwidth_ppm <- wshw
    hwidth_dp <- hwidth_ppm / x$ppm_nstep
    center_dp <- x$n / 2
    right_dp <- center_dp + hwidth_dp
    left_dp <- center_dp - hwidth_dp
    center_ppm <- x$ppm[center_dp]
    right_ppm <- x$ppm[right_dp]
    left_ppm <- x$ppm[left_dp]
    if (left_dp <= 1 || right_dp >= x$n) stop("WSR is out of range")
    named(
        left_ppm, right_ppm, center_ppm, hwidth_ppm,
        left_dp, right_dp, center_dp, hwidth_dp
    )
}

# Helpers for find_peak #####

calc_second_derivative <- function(y) {
    n <- length(y)
    x <- c(NA, y[-n]) # x[i] == y[i-1]
    z <- c(y[-1], NA) # z[i] == y[i+1]
    d <- x + z - 2 * y
    d
}

get_right_border <- function(j, d, m) {
    r <- j + 1
    while (r < m) { # use r<m instead of r<=m because c4 requires d[r+1]
        c1 <- d[r] > d[r - 1]
        c2 <- d[r] >= d[r + 1]
        c3 <- d[r] < 0
        c4 <- d[r + 1] >= 0
        is_right_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_right_border)) {
            return(r)
        }
        r <- r + 1
    }
    return(NA)
}

get_left_border <- function(j, d) {
    l <- j - 1
    while (l > 1) { # use l>1 instead of l>=1 because c4 requires d[l-1]
        c1 <- d[l] > d[l + 1]
        c2 <- d[l] >= d[l - 1]
        c3 <- d[l] < 0
        c4 <- d[l - 1] >= 0
        is_left_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_left_border)) {
            return(l)
        }
        l <- l - 1
    }
    return(NA)
}

#' @noRd
#'
#' @title Get Peak Score
#'
#' @description
#' Calculate the score of a peak based on the sum of absolute second derivative
#' values of its datapoints.
#'
#' @param j <- Index of the peak center
#' @param l <- Index of the left border
#' @param r <- Index of the right border
#' @param a <- Absolute values of the second derivative for all data points
#'
#' @return The score of the peak.
#'
#' @examples
#'
#' #      ____________________________________________________
#' #     |____2________5___________9_______12____14____16_____|
#' #     |             x                                      |
#' #     |          x  x  x  x                    x           |
#' #     |       x  x  x  x  x  x              x  x           |
#' #     |_.__x__x__x__x__x__x__x__x__.__.__x__x__x__x__x__.__|
#' #          |----2---|-----4-----|        |--3--|--6--|
#'
#' y <- c( 0, 1, 2, 3, 4, 3, 3, 2, 1, 0, 0, 1, 2, 3, 1, 1, 0  )
#' a <- c(NA, 0, 0, 0, 2, 1, 1, 0, 0, 1, 1, 0, 0, 3, 2, 1, NA )
#' all.equal(a, abs(calc_second_derivative(y)))
#'
#' s1 <- get_peak_score( 5, 2,   9, a)
#' s2 <- get_peak_score(14, 12, 16, a)
#' stopifnot(s1 == min(sum(a[2:5]), sum(a[5:9])))
#' stopifnot(s2 == min(sum(a[12:14]), sum(a[14:16])))
get_peak_score <- function(j, l, r, a) {
    if (any(is.na(a[c(l, j, r)]))) {
        0
    } else {
        min(sum(a[l:j]), sum(a[j:r]))
    }
}

# Helpers for fit_lorentz_curves #####

#' @title Initialize Lorentz Curve Parameters
#' @param spec List with elements: `x`, `y`, `peak` where `peak` is a list with
#' elements `center`, `left`, `right` and `high`.
#' @noRd
init_lc <- function(spec, verbose = TRUE) {

    # Init values
    p <- spec$peak
    ir <- p$right[p$high]  #
    ic <- p$center[p$high] # Index of each peak triplet position (PTP)
    il <- p$left[p$high]   #
    lmr <- sort(unique(c(il, ic, ir))) # Combined PTP indices
    rr <- match(ir, lmr) #
    rc <- match(ic, lmr) # Rank of each PTP
    rl <- match(il, lmr) #
    x <- spec$sdp; y <- spec$y_smooth # X and Y value for each data point
    xlmr <- x[lmr]; # X value for each PTP
    yr <- y[ir]; yc <- y[ic]; yl <- y[il]; # Intensity of each PTP
    xr <- x[ir]; xc <- x[ic]; xl <- x[il]; # Position of each PTP

    # Replace shoulders
    as <- (yr > yc) & (yc > yl) # Ascending shoulders (AS)
    ds <- (yr < yc) & (yc < yl) # Descending shoulders (DS)
    xl[ds] <- 2 * xc[ds] - xr[ds] # Replace DS with mirrored points (MP)
    xr[as] <- 2 * xc[as] - xl[as] # Replace AS with MP
    yl[ds] <- yr[ds] # Replace DS with MP
    yr[as] <- yl[as] # Replace AS with MP

    # Calculate distances
    wr  <- xr - xr #
    wc  <- xc - xr # Express positions wr/wc/wl as "distance to right border"
    wl  <- xl - xr #
    wrc <- wr - wc; wrl <- wr - wl; wcl <- wc - wl # x - distance between PTPs
    yrc <- yr - yc; yrl <- yr - yl; ycl <- yc - yl # y - distance between PTPs

    # Estimate parameters
    w <- calc_w(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr)
    lambda <- calc_lambda(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr)
    A <- calc_A(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, lambda, w, xr)

    # Calculate contribution of each lorentz curve to each PTP data point
    Z <- matrix(0, nrow = length(lmr), ncol = length(ic)) # 3614 x 1227 urine_1
    idx_A_non_zero <- which(A != 0)
    for (j in idx_A_non_zero) {
        Z[, j] <- abs(A[j] * (lambda[j] / (lambda[j]^2 + (xlmr - w[j])^2)))
    }

    # Print MSE
    mse <- mean((y[lmr] - rowSums(Z))^2)
    if (verbose) logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, as, ds)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

refine_lc_v14 <- function(spec, Z) {

    # Init x and y values
    x <- spec$sdp; y <- spec$y_smooth # x and y value for each data point

    # Init peak related variables
    p <- spec$peak
    ir <- p$right[p$high]; ic <- p$center[p$high]; il <- p$left[p$high] # index of each peak triplet position (PTP)
    lmr <- sort(unique(c(il, ic, ir))) # combined PTP indices
    rr <- match(ir, lmr);  rc <- match(ic, lmr);   rl <- match(il, lmr) # rank  of each PTP
    xlmr <- x[lmr]; # x value for each PTP
    yr <- y[ir]; yc <- y[ic]; yl <- y[il]; # intensity of each PTP
    xr <- x[ir]; xc <- x[ic]; xl <- x[il]; # position of each PTP
    sl <- sc <- sr <- numeric(length(ic)); # sum of lorentz curves (SLC) at each PTP
    ql <- qc <- qr <- numeric(length(ic)); # ratio (SLC / original spectrum) at each PTP

    # Init distance related variables
    wr  <- wc  <- wl  <- numeric(length(ic));
    wrc <- wrl <- wcl <- numeric(length(ic));
    yrc <- yrl <- ycl <- numeric(length(ic));

    # Init lorentz curves parameters and matrices
    A <- lambda <- w <- numeric(length(ic))

    for (i in seq_along(il)) {

        # Calculate the sum of all lorentz curves for each data point
        sl[i] <- sum(Z[rl[i], ])
        sc[i] <- sum(Z[rc[i], ])
        sr[i] <- sum(Z[rr[i], ])

        # Calculate the proportion between original spectrum an the sum of the
        # lorentz curves for each peak triplets position
        ql[i] <- yl[i] / sl[i]
        qc[i] <- yc[i] / sc[i]
        qr[i] <- yr[i] / sr[i]

        # Calculate the new heights of the peak triplets
        yl[i] <- Z[rl[i], i] * ql[i]
        yc[i] <- Z[rc[i], i] * qc[i]
        yr[i] <- Z[rr[i], i] * qr[i]

        # Calculate mirrored points for ascending and descending shoulders
        if ((yl[i] < yc[i]) && (yc[i] < yr[i])) { # Ascending shoulder
            xr[i] <- 2 * xc[i] - xl[i]
            yr[i] <- yl[i]
        }
        if ((yl[i] > yc[i]) && (yc[i] > yr[i])) { # Descending shoulder
            xl[i] <- 2 * xc[i] - xr[i]
            yl[i] <- yr[i]
        }

        # Calculate distances between peak triplet positions and intensities
        wr[i] <- xr[i] - xr[i]
        wc[i] <- xc[i] - xr[i]
        wl[i] <- xl[i] - xr[i]
        wrc[i] <- wr[i] - wc[i]
        wrl[i] <- wr[i] - wl[i]
        wcl[i] <- wc[i] - wl[i]
        yrc[i] <- yr[i] - yc[i]
        yrl[i] <- yr[i] - yl[i]
        ycl[i] <- yc[i] - yl[i]

        # Estimate parameters
        w[i] <- calc_w(
            wr[i], wc[i], wl[i],
            yr[i], yc[i], yl[i],
            wrc[i], wrl[i], wcl[i],
            yrc[i], yrl[i], ycl[i],
            xr[i]
        )
        lambda[i] <- calc_lambda(
            wr[i], wc[i], wl[i],
            yr[i], yc[i], yl[i],
            wrc[i], wrl[i], wcl[i],
            yrc[i], yrl[i], ycl[i],
            xr[i]
        )
        A[i] <- calc_A(
            wr[i], wc[i], wl[i],
            yr[i], yc[i], yl[i],
            wrc[i], wrl[i], wcl[i],
            yrc[i], yrl[i], ycl[i],
            lambda[i], w[i], xr[i]
        )

        # Calculate contribution of each lorentz curve to each data point
        cond <- (w[i] == 0) || (lambda[i] == 0) || (A[i] == 0)
        Z[, i] <- if (cond) 0 else abs(A[i] * (lambda[i] / (lambda[i]^2 + (xlmr - w[i])^2)))

        if (w[i] == 0 || lambda[i] == 0 || A[i] == 0) {
            Z[, i] <- 0
        } else {
            Z[, i] <- abs(A[i] * (lambda[i] / (lambda[i]^2 + (xlmr - w[i])^2)))
        }
    }

    # Print MSE
    mse <- mean((y[lmr] - rowSums(Z))^2)
    logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, sl, sc, sr, ql, qc, qr)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

# Taken from Appendix E of Koh et. al. 2009
calc_w <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr) {
    t1 <- wr^2 * yr * ycl
    t2 <- wl^2 * yl * yrc
    t3 <- wc^2 * yc * yrl
    t4 <- 2 * wrc * yr * yc
    t5 <- 2 * wcl * yc * yl
    t6 <- 2 * wrl * yr * yl
    w <- (t1 + t2 - t3) / (t4 + t5 - t6) + xr
    w[is.nan(w)] <- 0 # If (t4 + t5 - t6) is 0, then w is NaN. In this case we set w to 0.
    w
}

# Taken from Appendix E of Koh et. al. 2009
calc_lambda <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr) {
    lambda <- -((sqrt(abs((-wc^4 * yc^2 * yrl^2 - wr^4 * yr^2 * ycl^2 - wl^4 * yrc^2 * yl^2 + 4 * wc * wl^3 * yc * ((-yr) + yc) * yl^2 + 4 * wc^3 * wl * yc^2 * yl * ((-yr) + yl) + 4 * wr^3 * yr^2 * ycl * (wc * yc - wl * yl) + 4 * wr * yr * (wc^3 * yc^2 * yrl - wc * wl^2 * yc * (yr + yc - 2 * yl) * yl + wl^3 * yrc * yl^2 - wc^2 * wl * yc * yl * (yr - 2 * yc + yl)) + 2 * wc^2 * wl^2 * yc * yl * (yr^2 - 3 * yc * yl + yr * (yc + yl)) + 2 * wr^2 * yr * (-2 * wc * wl * yc * yl * (-2 * yr + yc + yl) + wl^2 * yl * (yr * (yc - 3 * yl) + yc * (yc + yl)) + wc^2 * yc * (yr * (-3 * yc + yl) + yl * (yc + yl)))))))) / (2 * sqrt((wr * yr * ycl + wl * yrc * yl + wc * yc * ((-yr) + yl))^2))
    lambda[is.nan(lambda)] <- 0
    lambda
}

# Taken from Appendix E of Koh et. al. 2009
calc_A <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, lambda, w, xr) {
    A <- (-4 * wrc * wrl * wcl * yr * yc * yl * (wr * yr * ycl + wl * yl * yrc + wc * yc * (-yrl)) * lambda) / (wrc^4 * yr^2 * yc^2 - 2 * wrc^2 * yr * yc * (wrl^2 * yr + wcl^2 * yc) * yl + (wrl^2 * yr - wcl^2 * yc)^2 * yl^2)
    A[is.nan(A)] <- 0
    A
}

# General helpers #####

#' @noRd
#' @title Calculate Lorentz Curve values
#'
#' @description
#' Calculates the values of a Lorentz Curve for a vector of input values `x`.
#' The Lorentz Curve is defined as \eqn{A \cdot \frac{\lambda}{\lambda^2 +
#' (x_i - x_0)^2}}.
#'
#' @param x Numeric vector of x values.
#' @param x0 Center of the Lorentz curve.
#' @param A Amplitude parameter of the Lorentz curve.
#' @param lambda Half width at half height of the Lorentz curve.
#'
#' @return Numeric vector of y values.
#'
#' @examples
#' x <- 1:10
#' x0 <- 5
#' A <- 10
#' lambda <- 2
#' y1 <- lorentz(x, x0, A, lambda)
#' y2 <- A * pi * dcauchy(x, location = x0, scale = lambda)
#' stopifnot(all.equal(y1, y2))
lorentz <- function(x, x0, A, lambda, lcpar = NULL) {
    # For details see [Wikipedia > Cauchy_distribution > Properties_of_PDF]
    # (https://en.wikipedia.org/wiki/Cauchy_distribution#Properties_of_PDF)
    # in particular the formula below sentence "In physics, a three-parameter
    # Lorentzian function is often used".
    if (!is.null(lcpar)) {
        nams <- names(lcpar)
        if ("A" %in% nams) A <- lcpar[["A"]]
        if ("lambda" %in% nams) lambda <- lcpar[["lambda"]]
        if ("x_0" %in% nams) x0 <- lcpar[["x_0"]]
        if ("x0" %in% nams) x0 <- lcpar[["x0"]]
        if ("w" %in% nams) x0 <- lcpar[["w"]]
    }
    A * (lambda / (lambda^2 + (x - x0)^2))
}

lorentz_sup <- function(x, x0, A, lambda, lcpar = NULL) {
    if (!is.null(lcpar)) {
        nams <- names(lcpar)
        if ("A" %in% nams) A <- lcpar[["A"]]
        if ("lambda" %in% nams) lambda <- lcpar[["lambda"]]
        if ("x_0" %in% nams) x0 <- lcpar[["x_0"]]
        if ("x0" %in% nams) x0 <- lcpar[["x0"]]
        if ("w" %in% nams) x0 <- lcpar[["w"]]
    }
    sapply(x, function(xi) {
        sum(abs(A * (lambda / (lambda^2 + (xi - x0)^2))))
    })
}

#' @noRd
#' @title Calculate Lorentz Curve Integrals
#' @description
#' Calculates the integral of a Lorentz curve for a vector of input values `x`.
lorentz_int <- function(x0, A, lambda, lcpar = NULL, limits = NULL) {
    if (is.list(lcpar)) {
        nams <- names(lcpar)
        if ("A" %in% nams) A <- lcpar$A
        if ("lambda" %in% nams) lambda <- lcpar$lambda
        if ("x_0" %in% nams) x0 <- lcpar$x_0
        if ("x0" %in% nams) x0 <- lcpar$x0
        if ("w" %in% nams) x0 <- lcpar$w
    }
    if (is.null(limits)) {
        A * pi
    } else {
        a <- min(limits)
        b <- max(limits)
        A * (atan((b - x0) / lambda) - atan((a - x0) / lambda))
    }
}

mse <- function(y, yhat, normed = FALSE) {
    if (normed) {
        mean(((y / sum(y)) - (yhat / sum(yhat)))^2)
    } else {
        mean((y - yhat)^2)
    }
}

#' @noRd
#' @description
#' Before version 1.2 of 'metabodecon', the deconvolution functions
#' `generate_lorentz_curves()` and `MetaboDecon1D()` wrote their output
#' partially as txt files to their input folder. The txt files were named
#' "SPEC_NAME parameter.txt" and "SPEC_NAME approximated_spectrum.txt". Since
#' version 1.2 these txt files are no longer created by default, to prevent
#' accidental modifications of the input folders. However, to stay backwards
#' compatible, functions that used to read "SPEC_NAME parameter.txt" and
#' "SPEC_NAME approximated_spectrum.txt" still accept them as input (e.g.
#' `gen_feat_mat()`). I.e., in order to test this functionality, we still need a
#' way to create the corresponding txt files (which is no longer done by
#' `generate_lorentz_curves()`). That's the purpose of this function: it takes
#' the output of `generate_lorentz_curves()` as input and creates the (now
#' deprecated) "SPEC_NAME parameter.txt" and "SPEC_NAME
#' approximated_spectrum.txt" in folder `outdir`.
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
