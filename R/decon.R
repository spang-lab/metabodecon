# Public API #####

#' @export
#'
#' @title Deconvolute one or more NMR spectra
#'
#' @description Deconvolutes NMR spectra by modeling each detected signal within
#' a spectrum as Lorentz Curve.
#'
#' @param ask Logical. Whether to ask for user input during the deconvolution
#' process. If FALSE, the provided default values will be used.
#'
#' @param x A `spectrum` or `spectra` object as described in [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
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
#' @param use_rust Logical. Whether to use the Rust backend for deconvolution.
#' Requires the [mdrb](https://github.com/spang-lab/mdrb) package. If TRUE and
#' mdrb is missing, an error is thrown. If FALSE, the R implementation is used.
#' If NULL, the Rust backend is used if available, otherwise the R implementation
#' is used.
#'
#' @return A 'decon2' object as described in [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
#' @details
#'
#' First, an automated curvature based signal selection is performed. Each
#' signal is represented by 3 data points to allow the determination of initial
#' Lorentz curves. These Lorentz curves are then iteratively adjusted to
#' optimally approximate the measured spectrum.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' ## Deconvolute a single spectrum
#' spectrum <- sim[1]
#' decon <- deconvolute(spectrum)
#'
#' ## Read multiple spectra from disk and deconvolute at once
#' spectra_dir <- metabodecon_file("sim_subset")
#' spectra <- read_spectra(spectra_dir)
#' decons <- deconvolute(spectra, sfr = c(3.55, 3.35))
deconvolute <- function(x,
    nfit=3,    smopts=c(2, 5), delta=6.4,    sfr=NULL,   wshw=0,
    ask=FALSE, force=FALSE,    verbose=TRUE, nworkers=1, use_rust=FALSE
) {
    # Check inputs
    stopifnot(
        is_spectrum_or_spectra(x),    is_int_or_null(nfit, 1),
        is_int_or_null(smopts, 2),    is_num_or_null(delta, 1),
        is_num_or_null(sfr, 2),       is_num_or_null(wshw, 1),
        is_bool(ask, 1),              is_bool(force, 1),
        is_bool(verbose, 1),          is_int(nworkers, 1),
        is_bool_or_null(use_rust, 1)
    )

    # Set suitable defaults
    sfr <- sfr %||% quantile(x$cs %||% x[[1]]$cs, c(0.9, 0.1))
    if (isTRUE(use_rust)) check_mdrb(stop_on_fail = TRUE)
    if (is.null(use_rust)) use_rust <- check_mdrb()

    # Perform deconvolution
    decons2 <- deconvolute_spectra(x,
        nfit, smopts, delta, sfr, wshw,
        ask, force, verbose, bwc=2,
        use_rust, nw=nworkers, igr=list(), rtyp="decon2"
    )

    # Convert and return
    if (length(decons2) == 1) decons2[[1]] else decons2
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
#' @author 2024-2025 Tobias Schmidt: initial version.
deconvolute_spectra <- function(x,
    nfit=3,        smopts=c(2, 5),    delta=6.4,      sfr=c(3.55, 3.35),
    wshw=0,        ask=FALSE,         force=FALSE,    verbose=TRUE,
    bwc=2,         use_rust=FALSE,    nw=1,           igr=list(),
    rtyp="idecon"
) {

    # Check inputs
    assert(
        is_spectrum(x) || is_spectra(x),
        is_int(nfit, 1), is_int(smopts, 2), is_num(delta, 1),
        is_bool(ask, 1), is_bool(force, 1), is_bool(verbose, 1),
        is_num(bwc, 1),
        is_num(sfr, 2)  || is_list_of_nums(sfr, length(x), 2),
        is_num(wshw, 1) || is_list_of_nums(wshw, length(x), 1),
        if (rtyp == "rdecon") isTRUE(use_rust) else is_bool(use_rust),
        is_char(rtyp, 1, "(decon[0-2]|idecon|rdecon)")
    )

    # Configure logging
    if (!verbose) local_options(toscutil.logf.file = nullfile())

    # Init locals
    spectra <- as_spectra(x)
    ns <- length(spectra)
    nc2 <- ceiling(detectCores() / 2)
    nw_apply <- if (use_rust) 1 else if (nw == "auto") min(nc2, ns) else nw
    nw_apply_str <- if (nw_apply == 1) "1 worker" else sprintf("%d workers", nw_apply)
    nw_deconv <- if (!use_rust) 1 else if (nw == "auto") min(nc2, ns) else nw
    ns_str <- if (ns == 1) "1 spectrum" else sprintf("%d spectra", ns)
    adjno <- get_adjno(spectra, ask)
    sfr_list <- get_sfr(spectra, sfr, ask, adjno)
    wshw_list <- get_wshw(spectra, wshw, ask, adjno)
    smopts_list <- get_smopts(spectra, smopts)
    igr_list <- list(list())

    # Deconvolute spectra
    logf("Starting deconvolution of %s using %s", ns_str, nw_apply_str)
    starttime <- Sys.time()
    decon_list <- mcmapply(nw_apply, deconvolute_spectrum,
        spectra,
        nfit, smopts_list, delta, sfr_list, wshw_list,
        ask, force, verbose, bwc,
        use_rust, nw_deconv, igr_list, rtyp
    )
    decons <- as_collection(decon_list, rtyp)
    duration <- format(round(Sys.time() - starttime, 3))
    logf("Finished deconvolution of %s in %s", ns_str, duration)

    # Return
    decons
}

#" EXAMPLES
#" urine_1_path <- metabodecon_file("urine_1")
#" urine_1 <- read_spectrum(urine_1_path)
#" system.time(deconvolute_spectrum(urine_1, use_rust = TRUE, nw = 4))
#' @noRd
#' @inheritParams deconvolute_spectra
#' @author 2024-2025 Tobias Schmidt: initial version.
deconvolute_spectrum <- function(x,
    nfit=3, smopts=c(2, 5), delta=6.4, sfr=c(3.55, 3.35), wshw=0,
    ask=FALSE, force=FALSE, verbose=TRUE, bwc=2,
    use_rust=FALSE, nw=1, igr=list(), rtyp="idecon"
) {
    # Check inputs
    assert(
        is_spectrum(x),
        is_int(nfit, 1),  is_int(smopts, 2),     is_num(delta, 1),
        is_num(sfr, 2),   is_num(wshw, 1),       is_bool(force, 1),
        is_num(bwc, 1),   is_bool(use_rust, 1),  is_int(nw, 1),
        is_list_of_nums(igr, nv=2),
        is_char(rtyp, 1, "(decon[0-2]|idecon|rdecon)")
    )

    # Init locals
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())
    name <- get_name(x)
    suffix <- if (use_rust) " using Rust backend" else ""
    args <- get_args(deconvolute_spectrum, ignore = "x")

    # Deconvolute
    logf("Starting deconvolution of %s%s", name, suffix)
    if (use_rust) {
        mdrb_spectrum <- mdrb::Spectrum$new(x$cs, x$si, sfr)
        mdrb_deconvr <- mdrb::Deconvoluter$new()
        mdrb_deconvr$set_moving_average_smoother(smopts[1], smopts[2])
        mdrb_deconvr$set_noise_score_selector(delta)
        mdrb_deconvr$set_analytical_fitter(nfit)
        mdrb_decon <- if (nw > 1) {
            mdrb_deconvr$set_threads(nw)
            mdrb_deconvr$par_deconvolute_spectrum(mdrb_spectrum)
        } else {
            mdrb_deconvr$deconvolute_spectrum(mdrb_spectrum)
        }
        decon <- new_rdecon(x, args, mdrb_spectrum, mdrb_deconvr, mdrb_decon)
    } else {
        # Sys.setenv(RAYON_NUM_THREADS=nw) # Must be set before R is started
        spectrum <- x
        sf <- c(1e3, 1e6)
        
        # Create scaled intensities
        y_raw <- spectrum$si # Raw signal intensities
        y_scaled <- y_raw / sf[2] # Scaled signal intensities
        
        # Process spectrum step by step
        y_nows <- rm_water_signal(spectrum, y_scaled, wshw, bwc, sf)
        y_pos <- rm_negative_signals(y_nows)
        smooth_results <- smooth_signals(y_pos, smopts[1], smopts[2], bwc)
        
        peak_results <- find_peaks(smooth_results$y_smooth)
        filtered_peak <- filter_peaks(spectrum, peak_results$peak, sfr, delta, force, bwc, sf)
        lc_results <- fit_lorentz_curves(spectrum, filtered_peak, smooth_results$y_smooth, nfit, bwc, sf)
        
        # Build idecon-like object for conversion
        cs <- spectrum$cs
        n <- length(cs)
        dp <- seq(n - 1, 0, -1)
        sdp <- seq((n - 1) / sf[1], 0, -1 / sf[1])
        hz <- spectrum$meta$fq
        ppm_range <- diff(range(cs))
        ppm_max <- max(cs)
        ppm_min <- min(cs)
        ppm_step <- ppm_range / (n - 1)
        ppm_nstep <- ppm_range / n
        name <- spectrum$meta$name
        meta <- spectrum$meta
        
        decon <- list(
            y_raw = y_raw, y_scaled = y_scaled, n = n, dp = dp, sdp = sdp, sf = sf,
            ppm = cs, hz = hz, ppm_range = ppm_range, ppm_max = ppm_max, ppm_min = ppm_min,
            ppm_step = ppm_step, ppm_nstep = ppm_nstep, name = name, meta = meta,
            args = args, y_nows = y_nows, y_pos = y_pos, Z = smooth_results$Z,
            y_smooth = smooth_results$y_smooth, d = peak_results$d, peak = filtered_peak,
            lci = lc_results$lci, lca = lc_results$lca, lcr = lc_results$lcr
        )
        class(decon) <- "idecon"
    }
    logf("Formatting return object as %s", rtyp)
    convert <- switch(rtyp,
        "decon0"=as_decon0, "decon1"=as_decon1, "decon2"=as_decon2,
        "idecon"=as_idecon, "rdecon"=as_rdecon
    )
    decon <- convert(decon)
    logf("Finished deconvolution of %s", name)

    decon
}

# Helpers for deconvolute_spectra #####

#' @noRd
#' @param x A `spectra` object or any other metabodecon collection object.
#' @description
#' Get number of spectrum that should be used to adjust all others. If ask is
#' FALSE or every spectrum should be adjusted individually, zero is returned.
#' @author 2024-2025 Tobias Schmidt: initial version.
get_adjno <- function(x, ask) {
    assert(is_spectra(x), is_bool(ask, 1))
    if (isFALSE(ask) ||
        length(x)==1 ||
        isFALSE(get_yn_input("Use same parameters for all spectra?"))) return(0)
    numbers <- seq_along(x)
    names <- get_names(x)
    names_str <- paste(numbers, names, sep = ": ", collapse = ", ")
    prompt <- sprintf("Number of spectrum for adjusting parameters? (%s)", names_str)
    adjno <- get_num_input(prompt, min = 1, max = length(x), int = TRUE)
    adjno
}

#' @noRd
#' @description
#' Converts one or more SFR vectors to list of vectors of the correct length. If
#' `ask` is TRUE, let user confirm entries.
#' @param x Any metabodecon collection object.
#' @param sfr SFR defaults. Can be a vector of length 2 or a list. If a list is
#' provided, it must have the same length as `x`.
#' @param ask Ask user to confirm suggested defaults?
#' @param adjno Number of spectrum to show when user is asked for confirmation.
#' If 0, the user is asked to confirm each entry individually.
#' @author 2024-2025 Tobias Schmidt: initial version.
get_sfr <- function(x, sfr, ask, adjno) {
    n <- length(x)
    if (is_num(sfr, 2)) sfr <- rep(list(sfr), n)
    if (!is_list_of_nums(sfr, n, 2)) stop("sfr should be a [list of] num(2)")
    if (ask && adjno == 0) sfr <- lapply(seq_len(n), function(i) confirm_sfr(x[[i]], sfr[[i]]))
    if (ask && adjno >= 1) sfr <- rep(list(confirm_sfr(x[[adjno]], sfr[[adjno]])), n)
    names(sfr) <- get_names(x)
    sfr
}

#' @noRd
#' @description Same as [get_sfr()], but for WSHW.
#' @author 2024-2025 Tobias Schmidt: initial version.
get_wshw <- function(x, wshw, ask, adjno) {
    n <- length(x)
    if (is_num(wshw, 1)) wshw <- rep(list(wshw), n)
    if (is_num(wshw, n)) wshw <- as.list(wshw)
    if (!is_list_of_nums(wshw, n, 1)) stop("wshw should be a [list of] num(1)")
    if (ask && adjno == 0) wshw <- mapply(confirm_wshw, x, wshw, SIMPLIFY = FALSE)
    if (ask && adjno >= 1) wshw <- rep(list(confirm_wshw(x[[adjno]], wshw[[adjno]])), n)
    names(wshw) <- names(x)
    wshw
}

#' @noRd
#' @title Get Smoothing Options
#'
#' @description
#' Convert one or more SMOPTS vectors to list of vectors of the correct length.
#'
#' @param x
#' Any metabodecon collection object.
#'
#' @param smopts
#' Default smoothing options. Can be a vector of length 2 or a list of such
#' vectors. If a list if provided, it must have the same length as `x`.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
get_smopts <- function(x, smopts) {
    n <- length(x)
    if (is_int(smopts, 2)) smopts <- rep(list(smopts), n)
    else if (is_list_of_nums(smopts, n, 2)) {}
    else stop("smopts should be a [list of] int(2)")
    names(smopts) <- names(x)
    smopts
}

# Helpers for deconvolute_spectrum #####

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D. Added code for bwc > 1.
rm_water_signal <- function(spectrum, y_scaled, wshw, bwc, sf = c(1e3, 1e6)) {
    assert(is_spectrum(spectrum), is_num(y_scaled), is_num(wshw, 1), is_num(bwc, 1))
    logf("Removing water signal")
    cs <- spectrum$cs
    if (bwc >= 1) {
        ppm_center <- (cs[1] + cs[length(cs)]) / 2
        idx_wsr <- which(cs > ppm_center - wshw & cs < ppm_center + wshw)
        y_nows <- y_scaled
        y_nows[idx_wsr] <- min(y_nows)
    } else {
        wsr <- enrich_wshw(wshw, spectrum, sf)
        left <- wsr$left_dp
        right <- wsr$right_dp
        idx_wsr <- right:left # (1)
        y_nows <- y_scaled
        y_nows[idx_wsr] <- 0.01 / sf[2]
        # (1) Order is important here, because right and left are floats. Example:
        # right <- 3.3; left <- 1.4
        # right:left == c(3.3, 2.3)
        # left:right == c(1.4, 2.4)
    }
    y_nows
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D. Added code for bwc > 1.
rm_negative_signals <- function(y_nows) {
    logf("Removing negative signals")
    if (is.null(y_nows)) stop("Water signal not removed yet. Call `rm_water_signal()` first.")
    y_pos <- abs(y_nows)
    y_pos
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
#'
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
smooth_signals <- function(y_pos, reps = 2, k = 5, verbose = TRUE) {
    if (verbose) logf("Smoothing signals")
    if (k %% 2 == 0) stop("k must be odd")
    n <- length(y_pos) # number of data points in total
    ws <- floor(k / 2) # window size left/right of center
    ct <- seq_len(n) # center positions
    lb <- pmax(ct - ws, 1) # left borders
    rb <- pmin(ct + ws, n) # right borders
    nw <- rb - lb + 1 # number of data points in window
    Z <- data.frame(y_pos)
    for (j in seq_len(reps)) {
        zj <- Z[[j]]
        zk <- sapply(seq_len(n), function(i) sum(zj[lb[i]:rb[i]]))
        zk <- (1 / nw) * zk
        Z[[j + 1]] <- zk
    }
    list(Z = Z[, -1], y_smooth = Z[, reps + 1])
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
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
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

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
find_peaks <- function(y_smooth) {
    logf("Starting peak selection")
    d <- calc_second_derivative(y = y_smooth)
    a <- abs(d)
    m <- length(d)
    dl <- c(NA, d[-m]) # dl[i] == d[i-1]
    dr <- c(d[-1], NA) # dr[i] == d[i+1]
    center <- which(d < 0 & d <= dl & d < dr)
    peak <- data.frame(left = NA, center = center, right = NA, score = NA)
    for (i in seq_along(center)) {
        j <- center[i]
        l <- peak$left[i] <- get_left_border(j, d)
        r <- peak$right[i] <- get_right_border(j, d, m)
        peak$score[i] <- get_peak_score(j, l, r, a)
    }
    logf("Detected %d peaks", length(center))
    list(d = d, peak = peak)
}

#' @noRd
#' @title WORK IN PROGRESS
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D. Added code for bwc > 1.
#'
#' @examples
#' peak <- list(score = c(1, 5, 2), center = c(1, 2, 3))
#' sdp <- c(3, 2, 1)
#' ispec <- named(peak, sdp)
#' sfr <- list(left_sdp = 2.8, right_sdp = 1.2)
#' rm3 <- filtered_ispec <- filter_peaks(ispec, sfr)
#' rm2 <- filtered_ispec <- filter_peaks(ispec, sfr, delta = 1)
filter_peaks <- function(spectrum, peak, sfr, delta = 6.4, force = FALSE, bwc = 1, sf = c(1e3, 1e6)) {
    assert(is_spectrum(spectrum))
    logf("Removing peaks with low pscores")
    cs <- spectrum$cs
    n <- length(cs)
    # Create scaled data point numbers for backwards compatibility
    sdp <- seq((n - 1) / sf[1], 0, -1 / sf[1])
    plb <- peak$left
    prb <- peak$right
    pct <- peak$center
    psc <- peak$score
    pok <- !is.na(plb) & !is.na(pct) & !is.na(prb)
    if (bwc < 1) {
        sfr_enriched <- enrich_sfr(sfr, spectrum, sf)
        in_left_sfr <- sdp[pct] >= sfr_enriched$left_sdp
        in_right_sfr <- sdp[pct] <= sfr_enriched$right_sdp
    } else {
        in_left_sfr <- cs[pct] >= max(sfr)
        in_right_sfr <- cs[pct] <= min(sfr)
    }
    in_sfr <- in_left_sfr | in_right_sfr
    if (sum(in_sfr) > 1) {
        mu <- mean(psc[in_sfr])
        sigma <- sd(psc[in_sfr])
    } else {
        if (!force) stop(paste(
            "Not enough signals found in signal free region.",
            "Please double check deconvolution parameters."
        ))
        logf(paste(
            "Not enough signals found in signal free region.",
            "This is a clear indication that the deconvolution parameters",
            "are not set correctly. Continuing anyways without dynamic peak",
            "filtering, because `force` is TRUE. Note that this might",
            "increase runtime drastically."
        ))
        mu <- 0
        sigma <- 0
    }
    peak$high <- pok & (psc > mu + delta * sigma)
    peak$region <- "norm"
    peak$region[in_left_sfr] <- "sfrl"
    peak$region[in_right_sfr] <- "sfrr"
    logf("Removed %d peaks", sum(!peak$high))
    peak
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D. Added code for bwc > 1.
fit_lorentz_curves <- function(spectrum, peak, y_smooth, nfit = 3, bwc = 1, sf = c(1e3, 1e6)) {
    logf("Initializing Lorentz curves")
    lci <- lc <- init_lc(spectrum, peak, y_smooth, sf) # Lorentz Curves Initialized
    lca <- vector("list", length = nfit) # Lorentz Curves Approximated
    logf("Refining Lorentz Curves")
    for (i in 1:nfit) lca[[i]] <- lc <- refine_lc_v14(spectrum, peak, y_smooth, lc$Z, sf)
    A <- lc$A
    lambda <- lc$lambda
    w <- lc$w
    cs <- spectrum$cs
    n <- length(cs)
    if (bwc < 1) {
        sdp <- seq((n - 1) / sf[1], 0, -1 / sf[1])
        limits <- c(0, max(sdp) + (1 / sf[1]))
        integrals <- lorentz_int(w, A, lambda, limits = limits)
    } else {
        integrals <- A * (- pi)
    }
    lcr <- list(A = A, lambda = lambda, w = w, integrals = integrals)
    list(lci = lci, lca = lca, lcr = lcr)
}

# Helpers for get_sfr and get_wshw #####

#' @noRd
#' @description Repeatedly ask the user to confirm/refine SFR borders.
#' @param x Any metabodecon object.
#' @author 2024-2025 Tobias Schmidt: initial version.
confirm_sfr <- function(x, sfr = c(11.44494, -1.8828)) {
    si <- x$y_scaled %||% x$si
    cs <- x$ppm %||% x$cs
    cs_min <- min(cs)
    cs_max <- max(cs)
    plot_sfr(cs, si, sfr)
    sfr_ok <- get_yn_input("Borders of signal free region (green) correctly selected?")
    while (!sfr_ok) {
        get_border <- function(msg) get_num_input(msg, cs_min, cs_max)
        sfr[1] <- get_border("Choose another left border: [e.g. 12]")
        sfr[2] <- get_border("Choose another right border: [e.g. -2]")
        plot_sfr(cs, si, sfr)
        sfr_ok <- get_yn_input("Borders of signal free region (green) correctly selected?")
    }
    sfr
}

#' @noRd
#' @description Repeatedly ask the user to confirm/refine the WSHW.
#' @author 2024-2025 Tobias Schmidt: initial version.
confirm_wshw <- function(x, wshw) {
    cs <- x$cs %||% x$ppm
    si <- x$si %||% x$y_scaled
    plot_ws(cs, si, wshw)
    ws_ok <- get_yn_input("Water artefact fully inside blue area?")
    while (!ws_ok) {
        wshw <- get_num_input("Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154]")
        plot_ws(cs, si, wshw)
        ws_ok <- get_yn_input("Water artefact fully inside blue area?")
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
#' @author 2024-2025 Tobias Schmidt: initial version.
enrich_sfr <- function(sfr, spectrum, sf = c(1e3, 1e6)) {
    assert(is_spectrum(spectrum))
    left_ppm <- sfr[1]
    right_ppm <- sfr[2]
    # Calculate values needed from spectrum
    cs <- spectrum$cs
    n <- length(cs)
    ppm_max <- max(cs)
    ppm_range <- diff(range(cs))
    ppm_nstep <- ppm_range / n # Backwards compatible calculation
    # Convert ppm to data point indices and scaled data point indices
    left_dp <- (n + 1) - (ppm_max - left_ppm) / ppm_nstep
    left_sdp <- left_dp / sf[1]
    right_dp <- (n + 1) - (ppm_max - right_ppm) / ppm_nstep
    right_sdp <- right_dp / sf[1]
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
#' @author 2024-2025 Tobias Schmidt: initial version.
enrich_wshw <- function(wshw, spectrum, sf = c(1e3, 1e6)) {
    assert(is_spectrum(spectrum))
    # Calculate values needed from spectrum
    cs <- spectrum$cs
    n <- length(cs)
    ppm_range <- diff(range(cs))
    ppm_nstep <- ppm_range / n # Backwards compatible calculation
    # Calculate water signal region
    hwidth_ppm <- wshw
    hwidth_dp <- hwidth_ppm / ppm_nstep
    center_dp <- n / 2
    right_dp <- center_dp + hwidth_dp
    left_dp <- center_dp - hwidth_dp
    center_ppm <- cs[center_dp]
    right_ppm <- cs[right_dp]
    left_ppm <- cs[left_dp]
    if (left_dp <= 1 || right_dp >= n) stop("WSR is out of range")
    named(
        left_ppm, right_ppm, center_ppm, hwidth_ppm,
        left_dp, right_dp, center_dp, hwidth_dp
    )
}

# Helpers for find_peak #####

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
calc_second_derivative <- function(y) {
    n <- length(y)
    x <- c(NA, y[-n]) # x[i] == y[i-1]
    z <- c(y[-1], NA) # z[i] == y[i+1]
    d <- x + z - 2 * y
    d
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
get_right_border <- function(j, d, m) {
    r <- j + 1
    while (r < m) { # use r<m instead of r<=m because c4 requires d[r+1]
        c1 <- d[r] > d[r - 1]
        c2 <- d[r] >= d[r + 1]
        c3 <- d[r] < 0
        c4 <- d[r + 1] >= 0
        is_right_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_right_border)) return(r)
        r <- r + 1
    }
    NA
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
get_left_border <- function(j, d) {
    l <- j - 1
    while (l > 1) { # use l>1 instead of l>=1 because c4 requires d[l-1]
        c1 <- d[l] > d[l + 1]
        c2 <- d[l] >= d[l - 1]
        c3 <- d[l] < 0
        c4 <- d[l - 1] >= 0
        is_left_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_left_border)) return(l)
        l <- l - 1
    }
    NA
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
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
#'
#' @examples
#'
#' #      ____________________________________________________
#' #     |____2________5___________9_______12____14____16_____|
#' #     |             x                                      |
#' #     |          x  x  x  x                    x           |
#' #     |       x  x  x  x  x  x              x  x           |
#' #     |_.__x__x__x__x__x__x__x__x__.__.__x__x__x__x__x__.__|
#'
#' y <- c( 0, 1, 2, 3, 4, 3, 3, 2, 1, 0, 0, 1, 2, 3, 1, 1, 0  )
#' a <- c(NA, 0, 0, 0, 2, 1, 1, 0, 0, 1, 1, 0, 0, 3, 2, 1, NA )
#' all.equal(a, abs(calc_second_derivative(y)))
#'
#' s1 <- get_peak_score( 5, 2,   9, a)
#' s2 <- get_peak_score(14, 12, 16, a)
#' stopifnot(s1 == min(sum(a[ 2:5]),  sum(a[ 5:9] )))
#' stopifnot(s2 == min(sum(a[12:14]), sum(a[14:16])))
get_peak_score <- function(j, l, r, a) {
    if (any(is.na(a[c(l, j, r)]))) {
        0
    } else {
        min(sum(a[l:j]), sum(a[j:r]))
    }
}

# Helpers for fit_lorentz_curves #####

#' @noRd
#' @title Initialize Lorentz Curve Parameters
#' @param spec
#' List with elements: `x`, `y`, `peak` where `peak` is a list with elements
#' `center`, `left`, `right` and `high`.
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
init_lc <- function(spectrum, peak, y_smooth, sf = c(1e3, 1e6), verbose = TRUE) {

    # Init values
    p <- peak
    ir <- p$right[p$high]  #
    ic <- p$center[p$high] # Index of each peak triplet position (PTP)
    il <- p$left[p$high]   #
    lmr <- sort(unique(c(il, ic, ir))) # Combined PTP indices
    rr <- match(ir, lmr) #
    rc <- match(ic, lmr) # Rank of each PTP
    rl <- match(il, lmr) #
    # Create scaled data point numbers for backwards compatibility
    cs <- spectrum$cs
    n <- length(cs)
    x <- sdp <- seq((n - 1) / sf[1], 0, -1 / sf[1])
    y <- y_smooth # X and Y value for each data point
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
    mse <- mse(y[lmr], rowSums(Z))
    if (verbose) logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, as, ds)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.\cr
refine_lc_v14 <- function(spectrum, peak, y_smooth, Z, sf = c(1e3, 1e6)) {

    # Init x and y values
    cs <- spectrum$cs
    n <- length(cs)
    x <- seq((n - 1) / sf[1], 0, -1 / sf[1])  # scaled data point numbers for backwards compatibility
    y <- y_smooth # x and y value for each data point

    # Init peak related variables
    p <- peak
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
    mse <- mse(y[lmr], rowSums(Z))
    logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, sl, sc, sr, ql, qc, qr)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D
#' based on Appendix E of Koh et. al. 2009.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
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

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D
#' based on Appendix E of Koh et. al. 2009.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
calc_lambda <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr) {
    lambda <- -((sqrt(abs((-wc^4 * yc^2 * yrl^2 - wr^4 * yr^2 * ycl^2 - wl^4 * yrc^2 * yl^2 + 4 * wc * wl^3 * yc * ((-yr) + yc) * yl^2 + 4 * wc^3 * wl * yc^2 * yl * ((-yr) + yl) + 4 * wr^3 * yr^2 * ycl * (wc * yc - wl * yl) + 4 * wr * yr * (wc^3 * yc^2 * yrl - wc * wl^2 * yc * (yr + yc - 2 * yl) * yl + wl^3 * yrc * yl^2 - wc^2 * wl * yc * yl * (yr - 2 * yc + yl)) + 2 * wc^2 * wl^2 * yc * yl * (yr^2 - 3 * yc * yl + yr * (yc + yl)) + 2 * wr^2 * yr * (-2 * wc * wl * yc * yl * (-2 * yr + yc + yl) + wl^2 * yl * (yr * (yc - 3 * yl) + yc * (yc + yl)) + wc^2 * yc * (yr * (-3 * yc + yl) + yl * (yc + yl)))))))) / (2 * sqrt((wr * yr * ycl + wl * yrc * yl + wc * yc * ((-yr) + yl))^2))
    lambda[is.nan(lambda)] <- 0
    lambda
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D
#' based on Appendix E of Koh et. al. 2009.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
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
#' @details
#' 1. The argument names as based on the names used by Koh et al. (2009).
#' 2. In Wikipedia, Lorentz Curves are described in article
#' [Cauchy_distribution]. The formula below sentence "In physics, a
#' three-parameter Lorentzian function is often used" (section
#' [Properties_of_PDF]) is equivalent to the one used by Koh. et al (2009),
#' although the variables have different names.
#'
#' [Cauchy_distribution]: https://en.wikipedia.org/wiki/Cauchy_distribution
#' [Properties_of_PDF]: https://en.wikipedia.org/wiki/Cauchy_distribution#Properties_of_PDF
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
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

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' @author 2024-2025 Tobias Schmidt: initial version.
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

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
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
#' @author 2024-2025 Tobias Schmidt: initial version.
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

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
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
