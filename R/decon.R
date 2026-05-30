# API #####

#' @export
#'
#' @title Deconvolute one or more NMR spectra
#'
#' @description Deconvolutes NMR spectra by modeling each detected signal within
#' a spectrum as Lorentz Curve.
#'
#' @param x A `spectrum` or `spectra` object as described in [metabodecon::metabodecon-classes].
#'
#' @param delta Threshold for peak filtering. Higher values result in more peaks
#' being filtered out. A peak is filtered if its score is below \eqn{\mu +
#' \sigma \cdot \delta}{mu + s * delta}, where \eqn{\mu}{mu} is the average
#' peak score in the signal-free region (SFR), and \eqn{\sigma}{s} is the
#' standard deviation of peak scores in the SFR. See 'Details'.
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
#' @param smit Integer. Number of smoothing iterations. See 'Details'.
#'
#' @param smws Integer. Smoothing window size (number of data points; must be
#' odd). See 'Details'.
#'
#' @param verbose Logical. Whether to print log messages during the
#' deconvolution process.
#'
#' @param use_rust Controls the deconvolution backend. `FALSE` or any numeric
#' value `< 1` (default) uses the R implementation. `TRUE` or any numeric
#' value `>= 1` uses the Rust backend via
#' [mdrb](https://github.com/spang-lab/mdrb). `NULL` auto-detects: uses Rust
#' if available, otherwise R. When set to `TRUE` / `>= 1` and mdrb is not
#' installed, an error is thrown.
#'
#' @param npmax Integer. Maximum number of peaks allowed in the result. If
#' `npmax >= 1`, the `nfit`, `smit`, `smws` and `delta` arguments are ignored
#' and a grid search over predefined parameter combinations is performed
#' instead. The combination with the smallest residual area ratio and fewer
#' than `npmax` peaks is selected. Grid search results are cached to disk
#' automatically.
#'
#' @param igrs Ignore regions. List of length-2 numeric vectors specifying the
#' start and endpoints of the chemical shift regions to ignore during
#' deconvolution. Peaks whose centers fall inside any ignore region are
#' excluded from fitting.
#'
#' @return A 'decon2' object as described in [metabodecon::metabodecon-classes].
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
#' spectrum <- sim[[1]]
#' decon <- deconvolute(spectrum)
#'
#' ## Read multiple spectra from disk and deconvolute at once
#' spectra_dir <- metabodecon_file("sim_subset")
#' spectra <- read_spectra(spectra_dir)
#' decons <- deconvolute(spectra, sfr = c(3.55,3.35))
deconvolute <- function(
    x,
    nfit=3, smit=2, smws=5, delta=6.4, npmax=0,
    sfr=NULL, igrs=list(),
    use_rust=FALSE, verbose=TRUE, nworkers=1
) {

    # Check inputs
    stopifnot(
        inherits(x, "spectrum") || inherits(x, "spectra"), is_int(nfit, 1),
        is_int_or_null(smit, 1),   is_int_or_null(smws, 1),
        is_num_or_null(delta, 1),  is_int(npmax, 1),
        is_num_or_null(sfr, 2),    is_list_of_nums(igrs, nv=2),
        is_bool_or_num(use_rust),  is_int(nworkers, 1),
        is_bool(verbose, 1)
    )
    if (use_rust >= 1) check_mdrb(stop_on_fail = TRUE)

    # Perform deconvolution
    decons2 <- deconvolute_spectra(
        x=x, nfit=nfit, smit=smit, smws=smws, delta=delta, npmax=npmax,
        sfr=sfr, igrs=igrs,
        use_rust=use_rust, verbose=verbose, nworkers=nworkers
    )

    # Convert and return
    if (length(decons2) == 1) decons2[[1]] else decons2
}

#' @export
#' @rdname deconvolute
#'
#' @title Default deconvolution-parameter grid
#'
#' @description
#' Returns the default grid of `(nfit, smit, smws, delta)` combinations used
#' by [metabodecon::deconvolute()] when `npmax >= 1`. Useful as the `deg`
#' argument to [metabodecon::fit_mdm()].
#'
#' @param conf Character string selecting a configuration. Currently only
#'   `"default"` is supported.
#'
#' @return A data frame with columns `nfit`, `smit`, `smws`, `delta`.
#'
#' @examples
#' get_deg()
get_deg <- function(conf="default") {
    expand.grid2(nfit=10, smit=1:3, smws=c(3,5,7,9), delta=(1:5)*1.6)
}

# Internal #####

#' @noRd
#' @author 2024-2026 Tobias Schmidt: initial version.
deconvolute_spectra <- function(
    x,
    nfit=3, smit=2, smws=5, delta=6.4, npmax=0,
    sfr=NULL, igrs=list(),
    use_rust=FALSE, verbose=TRUE, nworkers=1, # end of public args
    force=FALSE, full=TRUE
) {

    # Init locals
    if (!verbose) local_options(toscutil.logf.file = nullfile())
    x <- as_spectra(x)
    ns <- length(x)
    nw <- min(half_cores(), ns, nworkers)

    # Attach grid-search results to each spectrum (only when npmax > 0)
    if (npmax >= 1) x <- grid_deconvolute_spectra(
        x=x, sfr=sfr, verbose=verbose, nworkers=nw, use_rust=use_rust
    )

    # Deconvolute spectra. The shared chemical-shift grid (`cssh`) and
    # per-spectrum supsh used by alignment are not built here: a shared
    # grid is a property of *a collection being aligned* (clupa) rather
    # than of a single deconvolution, so clupa() owns that construction.
    logf("Starting deconvolution (spectra: %d, workers: %d)", ns, nw)
    starttime <- Sys.time()
    args <- get_args(deconvolute_spectrum, ignore=c("x"))
    decon_list <- mcmapply(nw, deconvolute_spectrum, x, MoreArgs = args)
    decons <- as_decons2(decon_list)
    duration <- format(round(Sys.time() - starttime, 3))
    logf("Finished deconvolution %s", duration)
    decons
}

#' @noRd
#' @author 2024-2026 Tobias Schmidt: initial version.
#' @examples
#' s <- deconvolute_spectrum(sap, smit=1, smws=3, delta=3, sfr=c(3.2,-3.2))
#' x <- read_spectrum(metabodecon_file("urine_1"))
#' u <- grid_deconvolute_spectrum(x, use_rust=TRUE)
#' d <- deconvolute_spectrum(x, npmax=1000)
deconvolute_spectrum <- function(
    x, nfit=3, smit=2, smws=5, delta=6.4, npmax=0,
    sfr=NULL, igrs=list(),
    use_rust=FALSE, verbose=TRUE, # end of public args
    force=FALSE, full=TRUE
) {

    # Init locals
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())
    sfr <- sfr %||% quantile(x$cs, c(0.9, 0.1))
    name <- get_name(x)
    backend <- if (use_rust >= 1) "Rust" else "R"
    suffix <- sprintf(" using %s backend", backend)

    # Pick best params from attached grid (when npmax > 0)
    if (npmax >= 1) {
        best <- pick_best_params(x, npmax, name)
        nfit <- best$nfit; smit <- best$smit
        smws <- best$smws; delta <- best$delta
    }
    args <- get_args(deconvolute_spectrum, ignore=c("x", "verbose", "use_rust", "force"))

    # Deconvolute with given/optimal parameters
    logf("Starting deconvolution of %s%s", name, suffix)
    decon <- if (use_rust >= 1) {
        deconvolute_spectrum_rust(x, args, sfr, igrs, nfit, smit, smws, delta, full)
    } else {
        deconvolute_spectrum_r(x, args, sfr, igrs, nfit, smit, smws, delta, force, full)
    }

    # Cache and return
    logf("Finished deconvolution of %s", name)
    decon
}

#' @noRd
#' @title Pick best params from attached `deg` grid for `npmax`
#' @description
#' Selects the parameter row with smallest area ratio (`ar`) among rows with
#' `np > 0` and `np < npmax`. Falls back to the row with the smallest `np`
#' when no row satisfies `np < npmax`.
#' @return A length-1 list with elements `nfit`, `smit`, `smws`, `delta`.
#' @author 2024-2026 Tobias Schmidt: initial version.
pick_best_params <- function(x, npmax, name) {
    if (is.null(x$deg)) stop(
        "deconvolute_spectrum() requires x$deg when npmax >= 1. ",
        "Call grid_deconvolute_spectra() first."
    )
    G <- x$deg
    G <- G[G$np > 0, ]
    # The Rust backend sometimes produces zero peaks if SFR and Delta are
    # too small. We ignore these and select only from the valid results.
    if (nrow(G) == 0) stop("All parameter sets produced zero peaks.")
    sub <- G[G$np < npmax, ]
    if (nrow(sub) == 0) {
        ln1 <- "No parameter set found with np < %d."
        ln2 <- "Using set with smallest np (%d) instead."
        logf(paste(ln1, ln2), npmax, min(G$np))
        best <- G[G$np == min(G$np), , drop=FALSE]
    } else {
        best <- sub[sub$ar == min(sub$ar), , drop=FALSE]
    }
    fmt <- "Best params for %s: nfit=%d, smit=%d, smws=%d, delta=%.2f"
    logf(fmt, name, best$nfit[1], best$smit[1], best$smws[1], best$delta[1])
    list(
        nfit=best$nfit[1], smit=best$smit[1],
        smws=best$smws[1], delta=best$delta[1]
    )
}

#' @noRd
#' @title Build a `decon2` object using the R deconvolution backend.
#' @author 2024-2026 Tobias Schmidt: initial version.
deconvolute_spectrum_r <- function(
    x, args, sfr, igrs, nfit, smit, smws, delta, force, full=TRUE
) {
    cs <- x$cs
    si <- x$si
    sfr_igr <- list(c(Inf, max(sfr)), c(min(sfr), -Inf))
    igrs <- c(sfr_igr, igrs)
    sm <- smooth_signals2(si, smit, smws)
    peaks <- find_peaks2(sm)
    peaks <- filter_peaks2(peaks, cs, sfr, delta, force, igrs)
    lcpar <- fit_lorentz_curves2(cs, si, peaks, nfit)
    sit <- if (full) {
        data.frame(sm=sm, sup=lorentz_sup(cs, lcpar=lcpar))
    } else {
        data.frame(sm=sm)
    }
    decon <- list(cs=cs, si=si, meta=x$meta, args=args,
                  sit=sit, peak=peaks, lcpar=lcpar)
    class(decon) <- c("decon2", "spectrum")
    decon
}

#' @noRd
#' @title Build a `decon2` object using the Rust deconvolution backend.
#' @author 2024-2025 Tobias Schmidt: initial version.
deconvolute_spectrum_rust <- function(
    x, args, sfr, igrs, nfit, smit, smws, delta, full=TRUE
) {
    mdrb_spectrum <- mdrb::Spectrum$new(x$cs, x$si, sfr)
    mdrb_deconvr <- mdrb::Deconvoluter$new()
    mdrb_deconvr$set_moving_average_smoother(smit, smws)
    mdrb_deconvr$set_noise_score_selector(delta)
    mdrb_deconvr$set_analytical_fitter(nfit)
    for (r in igrs) mdrb_deconvr$add_ignore_region(r[1], r[2])
    mdrb_decon <- mdrb_deconvr$deconvolute_spectrum(mdrb_spectrum)
    cs <- mdrb_spectrum$chemical_shifts()
    si <- mdrb_spectrum$intensities()
    lcpar <- as.data.frame(mdrb_decon$lorentzians())[, c("x0", "A", "lambda")]
    sm <- smooth_signals2(si, smit, smws)        # Rust does not return sm
    sit <- if (full) {
        data.frame(sm=sm, sup=mdrb_decon$superposition_vec(cs))
    } else {
        data.frame(sm=sm)
    }
    peak <- get_peak(lcpar$x0, cs)
    decon <- list(
        cs=cs, si=si, meta=x$meta, args=args, sit=sit, peak=peak,
        lcpar=lcpar
    )
    class(decon) <- c("decon2", "spectrum")
    decon
}

#' @noRd
#'
#' @title Attach grid-search performance tables to spectra
#'
#' @description
#' For each spectrum in `x`, runs `grid_deconvolute_spectrum()` and attaches
#' the resulting performance grid as element `$deg`. Idempotent: spectra
#' that already carry a `$deg` element are left untouched. The enriched
#' `spectra` object is returned so workers in `deconvolute_spectra()` can
#' look up best parameters locally.
#'
#' @param deg
#' Deconvolution-parameter grid: a data frame with columns `nfit`, `smit`,
#' `smws`, `delta`. The default is the 60-cell cartesian product
#' `expand.grid(nfit=10, smit=1:3, smws=c(3,5,7,9), delta=(1:5)*1.6)`.
#' May also be a model-fitting grid (`mog`) — i.e. a data frame that
#' additionally has an `npmax` column — in which case only unique
#' `(nfit, smit, smws, delta)` rows with `npmax > 0` are used.
grid_deconvolute_spectra <- function(
    x, deg=expand.grid(nfit=10, smit=1:3, smws=c(3,5,7,9), delta=(1:5)*1.6),
    sfr=NULL, igrs=list(), verbose=TRUE, nworkers=1, use_rust=FALSE
) {
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())
    if (!is.null(deg) && "npmax" %in% names(deg)) {
        cols <- c("nfit", "smit", "smws", "delta")
        deg <- unique(deg[deg$npmax > 0, cols, drop=FALSE])
        rownames(deg) <- NULL
    }

    seen <- sapply(x, function(s) !is.null(s$deg))
    logf("Reusing grids for %d/%d spectra", sum(seen), length(x))
    if (all(seen)) return(invisible(x))

    logf("Running grid deconvolution for remaining %d spectra", sum(!seen))
    idx <- which(!seen)
    nw <- min(nworkers, length(idx))
    enriched <- mcmapply(
        nw, grid_deconvolute_spectrum, x=x[idx],
        MoreArgs = list(deg=deg, sfr=sfr, igrs=igrs, verbose=verbose, use_rust=use_rust)
    )
    for (k in seq_along(idx)) x[[idx[k]]] <- enriched[[k]]
    invisible(x)
}

#' @noRd
#'
#' @title Deconvolute one spectrum using a grid of parameters
#'
#' @inheritParams deconvolute_spectrum
#'
#' @param deg
#' Deconvolution-parameter grid. See `grid_deconvolute_spectra()` for details.
#'
#' @return
#' The input spectrum with a `$deg` element attached: `deg` augmented with
#' columns `ar` (residual-area / spectra-area) and `np` (number of peaks in
#' the deconvolution result). Idempotent: if `x$deg` is already set, `x` is
#' returned unchanged.
#'
#' @examples
#'
#' x <- read_spectrum(metabodecon_file("urine_1"))
#' xR <- grid_deconvolute_spectrum(x, use_rust=FALSE)
#' xRust <- grid_deconvolute_spectrum(x, use_rust=TRUE)
#'
grid_deconvolute_spectrum <- function(
    x,
    deg=NULL,
    sfr=NULL, igrs=list(), verbose=TRUE, use_rust=FALSE
) {
    if (is.null(deg)) deg <- expand.grid2(
        nfit=10, smit=1:3, smws=c(3,5,7,9), delta=(1:5)*1.6
    )
    if (!is.null(x$deg)) return(x)
    if (!verbose) local_options(toscutil.logf.file = nullfile())

    cols <- c("smit", "smws", "delta", "nfit")
    stopifnot(all(cols %in% colnames(deg)))
    deg$ar <- NA_real_; deg$np <- NA_integer_
    truepar <- x$meta$simpar
    if (!is.null(truepar)) deg$prarpx <- NA_real_
    ds_args <- modifyList(
        as.list(formals(deconvolute_spectrum)),
        list(x=x, sfr=sfr, igrs=igrs, verbose=FALSE, use_rust=use_rust, npmax=0)
    )

    logf("Grid deconvoluting %s", specname <- get_name(x))
    for (i in seq_len(nrow(deg))) {
        logf("Grid search iteration %d/%d", i, nrow(deg))
        ds_args[cols] <- deg[i, cols]
        d <- do.call(deconvolute_spectrum, ds_args)
        deg[i, "ar"] <- sum(abs(d$sit$sup - d$si)) / sum(abs(d$si))
        deg[i, "np"] <- nrow(d$lcpar)
        if (!is.null(truepar)) {
            deg[i, "prarpx"] <- calc_prarp(d, truepar=truepar)$prarpx
        }
    }
    logf("Finished grid deconvolution of %s", specname)

    x$deg <- deg
    x
}

# Helpers for deconvolute_spectrum #####

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
#' The input `spec` list with two additional entries:
#' - `spec$Z`: A data frame containing the intermediate smoothed values after each iteration.
#' - `spec$y_smooth`: A numeric vector of the smoothed values after
#'
#' @details
#' Applies a centered moving average. Boundary values are filled with
#' the mean of available neighbors to maintain vector length.
smooth_signals2 <- function(y, reps = 2, k = 5) {
    if (k %% 2 == 0) stop("k must be odd")
    n <- length(y)
    for (i in seq_len(reps)) {
        filter <- rep(1 / k, k)
        z <- stats::filter(y, filter, sides = 2) # (1)
        q <- (k - 1) / 2 # (2)
        for (j in seq_len(q)) {
            z[j] <- mean(y[1:(q + j)]) # (3)
            z[n - j + 1] <- mean(y[(n - q - j + 1):n]) # (4)
        }
        y <- as.numeric(z)
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
    y
}

find_peaks2 <- function(y) {
    logf("Starting peak selection")
    P <- .Call(find_peaks_c, as.double(y))
    logf("Detected %d peaks", nrow(P))
    P
}

#' @noRd
#' @title Filter Peaks with Low Scores
#' @description
#' Filters peaks by score. Peaks whose center ppm falls outside
#' `sfr` are used to estimate noise; signal-region peaks with scores below `mean
#' + delta * sd` are removed. Peaks inside any `igr` region are also removed.
#' @param peaks Data frame with columns `left`, `center`, `right`, `score`.
#' @param cs Chemical shifts (ppm), same length as the spectrum.
#' @param sfr Length-2 numeric: signal-free region boundaries in ppm.
#' @param delta Filtering threshold in standard deviations.
#' @param force If TRUE, proceed even if no SFR peaks are found.
#' @param igr List of length-2 numeric vectors (ignore regions in ppm).
#' @return Filtered data frame (rows with `high == TRUE` only).
#' @author 2026 Tobias Schmidt: initial version.
filter_peaks2 <- function(peaks, cs, sfr, delta = 6.4, force = FALSE,
                          igr = list()) {
    logf("Removing peaks with low scores")
    ppm_ct <- cs[peaks$center]
    in_sfr <- ppm_ct >= max(sfr) | ppm_ct <= min(sfr)
    # Compute noise statistics from SFR peaks before any removal
    if (sum(in_sfr) > 1) {
        mu <- mean(peaks$score[in_sfr])
        sigma <- sd(peaks$score[in_sfr])
    } else {
        if (!force) stop(
            "Not enough signals found in signal free region. ",
            "Please double check deconvolution parameters."
        )
        mu <- 0; sigma <- 0
    }
    high <- peaks$score > mu + delta * sigma
    if (length(igr) > 0) {
        in_igr <- vapply(ppm_ct, function(c) {
            any(vapply(igr, function(r) c >= min(r) & c <= max(r), logical(1)))
        }, logical(1))
        high <- high & !in_igr
    }
    out <- peaks[high, ]
    logf("Removed %d peaks", nrow(peaks) - nrow(out))
    out
}

#' @noRd
#' @title Fit Lorentz Curves (v2)
#' @description
#' Modular replacement for [metabodecon::fit_lorentz_curves()]. Works directly
#' in ppm and uses the same algorithm as the Rust backend: 3-point peak stencil
#' with iterative refinement. Returns a data frame with columns `x0`, `A`,
#' `lambda`.
#' @param cs Chemical shifts in ppm.
#' @param si Signal intensities (raw, unsmoothed).
#' @param peaks Data frame with columns `left`, `center`, `right` (indices).
#' @param nfit Number of refinement iterations.
#' @return Data frame with columns `x0` (ppm), `A`, `lambda` (ppm).
#' @author 2026 Tobias Schmidt: initial version.
fit_lorentz_curves2 <- function(cs, si, peaks, nfit = 3) {
    logf("Fitting Lorentz curves (%d iterations)", nfit)
    il <- peaks$left; ic <- peaks$center; ir <- peaks$right
    np <- length(ic)

    # Build 3-point stencils (x1=left, x2=center, x3=right)
    x1 <- cs[il]; x2 <- cs[ic]; x3 <- cs[ir]
    y1 <- si[il]; y2 <- si[ic]; y3 <- si[ir]

    # Mirror shoulders
    mirror <- function(x1, x2, x3, y1, y2, y3) {
        inc <- y1 <= y2 & y2 <= y3  # ascending
        dec <- y1 >= y2 & y2 >= y3  # descending
        x3[inc] <- 2 * x2[inc] - x1[inc]
        y3[inc] <- y1[inc]
        x1[dec] <- 2 * x2[dec] - x3[dec]
        y1[dec] <- y3[dec]
        list(x1 = x1, x3 = x3, y1 = y1, y3 = y3)
    }
    m <- mirror(x1, x2, x3, y1, y2, y3)
    x1 <- m$x1; x3 <- m$x3; y1 <- m$y1; y3 <- m$y3

    # Solve 3-equation system
    solve_params <- function(x1, x2, x3, y1, y2, y3) {
        # maximum_position (x0)
        num <- x1^2 * y1 * (y2 - y3) + x2^2 * y2 * (y3 - y1) + x3^2 * y3 * (y1 - y2)
        den <- 2 * ((x1 - x2) * y1 * y2 + (x2 - x3) * y2 * y3 + (x3 - x1) * y3 * y1)
        maxp <- num / den
        maxp[!is.finite(maxp)] <- 0
        # half_width2 (hw2 = lambda^2)
        left <- (y1 * (x1 - maxp)^2 - y2 * (x2 - maxp)^2) / (y2 - y1)
        right <- (y2 * (x2 - maxp)^2 - y3 * (x3 - maxp)^2) / (y3 - y2)
        hw2 <- pmax((left + right) / 2, .Machine$double.eps)
        hw2[!is.finite(hw2)] <- .Machine$double.eps
        # scale_factor_half_width (sfhw = A * lambda)
        sfhw <- y2 * (hw2 + (x2 - maxp)^2)
        sfhw[!is.finite(sfhw)] <- 0
        list(sfhw = sfhw, hw2 = hw2, maxp = maxp)
    }
    p <- solve_params(x1, x2, x3, y1, y2, y3)

    # Build reduced spectrum (3 points per peak, flattened, ORIGINAL values)
    rs_x <- as.numeric(rbind(cs[il], cs[ic], cs[ir]))
    rs_y <- as.numeric(rbind(si[il], si[ic], si[ir]))

    # Iterative refinement
    for (iter in seq_len(nfit)) {
        # Superposition at each reduced spectrum point
        sup <- lorentz_sup(rs_x, x0 = p$maxp, Al = p$sfhw, l2 = p$hw2)
        # Ratio: original / superposition
        ratio <- rs_y / sup
        ratio[!is.finite(ratio)] <- 1
        # Update stencil intensities
        rm <- matrix(ratio, nrow = 3)
        y1 <- y1 * rm[1, ]; y2 <- y2 * rm[2, ]; y3 <- y3 * rm[3, ]
        # Re-mirror shoulders
        m <- mirror(x1, x2, x3, y1, y2, y3)
        x1 <- m$x1; x3 <- m$x3; y1 <- m$y1; y3 <- m$y3
        # Re-solve
        p <- solve_params(x1, x2, x3, y1, y2, y3)
    }

    # Filter degenerate peaks (match Rust CHECK_PRECISION = 1e6 * eps)
    eps <- 1e6 * .Machine$double.eps
    ok <- p$sfhw > eps & p$hw2 > eps
    sfhw <- p$sfhw[ok]
    hw2 <- p$hw2[ok]
    maxp <- p$maxp[ok]

    # Convert (sfhw, hw2, maxp) → (x0, A, lambda) for decon2 format
    lambda <- sqrt(hw2)
    A <- sfhw / lambda
    data.frame(x0 = maxp, A = A, lambda = lambda)
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
#' 1. The argument names are based on the names used by Koh et al. (2009).
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
#' @author
#' 2024-2026 Tobias Schmidt: initial versions (v1-v3); 2026 compiled-C version.
#' @param x Positions at which to evaluate (numeric vector).
#' @param x0 Peak centres (numeric vector, length np).
#' @param A Peak amplitudes (length np). Only needed when Al is NULL.
#' @param lambda Peak half-widths (length np). Only needed when l2 is NULL.
#' @param lcpar Optional named list/data frame with fields A, lambda, x0/x_0/w.
#' @param Al Pre-computed |A * lambda| (length np). Avoids recomputation.
#' @param l2 Pre-computed lambda^2 (length np). Avoids recomputation.
lorentz_sup <- function(
    x, x0, A = NULL, lambda = NULL, lcpar = NULL, Al = NULL, l2 = NULL
) {
    if (!is.null(lcpar)) {
        nams <- names(lcpar)
        if ("A" %in% nams) A <- lcpar[["A"]]
        if ("lambda" %in% nams) lambda <- lcpar[["lambda"]]
        if ("x_0" %in% nams) x0 <- lcpar[["x_0"]]
        if ("x0" %in% nams) x0 <- lcpar[["x0"]]
        if ("w" %in% nams) x0 <- lcpar[["w"]]
    }
    if (is.null(Al)) Al <- abs(A * lambda)
    if (is.null(l2)) l2 <- lambda^2
    if (length(x0) != length(Al) || length(x0) != length(l2)) {
        stop("x0, Al, and l2 must have the same length.", call. = FALSE)
    }
    .Call(
        lorentz_sup_c, as.double(x), as.double(x0), as.double(Al), as.double(l2)
    )
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
#'
#' @title Calculate the PRARP Score
#'
#' @description
#' Calculates the PRARP score for a deconvolution. The PRARP score is the
#' product of the peak ratio and the area ratio and can be used to assess the
#' quality of a deconvolution. See 'Details' for more information on how the
#' score is calculated.
#'
#' @param decon A list containing the deconvolution results, as returned by
#' [metabodecon::generate_lorentz_curves()].
#'
#' @param lcpar A data frame containing the true parameters of the peaks.
#'
#' @return The PRARP score as numeric scalar. In addition, a plot is created to
#' visualize the deconvolution results.
#'
#' @details
#' The PRARP score is calculated as follows:
#'
#' peak_ratio = min(peaks_true, peaks_found) / max(peaks_true, peaks_found)
#' area_ratio = min(area_true,  area_found)  / max(area_true,  area_found)
#' prarp      = peak_ratio * area_ratio
#'
#' I.e., the score is close to 1 if the number of peaks and the area of the
#' peaks are similar in the true and found spectra and the score is close to 0
#' if the number of peaks and/or the area of the peaks are very different.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' ## Bad deconvolution (PRARP ~= 0.2)
#' decon <- generate_lorentz_curves_sim(sim[[1]], delta = 6.4)
#' truepar <- sim[[1]]$meta$simpar[c("A", "x0", "lambda")]
#' calc_prarp(decon, truepar)
#' plot_prarp(decon, truepar)
#'
#' ## Good deconvolution (PRARP ~= 0.64)
#' decon <- generate_lorentz_curves_sim(sim[[1]], delta = 0)
#' truepar <- sim[[1]]$meta$simpar[c("A", "x0", "lambda")]
#' calc_prarp(decon, truepar)
#' plot_prarp(decon, truepar)
#'
calc_prarp <- function(x, truepar = NULL, ...) {

    obj <- as_decon2(x)
    truepar <- truepar %||% obj$meta$simpar

    x0_true <- truepar$x0
    x0_found <- obj$lcpar$x0
    idx_closest_true_peak <- sapply(x0_found, function(x0) which.min(abs(x0_true - x0)))

    np_true <- length(truepar$x0)
    np_found <- length(idx_closest_true_peak)
    np_correct <- length(unique(idx_closest_true_peak))
    np_wrong <- np_found - np_correct
    peak_ratio   <- min(np_found, np_true) / max(np_found, np_true)
    peak_ratio_x <- np_correct / (np_true + np_wrong)

    area_spectrum <- sum(abs(obj$si))
    area_residuals <- sum(abs(obj$sit$sup - obj$si))
    area_ratio <- area_residuals / area_spectrum

    prarp <- peak_ratio * (1 - area_ratio)
    prarpx <- peak_ratio_x * (1 - area_ratio)

    named(
        prarpx, prarp, peak_ratio_x, peak_ratio,
        np_true, np_found, np_correct, np_wrong,
        area_ratio, area_spectrum, area_residuals
    )
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
calc_prarpx <- function(x, truepar = NULL, ...) {
    calc_prarp(x, truepar)$prarpx
}
