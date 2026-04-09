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
#' @param npmax Integer. Maximum number of peaks allowed in the result. If
#' `npmax >= 1`, the `nfit`, `smopts` and `delta` arguments are ignored and a
#' grid search over predefined parameter combinations is performed instead. The
#' combination with the smallest residual area ratio and fewer than `npmax`
#' peaks is selected. Grid search results are cached to disk by default; see
#' `cadir`.
#'
#' @param igrs Ignore regions. List of length-2 numeric vectors specifying the
#' start and endpoints of the chemical shift regions to ignore during
#' deconvolution. Peaks whose centers fall inside any ignore region are
#' excluded from fitting.
#'
#' @param cadir Directory for caching grid search results and deconvolution
#' results for `npmax >= 1`. Defaults to [metabodecon::decon_cachedir()]. If
#' `NULL`, caching is disabled. Pass a custom path to use a different cache
#' directory.
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
#' spectrum <- sim[[1]]
#' decon <- deconvolute(spectrum)
#'
#' ## Read multiple spectra from disk and deconvolute at once
#' spectra_dir <- metabodecon_file("sim_subset")
#' spectra <- read_spectra(spectra_dir)
#' decons <- deconvolute(spectra, sfr = c(3.55, 3.35))
deconvolute <- function(x,
    nfit=3,    smopts=c(2, 5), delta=6.4,    sfr=NULL,      wshw=0,
    ask=FALSE, force=FALSE,    verbose=TRUE, nworkers=1,    use_rust=FALSE,
    npmax=0,   igrs=list(),    cadir=decon_cachedir()
) {

    # Check inputs
    stopifnot(
        is_spectrum_or_spectra(x),    is_int_or_null(nfit, 1),
        is_int_or_null(smopts, 2),    is_num_or_null(delta, 1),
        is_num_or_null(sfr, 2),       is_num_or_null(wshw, 1),
        is_bool(ask, 1),              is_bool(force, 1),
        is_bool(verbose, 1),          is_int(nworkers, 1),
        is_bool_or_null(use_rust, 1), is_int(npmax, 1),
        is_list_of_nums(igrs, nv=2),  is_str_or_null(cadir)
    )

    # Set suitable defaults
    sfr <- sfr %||% quantile(x$cs %||% x[[1]]$cs, c(0.9, 0.1))
    if (isTRUE(use_rust)) check_mdrb(stop_on_fail = TRUE)
    if (is.null(use_rust)) use_rust <- check_mdrb()

    # Perform deconvolution
    decons2 <- deconvolute_spectra(x,
        nfit, smopts, delta, sfr, wshw,
        ask, force, verbose, bwc=2,
        use_rust, nw=nworkers, igr=igrs, rtyp="decon2",
        npmax=npmax, cadir=cadir
    )

    # Convert and return
    if (length(decons2) == 1) decons2[[1]] else decons2
}

# Internal main functions #####

#' @noRd
#' @inheritParams deconvolute
#' @param bwc Whether to produce results backwards compatible with
#' [metabodecon::MetaboDecon1D()]. If `bwc == 0`, bug fixes introduced after version 0.2.2 of
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
    rtyp="idecon", npmax=0,           cadir=decon_cachedir()
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
        is_int(nw, 1), is.list(igr),
        is_char(rtyp, 1, "(decon[0-2]|idecon|rdecon)"),
        is_int(npmax, 1), is_str_or_null(cadir)
    )

    # Configure logging
    if (!verbose) local_options(toscutil.logf.file = nullfile())

    # Init locals
    spectra <- as_spectra(x)
    ns <- length(spectra)
    nc2 <- ceiling(detectCores() / 2)
    nw_apply <- if (nw == "auto") min(nc2, ns) else nw
    nw_apply_str <- if (nw_apply == 1) "1 worker" else sprintf("%d workers", nw_apply)
    nw_deconv <- 1
    ns_str <- if (ns == 1) "1 spectrum" else sprintf("%d spectra", ns)
    adjno <- get_adjno(spectra, ask)
    sfr_list <- get_sfr(spectra, sfr, ask, adjno)
    wshw_list <- get_wshw(spectra, wshw, ask, adjno)
    smopts_list <- get_smopts(spectra, smopts)
    igr_list <- list(igr)
    cadir_list <- list(cadir)

    # Deconvolute spectra
    logf("Starting deconvolution of %s using %s", ns_str, nw_apply_str)
    starttime <- Sys.time()
    decon_list <- mcmapply(nw_apply, deconvolute_spectrum,
        spectra,
        nfit, smopts_list, delta, sfr_list, wshw_list,
        ask, force, verbose, bwc,
        use_rust, nw_deconv, igr_list, rtyp,
        npmax, cadir_list
    )
    decons <- as_collection(decon_list, rtyp)
    duration <- format(round(Sys.time() - starttime, 3))
    logf("Finished deconvolution of %s in %s", ns_str, duration)

    # Return
    decons
}

#' @noRd
#' @inheritParams deconvolute_spectra
#' @author 2024-2025 Tobias Schmidt: initial version.
#' @examples
#' x <- sap[[1]];
#' nfit <- 3; smopts <- c(1,3); delta <- 3; sfr <- c(3.2,-3.2); wshw <- 0;
#' ask <- FALSE; force <- FALSE; verbose <- FALSE; bwc <- 2;
#' use_rust <- FALSE; nw <- 1; igr <- list(); rtyp <- "idecon"
#' idecon <- deconvolute_spectrum(
#'      x, nfit, smopts, delta, sfr, wshw,
#'      ask, force, verbose, bwc,
#'      use_rust, nw, igr, rtyp
#' )
#'
#' x <- read_spectrum(metabodecon_file("urine_1"))
#' nfit <- 3; smopts <- c(2,5); delta <- 6;
#' sfr <- quantile(x$cs, c(0.9, 0.1)); wshw <- 0;
#' ask <- FALSE; force <- FALSE; verbose <- TRUE; bwc <- 2;
#' use_rust <- TRUE; nw <- 1; igr <- list(); rtyp <- "decon2"
#' npmax=1000
#' decon1 <- deconvolute_spectrum(
#'      x, nfit, smopts, delta, sfr, wshw,
#'      ask, force, verbose, bwc, use_rust, nw, igr, rtyp,
#'      npmax
#' )
deconvolute_spectrum <- function(x,
    nfit=3, smopts=c(2, 5), delta=6.4, sfr=c(3.55, 3.35), wshw=0,
    ask=FALSE, force=FALSE, verbose=TRUE, bwc=2,
    use_rust=FALSE, nw=1, igr=list(), rtyp="idecon",
    npmax=0, cadir=decon_cachedir()
) {

    # Check inputs
    assert(
        is_spectrum(x),
        is_int(nfit, 1),  is_int(smopts, 2),     is_num(delta, 1),
        is_num(sfr, 2),   is_num(wshw, 1),       is_bool(force, 1),
        is_num(bwc, 1),   is_bool(use_rust, 1),  is_int(nw, 1),
        is_list_of_nums(igr, nv=2),
        is_char(rtyp, 1, "(decon[0-2]|idecon|rdecon)"),
        is_int(npmax, 1), is_str_or_null(cadir)
    )

    # Init locals
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())
    name <- get_name(x)
    suffix <- if (use_rust) " using Rust backend" else ""

    # Stop early if deconvolution is cached
    if (npmax >= 1) {
        cache <- disk_cache(cadir)
        args <- list("ds", x, sfr, wshw, force, bwc, use_rust, npmax, rtyp)
        cachehash <- rlang::hash(args)
        if (cache$exists(cachehash)) {
            logf("Cache hit for %s (npmax=%d)", name, npmax)
            return(cache$get(cachehash))
        }
    }

    # Perform grid search (to replace given nfit, smopts and delta)
    if (npmax >= 1) {
        fmt <- "Starting grid deconvolution of %s using %s backend"
        logf(fmt, name, if (use_rust) "Rust" else "R")
        G <- grid_deconvolute_spectrum(x, sfr, verbose, use_rust, cadir)
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
        nfit <- best$nfit[1]
        smopts <- c(best$smit[1], best$smws[1])
        delta <- best$delta[1]
        fmt <- "Done. Best params: nfit=%d, smopts=c(%d, %d), delta=%.2f"
        logf(fmt, nfit, smopts[1], smopts[2], delta)
    }
    args <- get_args(deconvolute_spectrum, ignore = "x")

    # Deconvolute with given/optimal parameters
    logf("Starting deconvolution of %s%s", name, suffix)
    if (use_rust) {
        mdrb_spectrum <- mdrb::Spectrum$new(x$cs, x$si, sfr)
        mdrb_deconvr <- mdrb::Deconvoluter$new()
        mdrb_deconvr$set_moving_average_smoother(smopts[1], smopts[2])
        mdrb_deconvr$set_noise_score_selector(delta)
        mdrb_deconvr$set_analytical_fitter(nfit)
        for (r in igr) mdrb_deconvr$add_ignore_region(r[1], r[2])
        mdrb_decon <- if (nw > 1) {
            mdrb_deconvr$set_threads(nw)
            mdrb_deconvr$par_deconvolute_spectrum(mdrb_spectrum)
        } else {
            mdrb_deconvr$deconvolute_spectrum(mdrb_spectrum)
        }
        decon <- new_rdecon(x, args, mdrb_spectrum, mdrb_deconvr, mdrb_decon)
    } else {
        ispec <- as_ispec(x)
        ispec <- set(ispec, args=args)
        ispec <- rm_water_signal(ispec, wshw, bwc)
        ispec <- rm_negative_signals(ispec)
        ispec <- smooth_signals(ispec, smopts[1], smopts[2], bwc)
        ispec <- find_peaks(ispec)
        ispec <- filter_peaks(ispec, sfr, delta, force, bwc, igr)
        ispec <- fit_lorentz_curves(ispec, nfit, bwc)
        decon <- as_idecon(ispec)
    }

    # Format, cache and return
    logf("Formatting return object as %s", rtyp)
    convert <- switch(rtyp,
        "decon0"=as_decon0, "decon1"=as_decon1, "decon2"=as_decon2,
        "idecon"=as_idecon, "rdecon"=as_rdecon
    )
    decon <- convert(decon)
    if (npmax >= 1) cache$set(cachehash, decon)
    logf("Finished deconvolution of %s", name)
    decon
}

#' @noRd
#' @description
#' WORK IN PROGRESS.
#'
#' Planned replacement for `deconvolute_spectrum()`. Should directly produce a
#' decon2 object, so we can remove the `ispec` and `idecon` classes.
#'
#' @inheritParams deconvolute_spectra
#'
#' @examples
#'
#' x <- sap[[1]];
#' nfit <- 3; smopts <- c(1,3); delta <- 3; sfr <- c(3.2,-3.2); wshw <- 0;
#' ask <- FALSE; force <- FALSE; verbose <- FALSE; bwc <- 2;
#' use_rust <- FALSE; nw <- 1; igr <- list(); rtyp <- "idecon"
#'
#' x <- read_spectrum(metabodecon_file("urine_1"))
#' nfit <- 3; smopts <- c(2,5); delta <- 6.4; sfr <- c(3.55,3.35); wshw <- 0;
#' ask <- FALSE; force <- FALSE; verbose <- FALSE; bwc <- 2;
#' use_rust <- FALSE; nw <- 1; igr <- list(); rtyp <- "idecon"
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
deconvolute_spectrum2 <- function(x,
    nfit=3, smopts=c(2, 5), delta=6.4, sfr=NULL, wshw=0,
    ask=FALSE, force=FALSE, verbose=TRUE, bwc=2,
    use_rust=FALSE, nw=1, igr=NULL, rtyp="decon2",
    abs_neg=TRUE
) {
    # Check inputs
    sfr <- sfr %||% quantile(x$cs, c(0.9, 0.1))
    assert(
        is_spectrum(x),
        is_int(nfit, 1),  is_int(smopts, 2),     is_num(delta, 1),
        is_num(sfr, 2),   is_num(wshw, 1),       is_bool(force, 1),
        is_num(bwc, 1),   is_bool(use_rust, 1),  is_int(nw, 1),
        is.null(igr) || is_list_of_nums(igr, nv=2),
        is_char(rtyp, 1, "(decon[0-2]|idecon|rdecon)"),
        is_bool(abs_neg, 1)
    )
    # NULL igr = remove SFR peaks (Rust does this natively via
    # signal_boundaries; for R we add SFR as ignore regions).
    # list() = keep SFR peaks (old R behavior).
    rm_sfr <- is.null(igr)
    if (is.null(igr)) igr <- list()

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
        for (r in igr) mdrb_deconvr$add_ignore_region(r[1], r[2])
        mdrb_decon <- if (nw > 1) {
            mdrb_deconvr$set_threads(nw)
            mdrb_deconvr$par_deconvolute_spectrum(mdrb_spectrum)
        } else {
            mdrb_deconvr$deconvolute_spectrum(mdrb_spectrum)
        }
        decon <- new_rdecon(x, args, mdrb_spectrum, mdrb_deconvr, mdrb_decon)
    } else {
        cs <- x$cs; si <- x$si
        r_igr <- igr
        if (rm_sfr) {
            sfr_igr <- list(c(Inf, max(sfr)), c(min(sfr), -Inf))
            r_igr <- c(sfr_igr, r_igr)
        }
        wsrm <- rm_water_signal2(cs, si, wshw)
        nvrm <- if (abs_neg) rm_negative_signals2(wsrm) else wsrm
        sm <- smooth_signals2(nvrm, smopts[1], smopts[2])
        peaks <- find_peaks2(sm)
        peaks <- filter_peaks2(peaks, cs, sfr, delta, force, r_igr)
        lcpar <- fit_lorentz_curves2(cs, si, peaks, nfit)
        sup <- lorentz_sup(cs, lcpar = lcpar)
        sit <- data.frame(
            wsrm = wsrm * 1e6,
            nvrm = nvrm * 1e6,
            sm = sm * 1e6,
            sup = sup
        )
        mse_list <- list(
            raw = mse(si, sup, normed = FALSE),
            norm = mse(si, sup, normed = TRUE),
            sm = mse(sit$sm, sup, normed = FALSE),
            smnorm = mse(sit$sm, sup, normed = TRUE)
        )
        decon <- list(cs = cs, si = si, meta = x$meta, args = args,
            sit = sit, peak = peaks[, c("left", "center", "right")],
            lcpar = lcpar, mse = mse_list)
        class(decon) <- "decon2"
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

grid_deconvolute_spectra <- function(
    x,
    sfr=NULL,
    verbose=TRUE,
    nw=1,
    use_rust=FALSE,
    cadir = decon_cachedir() # (1)
    # (1) We share the cache between `grid_deconvolute_spectra()` and
    # `grid_deconvolute_spectrum()`, so we can easily check whether the grid
    # search for an individual spectrum has already been performed. That means
    # we only have to perform the grid search with multiple workers for the
    # spectra that haven't been seen yet.
) {

    assert(
        is_spectra(x),
        is_num_or_null(sfr, 2)  || is_list_of_nums(sfr, length(x), 2),
        is_bool(verbose, 1),
        is_int(nw, 1),
        is_bool(use_rust),
        is_str_or_null(cadir)
    )
    sfr <- sfr %||% quantile(x[[1]]$cs, c(0.9, 0.1))
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())

    cache <- disk_cache(cadir)
    hashes <- lapply(x, function(s) rlang::hash(list("gds", s, sfr, use_rust)))
    gridlist <- vector("list", length(x))
    logf("Checking for cached results")
    seen <- sapply(hashes, cache$exists)
    logf("Reading %d/%d results from cache", sum(seen), length(x))
    gridlist[seen] <- lapply(hashes[seen], function(h) cache$get(h))
    if (sum(seen) == length(x)) return(invisible(gridlist))

    logf("Running grid deconvolution for remaining %d spectra", sum(!seen))
    u <- x[!seen] # u == unseen spectra
    nw <- min(nw, length(u))
    gridlist[!seen] <- mcmapply(nw, grid_deconvolute_spectrum,
        x=u, sfr=list(sfr), verbose=verbose, use_rust=use_rust,
        cadir=list(cadir)
    )
    # Note: grid_deconvolute_spectrum() caches its own result, so we
    # don't need to write to the cache here.
    invisible(gridlist)
}

#' @noRd
#'
#' @title Deconvolute one spectrum using a grid of parameters
#'
#' @inheritParams deconvolute_spectrum
#'
#' @return
#' A data frame with columns `smit`, `smws`, `delta`, `nfit`, `ar` and `np`,
#' where `smit`, `smws`, `delta` and `nfit` are the parameters used for
#' deconvolution, `ar` is the area ratio (residual-area/spectra-area) and `np`
#' is the number of peaks in the deconvolution result.
#'
#' @examples
#'
#' x <- read_spectrum(metabodecon_file("urine_1"))
#' cd <- decon_cachedir()
#'
#' a <- Sys.time()
#' use_rust <- FALSE
#' gridR <- grid_deconvolute_spectrum(x, use_rust=use_rust, cadir=cd)
#' runtimeR <- Sys.time() - a
#'
#' b <- Sys.time()
#' use_rust <- TRUE
#' gridRust <- grid_deconvolute_spectrum(x, use_rust=use_rust, cadir=cd)
#' runtimeRust <- Sys.time() - b
#'
grid_deconvolute_spectrum <- function(
    x,
    sfr=NULL,
    verbose=TRUE,
    use_rust=FALSE,
    cadir=decon_cachedir()
) {

    assert(
        is_spectrum(x), is_num_or_null(sfr, 2), is_bool(verbose, 1),
        is_bool(use_rust, 1), is_str_or_null(cadir)
    )
    if (!verbose) local_options(toscutil.logf.file = nullfile())

    # Resolve sfr before hashing so that NULL and its equivalent numeric
    # value produce the same cache key.
    sfr <- sfr %||% quantile(x$cs %||% x[[1]]$cs, c(0.9, 0.1))

    # Read from cache if available
    hash <- rlang::hash(list("gds", x, sfr, use_rust))
    cache <- disk_cache(cadir)
    if (cache$exists(hash)) {
        logf("Cache hit for %s", get_name(x))
        return(cache$get(hash))
    }
    backend <- if (use_rust) "Rust" else "R"
    specname <- get_name(x)

    logf("Grid deconvoluting %s using %s", specname, backend)
    grid <- expand.grid(
        smit = c(2, 3),
        smws = c(3, 5, 7, 9),
        delta = 2:8,
        nfit = c(3, 4, 5)
    )
    default_args <- as.list(formals(deconvolute_spectrum))
    call_args <- list(
        x = x, sfr = sfr, verbose = FALSE,
        use_rust = use_rust, rtyp="decon2",
        npmax = 0
    )
    args <- modifyList(default_args, call_args)
    for (i in seq_len(nrow(grid))) {
        logf("Grid search iteration %d/%d", i, nrow(grid))
        row <- grid[i, ]
        args$nfit <- row[["nfit"]]
        args$smopts <- c(row[["smit"]], row[["smws"]])
        args$delta <- row[["delta"]]
        d <- do.call(deconvolute_spectrum, args)
        grid[i, "ar"] <- sum(abs(d$sit$sup - d$si)) / sum(abs(d$si))
        grid[i, "np"] <- nrow(d$lcpar)
    }
    logf("Finished grid deconvolution of %s", specname)
    cache$set(hash, grid)
    grid
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
#' @description Same as [metabodecon::get_sfr()], but for WSHW.
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
#' Convert one or more SMOPTS vectors to a list of vectors of the correct length.
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
rm_water_signal <- function(x, wshw, bwc) {
    assert(is_ispec(x), is_num(wshw, 1), is_num(bwc, 1))
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

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl:
#' Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt:
#' Extracted and refactored corresponding code from MetaboDecon1D.
#' Added code for bwc > 1.
#' Dropped support for bwc == 0.
#' @param x Chemical Shifts in parts per million (ppm)
#' @param y Signal Intensities in arbtary units (au)
#' @param wshw Half-width of the water artifact in ppm
#' @return
#' A numeric vector of the signal intensities with the water signal removed.
rm_water_signal2 <- function(x, y, wshw) {
    assert(is_num(x), is_num(y), is_num(wshw, 1))
    logf("Removing water signal")
    cntr <- (x[1] + x[length(x)]) / 2
    y[which(x > cntr - wshw & x < cntr + wshw)] <- min(y)
    y
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D. Added code for bwc > 1.
rm_negative_signals <- function(spec) {
    logf("Removing negative signals")
    errmsg <- "Water signal not removed yet. Call `rm_water_signal()` first."
    if (is.null(spec$y_nows)) stop(errmsg)
    spec$y_pos <- abs(spec$y_nows)
    spec
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt:
#' Extracted and refactored corresponding code from MetaboDecon1D.
#' Added code for bwc > 1.
#' Dropped support for bwc == 0.
rm_negative_signals2 <- function(y) {
    logf("Removing negative signals")
    abs(y)
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
#' The input `spec` list with two additional entries:
#' - `spec$Z`: A data frame containing the intermediate smoothed values after each iteration.
#' - `spec$y_smooth`: A numeric vector of the smoothed values after
#'
#' @details
#' Old and slow version producing the same results as the
#' implementation within `deconvolution` from `MetaboDecon1D_deconvolution.R`.
#'
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
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
#'
#' @title Smooth Signal Intensities using a Moving Average
#'
#' @description
#' Smoothes signal intensities by applying a [moving average](
#' https://en.wikipedia.org/wiki/Moving_average) filter with a window size of k.
#'
#' @param y Signal intensities in arbtary units (au)
#'
#' @param reps The number of times to apply the moving average.
#'
#' @param k The number of points within the moving average window. Must be odd,
#' so the smoothed point is in the middle of the window.
#'
#' @return
#' A numeric vector of the smoothed values.
#'
#' @author
#' 2020-2021 Martina Haeckl:
#' Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt:
#' Extracted and refactored corresponding code from MetaboDecon1D.
#' Added code for bwc > 1.
#' Dropped support for bwc == 0.
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

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
find_peaks <- function(spec) {
    spec$d <- calc_second_derivative(spec$y_smooth)
    spec$peak <- find_peaks2(spec$y_smooth)
    spec
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
find_peaks2 <- function(y) {
    logf("Starting peak selection")
    d <- calc_second_derivative(y)
    pc <- get_peak_centers_fast(d)
    rb <- get_right_borders_fast(d, pc)
    lb <- get_left_borders_fast(d, pc)
    sc <- get_peak_scores_fast(d, pc, lb, rb)
    P <- data.frame(left = lb, center = pc, right = rb, score = sc)
    logf("Detected %d peaks", length(pc))
    P
}

#' @noRd
#' @title Filter Peaks with Low Scores
#' @description
#' Modular replacement for [filter_peaks()]. Takes atomic vectors instead of
#' the monolithic `ispec` list. Peaks whose center ppm falls outside `sfr` are
#' used to estimate noise; signal-region peaks with scores below
#' `mean + delta * sd` are removed. Peaks inside any `igr` region are also
#' removed.
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
            any(vapply(igr, function(r) {
                c >= min(r) & c <= max(r)
            }, logical(1)))
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
#' Modular replacement for [fit_lorentz_curves()]. Works directly in ppm and
#' uses the same algorithm as the Rust backend: 3-point peak stencil with
#' iterative refinement. Returns a data frame with columns `x0`, `A`,
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

    # Mirror shoulders (match Rust PeakStencil::mirror_shoulder)
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

    # Solve 3-equation system (match Rust FitterAnalytical)
    solve_params <- function(x1, x2, x3, y1, y2, y3) {
        # maximum_position (x0)
        num <- x1^2 * y1 * (y2 - y3) + x2^2 * y2 * (y3 - y1) +
            x3^2 * y3 * (y1 - y2)
        den <- 2 * ((x1 - x2) * y1 * y2 + (x2 - x3) * y2 * y3 +
            (x3 - x1) * y3 * y1)
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
        sup <- lorentz_sup_raw(rs_x, p$sfhw, p$hw2, p$maxp)
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
    sfhw <- p$sfhw[ok]; hw2 <- p$hw2[ok]; maxp <- p$maxp[ok]

    # Convert (sfhw, hw2, maxp) → (x0, A, lambda) for decon2 format
    lambda <- sqrt(hw2)
    A <- sfhw / lambda
    data.frame(x0 = maxp, A = A, lambda = lambda)
}

#' @noRd
#' @title Superposition using (sfhw, hw2, maxp) parameterization
#' @description
#' Computes the superposition of Lorentzians at positions `x` using the
#' transformed parameters (sfhw, hw2, maxp) as used during fitting. This
#' avoids unnecessary sqrt/division during iteration.
#' @return Numeric vector of superposition values.
lorentz_sup_raw <- function(x, sfhw, hw2, maxp) {
    result <- numeric(length(x))
    for (j in seq_along(maxp)) {
        result <- result + sfhw[j] / (hw2[j] + (x - maxp[j])^2)
    }
    result
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
#' @param igr List of length-2 numeric vectors specifying chemical shift
#' regions to ignore. Peaks whose centers fall inside any ignore region are
#' excluded.
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
filter_peaks <- function(ispec, sfr, delta = 6.4, force = FALSE, bwc = 1,
                         igr = list()) {
    assert(is_ispec(ispec))
    logf("Removing peaks with low scores")
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
    ispec$peak$high <- pok & (psc > mu + delta * sigma)
    if (length(igr) > 0) {
        cp <- if (bwc < 1) sdp[pct] else ppm[pct]
        in_igr <- vapply(cp, function(c) {
            any(vapply(igr, function(r) {
                c >= min(r) & c <= max(r)
            }, logical(1)))
        }, logical(1))
        ispec$peak$high <- ispec$peak$high & !in_igr
    }
    ispec$peak$region <- "norm"
    ispec$peak$region[in_left_sfr] <- "sfrl"
    ispec$peak$region[in_right_sfr] <- "sfrr"
    logf("Removed %d peaks", sum(!ispec$peak$high))
    ispec
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D. Added code for bwc > 1.
fit_lorentz_curves <- function(spec, nfit = 3, bwc = 1) {
    logf("Initializing Lorentz curves")
    spec$lci <- lc <- init_lc(spec) # Lorentz Curves Initialized
    spec$lca <- vector("list", length = nfit) # Lorentz Curves Approximated
    logf("Refining Lorentz Curves")
    for (i in 1:nfit) spec$lca[[i]] <- lc <- refine_lc_v14(spec, lc$Z)
    A <- lc$A
    lambda <- lc$lambda
    w <- lc$w
    if (bwc < 1) {
        limits <- c(0, max(spec$sdp) + (1 / spec$sf[1]))
        integrals <- lorentz_int(w, A, lambda, limits = limits)
    } else {
        integrals <- A * (- pi)
    }
    spec$lcr <- list(A = A, lambda = lambda, w = w, integrals = integrals)
    spec
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl:
#' Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt:
#' Extracted and refactored corresponding code from MetaboDecon1D.
#' Added code for bwc > 1.
fit_lorentz_curves <- function(spec, nfit = 3, bwc = 1) {
    logf("Initializing Lorentz curves")
    spec$lci <- lc <- init_lc(spec) # Lorentz Curves Initialized
    spec$lca <- vector("list", length = nfit) # Lorentz Curves Approximated
    logf("Refining Lorentz Curves")
    for (i in seq_len(nfit)) spec$lca[[i]] <- lc <- refine_lc_v14(spec, lc$Z)
    A <- lc$A
    lambda <- lc$lambda
    w <- lc$w
    if (bwc < 1) {
        limits <- c(0, max(spec$sdp) + (1 / spec$sf[1]))
        integrals <- lorentz_int(w, A, lambda, limits = limits)
    } else {
        integrals <- A * (- pi)
    }
    spec$lcr <- list(A = A, lambda = lambda, w = w, integrals = integrals)
    spec
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
#' [metabodecon::filter_peaks()]. For details see `CHECK-2: signal free region calculation`
#' in `TODOS.md`. (Update 2025-09-14: TODOS are no longer tracked in a separate
#' file, but outside of the repository. To retrieve the last actively maintained
#' version of `TODOS.md`, checkout commit 8b1f61b, i.e., v1.5.0.)
#' @author 2024-2025 Tobias Schmidt: initial version.
enrich_sfr <- function(sfr, x) {
    assert(is_ispec(x) || is_idecon(x))
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
#' Same as [metabodecon::enrich_sfr()], but only requires the chemical shift vector instead
#' of the full `spectrum` or `idecon` object. That's easier to test and
#' maintain.
#' @author 2025 Tobias Schmidt: initial version.
enrich_sfr2 <- function(sfr, cs) {
    n <- length(cs)
    ppm_max <- max(cs)
    ppm_nstep <- diff(range(cs)) / n
    left_ppm <- sfr[1]
    right_ppm <- sfr[2]
    left_dp <- (n + 1) - (ppm_max - left_ppm) / ppm_nstep
    left_sdp <- left_dp / 1000
    right_dp <- (n + 1) - (ppm_max - right_ppm) / ppm_nstep
    right_sdp <- right_dp / 1000
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
#' [metabodecon::rm_water_signal()]. For details see `CHECK-3: water signal calculation` in
#' `TODOS.md`. (Update 2025-09-14: TODOS are no longer tracked in a separate
#' file, but outside of the repository. To retrieve the last actively maintained
#' version of `TODOS.md`, checkout commit 8b1f61b, i.e., v1.5.0.)
#' @author 2024-2025 Tobias Schmidt: initial version.
enrich_wshw <- function(wshw, x) {
    assert(is_ispec(x) || is_idecon(x))
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

#' @noRd
#' @description
#' Same as [metabodecon::enrich_wshw()], but only requires the chemical shift vector instead
#' of the full `spectrum` or `idecon` object. That's easier to test and
#' maintain.
#' @author 2025 Tobias Schmidt: initial version.
enrich_wshw2 <- function(wshw, cs) {
    n <- length(cs)
    ppm_nstep <- diff(range(cs)) / n
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
        is_right_root <- (d[r] < 0 && d[r + 1] >= 0)
        is_right_maximum <- (d[r] > d[r - 1] && d[r] >= d[r + 1])
        is_right_border <- is_right_root || is_right_maximum
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
        is_left_root <- (d[l] < 0 && d[l - 1] >= 0)
        is_left_maximum <- (d[l] > d[l + 1] && d[l] >= d[l - 1])
        is_left_border <- is_left_root || is_left_maximum
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

get_peak_centers_fast <- function(d) {
    dl <- c(NA, d[-length(d)])
    dr <- c(d[-1], NA)
    pc <- which(d < 0 & d <= dl & d < dr)
    pc
}

#' @noRd
#' @param d Vector of second derivative (of the smoothed signal intensities).
#' @param pc Peak center indices.
#' @author 2025 Tobias Schmidt: initial version.
get_right_borders_fast <- function(d, pc, bwc = 2) {
    # Example Inputs:
    # d <- c(-2, -4, -2, 1, 2, 1, -1, 1, -2, -1) # Second Derivative
    # pc <- c(2, 7, 9) # Corresponding Peak Center Indices

    dl <- c(NA, d[-length(d)])
    # Second Derivative shifted left by 1 position. Example:
    # dl == c(NA, -2, -4, -2, 1, 2, 1, -1, 1, -2)

    dr <- c(d[-1], NA)
    # Second Derivative shifted right by 1 position. Example:
    # dr == c(-4, -2, 1, 2, 1, -1, 1, -2, -1, NA)

    is_right_root <- (dr >= 0) & (d < 0)
    is_local_maximum <- (dl < d) & (dr <= d)
    isrbc <- is_right_root | is_local_maximum
    isrbc[is.na(isrbc)] <- FALSE
    # Is Right Border Candidate? Vector of booleans describing whether a certain
    # position fulfills the criteria to be a right border. Example:
    # isrbc == c(F, F, TRUE, F, TRUE, F, TRUE, TRUE, F, F)

    rbc <- c(which(isrbc), NA)
    # Right Border Candidates. Example:
    # rbc == c(3, 5, 7, 8, NA)

    rbci <- cumsum(isrbc) + 1
    # Right Border Candidate Indices. I.e., for a peak at position i, rbci[i]
    # gives the position of the nearest right border candidate in rbc. Example:
    #
    # rbci == c(1, 1, 2, 2, 3, 3, 4, 5, 5, 5)
    #
    # I.e., the first two positions (1:2) have their nearest right border
    # candidate at rbc[1] (which points to position 3), the next two positions
    # (3:4) have their nearest right border candidate at rbc[2] (== 5) and the
    # next two positions (5:6) have their nearest right border candidate at
    # rbc[3] (== 7). The next position (7) has its nearest right border
    # candidate at rbc[4] (== 8) and the last three positions have their nearest
    # right border candidate at rbc[5], which is NA, as there is no right border
    # candidate after position 8.

    rbc[rbci[pc]] # Right Borders.
}

get_left_borders_fast <- function(d, pc) {
    nd <- length(d)
    pcrev <- rev(length(d) - pc + 1)
    drev <- rev(d)
    lbrev <- get_right_borders_fast(drev, pcrev)
    lb <- rev(nd - lbrev + 1)
    lb
}

get_peak_scores_fast <- function(d, pc, lb, rb) {
    # Calculate interval scores as minimum of left-to-center and center-to-right
    # interval sums. Use cumsum for maximum speed. Prepend a zero so we can
    # calculate interval sums as cumsum[pc+1]-cumsum[lb] instead of
    # cumsum[pc]-cumsum[lb-1]. This way we don't need to handle the lb==1 case
    # separately. Set scores to zero if lb or rb is NA.
    a <- abs(d)
    cumsum0 <- c(0, cumsum(replace(a, is.na(a), 0)))
    sum_lj <- cumsum0[pc + 1] - cumsum0[lb]
    sum_jr <- cumsum0[rb + 1] - cumsum0[pc]
    score <- pmin(sum_lj, sum_jr)
    score[is.na(score)] <- 0
    score
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
    mse <- mse(y[lmr], rowSums(Z))
    if (verbose) logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, as, ds)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

#' @noRd
#' @title Initialize Lorentz Curve Parameters
#' @param p
#' data.frame with columns `center`, `left`, `right` and `high`, giving the
#' indices of the datapoints identified as peak triplet positions (PTPs) and
#' whether the peak is above the score threshold (`high == TRUE`) or not.
#' @param x Vector of scaled data points (SDP).
#' @param y Smoothed intensity values.
#' @return
#' A list with elements:
#' * `A`: Area parameter of each Lorentz curve.
#' * `lambda`: Width parameter of each Lorentz curve.
#' * `w`: Position parameter of each Lorentz curve.
#' * `Z`: Matrix with height of each Lorentz curve (col) at each PTP (row).
#' * `D`: Data frame with distances between PTPs.
#' * `P`: Data frame with information about each PTP.
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.
init_lc2 <- function(x, y, p, verbose = TRUE) {

    # Init values
    ir <- p$right[p$high]  #
    ic <- p$center[p$high] # Index of each peak triplet position (PTP)
    il <- p$left[p$high]   #
    lmr <- sort(unique(c(il, ic, ir))) # Combined PTP indices
    rr <- match(ir, lmr) #
    rc <- match(ic, lmr) # Rank of each PTP
    rl <- match(il, lmr) #
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
    mse <- mse(y[lmr], rowSums(Z))
    logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, sl, sc, sr, ql, qc, qr)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

#' @noRd
#' @author
#' 2020-2021 Martina Haeckl: Wrote initial version as part of MetaboDecon1D.\cr
#' 2024-2025 Tobias Schmidt: Extracted and refactored corresponding code from
#' MetaboDecon1D.\cr
refine_lc2 <- function(lc) {

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
    num <- -wc^4 * yc^2 * yrl^2 - wr^4 * yr^2 * ycl^2 - wl^4 * yrc^2 * yl^2
    num <- num + 4 * wc * wl^3 * yc * ((-yr) + yc) * yl^2
    num <- num + 4 * wc^3 * wl * yc^2 * yl * ((-yr) + yl)
    num <- num + 4 * wr^3 * yr^2 * ycl * (wc * yc - wl * yl)
    num <- num + 4 * wr * yr * (
        wc^3 * yc^2 * yrl
        - wc * wl^2 * yc * (yr + yc - 2 * yl) * yl
        + wl^3 * yrc * yl^2
        - wc^2 * wl * yc * yl * (yr - 2 * yc + yl)
    )
    num <- num + 2 * wc^2 * wl^2 * yc * yl * (yr^2 - 3 * yc * yl + yr * (yc + yl))
    num <- num + 2 * wr^2 * yr * (
        -2 * wc * wl * yc * yl * (-2 * yr + yc + yl)
        + wl^2 * yl * (yr * (yc - 3 * yl) + yc * (yc + yl))
        + wc^2 * yc * (yr * (-3 * yc + yl) + yl * (yc + yl))
    )
    den <- 2 * sqrt((wr * yr * ycl + wl * yrc * yl + wc * yc * ((-yr) + yl))^2)
    lambda <- -(sqrt(abs(num)) / den)
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
    num <- -4 * wrc * wrl * wcl * yr * yc * yl
    num <- num * (wr * yr * ycl + wl * yl * yrc + wc * yc * (-yrl))
    num <- num * lambda
    den <- wrc^4 * yr^2 * yc^2
    den <- den - 2 * wrc^2 * yr * yc * (wrl^2 * yr + wcl^2 * yc) * yl
    den <- den + (wrl^2 * yr - wcl^2 * yc)^2 * yl^2
    A <- num / den
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
#' 2024-2026 Tobias Schmidt: initial versions (v1, v2, v3).
lorentz_sup <- function(x, x0, A, lambda, lcpar = NULL) {
    if (!is.null(lcpar)) {
        nams <- names(lcpar)
        if ("A" %in% nams) A <- lcpar[["A"]]
        if ("lambda" %in% nams) lambda <- lcpar[["lambda"]]
        if ("x_0" %in% nams) x0 <- lcpar[["x_0"]]
        if ("x0" %in% nams) x0 <- lcpar[["x0"]]
        if ("w" %in% nams) x0 <- lcpar[["w"]]
    }
    v <- getOption("metabodecon.lorentz_sup_version", NULL)
    if (is.null(v)) {
        cfun <- get_lorentz_sup_c()
        v <- if (is.function(cfun)) 3 else 2
    }
    if (v == 1) {
        sapply(x, function(xi) {
            sum(abs(A * (lambda / (lambda^2 + (xi - x0)^2))))
        })
    } else if (v == 2) {
        Al <- abs(A * lambda)
        l2 <- lambda^2
        result <- numeric(length(x))
        for (j in seq_along(x0)) {
            result <- result + Al[j] / (l2[j] + (x - x0[j])^2)
        }
        result
    } else if (v == 3) {
        cfun <- get_lorentz_sup_c()
        cfun(x, x0, A, lambda)
    } else {
        stop("Invalid lorentz_sup_version: ", v)
    }
}

#' @noRd
#' @description
#' Returns the compiled C implementation of lorentz_sup, or FALSE if compilation
#' is not possible. On first call, attempts to compile via [inline::cfunction()]
#' and caches the result (function or FALSE) in
#' `options("metabodecon.lorentz_sup_c")` so subsequent calls are instant.
#' @author 2026 Tobias Schmidt: initial version.
get_lorentz_sup_c <- function() {
    cfun <- getOption("metabodecon.lorentz_sup_c")
    if (!is.null(cfun)) return(cfun)
    cfun <- tryCatch(
        compile_lorentz_sup_c(),
        error = function(e) FALSE
    )
    options(metabodecon.lorentz_sup_c = cfun)
    cfun
}

#' @noRd
#' @description
#' Compiles the C implementation of lorentz_sup via inline::cfunction.
#' Throws an error if 'inline' is not installed or compilation fails.
#' @author 2026 Tobias Schmidt: initial version.
compile_lorentz_sup_c <- function() {
    if (!requireNamespace("inline", quietly = TRUE)) {
        stop("Package 'inline' required for C version")
    }
    inline::cfunction(
        sig = c(
            x = "numeric", x0 = "numeric",
            A = "numeric", lambda = "numeric"
        ),
        body = '
            int nx = Rf_length(x), np = Rf_length(x0);
            double *px  = REAL(x),  *px0 = REAL(x0);
            double *pA  = REAL(A),  *plam = REAL(lambda);
            SEXP res = Rf_protect(Rf_allocVector(REALSXP, nx));
            double *pr = REAL(res);
            memset(pr, 0, nx * sizeof(double));
            double *Al = (double *)R_alloc(np, sizeof(double));
            double *l2 = (double *)R_alloc(np, sizeof(double));
            for (int j = 0; j < np; j++) {
                Al[j] = fabs(pA[j] * plam[j]);
                l2[j] = plam[j] * plam[j];
            }
            for (int i = 0; i < nx; i++) {
                double s = 0.0, xi = px[i];
                for (int j = 0; j < np; j++) {
                    double d = xi - px0[j];
                    s += Al[j] / (l2[j] + d * d);
                }
                pr[i] = s;
            }
            Rf_unprotect(1);
            return res;
        ',
        includes = '#include <math.h>\n#include <string.h>'
    )
}

#' @noRd
benchmark_lorentz_sup <- function() {

    x <- read_spectrum(metabodecon_file("urine_1"))
    rt <- function(expr) system.time(expr)[["elapsed"]]
    opts <- options()
    on.exit(options(opts))

    logf("Get compile time for C version")
    options(metabodecon.lorentz_sup_c = NULL)
    tc <- rt(cfun <- get_lorentz_sup_c())

    logf("Version NULL: auto-select (C if available, else R v2)")
    options(metabodecon.lorentz_sup_version = NULL)
    td0 <- rt(d0 <- deconvolute(x, verbose = FALSE, use_rust = FALSE))
    cs <- d0$cs
    lcpar <- d0$lcpar
    ts0 <- rt(s0 <- lorentz_sup(cs, lcpar = lcpar))

    logf("Version 1: loop over data points")
    options(metabodecon.lorentz_sup_version = 1)
    td1 <- rt(d1 <- deconvolute(x, verbose = FALSE, use_rust = FALSE))
    ts1 <- rt(s1 <- lorentz_sup(cs, lcpar = lcpar))

    logf("Version 2: loop over peaks")
    options(metabodecon.lorentz_sup_version = 2)
    td2 <- rt(d2 <- deconvolute(x, verbose = FALSE, use_rust = FALSE))
    ts2 <- rt(s2 <- lorentz_sup(cs, lcpar = lcpar))

    logf("Version 3: loop in C")
    options(metabodecon.lorentz_sup_version = 3)
    td3 <- rt(d3 <- deconvolute(x, verbose = FALSE, use_rust = FALSE))
    ts3 <- rt(s3 <- lorentz_sup(cs, lcpar = lcpar))

    logf("Version 1 + Rust backend")
    options(metabodecon.lorentz_sup_version = 1)
    tdr1 <- rt(dr1 <- deconvolute(x, verbose = FALSE, use_rust = TRUE))
    tdr1b <- rt(dr1b <- deconvolute(x, verbose = FALSE, use_rust = TRUE))
    tdr1c <- rt(dr1c <- deconvolute(x, verbose = FALSE, use_rust = TRUE))

    logf("Version 2 + Rust backend")
    options(metabodecon.lorentz_sup_version = 2)
    tdr2 <- rt(dr2 <- deconvolute(x, verbose = FALSE, use_rust = TRUE))

    logf("Version 3 + Rust backend")
    options(metabodecon.lorentz_sup_version = 3)
    tdr3 <- rt(dr3 <- deconvolute(x, verbose = FALSE, use_rust = TRUE))

    # Print results
    logf("")
    logf("Urine1: %d datapoints, %d peaks",
        length(cs), length(lcpar$A))
    logf("")
    logf("Correctness checks:")
    logf("  all.equal(d1, d0):  %s", isTRUE(all.equal(d1, d0)))
    logf("  all.equal(d2, d0):  %s", isTRUE(all.equal(d2, d0)))
    logf("  all.equal(d3, d0):  %s", isTRUE(all.equal(d3, d0)))
    logf("")
    logf("  all.equal(dr1, d0):  %s", isTRUE(all.equal(dr1, d0)))
    logf("  all.equal(dr2, d0):  %s", isTRUE(all.equal(dr2, d0)))
    logf("  all.equal(dr3, d0):  %s", isTRUE(all.equal(dr3, d0)))
    logf("")
    logf("  all.equal(dr2, dr1):  %s", isTRUE(all.equal(dr2, dr1)))
    logf("  all.equal(dr3, dr1):  %s", isTRUE(all.equal(dr3, dr1)))
    logf("")
    logf("  all.equal(s1, s0):  %s", isTRUE(all.equal(s1, s0)))
    logf("  all.equal(s2, s0):  %s", isTRUE(all.equal(s2, s0)))
    logf("  all.equal(s3, s0):  %s", isTRUE(all.equal(s3, s0)))
    logf("")
    logf("lorentz_sup timings:")
    logf("  v0 (auto):       %7.3f s  (%.1fx)", ts0, ts1 / ts0)
    logf("  v1 (dp-loop):    %7.3f s", ts1)
    logf("  v2 (pk-loop):    %7.3f s  (%.1fx)", ts2, ts1 / ts2)
    logf("  v3 (C-loop):     %7.3f s  (%.1fx)", ts3, ts1 / ts3)
    logf("  compile time v3: %7.3f s", tc)
    logf("")
    logf("deconvolute timings (R backend):")
    logf("  v0 (auto):       %7.3f s  (%.1fx)", td0, td1 / td0)
    logf("  v1 (dp-loop):    %7.3f s", td1)
    logf("  v2 (pk-loop):    %7.3f s  (%.1fx)", td2, td1 / td2)
    logf("  v3 (C-loop):     %7.3f s  (%.1fx)", td3, td1 / td3)
    logf("")
    logf("deconvolute timings (Rust backend):")
    logf("  v1 (dp-loop):    %7.3f s  (%.1fx)", tdr1, td1 / tdr1)
    logf("  v2 (pk-loop):    %7.3f s  (%.1fx)", tdr2, td1 / tdr2)
    logf("  v3 (C-loop):     %7.3f s  (%.1fx)", tdr3, td1 / tdr3)
    logf("")
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
            logf(
                "Skipping RDS save: confirmation required but not in",
                "interactive mode. For details see",
                "`help('generate_lorentz_curves')`."
            )
        }
    }
}
