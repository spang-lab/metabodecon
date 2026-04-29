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
#' Deprecated: setting a non-zero value will produce a deprecation warning.
#' Future versions will always use `wshw=0`.
#'
#' @param use_rust Controls the deconvolution backend. `FALSE` or any numeric
#' value `< 1` (default) uses the R implementation. `TRUE` or any numeric
#' value `>= 1` uses the Rust backend via
#' [mdrb](https://github.com/spang-lab/mdrb). `NULL` auto-detects: uses Rust
#' if available, otherwise R. When set to `TRUE` / `>= 1` and mdrb is not
#' installed, an error is thrown.
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
    nfit=3,     smopts=c(2,5),  delta=6.4,     sfr=NULL,    wshw=NULL,
    ask=NULL,   force=FALSE,    verbose=TRUE,  nworkers=1,  use_rust=FALSE,
    npmax=0,    igrs=list(),    cadir=decon_cachedir()
) {

    # Deprecation warnings for removed parameters
    if (!is.null(wshw)) lifecycle::deprecate_warn(
        "2.0.0", "deconvolute(wshw=)",
        details = "wshw will be removed. Future versions always use wshw=0."
    )
    if (!is.null(ask)) lifecycle::deprecate_warn(
        "2.0.0", "deconvolute(ask=)",
        details = "ask will be removed. Future versions always use ask=FALSE."
    )
    wshw <- wshw %||% 0
    ask  <- ask  %||% FALSE

    # Check inputs
    stopifnot(
        is_spectrum_or_spectra(x),    is_int_or_null(nfit, 1),
        is_int_or_null(smopts, 2),    is_num_or_null(delta, 1),
        is_num_or_null(sfr, 2),       is_num_or_null(wshw, 1),
        is_bool(ask, 1),              is_bool(force, 1),
        is_bool(verbose, 1),          is_int(nworkers, 1),
        is_bool_or_num(use_rust),     is_int(npmax, 1),
        is_list_of_nums(igrs, nv=2),  is_str_or_null(cadir)
    )

    # Set suitable defaults
    sfr <- sfr %||% quantile(x$cs %||% x[[1]]$cs, c(0.9, 0.1))
    if (use_rust == 1) check_mdrb(stop_on_fail = TRUE)

    # Perform deconvolution
    decons2 <- deconvolute_spectra(x,
        nfit, smopts, delta, sfr, wshw,
        ask, force, verbose,
        use_rust, nworkers=nworkers, igrs=igrs, rtyp="decon2",
        npmax=npmax, cadir=cadir
    )

    # Convert and return
    if (length(decons2) == 1) decons2[[1]] else decons2
}

# Internal main functions #####

#' @noRd
#' @inheritParams deconvolute
#' @author 2024-2025 Tobias Schmidt: initial version.
deconvolute_spectra <- function(
    x,                  nfit=3,          smopts=c(2, 5),          delta=6.4,
    sfr=c(3.55, 3.35),  wshw=0,          ask=FALSE,               force=FALSE,
    verbose=TRUE,       use_rust=FALSE,  nworkers=1,              igrs=list(),
    rtyp="decon2",      npmax=0,         cadir=decon_cachedir(),  hash=NULL
) {

    # Check inputs
    assert(
        is_spectrum(x) || is_spectra(x),
        is_int(nfit, 1), is_int(smopts, 2), is_num(delta, 1),
        is_bool(ask, 1), is_bool(force, 1), is_bool(verbose, 1),
        is_num(sfr, 2)  || is_list_of_nums(sfr, length(x), 2),
        is_num(wshw, 1) || is_list_of_nums(wshw, length(x), 1),
        is_bool_or_null(use_rust),
        is_int(nworkers, 1), is.list(igrs),
        is_char(rtyp, 1, "(decon[0-2]|rdecon)"),
        is_int(npmax, 1), is_str_or_null(cadir)
    )

    # Return cached result when hash matches
    if (!is.null(hash)) {
        dc <- getOption("metabodecon.decon_cache")
        if (!is.null(dc$key) && dc$key == hash) {
            logf("Skipping deconvolution (cache hit)")
            return(dc$decons)
        }
    }

    # Configure logging
    if (!verbose) local_options(toscutil.logf.file = nullfile())

    # Init locals
    spectra <- as_spectra(x)
    ns <- length(spectra)
    nc2 <- ceiling(detectCores() / 2)
    nw_apply <- if (nworkers == "auto") min(nc2, ns) else nworkers
    nw_apply_str <- if (nw_apply == 1) "1 worker" else sprintf("%d workers", nw_apply)
    nw_deconv <- 1
    ns_str <- if (ns == 1) "1 spectrum" else sprintf("%d spectra", ns)
    adjno <- get_adjno(spectra, ask)
    sfr_list <- get_sfr(spectra, sfr, ask, adjno)
    wshw_list <- get_wshw(spectra, wshw, ask, adjno)
    smopts_list <- get_smopts(spectra, smopts)
    igrs_list <- list(igrs)
    cadir_list <- list(cadir)

    # Deconvolute spectra
    logf("Starting deconvolution of %s using %s", ns_str, nw_apply_str)
    starttime <- Sys.time()
    decon_list <- mcmapply(nw_apply, deconvolute_spectrum,
        spectra,
        nfit, smopts_list, delta, sfr_list, wshw_list,
        ask, force, verbose,
        use_rust, nw_deconv, igrs_list, rtyp,
        npmax, cadir_list
    )
    decons <- as_collection(decon_list, rtyp)
    duration <- format(round(Sys.time() - starttime, 3))
    logf("Finished deconvolution of %s in %s", ns_str, duration)

    # Store in cache if hash was provided
    if (!is.null(hash)) {
        options(metabodecon.decon_cache = list(key = hash, decons = decons))
    }

    # Return
    decons
}

#' @noRd
#' @inheritParams deconvolute_spectra
#' @author 2024-2025 Tobias Schmidt: initial version.
#' @examples
#' x <- sap[[1]];
#' nfit <- 3; smopts <- c(1,3); delta <- 3; sfr <- c(3.2,-3.2); wshw <- 0;
#' ask <- FALSE; force <- FALSE; verbose <- FALSE;
#' use_rust <- FALSE; nworkers <- 1; igrs <- list(); rtyp <- "decon2"
#' decon <- deconvolute_spectrum(
#'      x, nfit, smopts, delta, sfr, wshw,
#'      ask, force, verbose,
#'      use_rust, nworkers, igrs, rtyp
#' )
#'
#' x <- read_spectrum(metabodecon_file("urine_1"))
#' nfit <- 3; smopts <- c(2,5); delta <- 6;
#' sfr <- quantile(x$cs, c(0.9, 0.1)); wshw <- 0;
#' ask <- FALSE; force <- FALSE; verbose <- TRUE;
#' use_rust <- TRUE; nworkers <- 1; igrs <- list(); rtyp <- "decon2"
#' npmax=1000
#' decon1 <- deconvolute_spectrum(
#'      x, nfit, smopts, delta, sfr, wshw,
#'      ask, force, verbose, use_rust, nworkers, igrs, rtyp,
#'      npmax
#' )
deconvolute_spectrum <- function(x,
    nfit=3, smopts=c(2, 5), delta=6.4, sfr=c(3.55, 3.35), wshw=0,
    ask=FALSE, force=FALSE, verbose=TRUE,
    use_rust=FALSE, nworkers=1, igrs=list(), rtyp="decon2",
    npmax=0, cadir=decon_cachedir()
) {

    # Check inputs
    assert(
        is_spectrum(x),
        is_int(nfit, 1),  is_int(smopts, 2),      is_num(delta, 1),
        is_num(sfr, 2),   is_num(wshw, 1),        is_bool(force, 1),
        is_bool_or_num(use_rust),  is_int(nworkers, 1),
        is_list_of_nums(igrs, nv=2),
        is_char(rtyp, 1, "(decon[0-2]|rdecon)"),
        is_int(npmax, 1), is_str_or_null(cadir)
    )

    # Init locals
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())
    name <- get_name(x)
    backend <- if (use_rust == 1) "Rust" else "R"
    suffix <- sprintf(" using %s backend", backend)

    # Stop early if deconvolution is cached
    if (npmax >= 1) {
        cache <- disk_cache(cadir)
        args <- list("ds", x, sfr, wshw, force, use_rust, npmax, rtyp)
        cachehash <- rlang::hash(args)
        if (cache$exists(cachehash)) {
            logf("Cache hit for %s (npmax=%d)", name, npmax)
            return(cache$get(cachehash))
        }
    }

    # Perform grid search (to replace given nfit, smopts and delta)
    if (npmax >= 1) {
        fmt <- "Starting grid deconvolution of %s using %s backend"
        logf(fmt, name, backend)
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
    decon <- if (use_rust == 1) {
        deconvolute_rust(x, sfr, smopts, delta, nfit, igrs, nworkers, args)
    } else {
        deconvolute_r(x, sfr, smopts, delta, nfit, force, igrs, args)
    }

    # Format, cache and return
    logf("Formatting return object as %s", rtyp)
    convert <- switch(rtyp,
        "decon0"=as_decon0, "decon1"=as_decon1, "decon2"=as_decon2,
        "rdecon"=as_rdecon
    )
    decon <- convert(decon)
    if (npmax >= 1) cache$set(cachehash, decon)
    logf("Finished deconvolution of %s", name)
    decon
}

#' @noRd
deconvolute_rust <- function(
    x, sfr, smopts, delta, nfit, igrs, nworkers, args
) {
    mdrb_spectrum <- mdrb::Spectrum$new(x$cs, x$si, sfr)
    mdrb_deconvr <- mdrb::Deconvoluter$new()
    mdrb_deconvr$set_moving_average_smoother(smopts[1], smopts[2])
    mdrb_deconvr$set_noise_score_selector(delta)
    mdrb_deconvr$set_analytical_fitter(nfit)
    for (r in igrs) mdrb_deconvr$add_ignore_region(r[1], r[2])
    mdrb_decon <- if (nworkers > 1) {
        mdrb_deconvr$set_threads(nworkers)
        mdrb_deconvr$par_deconvolute_spectrum(mdrb_spectrum)
    } else {
        mdrb_deconvr$deconvolute_spectrum(mdrb_spectrum)
    }
    new_rdecon(x, args, mdrb_spectrum, mdrb_deconvr, mdrb_decon)
}

#' @noRd
deconvolute_r <- function(
    x, sfr, smopts, delta, nfit, force, igrs, args
) {
    cs <- x$cs;
    si <- x$si
    sfr_igr <- list(c(Inf, max(sfr)), c(min(sfr), -Inf))
    igrs <- c(sfr_igr, igrs)
    sm <- smooth_signals2(si, smopts[1], smopts[2])
    peaks <- find_peaks2(sm)
    peaks <- filter_peaks2(peaks, cs, sfr, delta, force, igrs)
    lcpar <- fit_lorentz_curves2(cs, si, peaks, nfit)
    sup <- lorentz_sup(cs, lcpar = lcpar)
    sit <- data.frame(
        wsrm = si * 1e6, nvrm = si * 1e6,
        sm = sm * 1e6, sup = sup
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
    decon
}

grid_deconvolute_spectra <- function(
    x,
    sfr=NULL,
    verbose=TRUE,
    nworkers=1,
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
        is_int(nworkers, 1),
        is_bool_or_num(use_rust),
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
    nworkers <- min(nworkers, length(u))
    gridlist[!seen] <- mcmapply(nworkers, grid_deconvolute_spectrum,
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
    x, sfr=NULL, verbose=TRUE, use_rust=FALSE,
    smit = c(2,3), smws = c(3,5,7,9), delta = 2:8, nfit = c(3,4,5),
    cadir=decon_cachedir()
) {

    assert(
        is_spectrum(x), is_num_or_null(sfr, 2), is_bool(verbose, 1),
        is_bool_or_num(use_rust), is_str_or_null(cadir)
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
get_right_borders_fast <- function(d, pc) {
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
    .Call(lorentz_sup_c, as.double(x), as.double(x0),
          as.double(Al), as.double(l2))
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
