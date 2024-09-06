# Main ####

#' @export
#'
#' @title
#'  Generate Lorentz Curves from NMR Spectra
#'
#' @description
#'  Deconvolutes NMR spectra by modeling each detected signal within a spectrum as Lorentz Curve.
#'
#' @inheritParams read_spectrum
#' @param data_path
#'  Either the path to a directory containing measured NMR spectra, a dataframe as returned by [read_spectrum()], or a list of such dataframes.
#' @param make_rds
#'  Logical or character. If TRUE, stores results as an RDS file on disk. If a character string, saves the RDS file with the specified name. Should be set to TRUE if many spectra are evaluated to decrease computation time.
#' @param nfit
#'  Number of iterations for approximating the parameters for the Lorentz curves.
#' @param wshw
#'  Half-width of the water artifact in ppm.
#' @param sfr
#'  Numeric vector with two entries: the ppm positions for the left and right border of the signal-free region of the spectrum.
#' @param smopts
#'  Numeric vector with two entries: the number of smoothing iterations and the number of data points to use for smoothing (must be odd).
#' @param delta
#'  Threshold value to distinguish between signal and noise. The higher the value, the more peaks get filtered out. The exact definition is as follows: a peak `i` gets filtered out, if his score is lower than `mu + s * delta`, where `mu` is the average peak score within the signal free region (SFR) and `s` is the standard deviation of peak scores in the SFR.
#' @param delta
#'  \loadmathjax Threshold for peak filtering. Higher values result in more peaks being filtered out. A peak is filtered if its score is below \mjeqn{\mu + \sigma \cdot \delta}{mu + s * delta}, where \mjeqn{\mu}{mu} is the average peak score in the signal-free region (SFR), and \mjeqn{\sigma}{s} is the standard deviation of peak scores in the SFR.
#' @param sf
#'  Numeric vector with two entries: the factors to scale the x-axis and y-axis.
#' @param ask
#'  Logical. Whether to ask for user input during the deconvolution process. If FALSE, the provided default values will be used.
#' @param debug
#'  Logical. Whether to return additional intermediate results for debugging purposes.
#' @param nworkers
#'  Number of workers to use for parallel processing. If `"auto"`, the number of workers will be determined automatically. If a number greater than 1, it will be limited to the number of spectra.
#' @param share_stdout
#'  Whether to share the standard output (usually your terminal) of the main process with the worker processes. Only relevant if `nworkers` is greater than 1. Note that this can cause messages from different workers to get mixed up, making the output hard to follow.
#' @param force
#'  If FALSE, the function stops with an error message if no peaks are found in the signal free region (SFR), as these peaks are required as a reference for peak filtering. If TRUE, the function instead proceeds without peak filtering, potentially increasing runtime and memory usage significantly.
#' @param verbose
#'  Logical. Whether to print log messages during the deconvolution process.
#'
#' @return
#'  A 'GLCDecon' as described in [Metabodecon Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
#' @details
#'  First, an automated curvature based signal selection is performed. Each signal is represented by 3 data points to allow the determination of initial Lorentz curves. These Lorentz curves are then iteratively adjusted to optimally approximate the measured spectrum.
#'
#' @examples
#' # Define the paths to the example datasets we want to deconvolute:
#' # `sim_dir`: directory containing 16 simulated spectra
#' # `sim_01`: path to the first spectrum in the `sim` directory
#' # `sim_01_spec`: the first spectrum in the `sim` directory as a dataframe
#' sim_dir <- metabodecon_file("sim_subset")
#' sim_1_dir <- file.path(sim_dir, "sim_01")
#' sim_2_dir <- file.path(sim_dir, "sim_02")
#' sim_1_spec <- read_spectrum(sim_1_dir)
#' sim_2_spec <- read_spectrum(sim_2_dir)
#' sim_12_specs <- list(sim_1_spec, sim_2_spec)
#'
#' # Define a little wrapper function so we don't have to provide all parameters
#' # every time we want to start the deconvolution procedure:
#' glc2 <- function(data_path) {
#'     generate_lorentz_curves(
#'         data_path = data_path,
#'         ask = FALSE,
#'         sfr = c(3.42, 3.58),
#'         ws = 0,
#'         smopts = c(1, 5),
#'         delta = 0.1,
#'         nworkers = 2,
#'         verbose = FALSE
#'     )
#' }
#'
#' # Deconvolute each input:
#' decon_sim_1_dir <- glc2(sim_1_dir)
#' decon_sim_1_spec <- glc2(sim_1_spec)
#' decon_sim_12_specs <- glc2(sim_12_specs)
#' decon_sim_dir <- glc2(sim_dir)
#'
#' # Make sure the results for the first spectrum are the same:
#' compare <- function(decon1, decon2) {
#'     ignore <- which(names(decon1) %in% c("number_of_files", "filename"))
#'     equal <- all.equal(decon1[-ignore], decon2[-ignore])
#'     stopifnot(isTRUE(equal))
#' }
#' compare(decon_sim_1_dir, decon_sim_1_spec)
#' compare(decon_sim_1_dir, decon_sim_12_specs[[1]])
#' compare(decon_sim_1_dir, decon_sim_dir[[1]])
#'
#' # Below example uses data from a real NMR experiment, instead of (small)
#' # simulated datasets and therefor requires multiple seconds to run. Because
#' # `ask` is TRUE in this example (the default value), the user will be asked
#' # for input during the deconvolution. To avoid this, set `ask = FALSE`.
#' \dontrun{
#' example_datasets <- download_example_datasets()
#' urine_1 <- file.path(example_datasets, "bruker/urine/urine_1")
#' decon_urine_1 <- generate_lorentz_curves(urine_1)
#' }
generate_lorentz_curves <- function(data_path = metabodecon_file("urine_1"),
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
                                    nworkers = 1,
                                    share_stdout = FALSE,
                                    force = FALSE,
                                    verbose = TRUE) {
    opts <- if (!verbose) options(toscutil.logf.file = nullfile())
    on.exit(options(opts), add = TRUE)
    spectra <- as_spectra(data_path, file_format, expno, procno)
    gspecs <- as_gspecs(spectra, sf)
    adjno <- get_adjno(gspecs, sfr, wshw, ask)
    gspecs <- get_sfrs(gspecs, sfr, ask, adjno)
    gspecs <- get_wsrs(gspecs, wshw, ask, adjno)
    gdecons <- deconvolute_gspecs(gspecs, smopts, delta, nfit, nworkers, share_stdout, force)
    # CONTINUE HERE
    # TODO: Replace sapply call below with as_decons2(gdecons)
    decons <- sapply(gdecons, function(d) d$ret, simplify = FALSE)
    store_as_rds(decons, make_rds, data_path)
    decons <- if (debug) gdecons else decons
    if (nfiles == 1) decons[[1]] else decons
}

# Helpers ####

deconvolute_gspecs <- function(gspecs,
                               smopts,
                               delta,
                               nfit,
                               nworkers,
                               share_stdout,
                               force) {
    nfiles <- length(gspecs)
    if (nworkers == "auto") nworkers <- ceiling(parallel::detectCores() / 2)
    nworkers <- min(nworkers, length(gspecs))
    starttime <- Sys.time()
    if (nworkers == 1) {
        logf("Starting deconvolution of %d spectra using 1 worker", nfiles)
        gdecons <- lapply(seq_len(nfiles), function(i) {
            deconvolute_spectrum(gspecs[[i]], smopts, delta, nfit, nfiles, force)
        })
    } else {
        logf("Starting %d worker processes", nworkers)
        cl <- parallel::makeCluster(nworkers, outfile = if (share_stdout) "" else nullfile())
        on.exit(parallel::stopCluster(cl), add = TRUE)
        logf("Exporting required functions and data to workers")
        exportvars <- c("logf", "fg", "deconvolute_spectrum", "gspecs", "smopts", "delta", "nfit", "nfiles", "debug")
        parallel::clusterExport(cl, exportvars, envir = environment())
        logf("Starting deconvolution of %d spectra using %d workers", nfiles, nworkers)
        gdecons <- parallel::parLapply(cl, seq_len(nfiles), function(i) {
            opts <- options(toscutil.logf.sep1 = sprintf(" PID %d ", Sys.getpid()))
            on.exit(options(opts), add = TRUE)
            deconvolute_spectrum(gspecs[[i]], smopts, delta, nfit, nfiles, force)
        })
    }
    names(gdecons) <- names(gspecs)
    endtime <- Sys.time()
    duration <- endtime - starttime
    logf("Finished deconvolution of %d spectra in %s", nfiles, format(round(duration, 3)))
    gdecons
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

# GLC Classes ####

#' @noRd
#' @title Convert to 'gspec' or 'gspecs'
#' @description Convert an object to a `gspec` or `gspecs` object, which are objects used internally by `generate_lorentz_curves()` to represent one or more spectra. For details see [metabodecon_classes].
#' @param x An object of type `spectrum` or `spectra`.
#' @param sf Numeric vector with two entries: the factors to scale the x-axis and y-axis.
#' @return An object of type `gspec` or `gspecs`.
#' @examples
#' sim <- metabodecon_file("sim_subset")
#' sim_spectra <- read_spectra(sim)
#' sim_gspecs <- as_gspecs(sim_spectra)
#' sim_01_spectrum <- sim_spectra[[1]]
#' sim_01_gspec <- as_gspec(sim_01_spectrum)
as_gspec <- function(x, sf = c(1e3, 1e6)) {
    if (!is_spectrum(x)) stop("Input must be a spectrum object, not ", class(x))
    y_raw <- x$si # Raw signal intensities
    y_scaled <- y_raw / sf[2] # Scaled signal intensities
    n <- length(y_raw) # Number of data points
    dp <- seq(n - 1, 0, -1) # Data point numbers
    sdp <- seq((n - 1) / sf[1], 0, -1 / sf[1]) # Scaled data point numbers [^1]
    ppm <- x$cs # Parts per million
    hz <- x$fq # Frequency in Hz
    ppm_range <- diff(range(x$cs)) # Range of the chemical shifts in ppm.
    ppm_max <- max(x$cs) # Maximum chemical shift in ppm.
    ppm_min <- min(x$cs) # Minimum chemical shift in ppm.
    ppm_step <- ppm_range / (n - 1) # Step size calculated correctly.
    ppm_nstep <- ppm_range / n # Wrong, but backwards compatible [^2].
    name <- x$name # Name of the spectrum
    g <- locals(without = "x")
    structure(g, class = "gspec")
    # [^1]: Same as `dp / sf[1]`, but with slight numeric differences, so we stick with the old calculation method for backwards compatibility.
    # [^2]: Example: ppm = 11, 23, 35, 47 ==> ppm_step == 12, ppm_nstep ~= 10.6 (not really useful, but we need it for backwards compatibility with MetaboDecon1D results)
}

#' @noRd
#' @rdname as_gspec
as_gspecs <- function(x, sf = c(1e3, 1e6)) {
    if (!is_spectra(x)) stop("Input must be of class spectra, not ", class(x))
    gg <- structure(lapply(x, as_gspec, sf = sf), class = "gspecs")
    gg <- set_names(gg, get_names(x))
    gg
}

#' @export
print.gspec <- function(x, ...) {
    str(x, 1)
}

#' @export
print.gspecs <- function(x, ...) {
    str(x, 1)
}

#' @export
print.decon <- function(x, ...) {
    str(x, 1)
}

#' @export
print.decons <- function(x, ...) {
    str(x, 1)
}
