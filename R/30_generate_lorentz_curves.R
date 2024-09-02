# Exported #####

#' @export
#' @title Generate Lorentz Curves from NMR Spectra
#' @description
#' Deconvolutes NMR spectra by modeling each detected signal within a spectrum as Lorentz Curve.
#' @param data_path Either the path to a directory containing measured NMR spectra, a dataframe as returned by [read_spectrum()], or a list of such dataframes.
#' @param file_format Format of the spectra files. Either `"bruker"` or `"jcampdx"`. Only relevant if `data_path` is a directory.
#' @param make_rds Logical or character. If TRUE, stores results as an RDS file on disk. If a character string, saves the RDS file with the specified name. Should be set to TRUE if many spectra are evaluated to decrease computation time.
#' @param expno The experiment number for the spectra files, e.g., `"10"`. Only relevant if `data_path` is a directory and `file_format` is `"bruker"`.
#' @param procno The processing number for the spectra, e.g., `"10"`. Only relevant if `data_path` is a directory and `file_format` is `"bruker"`.
#' @param nfit Number of iterations for approximating the parameters for the Lorentz curves.
#' @param wshw Half-width of the water artifact in ppm.
#' @param sfr Numeric vector with two entries: the ppm positions for the left and right border of the signal-free region of the spectrum.
#' @param smopts Numeric vector with two entries: the number of smoothing iterations and the number of data points to use for smoothing (must be odd).
#' @param delta Threshold value to distinguish between signal and noise. The higher the value, the more peaks get filtered out. The exact definition is as follows: a peak `i` gets filtered out, if his score is lower than `mu + s * delta`, where `mu` is the average peak score within the signal free region (SFR) and `s` is the standard deviation of peak scores in the SFR.
#' @param delta \loadmathjax Threshold for peak filtering. Higher values result in more peaks being filtered out. A peak is filtered if its score is below \mjeqn{\mu + \sigma \cdot \delta}{mu + s * delta}, where \mjeqn{\mu}{mu} is the average peak score in the signal-free region (SFR), and \mjeqn{\sigma}{s} is the standard deviation of peak scores in the SFR.
#' @param sf Numeric vector with two entries: the factors to scale the x-axis and y-axis.
#' @param ask Logical. Whether to ask for user input during the deconvolution process. If FALSE, the provided default values will be used.
#' @param debug Logical. Whether to return additional intermediate results for debugging purposes.
#' @param nworkers Number of workers to use for parallel processing. If `"auto"`, the number of workers will be determined automatically. If a number greater than 1, it will be limited to the number of spectra.
#' @param share_stdout Whether to share the standard output (usually your terminal) of the main process with the worker processes. Only relevant if `nworkers` is greater than 1. Note that this can cause messages from different workers to get mixed up, making the output hard to follow.
#' @param force If FALSE, the function stops with an error message if no peaks are found in the signal free region (SFR), as these peaks are required as a reference for peak filtering. If TRUE, the function instead proceeds without peak filtering, potentially increasing runtime and memory usage significantly.
#' @param verbose Logical. Whether to print log messages during the deconvolution process.
#' @return \loadmathjax
#' A list of deconvoluted spectra. Each deconvoluted spectrum is a list with the following elements:
#' * `number_of_files`: Number of deconvoluted spectra.
#' * `filename`: Name of the analyzed spectrum.
#' * `x_values`: Scaled datapoint numbers (SDP). Datapoints are numbered in descending order going from N - 1 to 0, where N equals the number of datapoints. Scaled data point numbers are obtained by dividing these numbers by the x-axis scale factor `sf[1]`. I.e., for a spectrum with 131072 datapoints and a scale factor of 1000, the first scale datapoint has value 131.071 and the last one has value 0.
#' * `x_values_ppm`: The chemical shift of each datapoint in ppm (parts per million).
#' * `y_values`: The scaled signal intensity (SSI) of each datapoint. Obtained by reading the raw intensity values from the provided `data_path` as integers and dividing them by the y-axis scale factor `sf[2]`.
#' * `spectrum_superposition`: Scaled signal intensity obtained by calculating the sum of all estimated Lorentz curves for each data point.
#' * `mse_normed`: Normalized mean squared error. Calculated as \mjeqn{\frac{1}{n} \sum_{i=1}^{n} (z_i - \hat{z}_i)^2}{1/n * sum((z_i - zhat_i)^2)} where \mjeqn{z_i}{z_i} is the normalized, smoothed signal intensity of data point i and \mjeqn{\hat{z}_i}{zhat_i} is the normalized superposition of Lorentz curves at data point i. Normalized in this context means that the vectors were scaled so the sum over all data points equals 1.
#' * `index_peak_triplets_middle`: Datapoint numbers of peak centers.
#' * `index_peak_triplets_left`: Datapoint numbers of left peak borders.
#' * `index_peak_triplets_right`: Datapoint numbers of right peak borders.
#' * `peak_triplets_middle`: Chemical shift of peak centers in ppm .
#' * `peak_triplets_left`: Chemical shift of left peak borders in ppm .
#' * `peak_triplets_right`: Chemical shift of right peak borders in ppm .
#' * `integrals`: Integrals of the Lorentz curves.
#' * `signal_free_region`: Borders of the signal free region of the spectrum in scaled datapoint numbers. Left of the first element and right of the second element no signals are expected.
#' * `range_water_signal_ppm`: Half width of the water signal in ppm. Potential signals in this region are ignored.
#' * `A`: Amplitude parameter of the Lorentz curves. Provided as negative number to maintain backwards compatibility with MetaboDecon1D. The area under the Lorentz curve is calculated as \mjeqn{A \cdot \pi}{A * pi}.
#' * `lambda`: Half width of the Lorentz curves in scaled data points. Provided as negative value to maintain backwards compatibility with MetaboDecon1D. Example: a value of -0.00525 corresponds to 5.25 data points. With a spectral width of 12019 Hz and 131072 data points this corresponds to a halfwidth at half height of approx. 0.48 Hz. The corresponding calculation is: (12019 Hz / 131071 dp) * 5.25 dp.
#' * `x_0`: Center of the Lorentz curves in scaled data points.
#' @details First, an automated curvature based signal selection is performed. Each signal is represented by 3 data points to allow the determination of initial Lorentz curves. These Lorentz curves are then iteratively adjusted to optimally approximate the measured spectrum.
#' @examples
#' # Define the paths to the example datasets we want to deconvolute:
#' # `sim_dir`: directory containing 16 simulated spectra
#' # `sim_01`: path to the first spectrum in the `sim` directory
#' # `sim_01_spec`: the first spectrum in the `sim` directory as a dataframe
#' sim_dir <- system.file("example_datasets/bruker/sim", package = "metabodecon")
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

    if (isFALSE(verbose)) {
        opts <- options(toscutil.logf.file = nullfile())
        on.exit(options(opts), add = TRUE)
    }
    # Read spectra and ask user for parameters
    if (is.character(data_path)) {
        spectra_ds <- read_spectra(data_path, file_format, expno, procno, raw = TRUE) # spectra in [deconvolute_spectra()] (DS) format
    } else if (is.data.frame(data_path)) {
        if (!all(c("cs", "si", "fq") %in% colnames(data_path))) stop("'data_path' must have columns 'si', 'cs' and 'fq'")
        spectra_ds <- list(data_path)
    } else if (is.list(data_path)) {
        if (!all(sapply(data_path, function(df) is.data.frame(df)))) stop("Invalid 'data_path' format.")
        dfs_ok <- sapply(data_path, function(df) all(c("cs", "si", "fq") %in% colnames(df)))
        if (!all(dfs_ok)) stop("Each spectrum in 'data_path' must have columns 'si', 'cs' and 'fq'.")
        spectra_ds <- data_path
    } else {
        stop("'data_path' must be a directory, a dataframe or a list of dataframes")
    }
    if (is.null(names(spectra_ds))) names(spectra_ds) <- paste0("spectrum_", seq_along(spectra_ds))
    spectra <- lapply(spectra_ds, as.glc_spectrum, sfx = sf[1], sfy = sf[2])
    adjno <- get_adjno(spectra, sfr, wshw, ask)
    spectra <- get_sfrs(spectra, sfr, ask, adjno)
    spectra <- get_wsrs(spectra, wshw, ask, adjno)

    # Deconvolute spectra
    nfiles <- length(spectra)
    nams <- names(spectra)
    if (nworkers == "auto") nworkers <- ceiling(parallel::detectCores() / 2)
    nworkers <- min(nworkers, length(spectra))
    starttime <- Sys.time()
    if (nworkers == 1) {
        logf("Starting deconvolution of %d spectra using 1 worker", nfiles)
        spectra <- lapply(seq_len(nfiles), function(i) {
            deconvolute_spectrum(nams[[i]], spectra[[i]], smopts, delta, nfit, nfiles, debug, force)
        })
    } else {
        logf("Starting %d worker processes", nworkers)
        cl <- parallel::makeCluster(nworkers, outfile = if (share_stdout) "" else nullfile())
        on.exit(parallel::stopCluster(cl), add = TRUE)
        logf("Exporting required functions and data to workers")
        exportvars <- c("logf", "fg", "deconvolute_spectrum", "nams", "spectra", "smopts", "delta", "nfit", "nfiles", "debug")
        parallel::clusterExport(cl, exportvars, envir = environment())
        logf("Starting deconvolution of %d spectra using %d workers", nfiles, nworkers)
        spectra <- parallel::parLapply(cl, seq_len(nfiles), function(i) {
            opts <- options(toscutil.logf.sep1 = sprintf(" PID %d ", Sys.getpid()))
            on.exit(options(opts), add = TRUE)
            deconvolute_spectrum(nams[[i]], spectra[[i]], smopts, delta, nfit, nfiles, debug)
        })
    }
    endtime <- Sys.time()
    duration <- endtime - starttime
    logf("Finished deconvolution of %d spectra in %s", nfiles, format(round(duration, 3)))
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
    if (nfiles == 1) ret[[1]] else ret
}

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

deconvolute_spectrum <- function(nam, spec, smopts, delta, nfit, nfiles, debug, force) {
    logf("Starting deconvolution of %s", nam)
    spec <- rm_water_signal(spec)
    spec <- rm_negative_signals(spec)
    spec <- smooth_signals(spec, reps = smopts[1], k = smopts[2])
    spec <- find_peaks(spec)
    spec <- filter_peaks(spec, delta, force)
    spec <- fit_lorentz_curves(spec, nfit)
    spec <- add_return_list(spec, nfiles, nam, debug)
    logf("Finished deconvolution of %s", nam)
    spec
}

rm_water_signal <- function(spec) {
    logf("Removing water signal")
    y <- spec$y_scaled
    left <- spec$wsr$left_dp
    right <- spec$wsr$right_dp
    y[right:left] <- 0.01 / spec$sfy # Order, i.e. `right:left` instead of `left:right`, is important here, because `right` and `left` are floats. Example: `right <- 3.3; left <- 1.4` ==> `right:left == c(3.3, 2.3)` and `left:right == c(1.4, 2.4)`.
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

#' @noRd
#' @title create backwards compatible return list
#' @param spec Deconvoluted spectrum as returned by [refine_lorentz_curves_v12()].
#' @param nfiles Number of deconvoluted spectrum.
#' @param nam Name of current spectrum.
#' @param debug Add debug info to the return list
#' @return The input spectrum with an additional list element `ret` containing the deconvolution results in a backwards compatible format.
add_return_list <- function(spec = glc(), nfiles = 1, nam = "urine_1", debug = TRUE) {
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
#' @title Calculate a superposition of Lorentz Curves
#' @description Calculates the superposition of Lorentz Curves for a set of x values and a set of parameter vectors (x0, A, lambda). The Lorentz Curve is defined as \mjeqn{A \cdot \frac{\lambda}{\lambda^2 + (x_i - x_0)^2}}.
#' @param x Numeric vector of x values.
#' @param x0 Centers of the Lorentz curves.
#' @param A Amplitudes of the Lorentz curves.
#' @param lambda Half widths at half height of the Lorentz curves.
#' @param iterover Character string specifying whether to iterate over the parameters or the x values. Must be either `"params"` or `"x"`.
#' @return Numeric vector of y values.
#' @examples lc(1:10, 4:5, 9:10, 2)
lcsp <- function(x, x0, A, lambda, iterover = "params") {
    message("TODO")
}

#' @noRd
#' @title Calculate Lorentz Curve values
#' @description Calculates the values of a Lorentz Curve for given x values. The Lorentz Curve is defined as \mjeqn{A \cdot \frac{\lambda}{\lambda^2 + (x_i - x_0)^2}}.
#' @param x Numeric vector of x values.
#' @param x0 Center of the Lorentz curve.
#' @param A Amplitude parameter of the Lorentz curve.
#' @param lambda Half width at half height of the Lorentz curve.
#' @return Numeric vector of y values.
#' @examples lc(1:10, 5, 10, 2)
lc <- function(x, x0, A, lambda) {
    # For details see formula below sentence "In physics, a three-parameter Lorentzian function is often used:" at [Wikipedia > Cauchy_distribution > #Properties_of_PDF](https://en.wikipedia.org/wiki/Cauchy_distribution#Properties_of_PDF).
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

# Private Deprecated #####

#' @title Generate Lorentz Curves from NMR Spectra
#' @description Deconvolutes NMR spectra and generates a Lorentz curve for each detected signal within a spectra.
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
    spectra <- lapply(spectra_ds, as.glc_spectrum, sfx = sf[1], sfy = sf[2])
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
        spec <- rm_water_signal(spec)
        spec <- rm_negative_signals(spec)
        spec <- smooth_signals(spec, reps = smopts[1], k = smopts[2])
        spec <- find_peaks(spec)
        spec <- filter_peaks(spec, delta)
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

#' @noRd
#' @title create backwards compatible return list
#' @param spec Deconvoluted spectrum as returned by [refine_lorentz_curves_v12()].
#' @param n Number of deconvoluted spectrum.
#' @param nam Name of current spectrum.
#' @param debug Add debug info to the return list
add_return_list_v12 <- function(spec = glc()$rv, n = 1, nam = "urine_1", debug = TRUE) {
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

# Prviate Experimental #####

filter_peaks_v13 <- function(ppm, # x values in ppm
                             pc, # peak center indices
                             ps, # peak scores
                             sfrl, # signal free region left in ppm
                             sfrr, # signal free region right in ppm
                             delta = 6.4) { # threshold parameter to distinguish between "real" peaks from noise
    if (any(is.na(ps))) stop("Peak scores must never be NA")
    logf("Removing peaks with low scores")
    in_sfr <- which(ppm[pc] >= ppm || ppm[pc] <= ppm)
    if (any(in_sfr)) {
        mu <- mean(ps[in_sfr])
        sigma <- sd(ps[in_sfr])
    } else {
        stop("No signals found in signal free region. Please double check deconvolution parameters.")
    }
    mu <- mean(ps[in_sfr]) # mean (greek letter mu)
    sigma <- sd(ps[in_sfr]) # standard deviation (greek letter sigma)
    gt_tau <- (ps >= mu + delta * sigma) # greater than threshold value (greek letter tau)
    logf("Removed %d peaks", sum(!gt_tau))
    list(in_sfr, gt_tau)
}
