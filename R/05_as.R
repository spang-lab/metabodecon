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

#' @noRd
#' @title create backwards compatible return list
#' @param spec Deconvoluted spectrum as returned by [refine_lorentz_curves_v12()].
#' @param nfiles Number of deconvoluted spectrum.
#' @param nam Name of current spectrum.
#' @return The input spectrum with an additional list element `ret` containing the deconvolution results in a backwards compatible format.
as_decon2 <- function(x) {
    if (!is_gdecon(x)) stop("Input must be a gdecon object, not ", class(x))

    # Prepare shortcuts to access the data
    sdp <- x$sdp
    ppm <- x$ppm
    hz <- x$hz
    dp <- x$dp
    y_raw <- x$y_raw
    y_smooth <- x$y_smooth
    A <- x$lcr$A
    lambda <- x$lcr$lambda
    x_0 <- x$lcr$w

    # Calculate spectrum superposition (takes approx. 2.2 seconds for urine_1)
    s <- sapply(sdp, function(x_i) sum(abs(A * (lambda / (lambda^2 + (x_i - x_0)^2)))))
    s_normed <- s / sum(s)

    # Calculate MSE_normed and MSE_normed_raw
    y_normed <- y_smooth / sum(y_smooth)
    y_raw_normed <- y_raw / sum(y_raw)
    mse_normed <- mean((y_normed - s_normed)^2)
    mse_normed_raw <- mean((y_raw_normed - s_normed)^2)

    # Create and return list
    structure(class = "decon2", .Data = list(
        number_of_files = 1,
        filename = x$name,
        x_values = x$sdp,
        x_values_ppm = x$ppm,
        y_values = x$y_smooth,
        spectrum_superposition = s,
        mse_normed = mse_normed,
        index_peak_triplets_middle = x$peak$center[x$peak$high],
        index_peak_triplets_left = x$peak$right[x$peak$high],
        index_peak_triplets_right = x$peak$left[x$peak$high],
        peak_triplets_middle = x$ppm[x$peak$center[x$peak$high]],
        peak_triplets_left = x$ppm[x$peak$right[x$peak$high]],
        peak_triplets_right = x$ppm[x$peak$left[x$peak$high]],
        integrals = x$lcr$integrals,
        signal_free_region = c(x$sfr$left_sdp, x$sfr$right_sdp),
        range_water_signal_ppm = x$wsr$hwidth_ppm,
        A = x$lcr$A,
        lambda = x$lcr$lambda,
        x_0 = x$lcr$w,
        # Fields only available in `decon2`, but not `decon1` (since v1.2.0)
        y_values_raw = x$y_raw,
        x_values_hz = x$hz,
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
    ))
}

as_decons2 <- function(x) {
    if (!is_gdecons(x)) stop("Input must be a gdecons object, not ", class(x))
    decons2 <- lapply(gdecons, as_decon2)
    names(decons2) <- names(x)
    class(decons2) <- "decons2"
    decons2
}
