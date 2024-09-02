as.glc_spectrum <- function(X, ...) {
    UseMethod("as.glc_spectrum", X)
}

#' @noRd
#' @title Convert normal Spectrum to GLC Spectrum
#' @description Takes a normal spectrum as returned by [read_spectrum()], [make_spectrum] or [simulate_spectrum] and converts it to a GLC spectrum, as expected by [generate_lorentz_curves()].
#' @param X A normal spectrum as returned by [read_spectrum()], [make_spectrum] or [simulate_spectrum].
#' @param sfx The scaling factor for the x-axis.
#' @param sfy The scaling factor for the y-axis.
#' @return A named list containing the following elements:
#' - `y_raw`: The raw signal intensities.
#' - `y_scaled`: The scaled signal intensities.
#' - `n`: The number of data points.
#' - `sfx`: The scaling factor for the x-axis.
#' - `sfy`: The scaling factor for the y-axis.
#' - `dp`: The data point numbers.
#' - `sdp`: The scaled data point numbers.
#' - `ppm`: The chemical shifts in ppm.
#' - `fq`: The frequencies in Hz.
#' - `ppm_min`: The minimum chemical shift in ppm.
#' - `ppm_max`: The maximum chemical shift in ppm.
#' - `ppm_range`: The range of the chemical shifts in ppm.
#' - `ppm_step`: The step size of the chemical shifts in ppm.
#' - `ppm_nstep`: The step size of the chemical shifts in ppm, calculated as `ppm_range / n`.
#' @examples
#' nrm_spec <- read_spectrum()
#' glc_spec <- as.glc_spectrum.spectrum(nrm_spec, 1e3, 1e6)
as.glc_spectrum.spectrum <- function(X, sfx, sfy) {
    y_raw <- X$si
    y_scaled <- y_raw / sfy
    n <- length(y_raw)
    ppm_range <- diff(range(X$cs))
    ppm_max <- max(X$cs)
    ppm_min <- min(X$cs)
    ppm_step <- ppm_range / (n - 1)
    ppm_nstep <- ppm_range / n
    # Example: data points in ppm = 1.1, 2.3, 3.5, 4.7
    # ==> ppm_step == 1.2
    # ==> ppm_nstep ~= 1.06 (not really useful, but we need it for backwards compatibility with MetaboDecon1D results)
    ppm <- X$cs # Parts per million
    hz <- X$fq # Frequency in Hz
    dp <- seq(n - 1, 0, -1) # Data point numbers
    sdp <- seq((n - 1) / sfx, 0, -1 / sfx) # Scaled data point numbers. (Same as `dp / sfx`, but with slight numeric differences, so we stick with the old calculation method for backwards compatibility)
    named(
        y_raw, y_scaled, # y-axis
        n, sfx, sfy, # misc
        dp, sdp, ppm, hz,# x-axis
        ppm_min, ppm_max, ppm_range, ppm_step, ppm_nstep # additional ppm info
    )
}
