
#' @rdname metabodecon1d_spectra
#'
#' @title MetaboDecon1D Spectra Objects
#'
#' @description
#' Each `metabodecon1d_spectrum` object represents a deconvoluted 1D NMR spectrum and can be created by calling [MetaboDecon1D()] with a specific `filepath`.
#' The structure of a `metabodecon1d_spectrum` object is described in the details section.
#' To test whether an object is a `metabodecon1d_spectrum` object, use [is_metabodecon1d_spectrum()].
#' Each `metabodecon1d_spectra` object is list of `metabodecon1d_spectrum` objects and can be created by calling [MetaboDecon1D()] without specifying `filepath`.
#' Testing for `metabodecon1d_spectra` objects can be done using [is_metabodecon1d_spectra()].
#'
#' @details
#' Each `metabodecon1d_spectrum` object is a list with the following elements: \loadmathjax
#'
#' - `filename`: Name of the analyzed spectrum.
#' - `x_values`: Scaled datapoint numbers (SDP). Datapoints are numbered in descending order going from N to 0, where N equals the . Scaled data point numbers are obtained by dividing these numbers by the scale factor of the x-axis. I.e., for a spectrum with 131072 datapoints and a scale factor of 1000, the first scale datapoint has value 131.071 and the last one has value 0.
#' - `x_values_ppm`: The chemical shift of each datapoint in ppm (parts per million).
#' - `y_values`: The scaled signal intensity (SSI) of each datapoint. Obtained by reading the raw intensity values from the provided `data_path` as integers and dividing them scale factor of the y-axis.
#' - `spectrum_superposition`: Scaled signal intensity obtained by calculating the sum of all estimated Lorentz curves for each data point.
#' - `mse_normed`: Normalized mean squared error. Calculated as \mjeqn{\frac{1}{n} \sum_{i=1}^{n} (z_i - \hat{z}_i)^2}{1/n * sum((z_i - zhat_i)^2)} where \mjeqn{z_i}{z_i} is the normalized, smoothed signal intensity of data point i and \mjeqn{\hat{z}_i}{zhat_i} is the normalized superposition of Lorentz curves at data point i. Normalized in this context means that the vectors were scaled so the sum over all data points equals 1.
#' - `peak_triplets_middle`: Chemical shift of peak centers in ppm.
#' - `peak_triplets_left`: Chemical shift of left peak borders in ppm.
#' - `peak_triplets_right`: Chemical shift of right peak borders in ppm.
#' - `index_peak_triplets_middle`: Datapoint numbers of peak centers.
#' - `index_peak_triplets_left`: Datapoint numbers of left peak borders.
#' - `index_peak_triplets_right`: Datapoint numbers of right peak borders.
#' - `integrals`: Integrals of the Lorentz curves.
#' - `signal_free_region`: Borders of the signal free region of the spectrum in scaled datapoint numbers. Left of the first element and right of the second element no signals are expected.
#' - `range_water_signal_ppm`: Half width of the water signal in ppm. Potential signals in this region are ignored.
#' - `A`: Amplitude parameter of the Lorentz curves as negated values. The area under the Lorentz curve is calculated as \mjeqn{A \cdot \pi}{A * pi}.
#' - `lambda`: Half width of the Lorentz curves in scaled data points as negated values. Example: a value of -0.00525 corresponds to 5.25 data points. With a spectral width of 12019 Hz and 131072 data points this corresponds to a halfwidth at half height of approx. 0.48 Hz. The corresponding calculation is: (12019 Hz / 131071 dp) * 5.25 dp.
#' - `x_0`: Center of the Lorentz curves in scaled data points.
#'
#' @examples
#' sim_subset_dir <- metabodecon_file("sim_subset")
#' ewobj <- evalwith(
#'      answers = c(10, 10, "y", 1, "n", 3.58, 3.42, "y", "n", 0, "y", "n"),
#'      expr = MetaboDecon1D(sim_subset_dir)
#' )
#' deconvs <- ewobj$rv
#' is_metabodecon1d_spectra(deconvs)
#' is_metabodecon1d_spectrum(deconvs[[1]])
is_metabodecon1d_spectrum <- function(x) {
    keys <- c(
        "number_of_files",
        "filename",
        "x_values",
        "x_values_ppm",
        "y_values",
        "spectrum_superposition",
        "mse_normed",
        "index_peak_triplets_middle",
        "index_peak_triplets_left",
        "index_peak_triplets_right",
        "peak_triplets_middle",
        "peak_triplets_left",
        "peak_triplets_right",
        "integrals",
        "signal_free_region",
        "range_water_signal_ppm",
        "A",
        "lambda",
        "x_0"
    )
    if (is.list(x) && all(keys %in% names(x))) TRUE else FALSE
}

#' @rdname metabodecon1d_spectra
#' @export
is_metabodecon1d_spectra <- function(x) {
    if (is.list(x) && all(sapply(x, is_metabodecon1d_spectrum))) TRUE else FALSE
}
