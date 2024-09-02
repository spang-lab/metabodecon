

#' @export
#' @title RawSpectra Class
#' @description
#' Each object of this class:
#'
#' - Represents one or more one dimensional NMR spectra.
#' - Is a list of [spectrum] objects.
#'
#' Objects of this class can be created using [read_spectra()].
#'
#' @examples
#' fq <- c(0.599771, 0.691471, 0.783171, 0.874871, 0.966571) + 600248506
#' cs <- c(7.164231, 7.164079, 7.163926, 7.163773, 7.163620)
#' si <- c(100861.5, 105161.5, 109285.3, 113199.5, 116719.7)
#' si1 <- si + rnorm(5, 0, 100)
#' si2 <- si + rnorm(5, 0, 100)
#' raw_spectrum1 <- spectrum(si1, cs, fq)
#' raw_spectrum2 <- spectrum(si2, cs, fq)
#' raw_spectra <- RawSpectra(raw_spectrum1, raw_spectrum2)
#' try(RawSpectra(raw_spectrum1, 2))
#' try(RawSpectra(raw_spectrum1, 2, raw_spectrum2, "a", "b"))
RawSpectra <- function(...) {
    raw_spectra <- list(...)
    class_list <- lapply(raw_spectra, class)
    is_ok <- sapply(class_list, function(x) "spectrum" %in% x)
    nok <- which(!is_ok)
    if (length(nok) == 1) stop("Argument ", nok, " is not of class spectrum")
    if (length(nok) >= 2) stop("Arguments ", collapse(nok, last = " and "), " are not of class spectrum")
    structure(raw_spectra, class = "RawSpectra")
}

#' @export
#' @name DeconvolutedSpectrum
#' @title DeconvolutedSpectrum Class
#' @description Each object of this class:
#' - Represents a single, one-dimensional, deconvoluted NMR spectrum.
#' - Is a list with elements `A`, `lambda`, `x_0`, `w`, `y_0`, `y`, `x`, `spectrum_superposition`.
#' Objects of this class can be created using [generate_lorentz_curves()].
#' @examples
#' TODO
#' @author Tobias Schmidt
#' @noRd
is_decon_obj <- function(x) {
    keys <- c(
        "number_of_files", "filename", "x_values", "x_values_ppm",
        "y_values", "spectrum_superposition", "mse_normed", "index_peak_triplets_middle",
        "index_peak_triplets_left", "index_peak_triplets_right", "peak_triplets_middle",
        "peak_triplets_left", "peak_triplets_right", "integrals", "signal_free_region",
        "range_water_signal_ppm", "A", "lambda", "x_0"
    )
    if (is.list(x) && all(keys %in% names(x))) TRUE else FALSE
}



#' @export
#' @name DeconvolutedSpectra
#' @title DeconvolutedSpectra Class
#' @description Each object of this class:
#' - Represents one or more one-dimensional, deconvoluted NMR spectra.
#' - Is a list of [DeconvolutedSpectrum] objects.
#' Objects of this class can be created using [generate_lorentz_curves()].
#' @examples
#' TODO
is_decon_list <- function(x) {
    if (is.list(x) && all(sapply(x, is_decon_obj))) TRUE else FALSE
}

#' @export
#' @name AlignedSpectra
#' @title AlignedSpectra Class
#' @description Each object of this class:
#' - Represents a set of one or more one-dimensional, deconvoluted and aligned NMR spectra.
#' - Is a [data.frame] where element i,j gives the integral value of the signal from spectrum i and chemical shift j.
#' Objects of this class can be created using [align_spectra()].
AlignedSpectra <- function(...) {
    stop("Not implemented")
}

# Required Methods #####

# plot
# print
# summary
# is.RawSpectrum
# is.RawSpectra
# is.DeconvolutedSpectrum
# is.DeconvolutedSpectra
# is.AlignedSpectra
# is.AlignedSpectrum
