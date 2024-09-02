#' @export
#' @title Create a Spectrum Object
#' @description Creates a spectrum object from the provided signal intensities, frequencies and chemical shifts.
#' @param si Numeric vector of signal intensities, ordered from highest to lowest corresponding chemical shift.
#' @param cs_max The highest chemical shift value in ppm, usually shown as left end of the spectrum.
#' @param cs_width The width of the spectrum in ppm.
#' @param fq_ref The reference frequency in Hz.
#' @param fq_width The width of the spectrum in Hz. Only used to check whether the values calculated from `cs_max`, `cs_width` and `fq_ref` match the provided value. If `NULL`, this check will be skipped.
#' @param force If `TRUE`, the function will not raise an error in case of discrepancies between the calculated and the provided spectrum width in Hz, but will print a info message instead. To hide this message as well, set `silent = TRUE`.
#' @param silent If `TRUE`, no output will be printed to the console.
#' @param name The name of the spectrum, e.g. "Blood 1" or "Urine Mouse X23D".
#' @param path The path to the spectrum file, e.g. "/example_datasets/bruker/urine/urine_1".
#' @param type The type of experiment, e.g. "H1 CPMG" or "H1 NOESY".
#' @param mfs The magnetic field strength in Tesla.
#' @return A list with class `spectrum` and the following elements:
#' - `si`: Signal intensities in Arbitrary Units.
#' - `cs`: Chemical shifts in ppm.
#' - `fq`: Frequencies in Hz.
#' - `name`: Name of the spectrum.
#' - `path`: Path to the spectrum file.
#' - `type`: Type of the spectrum.
#' - `mfs`: Magnetic field strength in Tesla.
#' @examples
#' si <- c(1, 1, 3, 7, 8, 3, 8, 5, 2, 1)
#' cs_max <- 14.8
#' cs_width <- 20.0
#' fq_ref <- 600.25 * 1e6
#' fq_width <- 12005
#' spectrum <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width)
#' spectrum2 <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width = 12010, force = FALSE)
make_spectrum <- function(si,
                          cs_max,
                          cs_width,
                          fq_ref,
                          fq_width = NULL,
                          force = FALSE,
                          silent = FALSE,
                          name = NULL,
                          path = NULL,
                          type = NULL,
                          mfs = NULL) {
    cs_min <- cs_max - cs_width # Lowest ppm value
    cs <- seq(cs_max, cs_max - cs_width, length.out = length(si)) # Chemical shift in parts per million
    fq_max <- fq_ref - (cs_min * 1e-6 * fq_ref)  # Highest frequency in Hz (corresponds to lowest ppm value)
    fq_min <- fq_ref - (cs_max * 1e-6 * fq_ref)  # Lowest frequency in Hz
    fq <- seq(fq_min, fq_max, length.out = length(si)) # Frequency in Hz
    fq_width_calc <- fq_max - fq_min
    if (!is.null(fq_width) && !isTRUE(all.equal(fq_width_calc, fq_width))) { # Check if calculated spectrum width in Hz matches the value provided by the user
        if (force) {
            stop(sprintf("Calculated spectrum width in Hz (%s) does not match the provided value (%s). Please read in the data manually or set `force = TRUE` to ignore this error. Please note that by doing so, all downstream calculations involving frequencies might be wrong, so be sure to double check the results.", round(fq_width_calc, 5), round(fq_width, 5)))
        } else if (!silent) {
            logf(sprintf("Calculated spectrum width in Hz (%s) does not match the provided value (%s). Continuing anyways, because `force` equals `TRUE`. Please note that all downstream calculations using frequencies might be wrong, so be sure to double check the results.", round(fq_width_calc, 5), round(fq_width, 5)))
        }
    }
    structure(named(si, cs, fq, name, path, type, mfs), class = "spectrum")
}

#' @export
#' @title Check whether an Object is a Spectrum
#' @description An object is considered a 'Spectrum', if it
#'
#' - is a list
#' - has all 'Mandatory Elements' listed below
#' - fulfills all 'Constraints' listed below
#'
#' # Mandatory Elements
#'
#' - `si`: Measured signal intensities in arbitrary units (au)
#' - `cs`: Corresponding "chemical shifts" in parts per pillion (ppm)
#' - `fq`: Corresponding frequencies in Hertz (Hz)
#'
#' # Optional Elements
#'
#' - `name`: Name of the spectrum, e.g. "Blood 1" or "Urine Mouse X23D"
#' - `path` are character vectors of length 1, e.g. "/example_datasets/bruker/urine/urine_1"
#' - `type`: Type of the spectrum, e.g. "H1 CPMG" or "H1 NOESY"
#' - `mfs`: Magnetic field strength in Tesla
#'
#' # Constraints
#'
#' - Elements `si`, `fq` and `cs` must have the same length.
#' - Each of the elements `name`, `type` and `path` must be a character vector of length 1 or NULL.
#' - Element `mfs` must be a numeric vector of length 1 or NULL.
#'
#' @examples
#' spectrum <- read_spectrum()
#' is.spectrum(spectrum)
is.spectrum <- function(xx, check.class = TRUE, check.contents = FALSE) {
    if (check.class) {
        if (!inherits(xx, "spectrum")) return(FALSE)
    }
    if (check.contents) {
        if (!is.list(xx)) return(FALSE)
        mandatory <- c("si", "fq", "cs")
        if (!all(mandatory %in% names(xx))) return(FALSE)
        lengths <- sapply(xx[mandatory], length)
        if (length(unique(lengths)) != 1) return(FALSE)
        optional <- c("name", "type", "path", "mfs")
        ok <- sapply(xx[optional], function(x) is.null(x) || (is.character(x) && length(x) == 1))
        if (!all(ok)) return(FALSE)
    }
    return(TRUE)
}

#' @export
#' @title Print a Spectrum Object
#' @description Prints the name, path, type, magnetic field strength, number of data points, chemical shifts and frequencies of a spectrum object.
#' @param x A spectrum object as returned by [make_spectrum()].
#' @examples
#' si <- c(1, 1, 3, 7, 8, 3, 8, 5, 2, 1)
#' cs_max <- 14.8
#' cs_width <- 20.0
#' fq_ref <- 600.25 * 1e6
#' fq_width <- 12005
#' spectrum <- read_spectrum()
#' print.spectrum(spectrum)
print.spectrum <- function(x, ...) {
    catf("Spectrum with %d data points\n", length(x$si))
    catf("- Signal Intensity Range: %.1f - %.1f\n", min(x$si), max(x$si))
    catf("- Frequency Range: %.1f - %.1f Hz\n", min(x$fq), max(x$fq))
    catf("- Chemical Shift Range: %.4f - %.4f ppm\n", max(x$cs), min(x$cs))
    catf("- Magnetic Field Stregth: %s\n", if (is.null(x$mfs)) "NULL" else paste(x$mfs, "T"))
    catf("- Name: %s\n", x$name %||% "NULL")
    catf("- Path: %s\n", x$path %||% "NULL")
    catf("- Type: %s\n", x$type %||% "NULL")
}



# Conversion functions
as.spectrum <- function(x) {
    UseMethod("as.my_class_a")
}
