# S3 Methods for Singletons #####

#' @name print_methods
#' @title S3 Methods for Printing Metabodecon Objects
#' @description S3 Methods for printing metabodecon objects as described in the [Metabodecon Classes](https://spang-lab.github.io/metabodecon/articles/).
#' @param x The object to print.
#' @param name Logical. If TRUE, the name of the object is printed before the object.
#' @param ... Not used. Only accepted to comply with generic [base::print()].
#' @examples
#' si <- c(1, 1, 3, 7, 8, 3, 8, 5, 2, 1)
#' cs_max <- 14.8
#' cs_width <- 20.0
#' fq_ref <- 600.25 * 1e6
#' fq_width <- 12005
#' spectrum <- read_spectrum()
#' print(spectrum)
NULL

#' @export
#' @rdname print_methods
print.spectrum <- function(x, name = FALSE, ...) {
    namestr <- if (name) paste0(x$meta$name %||% "NULL", ": ") else ""
    fmt <- "%sspectrum object (%d dp, %.1f to %.1f ppm)\n"
    catf(fmt, namestr, length(x$cs), max(x$cs), min(x$cs))
}

#' @export
#' @rdname print_methods
print.ispec <- function(x, name = FALSE, ...) {
    # fmt <- "%sispec object (%d dp, %.1f to %.1f ppm)\n"
    # namestr <- if (name) paste0(x$name %||% "NULL", ": ") else ""
    # catf(fmt, namestr, length(x$ppm), max(x$ppm), min(x$ppm))
    str(x, 1)
}

#' @export
#' @rdname print_methods
print.idecon <- function(x, name = FALSE, ...) {
    # ppm <- x$ppm
    # n <- length(ppm)
    # name <- if (name) paste0(x$name %||% "NULL", ": ") else ""
    # fmt <- "%sidecon object (%d dp, %.1f to %.1f ppm, %d peaks)\n"
    # catf(fmt, name, n, max(ppm), min(ppm), length(x$A))
    str(x, 1)
}

#' @export
#' @rdname print_methods
print.decon1 <- function(x, name = FALSE, ...) {
    ppm <- x$x_values_ppm
    n <- length(ppm)
    name <- if (name) paste0(x$filename %||% "NULL", ": ") else ""
    fmt <- "%sdecon1 object (%d dp, %.1f to %.1f ppm, %d peaks)\n"
    catf(fmt, name, n, max(ppm), min(ppm), length(x$A))
}

#' @export
#' @rdname print_methods
print.decon2 <- function(x, name = FALSE, ...) {
    str(x, 1)
}

# S3 Methods for Collections #####

#' @export
#' @rdname print_methods
print.spectra <- function(x, ...) {
    msg <- "spectra object consisting of %d spectrum objects:\n"
    catf(msg, length(x, ...))
    nams <- get_names(x, ...)
    msg <- "%s (%d datapoints from %.2f - %.2f ppm)\n"
    mapply(x, ..., nams, FUN = function(x, nam) {
        catf(msg, nam, length(x$si), min(x$cs), max(x$cs))
    })
    invisible(NULL)
}

#' @export
#' @rdname print_methods
print.ispecs <- function(x, ...) {
    # catf("ispecs object with %s ispec elements\n", length(x))
    # invisible(sapply(x, print, name = TRUE))
    str(x, 2, give.attr = FALSE)
}

#' @export
#' @rdname print_methods
print.idecons <- function(x, ...) {
    catf("idecons object with %s idecon elements\n", length(x))
    nams <- get_names(x)
    mapply(x, nams, FUN = function(xi, nam) {
        catf("%s: ", nam)
        print(xi, ...)
    })
}

#' @export
#' @rdname print_methods
print.decons1 <- function(x, ...) {
    catf("decons1 object with %s decon1 elements\n", length(x))
    invisible(sapply(x, print, name = TRUE))
}

#' @export
#' @rdname print_methods
print.decons2 <- function(x, ...) {
    catf("decons2 object with %s decon2 elements\n", length(x))
    invisible(sapply(x, print, name = TRUE))
}

# Internal Functions #####

capture.output2 <- function(..., collapse = "\n", trim = FALSE) {
    x <- utils::capture.output(...)
    if (trim) {
        x <- sapply(x, trimws)
    }
    if (!(identical(collapse, FALSE))) {
        x <- paste(x, collapse = collapse)
    }
    return(x)
}

dput2 <- function(..., collapse = " ", trim = TRUE) {
    x <- capture.output2(dput(...), collapse = collapse, trim = trim)
    return(x)
}

str2 <- function(...) {
    capture.output(str(...))
}

#' @noRd
#' @title Collapse a vector into a string
#' @description Collapses a vector into a single string, with elements separated by a specified separator. Essentially a shorthand for `paste(x, collapse = sep)`.
#' @param x A vector to collapse.
#' @param sep A string to use as the separator between elements. Default is ", ".
#' @return A single string with elements of x separated by sep.
#' @examples
#' collapse(c("a", "b", "c")) # "a, b, c"
#' collapse(1:5, sep = "-") # "1-2-3-4-5"
#' collapse(1:5, last = " and ") # "1, 2, 3, 4 and 5"
collapse <- function(x, sep = ", ", last = NULL) {
    if (is.null(last) || (n <- length(x)) == 1) {
        txt <- paste(x, collapse = sep)
    } else {
        txt <- paste(x[-n], collapse = sep)
        txt <- paste(txt, x[n], sep = last)
    }
    txt
}

#' @noRd
#' @description Fixed copy of [toscutil::logf()]. Can be replaced with original after issue [Fix: logf ignores file and append arguments](https://github.com/toscm/toscutil/issues/10) has been fixed.
logf <- function(fmt,
                 ...,
                 file = "",
                 append = TRUE,
                 prefix = function() now_ms(usetz = FALSE, color = "\033[1;30m"),
                 sep1 = " ",
                 sep2 = "",
                 end = "\n",
                 verbose = TRUE) {
    cat(prefix(), sep1, sprintf(fmt, ...), sep2, end, sep = "", file = file, append = append)
}

stopf <- function(fmt, ..., call. = TRUE, domain = NULL) {
    stop(sprintf(fmt, ...), call. = call., domain = domain)
}

get_logv <- function(verbose) {
    if (verbose) logf else function(...) NULL
}

human_readable <- function(x, unit, fmt = "%.1f") {
    m <- max(abs(x))
    # styler: off
    if      (m >= 1e+9) {prefix <- "G"; x <- x / 1e+9}
    else if (m >= 1e+6) {prefix <- "M"; x <- x / 1e+6}
    else if (m >= 1e+3) {prefix <- "k"; x <- x / 1e+3}
    else if (m <= 1e-3) {prefix <- "m"; x <- x / 1e-3}
    else if (m <= 1e-6) {prefix <- "u"; x <- x / 1e-6}
    else if (m <= 1e-9) {prefix <- "n"; x <- x / 1e-9}
    else                {prefix <- "";  x <- x / 1e+0}
    # styler: off
    xhr <- sprintf(fmt, x)
    fmtstr <- paste0(fmt, " %s%s")
    sprintf(fmtstr, x, prefix, unit)
}
