# General #####

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
                 file = .Options$toscutil.logf.file %||% "",
                 append = .Options$toscutil.logf.append %||% FALSE,
                 prefix = .Options$toscutil.logf.prefix %||% function() now_ms(usetz = FALSE, color = "\033[1;30m"),
                 sep1 = .Options$toscutil.logf.sep1 %||% " ",
                 sep2 = .Options$toscutil.logf.sep2 %||% "",
                 end = .Options$toscutil.logf.end %||% "\n") {
    cat(prefix(), sep1, sprintf(fmt, ...), sep2, end, sep = "", file = file, append = append)
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

# S3 Methods #####

#' @export
print.gspec <- function(x, ...) {
    catf("gspec object (%d dp, %.1f to %.1f ppm)\n", x$n, max(x$ppm), min(x$ppm))
}

#' @export
print.gdecon <- function(x, ...) {
    catf("gdecon object (%d dp, %.1f to %.1f ppm, %d peaks)\n", x$n, max(x$ppm), min(x$ppm), sum(x$peak$high))
}

#' @export
print.decon2 <- function(x, name = FALSE, ...) {
    ppm <- x$x_values_ppm
    n <- length(ppm)
    name <- if (name) paste0(x$filename %||% "NULL", ": ") else ""
    fmt <- "%sdecon2 object (%d dp, %.1f to %.1f ppm, %d peaks)\n"
    catf(fmt, name, n, max(ppm), min(ppm), length(x$A))
}

#' @export
print.gspecs <- function(x, ...) {
    catf("gspecs object with %s gspec elements\n", length(x))
    invisible(sapply(x, print, ...))
}

#' @export
print.gdecons <- function(x, ...) {
    catf("gdecons object with %s gdecon elements\n", length(x))
    nams <- get_names(x)
    mapply(x, nams, FUN = function(xi, nam) {
        catf("%s: ", nam)
        print(xi, ...)
    })
}

#' @export
print.decons2 <- function(x, ...) {
    catf("decons2 object with %s decon2 elements\n", length(x))
    invisible(sapply(x, print, name = TRUE, ...))
}
