#' @import grDevices
#' @import stats
#' @import utils
#' @import parallel
#' @import stats
#' @import graphics
#' @import grDevices

#' @import toscutil

# Package envs #####

tests <- list() # test cases (filled in later files directly after functions)
penv <- as.environment(list())

# File Handling #####

#' @noRd
#' @title Calculate a checksum for all files in a directory or a single file
#' @description Calculates a checksum for each file in a specified directory or a single file. If the input is a directory, the checksums are calculated recursively, meaning that it includes files in all subdirectories. The results are returned as a named vector, where the names are the relative file paths and the values are checksums.
#' @param path The directory or file to calculate checksums for.
#' @param method The method to use for calculating the checksum. Can be "size" (default) or "md5". If "size", the function returns the file sizes. If "md5", the function returns the MD5 hashes of the files.
#' @param ignore A character vector of regular expressions. Files matching any of these regular expressions will be ignored.
#' @return A named vector with file paths as names and hashes as values. If the input is a directory, the names will be the file paths relative to the directory. If the input is a file, the name will be the file name.
#' @details By default, the "checksum" calculated for each file is just its size. This method was chosen because it is the fastest available and typically sufficient for our needs. Traditional checksum methods, such as MD5, can present issues. For instance, PDF files may yield different checksums every time they are recreated, likely due to the inclusion of timestamps or similar metadata within the file.
#' @examples
#' checksum(pkg_file("R"))
#' checksum(pkg_file("R"), method = "md5")
checksum <- function(path, method = "size", ignore = c()) {
    calc_checksum <- switch(method,
        size = function(paths) file.info(paths)$size,
        md5 = function(paths) sapply(paths, digest::digest, algo = "md5")
    )
    if (isTRUE(file.info(path)$isdir)) {
        paths <- list.files(path, recursive = TRUE, full.names = TRUE)
        nams <- list.files(path, recursive = TRUE, full.names = FALSE)
        for (pattern in ignore) {
            paths <- paths[!grepl(pattern, paths)]
            nams <- nams[!grepl(pattern, nams)]
        }
    } else {
        paths <- path
        nams <- basename(path)
    }
    structure(calc_checksum(paths), names = nams)
}

#' @noRd
#' @description Recursively create a dirctory without warnings and return its path.
mkdirs <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, showWarnings = FALSE, recursive = TRUE)
    }
    path
}

clear <- function(dir) {
    unlink(dir, recursive = TRUE, force = TRUE)
    mkdirs(dir)
}

dir.size <- function(dir) {
    files <- list.files(dir, recursive = TRUE, full.names = TRUE)
    sum(file.info(files)$size)
}

normPath <- function(path, winslash = "/", mustWork = FALSE) {
    normalizePath(path, winslash = winslash, mustWork = mustWork)
}

#' @noRd
#' @title Return path to any file within this (installed) package
#' @param file (string) Relative path to file.
#' @param ... Arguments passed on to [system.file()].
#' @return Absolute path to `file` with '/' as file separator.
#' @examples
#' pkg_file("DESCRIPTION")
#' pkg_file()
pkg_file <- function(...) {
    system.file(..., package = "metabodecon")
}

# Input #####

readline <- function(...) {
    base::readline(...) # we must have our own copy of readline in the package namespace so we can mock it in tests
}

#' @noRd
#' @title Get numeric input from user
#' @description Prompts the user for input and checks if the input is a number between a minimum and maximum value. If the input is not valid, it keeps asking the user for input until they provide a valid response.
#' @param prompt The prompt to display to the user.
#' @param min The minimum valid value. Default is -Inf.
#' @param max The maximum valid value. Default is Inf.
#' @param int Whether the input should be an integer. Default is FALSE.
#' @return The user's input as a numeric value.
#' @examples
#' if (interactive()) {
#'     get_num_input("Enter a number between 1 and 10: ", min = 1, max = 10)
#' }
get_num_input <- function(prompt, min = -Inf, max = Inf, int = FALSE) {
    pat <- if (int) "^[+-]?[ ]*[0-9]+$" else "^[+-]?[ ]*[0-9]*\\.?[0-9]+$"
    typ <- if (int) "number" else "value"
    if (!endsWith(prompt, " ")) prompt <- paste0(prompt, " ")
    x <- trimws(readline(prompt = prompt))
    while (!(grepl(pat, x) && as.numeric(x) >= min && as.numeric(x) <= max)) {
        message("Error. Please enter a ", typ, " between ", min, " and ", max, ".")
        x <- trimws(readline(prompt = prompt))
    }
    x <- as.numeric(x)
    return(x)
}

#' @inherit get_num_input
#' @noRd
get_int_input <- function(prompt, min = -Inf, max = Inf) {
    x <- get_num_input(prompt, min, max, int = TRUE)
    return(x)
}

get_str_input <- function(prompt, valid) {
    x <- readline(prompt = prompt)
    n <- length(valid)
    valid_str <- if (n == 1) valid else paste(collapse(valid[-n]), "or", valid[n])
    while (!(x %in% valid)) {
        message("Error. Please type either ", valid_str, ".")
        x <- readline(prompt = prompt)
    }
    return(x)
}

#' @noRd
#' @title Get yes/no input from user
#' @description Prompts the user for input until they enter either 'y' or no 'n'. Returns TRUE if the user entered 'y' and FALSE if they entered 'n'.
#' @param prompt The prompt to display to the user.
#' @return TRUE if the user entered 'y' and FALSE if they entered 'n'.
#' @examples
#' if (interactive()) {
#'      show_dir <- get_yn_input("List dir content? (y/n) ")
#'      if (show_dir) print(dir())
#' }
get_yn_input <- function(prompt) {
    if (!grepl("(y/n)", prompt, fixed = TRUE)) {
        prompt <- paste0(prompt, " (y/n) ")
    }
    x <- get_str_input(prompt, c("y", "n"))
    y <- if (x == "y") TRUE else FALSE
    return(y)
}

# Print #####

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
collapse <- function(x, sep = ", ") {
    paste(x, collapse = sep)
}

#' @noRd
#' @title Concatenate and print with newline
#' @param ... Arguments to be concatenated and printed.
#' @examples
#' cat2("Hello", "world!")
cat2 <- function(..., col = NULL) {
    if (!is.null(col)) cat(col)
    cat(...)
    cat("\n")
    if (!is.null(col)) cat(RESET)
}

cat3 <- function(...,
                 sep = " ",
                 prefix = .Options$metabodecon.cat3_prefix %||% "",
                 suffix = "\n") {
    x <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS2 ")
    y <- paste(..., sep = sep)
    cat(x, prefix, y, suffix, sep = "")
}

catf <- function(fmt, ..., end = "", file = "", sep = " ", fill = FALSE, labels = NULL, append = FALSE) {
    cat(sprintf(fmt, ...), end = end, file = file, sep = sep, fill = fill, labels = labels, append = append)
}

msg <- function(..., sep = " ", appendLF = TRUE) {
    # message(paste(..., sep = sep), appendLF = appendLF)
    cat(format(Sys.time()), "")
    cat(..., sep = sep)
    if (appendLF) cat("\n")
}

msgf <- function(fmt, ..., appendLF = TRUE) {
    # message(sprintf(fmt, ...), appendLF = appendLF)
    cat(format(Sys.time()), "")
    cat(sprintf(fmt, ...))
    if (appendLF) cat("\n")
}

# Misc #####

named_list <- function(...) {
    structure(list(...), names = as.character(substitute(list(...)))[-1])
}

# Compare #####

#' @noRd
#' @title Check if strings represent integer values
#' @description Tests each element of a character vector to determine if it represents an integer value. It allows for optional leading and trailing whitespace, and optional leading + or - signs.
#' @param x A character vector where each element is tested to determine if it represents an integer.
#' @return A logical vector indicating TRUE for elements that represent integer values and FALSE otherwise.
#' @examples
#' fixed_point_notation <- c("2.0", "3.", ".4", "-5.", "-.6")
#' scientific_notation <- c("0.45e+04", "66e-05", "0.2e-3", "-33.e-1")
#' floats <- c(fixed_point_notation, scientific_notation)
#' ints <- c("5", "-5", "1234")
#' words <- c("Hello", "world", "!", "It was nice seeing you", ".")
#'
#' x <- is_float_str(floats)
#' y <- is_int_str(floats)
#' if (!all(x == TRUE) && all(y == FALSE)) stop("Test failed")
#'
#' x <- is_float_str(ints)
#' y <- is_int_str(ints)
#' if (!all(x == FALSE) && all(y = TRUE)) stop("Test failed")
#'
#' x <- is_float_str(words)
#' y <- is_int_str(words)
#' if (!all(x == FALSE) && all(y = FALSE)) stop("Test failed")
is_int_str <- function(x) {
    grepl("^\\s*[+-]?\\s*[0-9]+$", x, perl = TRUE)
}

#' @noRd
#' @inherit is_int_str
#' @title Check if strings represent floating-point numbers
#' @description Tests each element of a character vector to determine if it represents a floating-point number. It handles numbers in fixed-point notation (with or without a decimal point) and scientific notation. Allows for optional leading and trailing whitespace, and optional leading + or - signs.
is_float_str <- function(x) {
    grepl(
        paste0(
            "^\\s*[+-]?", # Optional leading spaces and sign at start of string
            "(\\d+\\.\\d*([eE][+-]?\\d+)?", # 1.e-3,  1.0e-3, 1.0e+3
            "|\\.\\d+([eE][+-]?\\d+)?",     # .1e-2, .1.0e-2,  .1e+4
            "|\\d+([eE][+-]?\\d+)",         # 1e-3,   1.0e-3,   1e+3
            ")$" # End of string
        ),
        x, perl = TRUE
    )
}

is_num <- function(x, n = NULL) {
    if (is.null(n)) {
        is.numeric(x)
    } else {
        is.numeric(x) && length(x) == n
    }
}

is_list_of_nums <- function(x, nl, nv) {
    if (is.null(nl) && is.null(nv)) {
        is.list(x) && all(sapply(x, is.numeric))
    } else if (is.null(nv)) {
        is.list(x) && length(x) == nl && all(sapply(x, is.numeric))
    } else {
        is.list(x) && length(x) == nl && all(sapply(x, is_num, n = nv))
    }
}

all_identical <- function(x) {
    all(sapply(x, identical, x[[1]]))
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

# Convert #####

#' @export
#' @title Calculate the Width of a Numeric Vector
#' @description This function calculates the width of a numeric vector by computing the difference between the maximum and minimum values in the vector.
#' @param x A numeric vector.
#' @return The width of the vector, calculated as the difference between its maximum and minimum values.
#' @examples
#' vec <- c(1, 3, 5, 7, 9)
#' width(vec)
width <- function(x) {
    diff(range(x))
}

#' @export
#' @title Convert values from unit A to unit B
#' @description Converts values from unit A to unit B using a conversion factor.
#' @param xa A numeric vector with values in unit A.
#' @param wa Width of a certain interval (e.g. the spectrum width) in unit A.
#' @param wb With of the same interval in unit B.
#' @return A numeric vector of values converted from unit A to unit B.
#' @examples
#' \donttest{
#' urine_1 <- pkg_file("example_datasets/bruker/urine/urine_1")
#' deconv <- generate_lorentz_curves(urine_1, ask = FALSE)[[1]]
#' sdp <- deconv$x_values
#' sdp_width <- diff(range(sdp))
#' ppm <- deconv$x_values_ppm
#' ppm_width <- diff(range(ppm))
#' lambda_sdp <- deconv$lambda
#' lambda_ppm <- convert_width(lambda_sdp, sdp_width, ppm_width)
#' }
convert_width <- function(xa, wa, wb) {
    xa * (wb / wa)
}

#' @export
#' @title Convert positions from unit A to unit B
#' @description Converts positions from unit A to unit B using a conversion factor.
#' @param xa A numeric vector with positions in unit A.
#' @param wa Width of a certain interval (e.g. the spectrum width) in unit A.
#' @param wb With of the same interval in unit B.
#' @param x0a The position of a reference point in unit A (e.g. the most right point of the spectrum).
#' @param x0b The position of the same reference point in unit B.
#' @return A numeric vector of positions converted from unit A to unit B.
#' @examples
#' \donttest{
#' urine_1 <- pkg_file("example_datasets/bruker/urine/urine_1")
#' deconv <- generate_lorentz_curves(urine_1, ask = FALSE)[[1]]
#' sdp <- deconv$x_values
#' ppm <- deconv$x_values_ppm
#' peak9_sdp <- deconv$x_0[1]
#' peak9_ppm <- convert_pos(peak9_sdp,
#'     width(sdp), width(ppm),  # width of the spectrum in both units
#'     sdp[3], ppm[3]           # position of any reference point in both units
#' )
#' stopifnot(all.equal(peak9_ppm, deconv$x_0_ppm[9]))
#' }
convert_pos <- function(xa, ya, yb) {
    wa <- ya[2] - ya[1]
    wb <- yb[2] - yb[1]
    y0a <- y0a[1]

    x0b + (xa - x0a) * (wb / wa)
}

#' @noRd
#' @title Calculate Magnetic Field Strength
#' @description Calculates the magnetic field strength based on the measured frequencies and chemical shifts of a spectrum as well as the gyromagnetic ratio for protons. For this to work, the following is assumed:
#' 1. The spectrum is a 1H spectrum
#' 2. The chemical shift of the reference is at 0.0 ppm
#' 3. The resonance frequency of the reference equals the resonance frequency of protons
#' @param X A data frame containing the spectrum data with at least two columns: `cs` for chemical shifts and `fq` for frequencies.
#' @return The magnetic field strength in Tesla.
#' @examples
#' X = read_spectrum(pkg_file("example_datasets/bruker/urine/urine_1"))
#' calc_B(X)
calc_B <- function(X = read_spectrum()) {
    # Example:
    # |---------------------------|----------------------|
    # 14.8 ppm (max)            0.0 (ref)       -5.2 (min)
    # 600.243 MHz (min)       600.252 (ref)  600.255 (max)
    cs_maxmin <- max(X$cs) - min(X$cs)
    cs_refmin <- 0.0000000 - min(X$cs)
    ratio <- cs_refmin / cs_maxmin
    fq_maxmin <- max(X$fq) - min(X$fq)
    fq_refmax <- ratio * fq_maxmin
    fq_ref <- max(X$fq) - fq_refmax # Frequency of the reference (which is equal the frequency of a proton)
    gamma <- 2.675e8 # Gyromagnetic ratio for protons
    B <- (2 * pi * fq_ref) / gamma
    B
}
