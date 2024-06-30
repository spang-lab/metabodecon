#' @import grDevices
#' @import stats
#' @import utils
#' @import parallel
#' @import stats
#' @import graphics
#' @import grDevices

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

is_num <- function(x, n = NULL) {
    if (is.null(n)) {
        is.numeric(x)
    } else {
        is.numeric(x) && length(x) == n
    }
}

is_list_of_nums <- function(x, nl, nv) {
    if (is.null(nv) && is.null(nv)) {
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

#' @noRd
#' @title Convert Parts per Million (ppm) to Data Points Numbers (dpn)
#' @description Converts parts per million (ppm) to data points numbers (dpn) for a given spectrum.
#' @param ppm A numeric vector of parts per million (ppm) values.
#' @param spectrum A list containing the spectrum data as returned by [load_bruker_spectrum()] or [read_spectrum()].
#' @param bwc Use the old, slightly incorrect method for conversion to maintain backwards compatibility with MetaboDecon1D results? For details see issue `Check: ppm to dpn conversion` in TODOS.md
ppm_to_dpn <- function(ppm, spectrum, bwc = TRUE) {
    dpn <- if (bwc) {
        # (ppm - spectrum$ppm_min) / spectrum$ppm_nstep + 1
        (spectrum$n + 1) - (spectrum$ppm_max - ppm) / spectrum$ppm_nstep
    } else {
        (ppm - spectrum$ppm_min) / spectrum$ppm_step
    }
    return(dpn)
}

#' @noRd
#' @title Convert Parts per Million (ppm) to Scaled Data Point Numbers (sdpn)
#' @description Converts parts per million (ppm) to scaled data point numbers (sdpn) for a given spectrum.
#' @inheritParams ppm_to_sdpn
ppm_to_sdpn <- function(ppm, spectrum, bwc = TRUE, sf = 1000) {
    dp <- ppm_to_dpn(ppm, spectrum, bwc)
    sdp <- dp / sf
    return(sdp)
}

#' @noRd
#' @title Convert Data Points (dp) to Parts per Million (ppm)
#' @description Converts data points (dp) to parts per million (ppm) for a given spectrum.
#' @param dp A numeric vector of data point (dp) values.
#' @param spectrum A list containing the spectrum data as returned by [load_bruker_spectrum()] or [read_spectrum()].
dp_to_ppm <- function(dp, spectrum) {
    ppm <- spectrum$ppm_min + dp * spectrum$ppm_step
    return(ppm)
}

#' @title Calculate Signal Width in Hz for deconvoluted NMR Spectra
#' @description Iterates over each deconvoluted spectrum in the provided list, calculates the full signal width at half height in Hz for each signal, and updates the spectrum object with the calculated width.
#' @param spectrum_data A list of spectrum data entries, where each entry is expected to have `x_values`, `lambda` and `peak_triplets_middle` among other properties.
#' @param sw_hz The spectral width in Hz.
#' @return The modified spectrum_data list where each element has an additional entry `lambda_hz`.
#' @examples
#' \donttest{
#' xp <- download_example_datasets()
#' dp <- file.path(xp, "bruker/urine/urine_1")
#' x <- generate_lorentz_curves(dp, ask = FALSE, nfit = 3, ncores = 1)
#' dspec <- x$urine_1
#' lambda_hz <- calc_signal_width_in_hz(dspec, sw_hz = 800)
#' }
#' @noRd
calc_signal_width_in_hz <- function(dspec, sw_hz) {
    nr_dpi <- dspec$x_values[1] * 1000 # number of data point intervals (131071 (2^17 - 1) for "128k spectra")
    nr_dp <- nr_dpi + 1 # number of data points
    hw_sdp <- dspec$lambda # half width in scaled data point intervals
    hw_dp <- hw_sdp * 1000 # half width in data point intervals
    w_dp <- 2 * hw_dp # width in data points intervals
    hzni <- sw_hz / nr_dp # hz n-interval (too small, because we divide by number of points instead of number of intervals)
    hzi <- sw_hz / nr_dpi # hz interval
    lambda_hz <- abs(hzni * w_dp)
}
