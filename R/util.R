#' @import grDevices
#' @import stats
#' @import utils

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

msg <- function(..., sep = " ", appendLF = TRUE) {
    message(paste(..., sep = sep), appendLF = appendLF)
}

msgf <- function(fmt, ..., appendLF = TRUE) {
    message(sprintf(fmt, ...), appendLF = appendLF)
}

normPath <- function(path, winslash = "/", mustWork = FALSE) {
    normalizePath(path, winslash = winslash, mustWork = mustWork)
}


str2 <- function(...) {
    capture.output(str(...))
}

#' @title Get numeric input from user
#' @description Prompts the user for input and checks if the input is a number between a minimum and maximum value. If the input is not valid, it keeps asking the user for input until they provide a valid response.
#' @param prompt The prompt to display to the user.
#' @param min The minimum valid value. Default is -Inf.
#' @param max The maximum valid value. Default is Inf.
#' @param int Whether the input should be an integer. Default is FALSE.
#' @return The user's input as a numeric value.
#' @examples \dontrun{
#' get_num_input("Enter a number between 1 and 10: ", min = 1, max = 10)
#' }
#' @noRd
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

#' @title Get string input from user
#' @description Prompts the user for input and checks if the input is in a list of valid responses. If the input is not valid, it keeps asking the user for input until they provide a valid response.
#' @param prompt The prompt to display to the user.
#' @param valid A vector of valid responses.
#' @return The user's input.
#' @examples \dontrun{
#' get_str_input("Enter a, b or c: ", c("a", "b", "c"))
#' }
#' @noRd
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

#' @title Get yes/no input from user
#' @description Prompts the user for input until they enter either 'y' or no 'n'. Returns TRUE if the user entered 'y' and FALSE if they entered 'n'.
#' @param prompt The prompt to display to the user.
#' @return TRUE if the user entered 'y' and FALSE if they entered 'n'.
#' @examples \dontrun{
#' show_dir <- get_yn_input("List dir content? (y/n) ")
#' if (show_dir) print(dir())
#' }
#' @noRd
get_yn_input <- function(prompt) {
    if (!grepl("(y/n)", prompt, fixed = TRUE)) {
        prompt <- paste0(prompt, " (y/n) ")
    }
    x <- get_str_input(prompt, c("y", "n"))
    y <- if (x == "y") TRUE else FALSE
    return(y)
}

#' @title Collapse a vector into a string
#' @description Collapses a vector into a single string, with elements separated by a specified separator. Essentially a shorthand for `paste(x, collapse = sep)`.
#' @param x A vector to collapse.
#' @param sep A string to use as the separator between elements. Default is ", ".
#' @return A single string with elements of x separated by sep.
#' @examples \dontrun{
#' collapse(c("a", "b", "c")) # "a, b, c"
#' collapse(1:5, sep = "-") # "1-2-3-4-5"
#' }
#' @noRd
collapse <- function(x, sep = ", ") {
    paste(x, collapse = sep)
}

#' @title Convert Parts per Million (ppm) to Data Points (dp)
#' @description This function converts parts per million (ppm) to data points (dp) for a given spectrum.
#' @param ppm A numeric vector of parts per million (ppm) values.
#' @param spectrum A list containing the spectrum data as returned by [load_bruker_spectrum()] or [load_jcampdx_spectrum()].
#' @param bwc Use the old, slightly incorrect method for conversion to maintain backwards compatibility with MetaboDecon1D results? For details see issue `Check: ppm to dp conversion` in TODOS.md
#' @noRd
ppm_to_dp <- function(ppm, spectrum, bwc = TRUE) {
    dp <- if (bwc) {
        # (ppm - spectrum$ppm_min) / spectrum$ppm_nstep + 1
        (spectrum$n + 1) - (spectrum$ppm_max - ppm) / spectrum$ppm_nstep
    } else {
        (ppm - spectrum$ppm_min) / spectrum$ppm_step
    }
    return(dp)
}

#' @title Convert Parts per Million (ppm) to Scaled Data Points (dp)
#' @description This function converts parts per million (ppm) to scaled data points (sdp) for a given spectrum.
#' @inheritParams ppm_to_dp
#' @noRd
ppm_to_sdp <- function(ppm, spectrum, bwc = TRUE, sf = 1000) {
    dp <- ppm_to_dp(ppm, spectrum, bwc)
    sdp <- dp / sf
    return(sdp)
}


#' @title Convert Data Points (dp) to Parts per Million (ppm)
#' @description This function converts data points (dp) to parts per million (ppm) for a given spectrum.
#' @param dp A numeric vector of data point (dp) values.
#' @param spectrum A list containing the spectrum data as returned by [load_bruker_spectrum()] or [load_jcampdx_spectrum()].
#' @noRd
dp_to_ppm <- function(dp, spectrum) {
    ppm <- spectrum$ppm_min + dp * spectrum$ppm_step
    return(ppm)
}

vcomp <- function(v1, v2) {
    if(!is.vector(v1) || is.list(v1)) stop("v1 must be a vector, but is:", class(v1))
    if(!is.vector(v2) || is.list(v2)) stop("v2 must be a vector, but is:", class(v2))
    if (identical(v1, v2)) {
        cat2("Identical")
    } else {
        capture.output(all_equal <- isTRUE(all.equal(v1, v2)))
        cat2(if (all_equal) "All equal" else "Different")
        x <- which(v1 != v2)
        d <- v1 - v2
        cat2("Objects:")
        cat2("length(v1):", length(v1))
        cat2("length(v2):", length(v2))
        cat2("head(v1):", head(v1))
        cat2("head(v2):", head(v2))
        cat2("tail(v1):", tail(v1))
        cat2("tail(v2):", tail(v2))
        cat2("class(v1):", class(v1))
        cat2("class(v2):", class(v2))
        cat2("attributes(v1):", attributes(v1))
        cat2("attributes(v2):", attributes(v2))
        cat2()
        cat2("Different Elements:")
        cat2("x <- which(v1 != v2)")
        cat2("length(x):", length(x))
        cat2("head(x)", head(x))
        cat2("tail(x)", tail(x))
        cat2()
        cat2("Differences")
        cat2("d <- v1 - v2")
        cat2("v1[head(x)]", v1[head(x)])
        cat2("v2[head(x)]", v2[head(x)])
        cat2("v1[tail(x)]", v1[tail(x)])
        cat2("v2[tail(x)]", v2[tail(x)])
        cat2("d[head(x)]", d[head(x)])
        cat2("d[tail(x)]", d[tail(x)])
    }
}

rescale <- function(x, range = c(0, 1)) {
    x_min <- min(x, na.rm = TRUE)
    x_max <- max(x, na.rm = TRUE)
    x_len <- x_max - x_min
    y_len <- max(range) - min(range)
    x_percent <- (x - x_min) / x_len
    y <- min(range) + x_percent * y_len
    y
}
