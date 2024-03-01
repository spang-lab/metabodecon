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
    pat <- if (int) "^[0-9]+$" else "^[0-9]*\\.[0-9]+$"
    typ <- if (int) "number" else "value"
    x <- trimws(readline(prompt = prompt))
    while (!(grepl(pat, x) && as.numeric(x) >= min && as.numeric(x) <= max)) {
        message("Error. Please enter a ", typ, " between ", min, " and ", max, ".")
        x <- trimws(readline(prompt = prompt))
    }
    x <- as.numeric(x)
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
