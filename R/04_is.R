# Type #####

#' @noRd
#' @examples
#' type(1:5) # "integer(5)"
#' type(NULL) # "NULL(0)"
#' type("Hello") # "character(1)"
#' type(1.0) # "numeric(1)"
#' type(NA) # "logical(1)"
#' type(c(NaN, Inf)) # "numeric(2)"
#' type(structure(list(a="a", b=1:5), class = "abc")) # "abc"
type <- function(x) {
    typ <- typeof(x)
    if (typ %in% c("logical", "integer", "double", "complex", "character", "raw")) {
        sprintf("%s(%d)", typ, length(x))
    } else {
        class(x)
    }
}

# String checking #####

is_str_or_null <- function(x) {
    is.null(x) || (is.character(x) && length(x) == 1)
}

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
            "|\\.\\d+([eE][+-]?\\d+)?", # .1e-2, .1.0e-2,  .1e+4
            "|\\d+([eE][+-]?\\d+)", # 1e-3,   1.0e-3,   1e+3
            ")$" # End of string
        ),
        x,
        perl = TRUE
    )
}

# Vector checking #####

is_num <- function(x, n = NULL) {
    if (is.null(n)) {
        is.numeric(x)
    } else {
        is.numeric(x) && length(x) == n
    }
}

is_int <- function(x, n = NULL) {
    if (is.null(x)) return(FALSE)
    x_is_int <- is.integer(x) || (is.numeric(x) && all(abs(x - round(x)) < sqrt(.Machine$double.eps)))
    has_correct_length <- is.null(n) || length(x) == n
    x_is_int && has_correct_length
}

is_char <- function(x, n = NULL, pattern = NULL) {
    !is.null(x) && is.character(x) && (is.null(n) || length(x) == n ) && (is.null(pattern) || all(grepl(pattern, x)))
}

is_bool <- function(x, n = NULL) {
    !is.null(x) && is.logical(x) && (is.null(n) || length(x) == n)
}

# List checking #####

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

# S3 Classes #####

is_gspec <- function(x) inherits(x, "gspec")

is_gdecon <- function(x) inherits(x, "gdecon")

is_gdecons <- function(x) inherits(x, "gdecons")
