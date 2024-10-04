# Any type #####

#' @name type_checking
#' @title Check the type of an object
#' @description Functions starting with `is_*` check objects for a specific type. Function `type()` is similar to `typeof` but also returns the length of the object for logical, integer, double, complex, character, and raw vectors. For other objects, it returns the class of the object.
#' @param x The object to check.
#' @param check_class Logical indicating whether to check the class of the object.
#' @param check_contents Logical indicating whether to check the contents of the object.
#' @param check_child_classes Logical indicating whether to check the class of each element of the object.
#' @return For `type()`, a character string indicating the type of the object. For `is_*` functions, a logical indicating whether the object is of the specified type.
#' @examples
#' obj <- structure(list(a="a", b=1:5), class = "abc")
#' type(obj)   # "abc"
#' type(1:5)   # "integer(5)"
#' type(NULL)  # "NULL(0)"
#' type("Hi")  # "character(1)"
#' type(1.0)   # "numeric(1)"
#' type(NA)    # "logical(1)"
#' type(NaN)   # "numeric(2)"
NULL

#' @export
#' @rdname type_checking
type <- function(x) {
    typ <- typeof(x)
    if (typ %in% c("logical", "integer", "double", "complex", "character", "raw")) {
        sprintf("%s(%d)", typ, length(x))
    } else {
        class(x)
    }
}

# S3 Singleton #####

#' @export
#' @rdname type_checking
is_spectrum <- function(x,
                        check_class = TRUE,
                        check_contents = FALSE) {
    # styler: off
    if (check_class && !inherits(x, "spectrum")) return(FALSE)
    if (!check_contents) return(TRUE)
    if (!is.list(x)) return(FALSE)
    mandatory <- c("si", "cs")
    if (!all(mandatory %in% names(x))) return(FALSE)
    # styler: on
    return(TRUE)
}

#' @export
#' @rdname type_checking
is_gspec <- function(x) inherits(x, "gspec")

#' @export
#' @rdname type_checking
is_gdecon <- function(x) inherits(x, "gdecon")

#' @export
#' @rdname type_checking
is_decon0 <- function(x) {
    is.list(x) && all(decon0_members %in% names(x))
}

#' @export
#' @rdname type_checking
is_decon1 <- function(x) inherits(x, "decon1")

#' @export
#' @rdname type_checking
is_decon2 <- function(x) inherits(x, "decon2")

# S3 Collection #####

#' @export
#' @rdname type_checking
is_spectra <- function(x,
                       check_class = TRUE,
                       check_contents = FALSE,
                       check_child_classes = FALSE) {
    # styler: off
    if (check_class && !inherits(x, "spectra")) return(FALSE)
    if (check_child_classes && !all(sapply(x, is_spectrum))) return(FALSE)
    if (!check_contents) return(TRUE)
    if (!is.list(x)) return(FALSE)
    if (!all(sapply(x, is_spectrum, check_contents = TRUE))) return(FALSE)
    # styler: on
    return(TRUE)
}

#' @export
#' @rdname type_checking
is_gspecs <- function(x) inherits(x, "gspecs")

#' @export
#' @rdname type_checking
is_gdecons <- function(x) inherits(x, "gdecons")

#' @export
#' @rdname type_checking
is_decons0 <- function(x) all(sapply(x, is_decon0))

#' @export
#' @rdname type_checking
is_decons1 <- function(x) inherits(x, "decons1")

#' @export
#' @rdname type_checking
is_decons2 <- function(x) inherits(x, "decons2")

# String #####

is_str <- function(x) {
    is.character(x) && length(x) == 1
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

is_existing_path <- function(x) {
    is_str(x) && file.exists(x)
}

# Vector #####

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

# List #####

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

is_equal <- function(x, y) {
    isTRUE(all.equal(x, y))
}