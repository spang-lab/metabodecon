# Public API #####

#' @export
#' @family {spectra functions}
#' @inherit read_spectrum
#' @title Test, Convert or Read Spectra Files from Disk
#' @description
#' `read_spectra()` reads spectra files from disk.
#' `is_spectra()` tests if an object is a `spectra` object.
#' `as_spectra()` convert an object to a `spectra` object.
#' For details about `spectra` objects, see [metabodecon_classes].
#' @return
#' `read_spectra()` and `as_spectra` return a `spectra` object.
#' `is_spectra()` returns a logical.
read_spectra <- function(data_path = pkg_file("example_datasets/bruker/urine"),
                         file_format = "bruker",
                         expno = 10,
                         procno = 10,
                         raw = FALSE,
                         silent = TRUE,
                         force = FALSE) {
    if (!file_format %in% c("bruker", "jcampdx")) {
        stop("Argument `file_format` should be either 'bruker' or 'jcampdx'")
    }
    dp <- normPath(data_path)
    jcampdx <- file_format == "jcampdx"
    bruker <- file_format == "bruker"
    r1_path <- file.path(dp, expno, "pdata", procno, "1r")
    r1_path_exists <- file.exists(r1_path)
    ends_with_dx <- grepl("\\.dx$", dp)
    if ((jcampdx && ends_with_dx) || (bruker && r1_path_exists)) {
        files <- basename(dp)
        paths <- dp
    } else if (jcampdx) {
        files <- dir(dp, pattern = "\\.dx$") # `.dx` files inside `path`
        paths <- dir(dp, pattern = "\\.dx$", full.names = TRUE)
    } else if (bruker) {
        files <- list.dirs(dp, recursive = FALSE, full.names = FALSE)
        paths <- list.dirs(dp, recursive = FALSE, full.names = TRUE) # folders inside `path`
        r1_paths <- file.path(paths, expno, "pdata", procno, "1r")
        r1_paths_exists <- file.exists(r1_paths)
        paths <- paths[r1_paths_exists]
        files <- files[r1_paths_exists]
    }
    if (length(files) == 0) {
        msg <- sprintf("No spectra found in directory '%s'.", data_path)
        if (file_format == "bruker") msg <- paste(msg, "Did you specify the correct `expno` and `procno`?")
        stop(msg)
    }
    spectra <- lapply(paths, function(path) {
        if (!silent) logf("Reading spectrum %s", path)
        read_spectrum(path, file_format, expno, procno, raw, silent, force)
    })
    names(spectra) <- files
    class(spectra) <- "spectra"
    invisible(spectra)
}

#' @export
#' @family {spectra functions}
#' @rdname read_spectra
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
#' @family {spectra functions}
#' @rdname read_spectra
#' @inheritParams read_spectra
#' @param ... Additional `spectrum` objects to put into the resulting `spectra` object..
#' @param nams Names to overwrite the individual spectrum names with. Must be NULL or a character vector of same length as `list(x, ...)`. If not provided, the names of the individual spectrum objects will be used, if available, else a default name of the form "spectrum_%d" where "%d" equals the spectrum number.
as_spectra <- function(x,
                       ...,
                       nams = NULL,
                       file_format = "bruker",
                       expno = 10,
                       procno = 10,
                       raw = FALSE,
                       silent = TRUE,
                       force = FALSE) {
    if (is_spectrum(x)) {
        xx <- structure(list(...), class = "spectra")
        ss <- set_names(xx, nams %||% get_names(xx))
    } else if (is.character(x)) {
        ss <- read_spectra(x, file_format, expno, procno, raw, silent, force)
    }
    ss
}

# Helpers #####

set_names <- function(x, nams) {
    if (!is.list(x)) stop("Input must be a list.")
    has_names <- all(sapply(x, function(e) "name" %in% names(e)))
    has_meta_names <- all(sapply(x, function(e) "name" %in% names(e$meta)))
    names(x) <- nams
    if (has_names) for (i in seq_along(x)) x[[i]]$name <- nams[[i]]
    if (has_meta_names) for (i in seq_along(x)) x[[i]]$meta$name <- nams[[i]]
    x
}

get_names <- function(x, default = "spectrum_%d") {
    stopifnot(class(x)[1] %in% c("gdecons", "gspecs", "spectra"))
    dn <- get_default_names(x, default)
    en <- names(x) # Element name
    sn <- sapply(x, function(s) s$meta$name %||% s$name) # Spectrum name
    sapply(seq_along(x), function(i) sn[i] %||% en[i] %||% dn[i])
}

get_default_names <- function(x, default) {
    if (length(default) == 1 && grepl("%d", default))
        return(sprintf(default, seq_along(x)))
    if (length(unique(default)) == length(x))
        return(default)
    stop("Default names must be a single string with a `%d` placeholder or a character vector of same length as the spectra object.")
}

#' @export
print.spectra <- function(xx) {
    msg <- "spectra object consisting of %d spectrum objects:\n"
    catf(msg, length(xx))
    nams <- get_names(xx)
    msg <- "%s (%d datapoints from %.2f - %.2f ppm)\n"
    mapply(xx, nams, FUN = function(x, nam) {
        catf(msg, nam, length(x$si), min(x$cs), max(x$cs))
    })
    invisible(NULL)
}
