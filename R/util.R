# Convert between Units (Public) ###############################################

#' @export
#' @title Convert from unit A to unit B
#'
#' @description
#' Converts positions/widths from unit A to unit B. If the direction of units  A
#' and B is reversed, the width's sign will be reversed as well. To keep  widths
#' strictly positive, wrap the result with `abs()`.
#'
#' @param xa
#' A numeric vector specifying widths/positions in unit A.
#'
#' @param ya,yb
#' A numeric vector specifying the positions of at least two points in unit A /
#' unit B.
#'
#' @return
#' A numeric vector of values converted from unit A to unit B.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' ya <- c(244, 246, 248, 250, 252)
#' yb <- c(15, 10, 5, 0, -5)
#' convert_width(c(2, 4, 8), ya, yb)
#' convert_pos(c(247, 249), ya, yb)
convert_pos <- function(xa, ya, yb) {
    # Example:
    #     |-------------w----------|
    #     |----d----|              |
    # xa  |        247             |
    # ya  244   246 | 248   250   252
    # yb  15    10  | 5     0     -5
    # ==>
    # wa  =  ya[n] - ya[1] = 252 - 244   = 8
    # wb  =  yb[n] - yb[1] = (-5) - 15   = 20
    # da  =  xa    - ya[1] = 247 - 244   = 3
    # db  =  da/wa * wb    = (3/8) * -20 = -7.5 (because da/wa == db/wb)
    # xb  =  yb[1] + db    = 15 + (-7.5) = 7.5
    if (length(ya) != length(yb)) stop("ya and yb must have the same length.")
    n <- length(ya)
    wa <- ya[n] - ya[1]
    wb <- yb[n] - yb[1]
    da <- xa - ya[1]
    db <- (da / wa) * wb
    xb <- yb[1] + db
    xb
}

#' @export
#' @rdname convert_pos
convert_width <- function(xa, ya, yb) {
    n <- length(ya)
    wa <- ya[n] - ya[1]
    wb <- yb[n] - yb[1]
    xa * (wb / wa)
}

#' @export
#' @title Calculate the Width of a Numeric Vector
#'
#' @description
#' Calculates the width of a numeric vector by computing the difference between
#' the maximum and minimum values in the vector.
#'
#' @param x
#' A numeric vector.
#'
#' @return
#' The width of the vector, calculated as the difference between its maximum and
#' minimum values.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' vec <- c(1, 3, 5, 7, 9)
#' width(vec)
width <- function(x) {
    diff(range(x))
}

#' @export
#' @title Show head and tail rows of a matrix-like object
#'
#' @description
#' Returns the first and last `n` rows of a matrix or data frame. If the input
#' has fewer than `2*n` rows, overlapping rows are returned only once.
#'
#' @param x
#' A matrix or data frame.
#'
#' @param n
#' Number of rows to take from the top and bottom.
#'
#' @return
#' A subset of `x` containing head and tail rows.
#'
#' @author 2024-2026 Tobias Schmidt: initial version.
#'
#' @examples
#' x <- matrix(seq_len(30), nrow = 10)
#' headtail(x, n = 2)
headtail <- function(x, n = 6) {
    ok <- is.matrix(x) || is.data.frame(x)
    if (!ok) stop("x must be a matrix or data frame.")
    if (!is_int(n) || length(n) != 1 || n < 0) {
        stop("n must be a single non-negative integer.")
    }

    nr <- nrow(x)
    if (nr == 0) return(x)

    n_show <- min(as.integer(n), nr)
    i_top <- seq_len(n_show)
    i_bot <- seq.int(from = nr - n_show + 1, to = nr)
    i <- unique(c(i_top, i_bot))
    x[i, , drop = FALSE]
}

# Convert between Units (Private) #############################################

#' @noRd
#' @description
#' Converts a vector of chemical shifts given in ppm to a vector of frequencies
#' in Hz.
#' @param cs Vector of chemical shifts in ppm.
#' @param fqref Frequency of the reference molecule in Hz.
#' @author 2024-2025 Tobias Schmidt: initial version.
in_hz <- function(cs, fqref) {
    fqmax <- fqref - (min(cs) * 1e-6 * fqref) # Highest frequency in Hz
    fqmin <- fqref - (max(cs) * 1e-6 * fqref) # Lowest frequency in Hz
    fq <- seq(fqmin, fqmax, length.out = length(cs))
    fq <- fq[rev(order(cs))] # (1)
    # (1) Sort frequencies in descending order of chemical shifts, as the
    #     highest chemical shift corresponds to the lowest frequency.
    fq
}

#' @noRd
#' @description
#' Convert the signal free region border (SFR) from scaled data point numbers
#' (SDP) to parts per million (PPM) in a backwards compatible way, i.e. this
#' function does the same mistakes and numerical errors as the original
#' conversion in MetaboDecon1D.
#' @param sfr_sdp SFR in SDP.
#' @param sdp All datapoints in SDP.
#' @param ppm All datapoints in ppm.
#' @details
#' See 'CHECK-2: Signal free region (SFR) calculation' in TODOS.md. (Update
#' 2025-09-14: TODOS are no longer tracked in TODOS.md, but outside of the
#' repository. To retrieve the last actively maintained version of TODOS.md,
#' checkout commit 8b1f61b, i.e., v1.5.0.)
#' @author 2024-2025 Tobias Schmidt: initial version.
sfr_in_ppm_bwc <- function(sfr_sdp, sdp, ppm) {
    assert(
        is.numeric(sfr_sdp), length(sfr_sdp) == 2,
        is.numeric(sdp), length(sdp) >= 5,
        is.numeric(ppm), length(ppm) == length(sdp)
    )
    n <- length(sdp)
    sdp_step <- diff(range(sdp)) / (n - 1)
    ppm_nstep <- diff(range(ppm)) / n
    sfr_dp <- sfr_sdp / sdp_step
    sfr_ppm <- max(ppm) - (n + 1 - sfr_dp) * ppm_nstep
    sfr_ppm
}

#' @noRd
#' @rdname sfr_in_ppm_bwc
#' @param sfr_ppm Signal free region borders (SFR) in parts per million (PPM).
#' @author 2024-2025 Tobias Schmidt: initial version.
sfr_in_sdp_bwc <- function(sfr_ppm, ppm, sf) {
    sfr_left  <- max(sfr_ppm)
    sfr_right <- min(sfr_ppm)
    ppm_left  <- max(ppm)
    ppm_right <- min(ppm)
    n <- length(ppm)
    ppm_nstep <- (ppm_left - ppm_right) / n
    dp_left   <- (n + 1) - (ppm_left - sfr_left)  / ppm_nstep
    dp_right  <- (n + 1) - (ppm_left - sfr_right) / ppm_nstep
    sdp_left  <- dp_left  / sf[1]
    sdp_right <- dp_right / sf[1]
    c(sdp_left, sdp_right)
}

# File Handling (Private) #####################################################

#' @noRd
#' @title Calculate a checksum for all files in a directory or a single file
#'
#' @description
#' Calculates a checksum for each file in a specified directory or a single
#' file. If the input is a directory, the checksums are calculated recursively,
#' meaning that it includes files in all subdirectories. The results are
#' returned as a named vector, where the names are the relative file paths and
#' the values are checksums.
#'
#' @param path The directory or file to calculate checksums for.
#'
#' @param method The method to use for calculating the checksum. Can be "size"
#' (default) or "md5". If "size", the function returns the file sizes. If "md5",
#' the function returns the MD5 hashes of the files.
#'
#' @param ignore A character vector of regular expressions. Files matching any
#' of these regular expressions will be ignored.
#'
#' @return A named vector with file paths as names and hashes as values. If the
#' input is a directory, the names will be the file paths relative to the
#' directory. If the input is a file, the name will be the file name.
#'
#' @details By default, the "checksum" calculated for each file is just its
#' size. This method was chosen because it is the fastest available and
#' typically sufficient for our needs. Traditional checksum methods, such as
#' MD5, can present issues. For instance, PDF files may yield different
#' checksums every time they are recreated, likely due to the inclusion of
#' timestamps or similar metadata within the file.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
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
#' @description
#' Recursively create a dirctory without warnings and return its path.
#' @author 2024-2025 Tobias Schmidt: initial version.
mkdirs <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, showWarnings = FALSE, recursive = TRUE)
    }
    path
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
clear <- function(dir) {
    unlink(dir, recursive = TRUE, force = TRUE)
    mkdirs(dir)
}

#' @noRd
#' @description
#' Return a no-op cache object with the same methods we use from
#' `cachem::cache_disk()`.
#' @author 2026 Tobias Schmidt: initial version.
null_cache <- function() {
    cache <- new.env(parent = emptyenv())
    cache$exists <- function(key) FALSE
    cache$get <- function(key) NULL
    cache$set <- function(key, value) invisible(NULL)
    cache
}

#' @noRd
#' @description
#' Return a disk cache unless `path` is `NULL`, in which case a no-op cache is
#' returned. The disk cache is configured without size, age, or entry-count
#' limits, and with an effectively disabled pruning throttle.
#' @author 2026 Tobias Schmidt: initial version.
disk_cache <- function(path) {
    if (is.null(path)) {
        null_cache()
    } else {
        cachem::cache_disk(
            dir = path,
            max_size = Inf,
            max_age = Inf,
            max_n = Inf,
            prune_rate = .Machine$integer.max
        )
    }
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
norm_path <- function(path, winslash = "/", mustWork = FALSE) {
    normalizePath(path, winslash = winslash, mustWork = mustWork)
}

#' @noRd
#' @title Return path to any file within this (installed) package
#' @param file (string) Relative path to file.
#' @param ... Arguments passed on to [system.file()].
#' @return Absolute path to `file` with '/' as file separator.
#' @author 2024-2025 Tobias Schmidt: initial version.
#' @examples
#' pkg_file("DESCRIPTION")
#' pkg_file()
pkg_file <- function(...) {
    system.file(..., package = "metabodecon")
}

#' @noRd
#' @title Store Object in File
#' @description Stores the object returned by the provided expression in the
#' provided path.
#' @param ... Additional arguments passed to the device function
#' @author 2024-2025 Tobias Schmidt: initial version.
#' @examples
#' tmp_png <- tempfile(fileext = ".png")
#' tmp_pdf <- tempfile(fileext = ".pdf")
#' store(plot(1:5), tmp_pdf)
#' store(plot(1:10), tmp_png, width = 10, units = "in", res= 72)
#' store(quote(plot(1:20)), quoted = TRUE)
store <- function(expr,
                  path = tempfile(fileext = ".pdf"),
                  verbose = TRUE,
                  quoted = FALSE,
                  format = tolower(tools::file_ext(path)),
                  ...) {
    path <- norm_path(path)
    devs <- c("png", "jpeg", "tiff", "bmp", "svg", "pdf")
    call <- if (quoted) eval else force
    if (verbose) logf("Writing %s", path)
    if (format == "jpg") format <- "jpeg"
    if (format == "rds") {
        saveRDS(call(expr), path)
    } else if (format %in% devs) {
        with_dev <- asNamespace("withr")[[paste0("with_", format)]]
        with_dev(path, call(expr), ...)
    } else {
        stop(sprintf("Unsupported format '%s'", format))
    }
}

# Reading User Input (Private) #####

#' @noRd
#' @description Copy of Readline
#' @details
#' We must have our own copy so we can replace if with mock functions during
#' testing. That's not (easily) possible for base functions.
#' @author 2024-2025 Tobias Schmidt: initial version.
readline <- function(...) {
    base::readline(...)
}

#' @noRd
#'
#' @title Get numeric input from user
#'
#' @description
#' Prompts the user for input and checks if the input is a number between a
#' minimum and maximum value. If the input is not valid, it keeps asking the
#' user for input until they provide a valid response.
#'
#' @param prompt The prompt to display to the user.
#' @param min The minimum valid value. Default is -Inf.
#' @param max The maximum valid value. Default is Inf.
#' @param int Whether the input should be an integer. Default is FALSE.
#'
#' @return The user's input as a numeric value.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' if (interactive()) {
#'      x <- get_num_input("Enter a number between 1 and 10: ", min = 1, max = 10)
#'      y <- get_num_input("Enter a number between 1 and 10: ", min = 1, max = 10)
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

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
get_int_input <- function(prompt, min = -Inf, max = Inf) {
    get_num_input(prompt, min, max, int = TRUE)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
get_str_input <- function(prompt, valid) {
    x <- readline(prompt = prompt)
    n <- length(valid)
    valid_str <- if (n == 1) valid else paste(collapse(valid[-n]), "or", valid[n])
    while (!(x %in% valid)) {
        message("Error. Please type either ", valid_str, ".")
        x <- readline(prompt = prompt)
    }
    x
}

#' @noRd
#'
#' @title Get yes/no input from user
#'
#' @description
#' Prompts the user for input until they enter either 'y' or no 'n'. Returns
#' TRUE if the user entered 'y' and FALSE if they entered 'n'.
#'
#' @param prompt The prompt to display to the user.
#'
#' @return
#' TRUE if the user entered 'y' and FALSE if they entered 'n'.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' if (interactive()) {
#'     show_dir <- get_yn_input("List dir content? (y/n) ")
#'     if (show_dir) print(dir())
#' }
get_yn_input <- function(prompt) {
    if (!grepl("(y/n)", prompt, fixed = TRUE)) {
        prompt <- paste0(prompt, " (y/n) ")
    }
    x <- get_str_input(prompt, c("y", "n"))
    if (x == "y") TRUE else FALSE
}

# Operators (Private) #########################################################

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
`%&&%` <- function(x, y) {
    if (is.null(x)) x else y
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
`%==%` <- function(x, y) {
    isTRUE(all.equal(x, y))
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
`%!=%` <- function(x, y) {
    !isTRUE(all.equal(x, y))
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
`%===%` <- function(x, y) {
    identical(x, y)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
`%!==%` <- function(x, y) {
    !identical(x, y)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
`%notin` <- function(x, y) {
    !(x %in% y)
}


# Print Functions (Private) #####

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
esc <- list(
    # https://en.wikipedia.org/wiki/ANSI_escape_code
    reset = "\033[0m", # Reset all formatting
    bold = "\033[1m", # Bold text
    dim = "\033[2m", # Dim text
    italic = "\033[3m", # Italic text
    underline = "\033[4m", # Underlined text
    slow_blink = "\033[5m", # Blinking text
    rapid_blink = "\033[6m", # Rapid blinking text
    reverse = "\033[7m", # Reverse video (swap foreground and background color)
    hidden = "\033[8m", # Hidden text (useful for passwords)
    strikthrough = "\033[9m", # Strikethrough text
    cursor_up = "\033[A",
    cursor_down = "\033[B",
    cursor_forward = "\033[C",
    cursor_back = "\033[D",
    cursor_next_line = "\033[E",
    cursor_prev_line = "\033[F",
    cursor_start_of_line = "\033[G",
    cursor_pos_save = "\033[s",
    cursor_pos_popgs = "\033[u",
    black = "\033[30m",
    red = "\033[31m",
    green = "\033[32m",
    yellow = "\033[33m",
    blue = "\033[34m",
    magenta = "\033[35m",
    cyan = "\033[36m",
    white = "\033[37m",
    reset = "\033[0m",
    bright_black = "\033[90m",
    bright_red = "\033[91m",
    bright_green = "\033[92m",
    bright_yellow = "\033[93m",
    bright_blue = "\033[94m",
    bright_magenta = "\033[95m",
    bright_cyan = "\033[96m",
    bright_white = "\033[97m",
    bg_black = "\033[40m",
    bg_red = "\033[41m",
    bg_green = "\033[42m",
    bg_yellow = "\033[43m",
    bg_blue = "\033[44m",
    bg_magenta = "\033[45m",
    bg_cyan = "\033[46m",
    bg_white = "\033[47m",
    bg_bright_black = "\033[100m",
    bg_bright_red = "\033[101m",
    bg_bright_green = "\033[102m",
    bg_bright_yellow = "\033[103m",
    bg_bright_blue = "\033[104m",
    bg_bright_magenta = "\033[105m",
    bg_bright_cyan = "\033[106m",
    bg_bright_white = "\033[107m"
)

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
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

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
dput2 <- function(..., collapse = " ", trim = TRUE) {
    x <- capture.output2(dput(...), collapse = collapse, trim = trim)
    return(x)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
str0 <- function(object,
                 max.level = 1,
                 give.attr = FALSE, ...) {
    str(
        object,
        max.level = max.level,
        give.attr = give.attr,
        ...
    )
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
str2 <- function(...) {
    capture.output2(str(...))
}

#' @noRd
#' @title Collapse a vector into a string
#' @description
#' Collapses a vector into a single string, with elements separated by a
#' specified separator. Essentially a shorthand for
#' `paste(x, collapse = sep)`.
#' @param x A vector to collapse.
#' @param sep A string to use as the separator between elements. Default is ", ".
#' @return A single string with elements of x separated by sep.
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' @description
#' Fixed copy of [toscutil::logf()]. Can be replaced with original after issue
#' [Fix: logf ignores file and append
#' arguments](https://github.com/toscm/toscutil/issues/10) has been fixed.
#' @author 2024-2025 Tobias Schmidt: initial version.
logf <- function(fmt,
                 ...,
                 file = getOption("toscutil.logf.file", ""),
                 append = TRUE,
                 prefix = function() now_ms(usetz = FALSE, color = "\033[1;30m"),
                 sep1 = " ",
                 sep2 = "",
                 end = "\n") {
    cat(prefix(), sep1, sprintf(fmt, ...), sep2, end, sep = "", file = file, append = append)
}

#' @noRd
catft <- function(fmt, ...) {

    now <- Sys.time()

    # Finish previous log by printing it's duration
    prev <- getOption("metabodecon.catft.time", NULL)
    dt <- if (is.null(prev)) "" else sprintf(" (%.1fs)\n", as.numeric(now - prev))
    istty <- isatty(stdout())
    if (istty) dt <- paste0(esc$bright_cyan, dt, esc$reset)
    cat(dt)

    # Start new log message with timestamp and formatted text
    txt <- sprintf(fmt, ...)
    ts <- format(now, "%Y-%m-%d %H:%M:%S")
    if (istty) txt <- paste0(esc$bright_black, txt, esc$reset)
    msg <- sprintf("%s %s", ts, txt)
    cat(msg)
    options(metabodecon.catft.time = now)
}

#' @noRd
#' @title Log formatted messages if verbosity is high enough
#' @description
#' Logs `fmt` via `logf()` if `verbosity >= 1` (`logv`) or `>= 2` (`logvv`).
#' The `verbosity` argument defaults to `caller_verbosity()`, so it is
#' automatically inherited from the calling environment when not supplied.
#' @param fmt A format string (character) or any object to print and log.
#' @param ... Values passed to `sprintf()` when `fmt` is character.
#' @param verbosity
#' Numeric verbosity level. Defaults to `caller_verbosity()`, which looks up
#' `verbosity` in the calling environment.
#' @return
#' The timestamp (POSIXct) of the log entry (like `logf()` does).
logv <- function(fmt, ..., verbosity = caller_verbosity()) {
    if (verbosity >= 1) logf(fmt, ...) else invisible(Sys.time())
}

logvv <- function(fmt, ..., verbosity = caller_verbosity()) {
    if (verbosity >= 2) logf(fmt, ...) else invisible(Sys.time())
}

#' @noRd
#' @title Inherit verbosity from the calling environment
#'
#' @description
#' Looks first for `verbose` and then for `verbosity` in the caller's
#' environment, walking the lexical parent chain via
#' `mget(..., inherits = TRUE)`. Intended to be used as a default argument:
#' `f <- function(verbosity = caller_verbosity()) { ... }`. This makes
#' `verbosity` behave like an inherited parameter — the caller's value flows
#' down automatically without having to be passed explicitly at every call
#' site.
#'
#' @details
#' **How the frame arithmetic works**
#'
#' The canonical use pattern is `f <- function(verbosity = caller_verbosity())`.
#' In R, default arguments are evaluated lazily inside the function's own
#' frame. So when `f` is called without `verbosity`, the promise
#' `caller_verbosity()` is forced inside `f`'s frame. At that point:
#'
#' - `parent.frame(1)` inside `caller_verbosity` is `f`'s frame.
#' - `parent.frame(2)` is `f`'s caller — exactly the environment from which
#'   `verbose` or `verbosity` should be inherited.
#'
#' The lookup then walks lexical parents via `inherits = TRUE`, so it also
#' resolves `verbose` or `verbosity` bound in enclosing closures (e.g.
#' anonymous functions passed to `lapply`).
#'
#' **Why implicit inheritance is acceptable here**
#'
#' As a general principle, relying on implicit parameter inheritance — where
#' a value is silently picked up from the calling scope rather than passed
#' explicitly — is a bad practice. It makes the call graph opaque: a reader
#' cannot see what value a function receives just by reading the call site,
#' which complicates both understanding and debugging.
#'
#' For `verbosity`, however, both of these concerns are greatly weakened:
#'
#' - *It does not affect correctness.* `verbosity` only controls how much
#'   is printed; it has no effect on computed results. A wrong value causes
#'   too much or too little output, never a wrong answer.
#' - *Explicit passing adds pervasive clutter.* `verbosity` must be threaded
#'   through every single function in the call chain. In a deep pipeline this
#'   means adding a `verbosity` argument and a `verbosity = verbosity`
#'   call-site token to dozens of functions, none of which add any
#'   informational value to the code as documentation.
#'
#' The combination — negligible risk from implicit inheritance, high
#' readability cost from explicit passing — makes `caller_verbosity()` an
#' acceptable exception to the general rule.
#'
#' @param default
#' Numeric value returned when neither `verbose` nor `verbosity` is found in
#' the calling environment or any of its lexical parents. Defaults to `1`.
#'
#' @return
#' A single numeric value: the resolved `verbose`/`verbosity` value, or
#' `default` if neither is found.
#'
#' @examples
#' d <- function(verbosity = caller_verbosity()) {
#'     logv("d: running")         # prints when verbosity >= 1
#'     logvv("d: verbose detail") # prints when verbosity >= 2
#' }
#' b <- function(verbosity = caller_verbosity()) d()
#' a <- function(verbosity = 1) b()
#'
#' a(verbosity = 0) # Expected: nothing printed
#' a(verbosity = 1) # Expected: "d: running" printed
#' a(verbosity = 2) # Expected: both lines printed
caller_verbosity <- function(default = 1) {
    env <- parent.frame(2) # skip caller_verbosity; land in the caller's caller
    mget("verbose", envir=env, inherits=TRUE, ifnotfound=list(NULL))[[1]] %||%
        mget("verbosity", envir=env, inherits=TRUE, ifnotfound=list(NULL))[[1]] %||%
        default
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
human_readable_difftime <- function(x) {
    sec <- as.numeric(x, units = "secs")
    out <- character(length(sec))
    for (i in seq_along(sec)) {
        s <- sec[i]
        if (is.na(s)) {
            out[i] <- NA_character_
            next
        }
        sign <- if (s < 0) "-" else ""
        s <- abs(s)
        if (s >= 3600) {
            h <- floor(s / 3600)
            m <- floor((s %% 3600) / 60)
            out[i] <- sprintf("%s%dh %dmin", sign, h, m)
        } else if (s >= 60) {
            m <- floor(s / 60)
            ss <- floor(s %% 60)
            out[i] <- sprintf("%s%dmin %ds", sign, m, ss)
        } else {
            out[i] <- sprintf("%s%.2fs", sign, s)
        }
    }
    out
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
human_readable <- function(x, unit, fmt = "%.1f") {
    if (inherits(x, "difftime")) return(human_readable_difftime(x))
    m <- max(abs(x))
    # styler: off
    if      (m >= 1e+9) { prefix <- "G"; x <- x / 1e+9 }
    else if (m >= 1e+6) { prefix <- "M"; x <- x / 1e+6 }
    else if (m >= 1e+3) { prefix <- "k"; x <- x / 1e+3 }
    else if (m <= 1e-3) { prefix <- "m"; x <- x / 1e-3 }
    else if (m <= 1e-6) { prefix <- "u"; x <- x / 1e-6 }
    else if (m <= 1e-9) { prefix <- "n"; x <- x / 1e-9 }
    else                { prefix <- "";  x <- x / 1e+0 }
    # styler: off
    fmtstr <- paste0(fmt, " %s%s")
    sprintf(fmtstr, x, prefix, unit)
}

# Type Checking (Private) #####

#' @noRd
#'
#' @name Type Checking
#'
#' @description
#' Returns a string giving the type and length of the input object for logical,
#' integer, double, complex, character, and raw vectors. For other objects, the
#' class of the object is returned.
#'
#' @param x
#' The object to check.
#'
#' @return
#' Aa character string indicating the type of the object.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' obj <- structure(list(a="a", b=1:5), class = "abc")
#' type(obj)   # "abc"
#' type(1:5)   # "integer(5)"
#' type(NULL)  # "NULL(0)"
#' type("Hi")  # "character(1)"
#' type(1.0)   # "numeric(1)"
#' type(NA)    # "logical(1)"
#' type(NaN)   # "numeric(2)"
type <- function(x) {
    typ <- typeof(x)
    if (typ %in% c("logical", "integer", "double", "complex", "character", "raw")) {
        sprintf("%s(%d)", typ, length(x))
    } else {
        class(x)
    }
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_num <- function(x, n = NULL) {
    is.numeric(x) && (is.null(n) || length(x) == n)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_int <- function(x, n = NULL) {
    x_is_int <- is.integer(x) || (is.numeric(x) && all(abs(x - round(x)) < sqrt(.Machine$double.eps)))
    has_correct_length <- is.null(n) || length(x) == n
    x_is_int && has_correct_length
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_char <- function(x, n = NULL, pattern = NULL) {
    is.character(x) &&
        (is.null(n) || length(x) == n ) &&
        (is.null(pattern) || all(grepl(pattern, x)))
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_bool <- function(x, n = NULL) {
    is.logical(x) && (is.null(n) || length(x) == n)
}

is_bool_or_num <- function(x, n = 1) {
    (is.logical(x) || (is.numeric(x)) && (is.null(n) || length(x) == n))
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_str <- function(x) {
    is.character(x) && length(x) == 1
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_num_or_null <- function(x, n = NULL) {
    is.null(x) || is_num(x, n)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_int_or_null <- function(x, n = NULL) {
    is.null(x) || is_int(x, n)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_char_or_null <- function(x, n = NULL, pattern = NULL) {
    is.null(x) || is_char(x, n, pattern)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_bool_or_null <- function(x, n = NULL) {
    is.null(x) || is_bool(x, n)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_str_or_null <- function(x) {
    is.null(x) || is_str(x)
}

#' @noRd
#'
#' @title Check if strings represent integer values
#'
#' @description
#' Tests each element of a character vector to determine if it represents an
#' integer value. It allows for optional leading and trailing whitespace, and
#' optional leading + or - signs.
#'
#' @param x
#' A character vector where each element is tested to determine if it represents
#' an integer.
#'
#' @return
#' A logical vector indicating TRUE for elements that represent integer values
#' and FALSE otherwise.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
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
#'
#' @title Check if strings represent floating-point numbers
#'
#' @description
#' Tests each element of a character vector to determine if it represents a
#' floating-point number. It handles numbers in fixed-point notation (with or
#' without a decimal point) and scientific notation. Allows for optional leading
#' and trailing whitespace, and optional leading + or - signs.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
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

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_existing_path <- function(x) {
    is_str(x) && file.exists(x)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_existing_dirpath <- function(x) {
    is_str(x) && dir.exists(x)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_equal <- function(x, y) {
    isTRUE(all.equal(x, y))
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
all_identical <- function(x) {
    all(sapply(x, identical, x[[1]]))
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_list_of_nums <- function(x, nl = NULL, nv = NULL) {
    if (!is.list(x)) return(FALSE)
    if (!(is.null(nl) || length(x) == nl)) return(FALSE)
    if (!(is.null(nv) || all(lengths(x) == nv))) return(FALSE)
    return(TRUE)
}

# Misc (Private) ##############################################################

cut2 <- function(x, breaks) {
    if (breaks == 1) {
        rep(1L, length(x))
    } else {
        cut(x, breaks, labels = FALSE)
    }
}

rbindlist <- function(x) {
    do.call(rbind, x)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
called_from_globalenv <- function() {
    identical(parent.frame(), .GlobalEnv)
}

#' @noRd
#' @title Get Named Function Arguments
#'
#' @description
#' Extracts the **named** arguments of a function as a named list. Variadic
#' arguments, i.e. `...`, are not included in the list. Missing values are
#' provided as "empty symbols".
#'
#' @param func
#' If provided, only arguments of this function are extracted from the
#' environment. See 'Details'.
#'
#' @param ignore
#' A character vector of argument names to ignore.
#'
#' @param env
#' The environment to extract the arguments from.
#'
#' @return
#' A named list of arguments.
#'
#' @details
#' Calling `args <- get_args(f)` as first statement in a function `f` produces
#' the same as `args <- as.list(environment())` (assuming no values were
#' provided via `...`). The advantage of using `get_args()` is, that it allows
#' to exclude certain arguments from the list and that it can be used
#' interactively from the global environment during function development.
#' (Calling `as.list(environment())` from the global environment would convert
#' the complete global environment into a list, meaning it can be extremely
#' slow.)
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' f <- function(a, b = 1, c = NULL, ...) {
#'      args <- get_args(f, ignore = c("a"))
#'      # do some calculations
#'      args
#' }
#' f(10, 20) # list(b = 20, c = NULL)
#'
#' g <- function(a, b = 1, c = NULL, ...) {
#'      get_args()
#' }
#' g(0, 1, 2, 3, 4) # list(a = 0, b = 1, c = 2   )
#' xx <- g()        # list(a = ,  b = 1, c = NULL)
#' is.symbol(xx$a)  # TRUE
get_args <- function(func = NULL, ignore = character(), env = parent.frame())
{
    assert(
        is.null(ignore) || is.character(ignore),
        is.environment(env),
        is.null(func) || is.function(func)
    )
    if (is.null(func)) {
        args <- as.list(env)
        args[ignore] <- NULL
    }  else {
        argnames <- names(formals(func))
        argnames <- argnames[!argnames %in% ignore]
        args <- sapply(argnames, function(name) env[[name]], simplify = FALSE)
    }
    args
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
empty_df <- function(names) {
    df <- data.frame(matrix(ncol = length(names), nrow = 0))
    colnames(df) <- names
    df
}

#' @noRd
#'
#' @title Calculate Magnetic Field Strength
#'
#' @description
#' Calculates the magnetic field strength based on the measured frequencies and
#' chemical shifts of a spectrum as well as the gyromagnetic ratio for protons.
#'
#' The function makes the following assumptions:
#'
#' 1. The spectrum is a 1H spectrum
#' 2. The chemical shift of the reference is at 0.0 ppm
#' 3. The resonance frequency of the reference equals the resonance frequency of
#'    protons
#'
#' @param x
#' A `spectrum` object as described in [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
#' @return
#' The magnetic field strength in Tesla.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' x <- read_spectrum(pkg_file("example_datasets/bruker/urine/urine_1"))
#' calc_B(x)
calc_B <- function(x = read_spectrum()) {
    # Example:
    # |---------------------------|----------------------|
    # 14.8 ppm (max)            0.0 (ref)       -5.2 (min)
    # 600.243 MHz (min)       600.252 (ref)  600.255 (max)
    cs_maxmin <- max(x$cs) - min(x$cs)
    cs_refmin <- 0.0000000 - min(x$cs)
    ratio <- cs_refmin / cs_maxmin
    fq_maxmin <- max(x$fq) - min(x$fq)
    fq_refmax <- ratio * fq_maxmin
    fq_ref <- max(x$fq) - fq_refmax # Frequency of the reference (which is equal the frequency of a proton)
    gamma <- 2.675e8 # Gyromagnetic ratio for protons
    B <- (2 * pi * fq_ref) / gamma
    B
}

#' @noRd
#'
#' @title Recursive Object Size Printer
#'
#' @description
#' Prints the size of an object and its subcomponents recursively, similar to
#' the `du` command in Unix. This function is useful for understanding the
#' memory footprint of complex objects in R.
#'
#' @param obj
#' The object to analyze.
#'
#' @param pname
#' Internal parameter for recursive calls, representing the parent name of the
#' current object. Should not be used directly.
#'
#' @param level
#' Internal parameter for recursive calls, indicating the current depth of
#' recursion. Should not be used directly.
#'
#' @param max_level
#' The maximum depth of recursion. Default is 1. Increase this to explore deeper
#' into nested structures.
#'
#' @param max_len
#' The maximum number of elements in a list to recurse into. Default is 50.
#' Lists with more elements will not be explored further.
#'
#' @param unit
#' The unit of measurement for object sizes. Can be "GB", "MB", "KB", or "B".
#' Default is "MB".
#'
#' @return
#' The total size of the object, including its subcomponents, as a character
#' string with the specified unit.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' obj <- list(
#'     a = rnorm(1000),
#'     b = list(
#'         c = rnorm(2000),
#'         d = list(
#'             e = rnorm(3000),
#'             f = rnorm(4000)
#'         )
#'     )
#' )
#' du(obj)
#' du(obj, max_level = Inf, unit = "KB")
du <- function(obj, pname = "", level = 0, max_level = 1, max_len = 50, unit = "MB") {
    match.arg(unit, c("GB", "MB", "KB", "B"))
    denom <- switch(unit,
        "GB" = 1e9,
        "MB" = 1e6,
        "KB" = 1e3,
        "B" = 1
    )
    size <- function(x) paste(round(object.size(x) / denom, 2), unit)
    for (name in names(obj)) {
        x <- obj[[name]]
        indent <- collapse(rep("    ", level), "")
        cat2(indent, name, " ", size(x), sep = "")
        if (is.list(x) && length(x) < max_len && level < max_level) {
            du(x, name, level + 1, max_level, unit = unit)
        }
    }
    if (level == 0) {
        total <- size(obj)
        cat2("Total:", total)
        invisible(total)
    }
}

#' @export
#'
#' @title Print the Structure of a Directory Tree
#'
#' @description
#' Prints the structure of a directory tree up to a specified maximum level of
#' depth. It lists all files and directories under the specified path,
#' displaying them in a tree-like structure.
#'
#' @param path The root path from which to start listing the directory
#' structure.
#'
#' @param max.level The maximum depth of directories to list.
#'
#' @param max.entries Maximum number of children to print per directory.
#' If a directory has more entries than this limit, the first
#' `ceiling(max.entries / 2)` and last `floor(max.entries / 2)` children are
#' shown, with `... [N]` in between, where `N` is the number of omitted
#' children.
#'
#' @param show.counts Logical. If `TRUE`, prints the number of direct children
#' (files + subdirectories) in brackets after each directory name. Disabled by
#' default.
#'
#' @param files.first Logical. If `TRUE`, files are listed before
#' subdirectories. If `FALSE` (default), subdirectories are listed first.
#'
#' @param level Internal parameter used for recursion, indicating the current
#' level of depth.
#'
#' @param prefix Internal parameter used for formatting the printed tree
#' structure.
#'
#' @return
#' NULL, called for its side effect of printing the directory structure.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' metabodecon_dir <- system.file(package = "metabodecon")
#' tree(metabodecon_dir, max.level = 1)
tree <- function(path,
                 max.level = 2,
                 max.entries = Inf,
                 show.counts = FALSE,
                 files.first = FALSE,
                 level = 0,
                 prefix = "") {

    if (level == max.level) return()

    entries <- list.files(path, full.names = TRUE)
    dirs <- entries[isdir <- file.info(entries)$isdir]
    files <- entries[!isdir]
    dirs <- sort(dirs)
    files <- sort(files)
    all_entries <- if (isTRUE(files.first)) c(files, dirs) else c(dirs, files)
    num_entries <- length(all_entries)

    if (level == 0) {
        if (isTRUE(show.counts)) {
            root_n_entries <- count_entries_in_dir(path)
            cat(path, " [", root_n_entries, "]", "\n", sep = "")
        } else {
            cat(path, "\n", sep = "")
        }
    }

    if (!is.numeric(max.entries) || length(max.entries) != 1 || is.na(max.entries)) {
        stop("max.entries must be a single non-missing numeric value.")
    }
    if (is.infinite(max.entries)) {
        n_show <- num_entries
    } else {
        n_show <- min(num_entries, max(0, floor(max.entries)))
    }
    n_show <- as.integer(n_show)
    truncated <- num_entries > n_show
    n_top <- if (truncated) ceiling(n_show / 2) else n_show
    n_bottom <- if (truncated) floor(n_show / 2) else 0L
    shown_top <- if (n_top > 0) all_entries[seq_len(n_top)] else character(0)
    shown_bottom <- if (n_bottom > 0) tail(all_entries, n_bottom) else character(0)

    printed_entries <- if (truncated) {
        c(shown_top, "...", shown_bottom)
    } else {
        shown_top
    }
    for (i in seq_along(printed_entries)) {
        entry <- printed_entries[i]
        is_ellipsis <- identical(entry, "...")
        is_last <- i == length(printed_entries)
        prefix2 <- if (is_last) "\u2514\u2500\u2500 " else "\u251C\u2500\u2500 "
        if (is_ellipsis) {
            n_omitted <- num_entries - n_show
            cat(prefix, prefix2, "... [", n_omitted, "]", "\n", sep = "")
            next
        }
        entry_is_dir <- isTRUE(file.info(entry)$isdir)
        if (entry_is_dir) {
            if (isTRUE(show.counts)) {
                n_entries <- count_entries_in_dir(entry)
                cat(prefix, prefix2, basename(entry), "/", " [", n_entries, "]", "\n", sep = "")
            } else {
                cat(prefix, prefix2, basename(entry), "/", "\n", sep = "")
            }
        } else {
            cat(prefix, prefix2, basename(entry), "\n", sep = "")
        }
        new_prefix <- if (is_last) paste0(prefix, "    ") else paste0(prefix, "\u2502   ")
        if (entry_is_dir) {
            tree(
                path = entry,
                max.level = max.level,
                max.entries = max.entries,
                show.counts = show.counts,
                files.first = files.first,
                level = level + 1,
                prefix = new_prefix
            )
        }
    }

    invisible(NULL)
}

#' @export
#' @rdname tree
tree_preview <- function(path,
                         max.level = 1,
                         max.entries = 9,
                         show.counts = TRUE,
                         files.first = TRUE) {
    tree(
        path = path,
        max.level = max.level,
        max.entries = max.entries,
        show.counts = show.counts,
        files.first = files.first
    )
}

#' @noRd
#' @description Count direct children (files and subdirectories) in a directory.
count_entries_in_dir <- function(dir) {
    entries <- list.files(dir, full.names = TRUE)
    as.integer(length(entries))
}

#' @noRd
#' @description Update values of a list and return the modified list.
#' @details
#' Setting a value to NULL will **not** remove the key from the list, but set
#' the value to NULL. I.e. setting values using `set` is closer to `[<-` than
#' to `[[<-`. Example:
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' # Inputs
#' xx <- list(a=1, b=2:5, d="a")
#'
#' # Examples
#' yy <- set(xx, a=10, d=4, f=sqrt)  # Set a and d. Add f.
#' zz <- set(xx, a=NULL)             # Set a to NULL.
#'
#' # Outputs
#' stopifnot(identical(yy, list(a=10, b=2:5, d=4, f=sqrt)))
#' stopifnot(identical(zz, list(a=NULL, b=2:5, d="a")))
#'
#' # For comparison
#' xx$a <- NULL                      # Removes a
#' xx[["a"]] <- NULL                 # Removes a
#' xx["a"] <- list(NULL)             # Set a to NULL (like `set`)
#'
set <- function(...) {
    args <- list(...)
    obj <- args[[1]]
    vals <- args[-1]
    keys <- names(vals)
    mapply(function(k, v) obj[k] <<- list(v), keys, vals)
    obj
}

#' @noRd
#' @title Pop named Element from List
#' @description
#' Returns the value of `obj[[key]]` and sets `obj[[key]]` to NULL. This
#' function modifies the input list in-place, which is different from most other
#' R functions.
#' @param obj A list.
#' @param obj A list.
#' @param key The key to pop from the list.
#' @param default The default value to return if the key does not exist in the list.
#' @author 2024-2025 Tobias Schmidt: initial version.
#' @examples
#' tmp <- list(a = NULL, b = 2); pop(tmp, "a")    # Returns NULL
#' tmp <- list(a = NULL, b = 2); pop(tmp, "a", 5) # Returns NULL
#' tmp <- list(a = NULL, b = 2); pop(tmp, "b")    # Returns 2
#' tmp <- list(a = NULL, b = 2); pop(tmp, "b", 5) # Returns 2
#' tmp <- list(a = NULL, b = 2); pop(tmp, "z")    # Returns NULL
#' tmp <- list(a = NULL, b = 2); pop(tmp, "z", 5) # Returns 5
pop <- function(obj, key, default = NULL, env = parent.frame()) {
    val <- obj[[key]]
    if (!is.null(default) && is.null(val) && !exists(key, obj)) {
        # If the value of 'key' is null, it can have two reasons:
        # 1. The key does not exist in 'obj'. In this case we want to return the 'default' value.
        # 2. The key exists in 'obj', but is explicitly set to NULL. In this case we want to return NULL.
        # If 'default' is NULL, we can skip this check as we would return NULL in any case.
        val <- default
    }
    expr <- substitute(obj[[key]] <- NULL)
    eval(expr, env) # Remove key from obj in calling env
    val
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
timestamp <- function() {
    format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
}

#' @export
#' @title Make transparent
#' @description Make a color transparent by adding an alpha channel.
#' @param col Character string specifying the color to make transparent.
#' @param alpha Numeric value between 0 and 1 specifying the transparency level.
#' @return A character string representing the color with an alpha channel.
#' @author 2024-2025 Tobias Schmidt: initial version.
#' @examples
#' transp("violet", 0.08)
#' transp("black", 0.5)
transp <- function(col = "violet", alpha = 0.08) {
    col <- col2rgb(col)[, 1] / 255
    rgb(col[1], col[2], col[3], alpha = alpha)
}

#' @noRd
#' @description
#' Returns half the available CPU cores (rounded up).
half_cores <- function(nmax = Inf) {
    min(ceiling(parallel::detectCores() / 2), nmax)
}

#' @noRd
#' @description
#' Multi core version of mapply with automatic logging of worker output.
#' For `nw == 1` normal `mapply` is used.
#' @author 2024-2026 Tobias Schmidt: initial version.
mcmapply <- function(
    nw, FUN, ...,
    MoreArgs = NULL, RECYCLE = TRUE, SIMPLIFY = FALSE,
    USE.NAMES = TRUE, .scheduling = c("static", "dynamic"),
    loadpkg = TRUE, log = TRUE
) {
    if (nw == 1) {
        mapply(
            FUN, ...,
            MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY,
            USE.NAMES = USE.NAMES
        )
    } else if (.Platform$OS.type == "windows") {
        .scheduling <- match.arg(.scheduling)
        cl <- get_worker_pool(nw, loadpkg = loadpkg, log = log)
        on.exit(stopCluster(cl), add = TRUE, after = FALSE)
        clusterMap(
            cl, FUN, ...,
            MoreArgs = MoreArgs, RECYCLE = RECYCLE, SIMPLIFY = SIMPLIFY,
            USE.NAMES = USE.NAMES, .scheduling = .scheduling
        )
    } else {
        parallel::mcmapply(
            FUN, ...,
            MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY,
            USE.NAMES = USE.NAMES, mc.cores = nw
        )
    }
}

#' @noRd
#' @description
#' Create a worker pool with package loading and per-worker logging.
#' @author 2026 Tobias Schmidt: initial version.
get_worker_pool <- function(nw, loadpkg = TRUE, log = TRUE) {
    logf("Creating pool of worker processes")
    cl <- makeCluster(nw)
    if (loadpkg) {
        expr <- parse(text = 'attach(asNamespace("metabodecon"))')
        clusterCall(cl, eval, expr)
    }
    if (log) {
        logfiles <- get_worker_logs(nw)
        logf("For progress output see:")
        sapply(logfiles, logf)
        clusterApply(cl, logfiles, sink, append = TRUE, type = "output")
    }
    cl
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
get_worker_logs <- function(nw, create = TRUE) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%OS3")
    timestamp <- gsub(".", "_", timestamp, fixed = TRUE)
    logdirrel <- file.path("logs", timestamp)
    logdir <- tmpdir(subdir = logdirrel, create = TRUE)
    logfiles <- paste0("worker_", seq_len(nw), ".log")
    logpaths <- file.path(logdir, logfiles)
    if (create) file.create(logpaths)
    logpaths
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
load_all <- function() {
    x <- Sys.time()
    logf("Calling: pkgload::load_all(reset = TRUE)")
    pkgload::load_all(reset = TRUE, quiet = TRUE)
    logf("Calling: pkgload_env$insert_global_shims(force = TRUE)")
    pkgload_env <- environment(pkgload::load_all)
    pkgload_env$insert_global_shims(force = TRUE)
    diff <- Sys.time() - x
    logf("Elapsed: %s", format(diff))
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
document <- function() {
    x <- Sys.time()
    logf("Calling: devtools::document(quiet = TRUE)")
    devtools::document(quiet = TRUE)
    logf("Calling: pkgload_env$insert_global_shims(force = TRUE)")
    pkgload_env <- environment(pkgload::load_all)
    pkgload_env$insert_global_shims(force = TRUE)
    diff <- Sys.time() - x
    logf("Elapsed: %s", format(diff))
}

#' @noRd
#' @author 2026 Tobias Schmidt: initial version.
style <- function() {
    stopifnot(loaded_via_devtools())
    style_func <- function(...) {
        transformers <- styler::tidyverse_style(indent_by = 4, ...)
        unindent_fun_dec <- transformers$indention$unindent_fun_dec
        formals(unindent_fun_dec)$indent_by <- 4L
        transformers$indention$unindent_fun_dec <- unindent_fun_dec
        transformers$space$set_no_space_around_eq_sub <- function(pd) {
            is_eq_sub <- pd$token == "EQ_SUB"
            is_before_eq_sub <- c(is_eq_sub[-1], FALSE)
            pd$spaces[is_eq_sub | is_before_eq_sub] <- 0L
            pd
        }
        transformers
    }
    styler::style_pkg(style = style_func)
}


# On Load (Private) #####

#' @noRd
#'
#' @description
#' Acts like [stopifnot()] during development but does nothing in production.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @details
#' In the package source code, this function is defined as a copy of
#' [stopifnot()]. However, during package loading, it is replaced with an empty
#' function unless the package is loaded via [devtools::load_all()]. The actual
#' replacement is implemented in [.onLoad()].
#'
#' The idea is that exported functions should use plain [stopifnot()] to
#' validate their inputs, whereas private functions should use
#' [metabodecon::assert()] instead. This approach allows us to use rigorous type
#' checking during development without impacting performance in production.
#'
#' If we need to keep assertions enabled in production, we can set the option
#' `metabodecon.assert` to `stopifnot` before loading the package.
#'
#' If we want to disable assertions during development, e.g. to get realistic
#' runtime estimates, we can set the option `metabodecon.assert` to
#' `function(...) {}` before calling `devtools::load_all()`.
#'
#' Example:
#'
#' ```r
#' # Steps:
#' # (1) Load metabodecon with assertions disabled
#  # (2) Unload metabodecon
#  # (3A) Configure stopifnot as the assertion function OR
#  # (3B) Configure empty function as assertion function
#  # (4) Reload metabodecon
#' library(metabodecon)                           # (1)
#' unloadNamespace("metabodecon")                 # (2)
#' options(metabodecon.assert = stopifnot)        # (3A)
#' options(metabodecon.assert = function(...) {}) # (3B)
#' library(metabodecon)                           # (4)
#' ```
assert <- stopifnot

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
.onLoad <- function(libname, pkgname) {
    pkgenv <- topenv()

    # Only enable assertions during development, as they are quite runtime heavy
    # and we want to catch problems early.
    if (!loaded_via_devtools()) pkgenv$assert <- function(...) {}

    # Allow enabling assertions even when loaded via library by setting the
    # option `metabodecon.assert` to `stopifnot` before loading the package.
    if (!is.null(x <- .Options$metabodecon.assert)) pkgenv$assert <- x

    # Create a folder 'cache' inside the persistent data directory if we are on
    # a development machine. This allows us to keep cache data across R sessions
    # during development, which can speed up development iterations a lot.
    if (loaded_via_devtools()) {
        cache_dir <- file.path(datadir_persistent(), "cache")
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        options(metabodecon.aki_cache = cache_dir)
    }
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
loaded_via_devtools <- function() {
    pkg_dir <- dirname(system.file("DESCRIPTION", package = "metabodecon"))
    dir.exists(file.path(pkg_dir, "inst"))
}
