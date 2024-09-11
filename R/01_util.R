#' @import grDevices
#' @import stats
#' @import utils
#' @import parallel
#' @import stats
#' @import graphics
#' @import grDevices
#' @import mathjaxr

#' @import toscutil

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
#'     show_dir <- get_yn_input("List dir content? (y/n) ")
#'     if (show_dir) print(dir())
#' }
get_yn_input <- function(prompt) {
    if (!grepl("(y/n)", prompt, fixed = TRUE)) {
        prompt <- paste0(prompt, " (y/n) ")
    }
    x <- get_str_input(prompt, c("y", "n"))
    y <- if (x == "y") TRUE else FALSE
    return(y)
}

# Interactive #####

#' @noRd
#' @title Recursive Object Size Printer
#' @description Prints the size of an object and its subcomponents recursively, similar to the `du` command in Unix. This function is useful for understanding the memory footprint of complex objects in R.
#' @param obj The object to analyze.
#' @param pname Internal parameter for recursive calls, representing the parent name of the current object. Should not be used directly.
#' @param level Internal parameter for recursive calls, indicating the current depth of recursion. Should not be used directly.
#' @param max_level The maximum depth of recursion. Default is 1. Increase this to explore deeper into nested structures.
#' @param max_len The maximum number of elements in a list to recurse into. Default is 50. Lists with more elements will not be explored further.
#' @param unit The unit of measurement for object sizes. Can be "GB", "MB", "KB", or "B". Default is "MB".
#' @return The total size of the object, including its subcomponents, as a character string with the specified unit.
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
#' @title Print the Structure of a Directory Tree
#' @description Prints the structure of a directory tree up to a specified maximum level of depth. It lists all files and directories under the specified path, displaying them in a tree-like structure.
#' @param path The root path from which to start listing the directory structure.
#' @param max.level The maximum depth of directories to list.
#' @param level Internal parameter used for recursion, indicating the current level of depth.
#' @param prefix Internal parameter used for formatting the printed tree structure.
#' @return NULL, called for its side effect of printing the directory structure.
#' @examples
#' metabodecon_dir <- system.file(package = "metabodecon")
#' tree(metabodecon_dir, max.level = 1)
tree <- function(path, max.level = 2, level = 0, prefix = "") {
    if (level == max.level) {
        return()
    }
    entries <- list.files(path, full.names = TRUE)
    dirs <- entries[isdir <- file.info(entries)$isdir]
    files <- entries[!isdir]
    all_entries <- sort(c(dirs, files))
    num_entries <- length(all_entries)
    if (level == 0) cat(path, "\n", sep = "")
    for (i in seq_along(all_entries)) {
        entry <- all_entries[i]
        is_last <- i == num_entries
        prefix2 <- if(is_last) "\u2514\u2500\u2500 " else "\u251C\u2500\u2500 "
        cat(prefix, prefix2, basename(entry), ifelse(file.info(entry)$isdir, "/", ""), "\n", sep = "")
        new_prefix <- if (is_last) paste0(prefix, "    ") else paste0(prefix, "\u2502   ")
        if (file.info(entry)$isdir) tree(entry, max.level, level + 1, new_prefix)
    }
    invisible(NULL)
}

# Operators #####

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

# Convert #####

#' @export
#' @title Calculate the Width of a Numeric Vector
#' @description Calculates the width of a numeric vector by computing the difference between the maximum and minimum values in the vector.
#' @param x A numeric vector.
#' @return The width of the vector, calculated as the difference between its maximum and minimum values.
#' @examples
#' vec <- c(1, 3, 5, 7, 9)
#' width(vec)
width <- function(x) {
    diff(range(x))
}

#' @export
#' @title Convert from unit A to unit B
#' @description `convert_pos` converts positions from unit A to unit B. `convert_width` converts widths from unit A to unit B. If the direction of units A and B is reversed, the width's sign will be reversed as well. To keep widths strictly positive, wrap the result with `abs()`.
#' @param xa For `convert_width` a numeric vector specifying widths in unit A. For `convert_pos` a numeric vector specifying positions in unit A.
#' @param ya A numeric vector giving the positions of at least two points in unit A.
#' @param yb A numeric vector giving the positions of the same points in unit B.
#' @return A numeric vector of values converted from unit A to unit B.
#' @examples
#' ya <- c(244, 246, 248, 250, 252)
#' yb <- c(15, 10, 5, 0, -5)
#' convert_width(c(2, 4, 8), ya, yb)
#' convert_pos(c(247, 249), ya, yb)
convert_width <- function(xa, ya, yb) {
    n <- length(ya)
    wa <- ya[n] - ya[1]
    wb <- yb[n] - yb[1]
    xa * (wb / wa)
}

#' @export
#' @rdname convert_width
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

#' @noRd
#' @title Calculate Magnetic Field Strength
#' @description Calculates the magnetic field strength based on the measured frequencies and chemical shifts of a spectrum as well as the gyromagnetic ratio for protons. For this to work, the following is assumed:
#' 1. The spectrum is a 1H spectrum
#' 2. The chemical shift of the reference is at 0.0 ppm
#' 3. The resonance frequency of the reference equals the resonance frequency of protons
#' @param x A spectrum object as described in [metabodecon_classes].
#' @return The magnetic field strength in Tesla.
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

# Misc #####

#' @noRd
#' @examples
#' xx <- list(a=1, b=2:5, d="a")
#' set(xx, a=10, d=4, f=sqrt)
set <- function(...) {
    args <- list(...)
    obj <- args[[1]]
    vals <- args[-1]
    keys <- names(vals)
    mapply(function(k, v) obj[[k]] <<- v, keys, vals)
    obj
}

# Doc Pages #####

#' @title Metabodecon Classes
#' @description
#'  Metabodecon introduces a set of classes to highlight the presence of certain elements in corresponding objects.
#'
#'  The order of elements may vary between different versions of Metabodecon, thus elements should always be accessed by name, for example, using `x$si` or `x[["cs"]]`.
#'  A short description of each class is given in the listing below.
#'
#'  -  `spectrum`: One NMR spectrum
#'  -  `decon0`: One deconvoluted NMR spectrum stored in [MetaboDecon1D()] format
#'  -  `decon1`: One deconvoluted NMR spectrum stored in [generate_lorentz_curves()] format
#'  -  `decon2`: One deconvoluted NMR spectrum stored in [deconvolute_gspec()] format
#'  -  `alignment`: One aligned NMR spectrum
#'  -  `spectra`: List of multiple NMR spectra
#'  -  `decons1`: List of multiple deconvoluted NMR stored in [MetaboDecon1D()] format
#'  -  `decons2`: List of multiple deconvoluted NMR stored in [generate_lorentz_curves()] format
#'  -  `decons3`: List of multiple deconvoluted NMR stored in [deconvolute_gspec()] format
#'  -  `alignments`: List of multiple aligned NMR spectra
#'
#'  More details can be found in Metabodecon's online documentation at [Metabodecon Classes](https://spang-lab.github.io/metabodecon/articles/Metabodecon-Classes.html).
metabodecon_classes <- NULL
