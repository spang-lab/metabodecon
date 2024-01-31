#' @name penv
#' @title Environment for storing global package data
#' @description `penv` is an environment that contains the following elements:
#' * `fn_backups`: used by [mock()] to store original functions so that [restore()] can restore them.
#' * `option_backups`: used by [push_option()] to store original options so that [restore()] can restore them.
#' * `open_conns`: used by [redirect()] to store open connections so that [restore()] can close them.
#' * `wds`: used by [pushd()] so that [popd()] and [restore()] can restore the working directory.
#' @noRd
penv <- as.environment(list(
    fn_backups = list(),
    option_backups = list(),
    open_conns = list(),
    wds = list()
))

#' @name pushd
#' @title Push current working directory and change to new directory
#' @description This function stores the current working directory in `penv$wds` and then changes to the specified new directory.
#' @param newdir The new directory to switch to.
#' @return Invisible NULL.
#' @noRd
pushd <- function(newdir) {
    olddir <- getwd()
    penv$wds <- c(olddir, penv$wds)
    setwd(newdir)
    invisible(NULL)
}

#' @name popd
#' @title Pop the last working directory and change to it
#' @description This function restores the last working directory stored in `penv$wds` and changes to it.
#' @param all Logical. If TRUE, all directories are popped. If FALSE, only the last directory is popped.
#' @return Invisible NULL.
#' @noRd
popd <- function(all = FALSE) {
    if (length(penv$wds) == 0) {
        return(invisible(NULL))
    } else if (all) {
        newdir <- penv$wds[[length(penv$wds)]]
        penv$wds <- NULL
    } else {
        newdir <- penv$wds[[1]]
        penv$wds <- penv$wds[-1]
    }
    setwd(newdir)
    return(invisible(NULL))
}

#' @title Do nothing
#' @param ... Not used
#' @return NULL
#' @noRd
pass <- function(...) {
    NULL
}

#' @title Concatenate and print with newline
#' @param ... Arguments to be concatenated and printed.
#' @examples \dontrun{
#' cat2("Hello, ", "world!")
#' }
#' @noRd
cat2 <- function(...) {
    cat(...)
    cat("\n")
}

#' @title Calculate a checksum for all files in a directory or a single file
#' @description This function calculates a checksum for each file in a specified directory or a single file. If the input is a directory, the checksums are calculated recursively, meaning that it includes files in all subdirectories. The results are returned as a named vector, where the names are the relative file paths and the values are checksums.
#' @param path The directory or file to calculate checksums for.
#' @param method The method to use for calculating the checksum. Can be "size" (default) or "md5". If "size", the function returns the file sizes. If "md5", the function returns the MD5 hashes of the files.
#' @param ignore A character vector of regular expressions. Files matching any of these regular expressions will be ignored.
#' @return A named vector with file paths as names and hashes as values. If the input is a directory, the names will be the file paths relative to the directory. If the input is a file, the name will be the file name.
#' @examples \dontrun{
#' checksum(system.file(package = "metabodecon")) # directory example
#' checksum(system.file("DESCRIPTION", package = "metabodecon"))
#' checksum(system.file("DESCRIPTION", package = "metabodecon"), method = "md5")
#' }
#' @details By default, the "checksum" calculated for each file is just its size. This method was chosen because it is the fastest available and typically sufficient for our needs. Traditional checksum methods, such as MD5, can present issues. For instance, PDF files may yield different checksums every time they are recreated, likely due to the inclusion of timestamps or similar metadata within the file.
#' @noRd
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
        checksums <- calc_checksum(paths)
        names(checksums) <- nams
    } else {
        checksums <- calc_checksum(path)
        names(checksums) <- basename(path)
    }
    checksums
}

#' @title Run tests with the option to skip slow tests
#' @description This function runs the tests in the current R package. It has the option to skip slow tests. The decision to skip slow tests is controlled by the SKIP_SLOW_TESTS environment variable. If 'all' is FALSE and 'skip_slow' is TRUE, then slow tests are skipped. The original value of the SKIP_SLOW_TESTS environment variable is restored after the tests are run.
#' @param all Logical. If TRUE, all tests are run. If FALSE, slow tests are skipped.
#' @return The result of devtools::test()
#' @noRd
test <- function(all=FALSE) {
    SKIP_SLOW_TESTS_OLD <- Sys.getenv("SKIP_SLOW_TESTS")
    Sys.setenv(SKIP_SLOW_TESTS = if (all) "FALSE" else "TRUE")
    on.exit(Sys.setenv(SKIP_SLOW_TESTS = SKIP_SLOW_TESTS_OLD), add = TRUE)
    devtools::test()
}

#' @title Prepares the output directory for a test case
#' @param fn The name of the function being tested.
#' @param tc The name of the test case.
#' @param inputs Paths to be copied to the output directory.
#' @param verbose Print info messages.
#' @return Path to the created output directory.
#' @examples \dontrun{
#' output_dir <- prepare_test_dir("myfunc", "1")
#' output_dir <- prepare_test_dir("myfunc", "2", inputs = pkg.file("misc/datasets/urine"))
#' }
#' @noRd
prepare_test_dir <- function(fn, tc, inputs = c(), verbose = TRUE) {
    msg <- if (verbose) cat2 else pass
    d <- tempdir()
    p <- normalizePath(file.path(d, "tests", fn, tc), winslash = "/", mustWork = FALSE)
    msg("Deleting and creating:", p)
    unlink(p, recursive = TRUE)
    dir.create(p, recursive = TRUE)
    if (length(inputs) > 0) {
        msg("Copying following paths:")
        msg(collapse(inputs, "\n"))
        file.copy(from = inputs, to = p, recursive = TRUE)
    }
    return(p)
}

#' @title Return path to any file within this (installed) package
#' @param file (string) Relative path to file.
#' @param ... Arguments passed on to [system.file()].
#' @return Absolute path to `file` with '/' as file separator.
#' @examples \dontrun{
#' pkg.file("DESCRIPTION")
#' pkg.file() # Path to the package root directory
#' }
#' @noRd
pkg.file <- function(...) {
    system.file(..., package = "metabodecon")
}

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
