# General #####

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
#' @return Path to old working directory.
#' @noRd
pushd <- function(newdir) {
    olddir <- getwd()
    penv$wds <- c(olddir, penv$wds)
    setwd(newdir)
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
#' @description This function runs the tests in the current R package. If `all` is TRUE, it well set environment variable `RUN_SLOW_TESTS` to "TRUE" so that all tests are run. If `all` is FALSE, it will set `RUN_SLOW_TESTS` to "FALSE" so that slow tests are skipped.
#' @param all Logical. If TRUE, all tests are run. If FALSE, slow tests are skipped.
#' @return The result of devtools::test()
#' @noRd
test <- function(all = FALSE) {
    RUN_SLOW_TESTS_OLD <- Sys.getenv("RUN_SLOW_TESTS")
    Sys.setenv(RUN_SLOW_TESTS = if (all) "TRUE" else "FALSE")
    on.exit(Sys.setenv(RUN_SLOW_TESTS = RUN_SLOW_TESTS_OLD), add = TRUE)
    devtools::test()
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

dir.size <- function(dir) {
    files <- list.files(dir, recursive = TRUE, full.names = TRUE)
    sum(file.info(files)$size)
}

read_spectrum_info <- function(path) {
    # Written by [MetaboDecon1D()] as follows:
    # utils::write.table(spectrum_info, path, sep=",", col.names=FALSE, append=FALSE)
    spectrum_info <- utils::read.table(path, sep = ",", header = FALSE, row.names = 1)
    spectrum_info <- `colnames<-`(spectrum_info, paste0("X", seq_along(colnames(spectrum_info))))
}

read_spectrum_output <- function(path) {
    # Written by [MetaboDecon1D()] as follows:
    # utils::write.table(spectrum_output, path, sep=",", col.names=FALSE, append=FALSE)
    spectrum_output <- utils::read.table(path, sep = ",", header = FALSE, row.names = 1)
    spectrum_output <- `colnames<-`(spectrum_output, paste0("X", seq_along(colnames(spectrum_output))))
}

#' @title Check if the size of each file in a directory is within a certain range
#' @description This function checks if the size of each file in a directory is within 90% to 110% of the expected size.
#' If a file size is not within this range, a message is printed and an error is thrown.
#' @param testdir A character string specifying the directory to check.
#' @param size_exp A named numeric vector where the names are filenames and the values are the expected file sizes.
#' @examples
#' \dontrun{
#' testdir <- tempdir()
#' file.create(file.path(testdir, "file1.txt"))
#' file.create(file.path(testdir, "file2.txt"))
#' size_exp <- c(file1.txt = 100, file2.txt = 200)
#' expect_file_size(testdir, size_exp)
#' }
#' @noRd
expect_file_size <- function(testdir, size_exp) {
    paths <- file.path(testdir, names(size_exp))
    size_obs <- file.info(paths)$size
    size_min <- size_exp * 0.9
    size_max <- size_exp * 1.1
    file_has_correct_size <- size_obs > size_exp * 0.9 & size_obs < size_exp * 1.1
    lapply(seq_along(size_exp), function(i) {
        if (!isTRUE(file_has_correct_size[i])) {
            message(sprintf("Size of %s is %s which is not between %s and %s", paths[i], size_obs[i], size_exp[i] * 0.9, size_exp[i] * 1.1))
        }
    })
    testthat::expect_true(all(file_has_correct_size))
}

#' Test if the structure of an object matches the expected string
#'
#' @param obj The object to test
#' @param expected_str The expected structure of the object as a string. Can be obtained by calling `dput(capture.output(str(obj)))`.
#' @return A logical value indicating whether the structure of the object matches the expected string
#' @examples
#' expect_str(list(a = 1, b = 2), c("List of 2", " $ a: num 1", " $ b: num 2"))
#' @noRd
expect_str <- function(obj, expected_str) {
    testthat::expect_identical(capture.output(str(obj)), expected_str)
}

# With #####

#' @description Run expression with predefined global state
#' @param expr Expression to be evaluated.
#' @param testdir ID of the test directory. E.g. `"MetaboDecon1D/2"`. Will be created and populated with `inputs`, but cleared.
#' @param answers Answers to be returned by readline().
#' @param output Outputs to be captured. Can be "stdout" and/or "stderr".
#' @param message Path to the file where stdout should be redirected to.
#' @param plots Path to the pdf file where plots should be saved to.
#' @param datadir_temp State of the mocked temporary data directory. See details section.
#' @param datadir_persistent State of the mocked persistent data directory. See details section.
#' @param inputs Paths to be copied to the test directory. See [fill_with_inputs()] for details.
#' @param debug Call `browser()` before evaluating `expr`?
#' @param ... Additional keyword params passed to [push_option()]
#' @details The `datadir_temp` and `datadir_persistent` arguments accept values "missing", "filled" and "empty".
#' Setting a value unequal NULL causes the functions [datadir_temp()] and/or [datadir_persistent()] to be replaced with mocks functions pointing to fake directories. Functions depending on these functions will then use the fake directories instead of the real ones.
#' When set to "missing" the returned mock directory does not exist.
#' When set to "empty" it exists and is guaranteed to be empty.
#' When set to "filled", it is populated with example datasets.
#' Attention: the mocked functions, i.e. [datadir_temp()] and [datadir_persistent()] cannot be used directly inside the expression when called via `devtools::test()`. I'm not sure why, but it seems as if devtools and/or testthat have their own copies of the functions which are used when the expression is evaluated.
#' @return A list with elements `rv`, 'runtime`, `output`, `message`, `plots`, `testdir` and `inputs`
#' Element `rv` contains the return value of the evaluated expression.
#' Element `runtime` contains the runtime of the evaluated expression.
#' Elements `output`, `message` and `plots` are environments containing the captured output, message and plots, respectively. For further details see [redirect()].
#' Element `testdir` contains the path to the test directory.
#' Element `inputs` equals the `inputs` argument.
#' @noRd
with <- function(expr,
                 testdir = NULL,
                 answers = NULL,
                 output = NULL,
                 message = NULL,
                 plots = NULL,
                 datadir_temp = NULL,
                 datadir_persistent = NULL,
                 inputs = character(),
                 debug = FALSE,
                 ...) {
    on.exit(restore(), add = TRUE)
    push_option(...)
    testdir <- push_testdir(testdir)
    inputs <- fill_with_inputs(testdir, inputs)
    mocked_datadir_temp <- get_datadir_mock(type = "temp", state = datadir_temp)
    mocked_datadir_persistent <- get_datadir_mock(type = "persistent", state = datadir_persistent)
    mocked_readline <- get_readline_mock(answers)
    redirect_list <- redirect(output = output, message = message, plots = plots)
    errors <- character()
    withCallingHandlers(
        {
            if (debug) browser()
            testthat::with_mocked_bindings(
                tryCatch(
                    {
                        rv_runtime_list <- measure_runtime(expr)
                    },
                    error = function(e) {
                        restore(streams = "message")
                        stop(e)
                    }
                ),
                datadir_temp = mocked_datadir_temp,
                datadir_persistent = mocked_datadir_persistent,
                readline = mocked_readline
            )
        },
        warning = print_warning_as_message
    )
    return(c(rv_runtime_list, redirect_list, list(testdir = testdir, inputs = inputs)))
}

#' @title Redirect output, message, and plot streams
#' @description Redirects the output, message, and plot streams in R to either a specified file or a captured character vector.
#' If the target is "captured", the stream is captured into a character vector.
#' If the target is NULL, the stream is not redirected.
#' Otherwise, the target is assumed to be a file path, and the stream is redirected to that file.
#' @param output The target for the output stream. Defaults to "captured".
#' @param message The target for the message stream. Defaults to NULL.
#' @param plots The target for the plot stream. Defaults to NULL. If a file path is provided, plots will be saved to a PDF at that location.
#' @return A list with 3 environments with name `output`, `message` and `plots` that each have the following elements:
#'
#' -  conn:   e.g. `'file' int 3` or `Named int 2` (for plots)
#' -  path:   e.g. "output.txt" or "plots.pdf"
#' -  penv:   e.g. <environment: 0x0000019fb9fc7318>
#' -  stream: either "output", "message" or "plots"
#' -  target: e.g. "output.txt" or "plots.pdf" or "captured"
#' -  text:   e.g. `c("Hello", "World")` or `NULL`
#'
#' @examples
#' # Capture output
#'
#' redirects <- redirect(output = "captured")
#' cat("Hello\n")
#' cat("World\n")
#' restore()
#' ls.str(redirects$output)
#'
#' # Capture messages, redirect output to a file, and save plots to a PDF
#'
#' redirects <- redirect(output = "output.txt", message = "captured", plots = "plots.pdf")
#' print("This goes to output.txt")
#' message("This is captured")
#' warning("This is captured as well")
#' plot(1:10) # This plot is saved to plots.pdf
#' restore()
#' cat(redirects$message$text)
#' @noRd
redirect <- function(output = "captured", message = NULL, plots = NULL) {
    streams <- c("output", "message", "plots")
    targets <- list(output = output, message = message, plots = plots)
    redirects <- lapply(streams, function(stream) {
        if (!is.null(penv$open_conns[[stream]])) {
            stop(sprintf("Stream '%s' is already redirected. Call `restore()` first.", stream))
        }
        target <- targets[[stream]]
        if (is.null(target)) {
            return(environment())
        }
        path <- if (target == "captured") NULL else target
        text <- if (target != "captured") NULL else vector("character")
        if (stream == "plots") {
            grDevices::pdf(target)
            conn <- grDevices::dev.cur()
        } else {
            conn <- if (target == "captured") textConnection("text", "wr", local = TRUE) else file(target, open = "wt")
            sink(conn, type = stream)
            conn
        }
        penv$open_conns[[stream]] <- conn
        return(environment())
    })
    redirects <- stats::setNames(redirects, streams)
    redirects
}

#' @name mock_readline
#' @title Mock the readline function
#' @description This function mocks the readline function to return predefined answers.
#' @param answers A list of answers that the mocked readline function will return.
#' @return Invisible NULL.
#' @noRd
mock_readline <- function(answers) {
    if (!is.null(answers)) {
        readline_mock <- get_readline_mock(answers)
        patch("readline", readline_mock)
    }
}

#' @name mock_datadir
#' @title Mock the datadir function
#' @description This function mocks the datadir function to return paths to fake directories for temporary and persistent data.
#' @param type The type of data directory to mock. Can be "temp" or "persistent".
#' @param state The state of the data directory to mock. Can be "missing", "empty", or "filled".
#' @return Invisible NULL.
#' @noRd
mock_datadir <- function(type = c("temp", "persistent"), state = c("missing", "empty", "filled")) {
    if (is.null(state) || is.null(type)) {
        return()
    }
    type <- match.arg(type)
    state <- match.arg(state)
    if (type == "persistent") {
        datadir_persistent_mock <- get_datadir_mock(type = "persistent", state = state)
        patch("datadir_persistent", datadir_persistent_mock)
    } else {
        datadir_temp_mock <- get_datadir_mock(type = "temp", state = state)
        patch("datadir_temp", datadir_temp_mock)
    }
}

#' @title Restore mocked functions, redirected streams and the working directory
#' @description Restores all mocked functions, redirected streams and the working directory to their original state.
#' Usually called after [redirect()], [mock_readline()] or [mock_datadir].
#' Also used internally by [with()].
#' @param fns A character vector of function names to restore. If NULL, all functions are restored.
#' @param streams A character vector of streams to restore. Can be "output", "message", and/or "plots".
#' @param wd Logical. If TRUE, the working directory is restored by calling [popd()] with option `all=TRUE`.
#' @return Invisible NULL.
#' @noRd
restore <- function(fns = NULL,
                    streams = c("output", "message", "plots"),
                    wd = TRUE,
                    opts = NULL) {
    fns <- if (is.null(fns)) names(penv$fn_backups) else fns
    opts <- if (is.null(opts)) names(penv$option_backupss) else opts
    lapply(fns, restore_fn)
    lapply(streams, restore_stream)
    lapply(opts, restore_option)
    if (wd) popd(all = TRUE)
    invisible()
}

# Helpers #####

#' @title Patch a single function inside metabodecons namespace
#' @description
#' Replace a function in the metabodecon namespace with a replacement function.
#' This can be useful for testing.
#' The original function can be restored using [restore()].
#' The original function is backed up in the `penv$fn_backups` environment.
#' To list all currently patched functions, use `ls.str(penv$fn_backups)`.
#' Used internally by [mock_datadir()] and [mock_readline()].
#' @param fn The name of the function to replace, as a string.
#' @param repl The replacement function.
#' @examples
#' # Replace the `foo` function with a function that always returns 1
#' patch("foo", function() 1)
#'
#' # Restore the original `foo` function
#' restore("foo")
#' @details This function only patches the object inside `namespace:metabodecon`, but not inside `package:metabodecon` (exported objects), i.e. calling `fn` from outside the package will still call the original function. But functions defined inside the package will use the patched function.
#' @noRd
patch <- function(fn, repl) {
    if (is.null(penv$fn_backups[[fn]])) {
        penv$fn_backups[[fn]] <- get(fn, envir = asNamespace("metabodecon"))
    } else {
        stop(sprintf("Function '%s' is already patched. Call `restore()` first.", fn))
    }
    assign(fn, repl, pos = "package:metabodecon") # package:metabodecon contains the exported functions from the package, this is what is called when you enter `fn()` in the console
    utils::assignInNamespace(fn, repl, ns = "metabodecon") # namespace:metabodecon contains all functions defined in the package, this is what is used by other function defined in the package
}

#' @title Restore a single mocked function
#' @description Restores a single mocked function to its original state. Used internally by [restore()].
#' @param fn The name of the function to restore.
#' @return Invisible NULL.
#' @noRd
restore_fn <- function(fn) {
    if (!is.null(penv$fn_backups[[fn]])) {
        assign(fn, penv$fn_backups[[fn]], pos = "package:metabodecon")
        assignInNamespace(fn, penv$fn_backups[[fn]], ns = "metabodecon")
        penv$fn_backups[[fn]] <- NULL
    }
}

#' @title Restore a single redirected stream
#' @description Restores a single redirected stream to its original state. Used internally by [restore()].
#' @param stream The name of the stream to restore. Can be "output", "message", or "plots".
#' @return Invisible NULL.
#' @noRd
restore_stream <- function(stream) {
    conn <- penv$open_conns[[stream]]
    if (!is.null(conn)) {
        if (stream == "plots") {
            dev.off(conn)
        } else {
            sink(NULL, type = stream)
            close(conn)
        }
        penv$open_conns[[stream]] <- NULL
    }
}

restore_option <- function(opt) {
    if (!is.null(penv$option_backupss[[opt]])) {
        options(opt = penv$option_backupss[[opt]])
        penv$option_backupss[[opt]] <- NULL
    }
}

restore_warnings <- function() {
    restore_option("warn")
    restore_option("warning.expression")
}

#' @title Creates a mock readline function for testing
#' @description Creates a mock readline function that returns the next element from a character vector each time it's called.
#' Used internally by [mock_readline()].
#' @param texts A character vector of responses to be returned by the readline function.
#' @return A function that mimics the readline function, returning the next element from `texts` each time it's called.
#' @examples \dontrun{
#' readline_mock <- get_readline_mock(c("yes", "no", "maybe"))
#' readline_mock("Continue? ") # Returns "yes"
#' readline_mock("Continue? ") # Returns "no"
#' readline_mock("Continue? ") # Returns "maybe"
#' }
#' @noRd
get_readline_mock <- function(texts, env = as.environment(list())) {
    if (is.null(texts)) {
        return(readline)
    }
    env$readline_called <- 0
    readline <- function(prompt = "") {
        env$readline_called <- env$readline_called + 1
        if (env$readline_called > length(texts)) {
            msg <- "readline called %s times, but only %s answers were provided."
            stop(sprintf(msg, env$readline_called, length(texts)))
        }
        message(prompt, appendLF = FALSE)
        message(texts[env$readline_called])
        return(texts[env$readline_called])
    }
}

#' @title Get a mock for the datadir functions
#' @description This function returns a function that when called, returns a path to a mock data directory.
#' The type and state of the mock data directory can be specified.
#' Used internally by [mock_datadir()].
#' @param type The type of data directory to mock. Can be "persistent" or "temp".
#' @param state The state of the data directory to mock. Can be "missing", "empty", or "filled".
#' @return A function that when called, returns a path to the mock data directory.
#' @examples
#' datadir_persistent_mock <- get_datadir_mock(type = "persistent", state = "missing")
#' datadir_temp_mock <- get_datadir_mock(type = "temp", state = "empty")
#' patch("datadir_persistent", datadir_persistent_mock)
#' patch("datadir_temp", datadir_temp_mock)
#' datadir_persistent()
#' datadir_temp()
#' restore()
#' @noRd
get_datadir_mock <- function(type = c("persistent", "temp"),
                             state = c("missing", "empty", "filled")) {
    if (is.null(penv$cached_zip)) {
        # Make sure we now the path of the cached zip file even after mocking the original functions.
        # Only call this on the first invocation, otherwise cache_example_datasets() might already use the mocked data dirs.
        penv$cached_zip <- cache_example_datasets()
    }
    if (is.null(state) && type == "persistent") {
        return(datadir_persistent)
    }
    if (is.null(state) && type == "temp") {
        return(datadir_temp)
    }
    type <- match.arg(type)
    state <- match.arg(state)
    p <- file.path(mockdir(), "datadir", type, state)
    p <- normalizePath(p, "/", mustWork = FALSE)
    switch(state,
        "missing" = unlink(p, recursive = TRUE, force = TRUE),
        "empty" = clear(p),
        "filled" = fill_with_example_datasets(dst_dir = p, src_zip = penv$cached_zip)
    )
    function() p
}

#' @name push_testdir
#' @title Push a test directory to the stack of working directories
#' @description This function creates a new test directory and sets it as working directory using [pushd()]. Used internally by [with()].
#' @param testdir The path of the test directory relative to [testdir()].
#' @return Path to old working directory or NULL if `testdir` is NULL.
#' @noRd
push_testdir <- function(testdir) {
    if (!is.null(testdir)) {
        wd <- file.path(testdir(), testdir)
        mkdirs(wd)
        pushd(wd)
        wd
    }
}

push_option <- function(...) {
    newopts <- list(...)
    oldopts <- options(newopts)
    penv$option_backups <- modifyList(penv$option_backups, oldopts)
    penv$option_backups
}

warnings_as_messages <- function() {
    warning.expression <- quote({
        message("Warning: ", conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    push_option(warn = 1)
    push_option(warning.expression = warning.expression)
}

fill_with_example_datasets <- function(dst_dir, src_zip) {
    dst_zip <- file.path(dst_dir, "example_datasets.zip")
    dst_zip_has_correct_size <- isTRUE(file.size(dst_zip) == xds$zip_size) # Implies existence.
    dst_zip_has_wrong_size <- !dst_zip_has_correct_size

    dst_subdir <- mkdirs(file.path(dst_dir, "example_datasets"))
    dst_subdir_has_enough_files <- isTRUE(length(dir(dst_dir, recursive = TRUE)) >= 1018)
    dst_subdir_is_missing_files <- !dst_subdir_has_enough_files
    # Dont check size because it takes too long (3.36s on first run on Win11 with fast SSD). Listing files can be done in 0.07s.

    if (dst_zip_has_wrong_size) {
        file.copy(src_zip, dst_zip, overwrite = TRUE)
    }
    if (dst_subdir_is_missing_files) {
        unlink(dst_subdir, recursive = TRUE, force = TRUE) # 0.45s
        utils::unzip(dst_zip, exdir = dst_dir) # 1.41s
    }
}

measure_runtime <- function(expr) {
    start_time <- Sys.time()
    rv <- expr
    end_time <- Sys.time()
    runtime <- end_time - start_time
    return(list(rv = rv, runtime = runtime))
}

print_warning_as_message <- function(w) {
    message("Warning: ", conditionMessage(w))
    invokeRestart("muffleWarning")
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

fill_with_inputs <- function(dstdir, inputs = c()) {
    if (!is.null(dstdir) && length(inputs) >= 0) {
        src_zip <- cache_example_datasets(extract = TRUE)
        src_dir <- gsub(".zip", "", src_zip, fixed = TRUE)
        bruker_dir <- file.path(src_dir, "bruker")
        jcampdx_dir <- file.path(src_dir, "jcampdx")
        inputs <- gsub("bruker", bruker_dir, inputs, fixed = TRUE)
        inputs <- gsub("jcampdx", jcampdx_dir, inputs, fixed = TRUE)
        file.copy(from = inputs, to = dstdir, recursive = TRUE)
    }
    inputs
}

readline <- function(...) {
    base::readline(...) # we must have our own copy of readline in the package namespace so we can mock it in tests
}

testdir <- function() {
    p <- file.path(tempdir(), "tests")
    normalizePath(p, "/", mustWork = FALSE)
}

mockdir <- function() {
    p <- file.path(tempdir(), "mocks")
    normalizePath(p, "/", mustWork = FALSE)
}
