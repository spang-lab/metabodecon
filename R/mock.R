# Main #####

#' @description Run expression with predefined global state
#' @param expr Expression to be evaluated.
#' @param testdir ID of the test directory.
#' @param answers Answers to be returned by readline().
#' @param output Outputs to be captured. Can be "stdout" and/or "stderr".
#' @param message Path to the file where stdout should be redirected to.
#' @param plots Path to the pdf file where plots should be saved to.
#' @param datadir_temp State of the mocked temporary data directory. See details section.
#' @param datadir_persistent State of the mocked persistent data directory. See details section.
#' @details The `datadir_temp` and `datadir_persistent` arguments accept values "missing", "filled" and "empty".
#' Setting a value unequal NULL causes the functions [datadir_temp()] and/or [datadir_persistent()] to be replaced with mocks functions pointing to fake directories. Functions depending on these functions will then use the fake directories instead of the real ones.
#' When set to "missing" the returned mock directory does not exist.
#' When set to "empty" it exists and is guaranteed to be empty.
#' When set to "filled", it is populated with example datasets.
#' Attention: the mocked functions, i.e. [datadir_temp()] and [datadir_persistent()] cannot be used directly inside the expression when called via `devtools::test()`. I'm not sure why, but it seems as if devtools and/or testthat have their own copies of the functions which are used when the expression is evaluated.
#' @return A list with elements `rv`, `output`, `message` and `plots`.
#' Element `rv` contains the return value of the evaluated expression.
#' Elements `output`, `message` and `plots` are environments containing the captured output, message and plots, respectively. For further details see [redirect()].
#' @noRd
with <- function(expr,
                 testdir = NULL,
                 answers = NULL,
                 output = NULL,
                 message = NULL,
                 plots = NULL,
                 datadir_temp = NULL,
                 datadir_persistent = NULL) {
    on.exit(restore(), add = TRUE)
    mock_readline(answers)
    mock_datadir(type = "temp", state = datadir_temp)
    mock_datadir(type = "persistent", state = datadir_persistent)
    push_testdir(testdir)
    retlist <- redirect(output = output, message = message, plots = plots)
    withCallingHandlers({
        retlist$rv <- expr
    }, warning = function(w) {
        message("Warning: ", conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    retlist
}

#' @title Redirect output, message, and plot streams
#' @description Redirects the output, message, and plot streams in R to either a specified file or a captured character vector.
#' If the target is "captured", the stream is captured into a character vector.
#' If the target is NULL, the stream is not redirected.
#' Otherwise, the target is assumed to be a file path, and the stream is redirected to that file.
#' @param output The target for the output stream. Defaults to "captured".
#' @param message The target for the message stream. Defaults to NULL.
#' @param plots The target for the plot stream. Defaults to NULL. If a file path is provided, plots will be saved to a PDF at that location.
#' @return A list of environments, one for each stream, containing the redirection settings.
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
            conn <- dev.cur()
        } else {
            conn <- if (target == "captured") textConnection("text", "wr", local = TRUE) else file(target, open = "wt")
            sink(conn, type = stream)
            conn
        }
        penv$open_conns[[stream]] <- conn
        return(environment())
    })
    redirects <- setNames(redirects, streams)
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
    opts <- if (is.null(opts)) names(penv$option_backups) else opts
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
    assignInNamespace(fn, repl, ns = "metabodecon") # namespace:metabodecon contains all functions defined in the package, this is what is used by other function defined in the package
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
    if (!is.null(penv$option_backups[[opt]])) {
        options(opt = penv$option_backups[[opt]])
        penv$option_backups[[opt]] <- NULL
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
    env$readline_called <- 0
    readline <- function(prompt = "") {
        env$readline_called <- env$readline_called + 1
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
#' @return Invisible NULL.
#' @noRd
push_testdir <- function(testdir) {
    if (!is.null(testdir)) {
        new_wd <- file.path(testdir(), testdir)
        mkdirs(new_wd)
        pushd(new_wd)
    }
}

push_option <- function(...) {
  options_list <- list(...)
  browser()
  for (option in names(options_list)) {
    penv$option_backups[[option]] <- getOption(option)
  }
  do.call(options, options_list)
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
    dst_subdir_has_enough_files <- isTRUE(length(dir(dst_dir, recursive=TRUE)) >= 1018)
    dst_subdir_is_missing_files <- !dst_subdir_has_enough_files
    # Dont check size because it takes too long (3.36s on first run on Win11 with fast SSD). Listing files can be done in 0.07s.

    if (dst_zip_has_wrong_size) {
        file.copy(src_zip, dst_zip, overwrite = TRUE)
    }
    if (dst_subdir_is_missing_files) {
        unlink(dst_subdir, recursive = TRUE, force = TRUE) # 0.45s
        unzip(dst_zip, exdir = dst_dir) # 1.41s
    }
}

# Deprecated #####

#' @title Evaluate expression with redirected output streams
#' @param stdout (string) path where stdout should be redirected to
#' @param stderr (string) path where stderr should be redirected to
#' @param plots (string) path where plots should be saved to
#' @param expr (expression) expression to be evaluated
#' @return The result of evaluating `expr`
#' @examples \dontrun{
#' with_redirects(stdout = "output.txt", plots = "plots.pdf", expr = {
#'     print("Starting plotting")
#'     plot(1:10)
#'     print("Done plotting")
#' })
#' }
#' @noRd
with_redirects <- function(stdout = NULL, stderr = NULL, plots = NULL, expr) {
    if (!is.null(stdout)) {
        sink(stdout, type = "output")
        on.exit(sink(file = NULL, type = "output"), add = TRUE)
    }
    if (!is.null(stderr)) {
        stderr_stream <- file(stderr, open = "wt")
        sink(stderr_stream, type = "message")
        on.exit(sink(file = NULL, type = "message"), add = TRUE)
        on.exit(close(stderr_stream), add = TRUE)
    }
    if (!is.null(plots)) {
        grDevices::pdf(plots)
        on.exit(grDevices::dev.off(), add = TRUE)
    }
    expr
}
