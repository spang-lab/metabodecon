# evalwith #####

#' @export
#' @title Evaluate an expression with predefined global state
#' @description
#' Evaluates an expression with a predefined global state, including the:
#'
#' - working directory (set via [setwd()])
#' - global options (set via [options()])
#' - graphical parameters (set via [par()])
#'
#' In addition to that, `evalwith` allows to:
#'
#' - Redirect or capture the output and/or message stream via [sink()]
#' - Measure the runtime of the evaluated expression via [system.time()]
#' - Creating a temporary test directory (inside [tmpdir()]) and populating it
#'   with input files according to `inputs`
#' - Predefine answers for calls to [readline()] happening during evaluation of
#'   `expr`
#' - Caching the result of the expression
#'
#' All changes to the global state are reverted after the expression has been
#' evaluated.
#'
#' @param expr Expression to be evaluated.
#'
#' @param testdir ID of the test directory. E.g. `"xyz/2"`. Will be created and
#' populated with `inputs`. To clear, use `clear(testdir("xyz/2"))`.
#'
#' @param answers Answers to be returned by readline().
#'
#' @param output Path to the file where output stream should be redirected to.
#' Use `"captured"` to capture the output.
#'
#' @param message Path to the file where message stream be redirected to. Use
#' `"captured"` to capture the messages.
#'
#' @param plot An expression opening a device, the string "captured" or a path
#' ending in ".pdf", ".svg", or ".png". Examples: `svg("tmp.svg")`,
#' `quote(pdf("tmp.pdf"))`, `"captured"`, `"tmp.png"`. Passing `"captured"` is
#' equivalent to passing `tempfile(fileext = ".png")`.
#'
#' @param datadir_temp State of the mocked temporary data directory. See details
#' section.
#'
#' @param datadir_persistent State of the mocked persistent data directory. See
#' details section.
#'
#' @param inputs Paths to be copied to the test directory before evaluating
#' `expr`.
#'
#' @param opts Named list of options to be set. See [options()].
#'
#' @param pars Named list of parameters to be set. See [par()].
#'
#' @param cache Logical indicating whether to cache the result of the
#' expression.
#'
#' @param overwrite Logical indicating whether to overwrite the cache file if it
#' already exists.
#'
#' @details
#' The `datadir_temp` and `datadir_persistent` arguments accept values
#' "missing", "filled" and "empty". Setting a value unequal NULL causes the
#' functions [datadir_temp()] and/or [datadir_persistent()] to be replaced with
#' mock functions pointing to fake directories. Functions depending on these
#' functions will then use the fake directories instead of the real ones. When
#' set to "missing" the returned mock directory does not exist. When set to
#' "empty" it exists and is guaranteed to be empty. When set to "filled", it is
#' populated with example datasets.
#'
#' Attention: the mocked functions, i.e. [datadir_temp()] and
#' [datadir_persistent()] cannot be used directly inside `expr` when called via
#' `devtools::test()`. I'm not sure why, but it seems as if devtools and/or
#' testthat have their own copies of the functions which are used when the
#' expression is evaluated.
#'
#' @return
#' A list containing with following elements:
#'
#' - `rv`: The return value of the expression.
#' - `runtime`: The "elapsed" runtime of the expression in seconds. Measured
#'   with [system.time()].
#' - `output`: The captured output.
#' - `message`: The captured messages.
#' - `plot`: The path to the saved plot.
#' - `testdir`: The path to the test directory.
#' - `inputs`: The paths to the copied input files.
#'
#' @examples
#' x1 <- evalwith(output = "captured", cat("Helloworld\n"))
#' str(x1)
#'
#' x2 <- evalwith(datadir_persistent = "missing", message = "captured", datadir())
#' str(x2)
#'
#' x3 <- evalwith(testdir = "dummy", inputs = "bruker/urine/urine_1", dir())
#' str(x3)
#'
#' x4 <- evalwith(Sys.sleep(0.02))
#' str(x4)
evalwith <- function(expr, # nolint: cyclocomp_linter.
                     testdir = NULL,
                     answers = NULL,
                     output = NULL,
                     message = NULL,
                     plot = NULL,
                     datadir_temp = c("default", "missing", "empty", "filled")[1],
                     datadir_persistent = c("default", "missing", "empty", "filled")[1],
                     inputs = character(),
                     opts = NULL,
                     pars = NULL,
                     cache = FALSE,
                     overwrite = FALSE) {
    if (isTRUE(cache)) {
        cachedir <- cachedir()
        cachefile <- file.path(cachedir, paste0(testdir, ".rds"))
        if (file.exists(cachefile) && isFALSE(overwrite)) {
            return(readRDS(cachefile))
        }
    }
    if (!is.null(testdir)) {
        testpath <- file.path(testdir(), testdir)
        mkdirs(testpath)
        owd <- setwd(testpath)
        on.exit(setwd(owd), add = TRUE)
        if (!is.null(inputs)) {
            pkg_inputpaths <- sapply(paste0("example_datasets/", inputs), pkg_file)
            if (any(pkg_inputpaths == "")) {
                xds_inputs <- inputs[pkg_inputpaths == ""]
                src_dir <- download_example_datasets()
                bruker_dir <- file.path(src_dir, "bruker")
                jcampdx_dir <- file.path(src_dir, "jcampdx")
                xds_inputpaths <- gsub("bruker", bruker_dir, xds_inputs, fixed = TRUE)
                xds_inputpaths <- gsub("jcampdx", jcampdx_dir, xds_inputpaths, fixed = TRUE)
            } else {
                xds_inputpaths <- c()
            }
            inputpaths <- c(pkg_inputpaths, xds_inputpaths)
            file.copy(from = inputpaths, to = testpath, recursive = TRUE)
        }
    }
    outvec <- vector("character")
    if (!is.null(output)) {
        outcon <- if (output == "captured") textConnection("outvec", "wr", local = TRUE) else file(output, open = "wt")
        sink(outcon, type = "output")
        on.exit(sink(NULL), add = TRUE)
        on.exit(close(outcon), add = TRUE)
    }
    msgvec <- vector("character")
    if (!is.null(message)) {
        msgcon <- if (message == "captured") textConnection("msgvec", "wr", local = TRUE) else file(message, open = "wt")
        sink(msgcon, type = "message")
        on.exit(sink(NULL, type = "message"), add = TRUE)
        on.exit(close(msgcon), add = TRUE)
    }
    dev_cur <- dev.cur()
    force(plot) # forces eval, useful if plot ~= `svg("abc.svg")`
    if (identical(dev_cur, dev.cur()) && !is.null(plot)) {
        if (is.expression(plot)) eval(plot) # plot == `quote(svg("abc.svg", width = 10))`
        if (identical(plot, "captured")) plot <- tempfile(fileext = ".png")
        if (grepl("\\.pdf$", plot)) {
            pdf(plot)
        } else if (grepl("\\.svg$", plot)) {
            svg(plot)
        } else if (grepl("\\.png$", plot)) {
            png(plot)
        } else {
            stop("plot must be an expression opening a device or a path ending in .pdf, .svg, or .png")
        }
    }
    if (!is.null(opts)) {
        oldopts <- options(opts)
        on.exit(options(oldopts), add = TRUE)
    }
    if (!is.null(pars)) {
        oldpars <- par(pars)
        on.exit(par(oldpars), add = TRUE)
    }
    if (!identical(dev_cur, dev.cur())) {
        # This must be done after both arguments `plot` and `pars` have been evaluated, as both can lead to a change in the graphical device.
        on.exit(while (!identical(dev_cur, dev.cur())) dev.off(), add = TRUE)
    }
    withCallingHandlers(
        testthat::with_mocked_bindings(
            code = {
                tryCatch(
                    {
                        runtime <- system.time(rv <- expr)[["elapsed"]]
                    },
                    error = function(e) {
                        sink(NULL, type = "message")
                        stop(e)
                    }
                )
            },
            datadir_temp = get_datadir_mock(type = "temp", state = datadir_temp),
            datadir_persistent = get_datadir_mock(type = "persistent", state = datadir_persistent),
            readline = get_readline_mock(answers),
            .package = if (loaded_via_devtools()) NULL else "metabodecon"
        ),
        warning = function(w) {
            message("Warning: ", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    retobj <- invisible(list(
        rv = rv, runtime = runtime,
        output = outvec, message = msgvec, plot = plot,
        testdir = testdir, inputs = inputs
    ))
    if (isTRUE(cache) && (!file.exists(cachefile) || isTRUE(overwrite))) saveRDS(retobj, cachefile)
    invisible(retobj)
}

# evalwith Helpers #####

#' @noRd
#' @title Creates a mock readline function for testing
#' @description Creates a mock readline function that returns the next element
#' from a character vector each time it's called.
#' Used internally by [mock_readline()].
#' @param texts A character vector of responses to be returned by the readline
#' function.
#' @return A function that mimics the readline function, returning the next
#' element from `texts` each time it's called.
#' @examples
#' readline_mock <- get_readline_mock(c("yes", "no", "maybe"))
#' readline_mock("Continue? ") # Returns "yes"
#' readline_mock("Continue? ") # Returns "no"
#' readline_mock("Continue? ") # Returns "maybe"
#' try(readline_mock("Continue? ")) # Throws error "readline called 4 times, but only 3 answers were provided"
get_readline_mock <- function(texts, env = as.environment(list())) {
    if (is.null(texts)) {
        return(readline)
    }
    env$readline_called <- 0
    readline <- function(prompt = "") {
        env$readline_called <- env$readline_called + 1
        message(prompt, appendLF = FALSE)
        if (env$readline_called > length(texts)) {
            msg <- "readline called %s times, but only %s answers were provided."
            stop(sprintf(msg, env$readline_called, length(texts)))
        }
        message(texts[env$readline_called])
        return(texts[env$readline_called])
    }
}

#' @noRd
#' @title Get a mock for the datadir functions
#' @description Returns a function that, when called, returns a path to a mock
#' data directory. The type and state of the mock data directory can be
#' specified. Used internally by [mock_datadir()].
#' @param type The type of data directory to mock. Can be "persistent" or
#' "temp".
#' @param state The state of the data directory to mock. Can be "missing",
#' "empty", or "filled".
#' @return A function that when called, returns a path to the mock data
#' directory.
#' @examples
#' datadir_persistent_mock <- get_datadir_mock(type = "persistent", state = "missing")
#' datadir_temp_mock <- get_datadir_mock(type = "temp", state = "empty")
#' datadir_persistent_mock()
#' datadir_temp_mock()
get_datadir_mock <- function(type = "temp", state = "default") {
    type <- match.arg(type, c("temp", "persistent"))
    state <- match.arg(state, c("default", "missing", "empty", "filled"))
    if (state == "default" && type == "persistent") {
        return(datadir_persistent)
    }
    if (state == "default" && type == "temp") {
        return(datadir_temp)
    }
    p <- norm_path(file.path(mockdir(), "datadir", type, state))
    if (state %in% c("missing", "empty")) unlink(p, recursive = TRUE, force = TRUE)
    if (state == "empty") mkdirs(p)
    if (state == "filled") download_example_datasets(dst_dir = p)
    function() p
}

loaded_via_devtools <- function() {
    pkg_dir <- dirname(system.file("DESCRIPTION", package = "metabodecon"))
    loaded_via_devtools <- dir.exists(file.path(pkg_dir, "inst"))
    return(loaded_via_devtools)
}

# testthat #####

#' @noRd
#' @title Run tests with the option to skip slow tests
#' @description Runs the tests in the current R package. If `all` is TRUE, it
#' well set environment variable `RUN_SLOW_TESTS` to "TRUE" so that all tests
#' are run. If `all` is FALSE, it will set `RUN_SLOW_TESTS` to "FALSE" so that
#' slow tests are skipped.
#' @param all Logical. If TRUE, all tests are run. If FALSE, slow tests are
#' skipped.
#' @return The result of devtools::test()
#' @examples
#' if (interactive()) {
#'     run_tests(all = FALSE)
#' }
run_tests <- function(all = FALSE) {
    RUN_SLOW_TESTS_OLD <- Sys.getenv("RUN_SLOW_TESTS")
    Sys.setenv(RUN_SLOW_TESTS = if (all) "TRUE" else "FALSE")
    on.exit(Sys.setenv(RUN_SLOW_TESTS = RUN_SLOW_TESTS_OLD), add = TRUE)
    devtools::test()
}

skip_if_slow_tests_disabled <- function() {
    if (!Sys.getenv("RUN_SLOW_TESTS") == "TRUE") {
        testthat::skip("Slow tests are disabled. Use `Sys.setenv(RUN_SLOW_TESTS=TRUE)` or `run_tests(all=TRUE)` to enable.")
    }
}

#' @title Check if the size of each file in a directory is within a certain
#' range
#' @description Check if the size of each file in a directory is within 90% to
#' 110% of the expected size.
#' If a file size is not within this range, a message is printed and an error is thrown.
#' @param testdir A character string specifying the directory to check.
#' @param size_exp A named numeric vector where the names are filenames and the
#' values are the expected file sizes.
#' @examples
#' \dontrun{
#' testdir <- tmpdir()
#' file.create(file.path(testdir, "file1.txt"))
#' file.create(file.path(testdir, "file2.txt"))
#' size_exp <- c(file1.txt = 100, file2.txt = 200)
#' expect_file_size(testdir, size_exp)
#' }
#' @noRd
expect_file_size <- function(testdir, size_exp) {
    paths <- file.path(testdir, names(size_exp))
    size_obs <- file.info(paths)$size
    file_has_correct_size <- isTRUE(size_obs > size_exp * 0.9 & size_obs < size_exp * 1.1)
    lapply(seq_along(size_exp), function(i) {
        if (!isTRUE(file_has_correct_size[i])) {
            message(sprintf("Size of %s is %s which is not between %s and %s", paths[i], size_obs[i], size_exp[i] * 0.9, size_exp[i] * 1.1))
        }
    })
    testthat::expect_true(all(file_has_correct_size))
}

#' @title Expect Structure
#' @description Tests if the structure of an object matches the expected string
#' @param obj The object to test
#' @param expected_str The expected structure of the object as a string. Can be
#' obtained by calling `dput(capture.output(str(obj)))`.
#' @return A logical value indicating whether the structure of the object
#' matches the expected string
#' @examples
#' expect_str(list(a = 1, b = 2), c("List of 2", " $ a: num 1", " $ b: num 2"))
#' @noRd
expect_str <- function(obj, expected_str) {
    testthat::expect_identical(capture.output(str(obj)), expected_str)
}
