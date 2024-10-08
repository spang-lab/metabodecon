# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Evalwith #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

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

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Evalwith Helpers #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @noRd
#' @title Creates a mock readline function for testing
#'
#' @description
#' Creates a mock readline function that returns the next element from a
#' character vector each time it's called. Used internally by [mock_readline()].
#'
#' @param texts A character vector of responses to be returned by the readline
#' function.
#'
#' @return
#' A function that mimics the readline function, returning the next element from
#' `texts` each time it's called.
#'
#' @examples
#' readline_mock <- get_readline_mock(c("yes", "no", "maybe"))
#' readline_mock("Continue? ") # Returns "yes"
#' readline_mock("Continue? ") # Returns "no"
#' readline_mock("Continue? ") # Returns "maybe"
#' try(readline_mock("Continue? ")) # Throws error
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

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Testthat #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

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

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Misc #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @noRd
#'
#' @title Calculate the PRARP Score
#'
#' @description
#' Calculates the PRARP score for a deconvolution. The PRARP score is the
#' product of the peak ratio and the area ratio and can be used to assess the
#' quality of a deconvolution. See 'Details' for more information on how the
#' score is calculated.
#'
#' @param decon A list containing the deconvolution results, as returned by
#' [generate_lorentz_curves()].
#'
#' @param lcpar A data frame containing the true parameters of the peaks.
#'
#' @return The PRARP score as numeric scalar. In addition, a plot is created to
#' visualize the deconvolution results.
#'
#' @details
#' The function first plots the deconvolution results for visual inspection and
#' then returns the PRARP score for the deconvolution.
#'
#' The plotting is done as follows:
#'
#' 1. Plot the deconvoluted spectrum using `plot_spectrum()`.
#' 2. Draw green circles around found peaks [^1].
#' 3. Draw red circles around missed peaks [^1].
#' 4. Draw red rectangles around falsely detected peaks. [^2]
#'
#' [^1]: we consider a peak as 'found' if there is at least one detected peak
#'       center within 0.001 ppm of the true peak position. If this is not the
#'       case, the peak is considered as 'missed'.
#' [^2]: we consider a peak as 'falsely detected' if there is no true peak
#'       center within 0.001 ppm of the detected peak position.
#'
#' In addition, a quality score is calculated as follows:
#'
#' quality    = peak_ratio * area_ratio
#' peak_ratio = min(peaks_true, peaks_found) / max(peaks_true, peaks_found)
#' area_ratio = min(area_true,  area_found)  / max(area_true,  area_found)
#'
#' I.e., the score is close to 1 if the number of peaks and the area of the
#' peaks are similar in the true and found spectra and the score is close to 0
#' if the number of peaks and/or the area of the peaks are very different.
#'
#' @examples
#' ## Bad deconvolution (PRARP ~= 0.2)
#' decon <- generate_lorentz_curves_sim(sim[[1]], delta = 6.4)
#' truepar <- sim[[1]]$meta$simpar[c("A", "x0", "lambda")]
#' calc_prarp(decon, truepar, show = TRUE)
#'
#' ## Good deconvolution (PRARP ~= 0.64)
#' decon <- generate_lorentz_curves_sim(sim[[1]], delta = 0)
#' truepar <- sim[[1]]$meta$simpar[c("A", "x0", "lambda")]
#' calc_prarp(decon, truepar, show = TRUE)
#'
calc_prarp <- function(decon, truepar, show = FALSE) {

    n_peaks_dcnv <- length(decon$A)
    n_peaks_true <- length(truepar$A)
    n_peaks_min <- min(n_peaks_dcnv, n_peaks_true)
    n_peaks_max <- max(n_peaks_dcnv, n_peaks_true)
    peak_ratio <- n_peaks_min / n_peaks_max

    area_spectrum <- sum(decon$y_values)
    area_residuals <- sum(abs(decon$spectrum_superposition - decon$y_values))
    area_min <- min(area_residuals, area_spectrum)
    area_max <- max(area_residuals, area_spectrum)
    area_ratio <- 1 - (area_min / area_max)

    prarp <- peak_ratio * area_ratio
    if (show) plot_prarp(decon, truepar, prarp, peak_ratio, area_ratio)
    prarp
}

plot_prarp <- function(decon, truepar, prarp, peak_ratio, area_ratio) {

    # Check which peaks are found correctly and which were missed
    d <- decon
    x <- decon$x_values_ppm
    y <- decon$y_values
    dcnvpar <- data.frame(A = d$A_ppm, x0 = d$x_0_ppm, lambda = d$lambda_ppm)
    dcnvpar$y0 <- calc_y0(x, y, x0 = dcnvpar$x0)
    truepar$y0 <- calc_y0(x, y, x0 = truepar$x0)
    for (i in seq_along(truepar$x0)) {
        j <- which.min(abs(dcnvpar$x0 - truepar$x0[i]))
        d <- abs(dcnvpar$x0[j] - truepar$x0[i])
        truepar$closest[i] <- j
        truepar$dist[i] <- d
        truepar$found[i] <- d < 0.001
    }
    for (i in seq_along(dcnvpar$x0)) {
        j <- which.min(abs(truepar$x0 - dcnvpar$x0[i]))
        d <- abs(truepar$x0[j] - dcnvpar$x0[i])
        dcnvpar$closest[i] <- j
        dcnvpar$dist[i] <- d
        dcnvpar$correct[i] <- d < 0.001
    }

    # Plot the deconvolution results
    plot_spectrum(decon, foc_rgn = c(0.25, 0.75), foc_only = TRUE)
    points(
        x = dcnvpar$x0,
        y = dcnvpar$y0,
        pch = 21,
        cex = 3,
        bg = "transparent",
        lwd = 2,
        col = ifelse(dcnvpar$correct, "green", "red")
    )
    points(
        x = truepar$x0[!truepar$found],
        y = truepar$y0[!truepar$found],
        pch = 21,
        cex = 3,
        bg = "transparent",
        lwd = 2,
        col = "orange"
    )
    legend(
        x = "right",
        legend = c("Found", "Missed", "False Finding"),
        pch = 21,
        col = c("green", "orange", "red"),
        pt.bg = "transparent"
    )
    format2 <- function(x) format(x, digits = 2, scientific = FALSE)
    legend(
        x = "topleft",
        legend = c(
            paste("MSE Normed:", format2(decon$mse_normed)),
            paste("MSE Normed Raw:", format2(decon$mse_normed_raw)),
            paste("Peak Ratio (PR):", format2(peak_ratio)),
            paste("Area Ratio (AR):", format2(area_ratio)),
            paste("PRARP:", format2(prarp))
        ),
        bty = "n"
    )
    named(prarp, peak_ratio, area_ratio)
}

calc_y0 <- function(x, y, x0) {
    i0 <- convert_pos(x0, x, 1:length(x))
    i0_floor <- floor(i0)
    i0_ceil <- ceiling(i0)
    i0_frac <- i0 - i0_floor
    y_floor <- y[i0_floor]
    y_ceil <- y[i0_ceil]
    y0 <- y_floor + (y_ceil - y_floor) * i0_frac
}
