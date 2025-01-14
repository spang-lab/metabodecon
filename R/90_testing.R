# Evalwith (Public) #####

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
        local_dir(testpath)
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
        on.exit(close(outcon), add = TRUE, after = FALSE)
        on.exit(sink(NULL), add = TRUE, after = FALSE)
    }
    msgvec <- vector("character")
    if (!is.null(message)) {
        msgcon <- if (message == "captured") textConnection("msgvec", "wr", local = TRUE) else file(message, open = "wt")
        sink(msgcon, type = "message")
        on.exit(close(msgcon), add = TRUE, after = FALSE)
        on.exit(sink(NULL, type = "message"), add = TRUE, after = FALSE)
    }
    dev_cur <- dev.cur()
    force(plot) # forces eval, useful if plot ~= `svg("abc.svg")`
    if (identical(dev_cur, dev.cur()) && !is.null(plot)) {
        if (is.expression(plot)) eval(plot) # plot == `quote(svg("abc.svg", width = 10))`
        if (identical(plot, "captured")) plot <- tempfile(fileext = ".png")
        if (grepl("\\.pdf$", plot)) pdf(plot)
        else if (grepl("\\.svg$", plot)) svg(plot)
        else if (grepl("\\.png$", plot)) png(plot)
        else stop("plot must be an expression opening a device or a path ending in .pdf, .svg, or .png")
    }
    if (!is.null(opts)) local_options(opts)
    if (!is.null(pars)) local_par(pars)
    if (!identical(dev_cur, dev.cur())) {
        # This must be done after both arguments `plot` and `pars` have been
        # evaluated, as both can lead to a change in the graphical device.
        on.exit(while (!identical(dev_cur, dev.cur())) dev.off(), add = TRUE, after = FALSE)
    }
    withCallingHandlers(
        testthat::with_mocked_bindings(
            code = tryCatch(
                expr = {
                    runtime <- system.time(rv <- expr)[["elapsed"]]
                },
                error = function(e) {
                    sink(NULL, type = "message")
                    stop(e)
                }
            ),
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

# Evalwith Helpers (Private) #####

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
    if (state == "filled") download_example_datasets(dst_dir = p, silent = TRUE)
    function() p
}

loaded_via_devtools <- function() {
    pkg_dir <- dirname(system.file("DESCRIPTION", package = "metabodecon"))
    loaded_via_devtools <- dir.exists(file.path(pkg_dir, "inst"))
    return(loaded_via_devtools)
}

# Testthat Helpers (private) #####

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
    on.exit(Sys.setenv(RUN_SLOW_TESTS = RUN_SLOW_TESTS_OLD), add = TRUE, after = FALSE)
    devtools::test()
}

skip_if_slow_tests_disabled <- function() {
    if (!Sys.getenv("RUN_SLOW_TESTS") == "TRUE") {
        testthat::skip("Slow tests (Use `Sys.setenv(RUN_SLOW_TESTS=TRUE)` or `run_tests(all=TRUE)` to enable).")
    }
}

skip_if_not_in_globenv <- function() {
    if (!identical(environment(), .GlobalEnv)) {
        testthat::skip("Manual tests (To run, open file and execute the code manually).")
    }
}

#' @noRd
#' @title Check if the size of each file in a directory is within a certain
#' range
#' @description Check if the size of each file in a directory is within 90% to
#' 110% of the expected size.
#' If a file size is not within this range, a message is printed and an error is thrown.
#' @param testdir A character string specifying the directory to check.
#' @param size_exp A named numeric vector where the names are filenames and the
#' values are the expected file sizes.
#' @examples
#' testdir <- tmpdir("examples/expect_file_size", create = TRUE)
#' cat("Helloworld\n", file = file.path(testdir, "file1.txt"))
#' cat("Goodbye\n", file = file.path(testdir, "file2.txt"))
#' size_exp <- c(file1.txt = 12, file2.txt = 9)
#' expect_file_size(testdir, size_exp)
expect_file_size <- function(testdir, size_exp) {
    paths <- file.path(testdir, names(size_exp))
    size_obs <- file.info(paths)$size
    file_has_correct_size <- size_obs > size_exp * 0.9 & size_obs < size_exp * 1.1
    lapply(seq_along(size_exp), function(i) {
        if (!isTRUE(file_has_correct_size[i])) {
            message(sprintf("Size of %s is %s which is not between %s and %s", paths[i], size_obs[i], size_exp[i] * 0.9, size_exp[i] * 1.1))
        }
    })
    testthat::expect_true(all(file_has_correct_size))
}

#' @noRd
#'
#' @title Expect Structure
#'
#' @description
#' Tests if the structure of an object matches the expected string
#'
#' @param
#' obj The object to test
#'
#' @param expected_str
#' The expected structure of the object as a string. Can be obtained by calling
#' `dput(capture.output(str(obj)))`.
#'
#' @return
#' A logical value indicating whether the structure of the object matches the
#' expected string.
#'
#' @examples
#' expect_str(list(a = 1, b = 2), c("List of 2", " $ a: num 1", " $ b: num 2"))
#'
expect_str <- function(obj, expected_str, ...) {
    testthat::expect_identical(capture.output(str(obj, ...)), expected_str)
}

# Misc (Private) #####

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
#' calc_prarp(decon, truepar)
#' plot_prarp(decon, truepar)
#'
#' ## Good deconvolution (PRARP ~= 0.64)
#' decon <- generate_lorentz_curves_sim(sim[[1]], delta = 0)
#' truepar <- sim[[1]]$meta$simpar[c("A", "x0", "lambda")]
#' calc_prarp(decon, truepar)
#' plot_prarp(decon, truepar)
#'
calc_prarp <- function(x, truepar = NULL, ...) {

    obj <- as_decon2(x, ...)
    truepar <- truepar %||% obj$meta$simpar

    x0_true <- truepar$x0
    x0_found <- obj$lcpar$x0
    idx_closest_true_peak <- sapply(x0_found, function(x0) which.min(abs(x0_true - x0)))

    np_true <- length(truepar$x0)
    np_found <- length(idx_closest_true_peak)
    np_correct <- length(unique(idx_closest_true_peak))
    np_wrong <- np_found - np_correct
    peak_ratio   <- min(np_found, np_true) / max(np_found, np_true)
    peak_ratio_x <- np_correct / (np_true + np_wrong)

    area_spectrum <- sum(obj$si)
    area_residuals <- sum(abs(obj$sit$sup - obj$si))
    area_ratio <- 1 - (min(area_residuals, area_spectrum) / max(area_residuals, area_spectrum))

    prarp <- peak_ratio * area_ratio
    prarpx <- peak_ratio_x * area_ratio

    named(prarpx, prarp, peak_ratio_x, peak_ratio, np_true, np_found, np_correct, np_wrong, area_ratio, area_spectrum, area_residuals)
}

plot_prarp <- function(decon, truepar) {

    # Calculate PRARP score
    obj <- calc_prarp(decon, truepar)
    prarp <- obj$prarp
    peak_ratio <- obj$peak_ratio
    area_ratio <- obj$area_ratio

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
    plot_spectrum(decon, foc_rgn = c(0, 1), foc_only = TRUE)
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

# Compare (Private) #####

pairwise_identical <- function(x) {
    lapply(seq_len(length(x) - 1), function(i) identical(x[[i]], x[[i + 1]]))
}

#' @noRd
#' @title Compare two vectors
#' @description Checks if `x` and `y` are identical, all.equal or different
#' vectors. If differences are greated than expected, the differing elements are
#' printed.
#' @param x First vector.
#' @param y Second vector.
#' @param xpct Expected result. 0==identical, 1==all.equal, 2==different,
#' 3==error.
#' @param silent Logical indicating whether to print the results.
#' @return 0==identical, 1==all.equal, 2==different, 3==error
#' @details The function compares the vectors `x` and `y` and prints the
#' results. If the vectors are different, the differing elements are printed as
#' follows:
#'
#' ```R
#'              x       y  i  a                b                     z
#' class  integer numeric  .  .                .                     .
#' length      10      10  2  2                2                     2
#' head        1L       1 5L 5L                5 -8.88178419700125e-16
#' tail       10L      10 9L 9L 8.99999999999999  1.06581410364015e-14
#' ```
#'
#' Where
#'
#' - `i` == the indices of the differing elements
#' - `a` == the differing elements from `x`
#' - `b` == the differing elements from `y`
#' - `z` == the difference between `a` and `b`
#'
#' @examples
#' x <- 1:10
#' y1 <- 1:10
#' y2 <- c(1:4, 5.000000000000001, 6:8, 8.99999999999999, 10)
#' y3 <- c(1:4, 5.1, 6:8, 8.9, 10)
#' y4 <- list("a", "b", "c")
#' vcomp(x, y1)
#' vcomp(x, y2)
#' vcomp(x, y3)
#' vcomp(x, y4)
vcomp <- function(x, y, xpct = 0, silent = FALSE) {
    isvec <- function(x) is.vector(x) && !is.list(x)
    callstr <- paste(deparse(sys.call()), collapse = "")
    callstr <- gsub("\\s+", " ", callstr)
    status <- { # styler: off
        if (identical(x, y)) 0
        else if (isTRUE(all.equal(x, y))) 1
        else if (isvec(x) && isvec(y)) 2
        else 3
    } # styler: on
    if (!silent) {
        msg <- c("identical", "all.equal", "different", "error")[status + 1]
        if (status == 3) msg <- paste0(msg, ": class(x)==", class(x), ", class(y)==", class(y))
        col <- c(esc$green, esc$yellow, esc$red, esc$red)[status + 1]
        cat2(callstr, " ", col, msg, esc$reset, sep = "")
        if (status %in% 1:2 && status != xpct) {
            line <- "----------------------------------------"
            cat2(col, line, esc$reset, sep = "")
            i <- which(x != y)
            a <- x[i]
            b <- y[i]
            z <- a - b # nolint: object_usage_linter
            rows <- c("class", "length", "head", "tail")
            cols <- c("x", "y", "i", "a", "b", "z")
            objs <- list(x = x, y = y, i = i, a = a, b = b, z = z)
            df <- as.data.frame(matrix(".", length(rows), length(cols), dimnames = list(rows, cols)))
            df["class", c("x", "y")] <- c(class(x), class(y))
            df["length", c("x", "y", "i")] <- c(length(x), length(y), length(i))
            df["length", c("a", "b", "z")] <- c(".", ".", ".")
            df["head", ] <- sapply(objs, function(obj) dput2(head(obj, 1)))
            df["tail", ] <- sapply(objs, function(obj) dput2(tail(obj, 1)))
            print(df)
            cat2(col, line, esc$reset, sep = "")
        }
    }
    invisible(status)
}

#' @noRd
#' @description Compares a spectrum deconvoluted with
#' [generate_lorentz_curves_v12()] with a spectrum deconvoluted with
#' [MetaboDecon1D()].
#' @param x Result of [generate_lorentz_curves_v12()].
#' @param y Result of [MetaboDecon1D()].
#' @examples
#' sim_01 <- metabodecon_file("sim_01")[1]
#' new <- generate_lorentz_curves_sim(sim_01)
#' old <- MetaboDecon1D(sim_01)
#' r <- compare_spectra(new, old, silent = FALSE)
compare_spectra <- function(new, old, silent = FALSE) { # styler: off

    # Check args
    if (!is_idecon(new)) stop("new must be an object of type idecon")
    if (!is_decon0(old)) stop("old must be an object of type decon0")
    msg <- "old$debuglist is missing (set debug = TRUE in MetaboDecon1D())"
    if (!"debuglist" %in% names(old)) stop(msg)
    dbg <- old$debuglist

    # Define comparison functions
    ident <- update_defaults(vcomp, xpct = 0, silent = silent)
    equal <- update_defaults(vcomp, xpct = 1, silent = silent)

    # Define result vector
    r <-  logical()

    # Compare values after spectrum has been read and scaled
    o2 <- old$debuglist$data
    r[1] <-  equal(new$y_raw, dbg$data$spectrum_y_raw)
    r[1] <-  ident(as.numeric(new$y_raw), as.numeric(dbg$data$spectrum_y_raw))
    r[2] <-  ident(new$y_scaled, dbg$data$spectrum_y)

    # Compare values after water signal removal
    o3 <- old$debuglist$wsrm
    new_sfr <- enrich_sfr(sfr = new$args$sfr, x = new)
    new_wsr <- enrich_wshw(new$args$wsr, new)
    r[3] <-  ident(new$n, o3$spectrum_length)
    r[4] <-  ident(new$sdp, o3$spectrum_x)
    r[5] <-  equal(new$ppm, o3$spectrum_x_ppm)
    r[6] <-  ident(new_sfr$left_sdp, o3$signal_free_region_left)
    r[7] <-  ident(new_sfr$right_sdp, o3$signal_free_region_right)
    r[8] <-  ident(new$wsr$left_dp, o3$water_signal_left)
    r[9] <-  ident(new$wsr$right_dp, o3$water_signal_right)
    r[10] <- ident(new$y_nows, o3$spectrum_y)

    # Compare values after negative value removal
    o4 <- old$debuglist$nvrm
    r[11] <- ident(new$y_pos, o4$spectrum_y)

    # Compare values after  smoothing
    o5 <- old$debuglist$smooth
    r[12] <- ident(new$y_smooth, o5$spectrum_y)

    # Compare values after peak selection
    o6 <- old$debuglist$peaksel
    r[14] <- equal(new$d, c(NA, o6$second_derivative[2, ], NA))
    r[16] <- equal(new$peak$center, o6$peaks_x + 1) # (1)
    r[17] <- equal(new$peak$right, o6$left_position[1, ] + 1) # (1)
    r[18] <- equal(new$peak$left, o6$right_position[1, ] + 1) # (1)
    # (1) MetaboDecon1D did not store NAs at the border, which is bad, because
    # you need to shift every index by one when you switch from
    # `second_derivative` to any other vector like `x_ppm` or `y_au`.


    # Filter peaks
    o7 <- old$debuglist$peakfilter
    o8 <- old$debuglist$peakscore
    border_is_na <- which(is.na(new$peak$left) | is.na(new$peak$right)) # (2)
    new_peak2    <- if (length(border_is_na) > 0) new$peak[-border_is_na, ] else new$peak
    left_pos     <- if (is.matrix(o7$left_position)) o7$left_position[1, ] else o7$left_position
    right_pos    <- if (is.matrix(o7$right_position)) o7$right_position[1, ] else o7$right_position
    new_idx_left   <- which(new_peak2$region == "sfrl")
    new_idx_right  <- which(new_peak2$region == "sfrr")
    new_idx_sfr    <- c(new_idx_left, new_idx_right)
    new_mean_score <- mean(c(new_peak2$score[new_idx_sfr]))
    new_mean_sd    <- sd(new_peak2$score[new_idx_sfr])
    r[19] <- equal(new_peak2$center,               o7$peaks_index + 1) # (1)
    r[20] <- equal(new_peak2$right,                left_pos + 1)
    r[21] <- equal(new_peak2$left,                 right_pos + 1)
    r[22] <- equal(new$sdp[new_peak2$center],      o7$peaks_x)
    r[25] <- ident(new_peak2$score,                o8$scores[1, ])
    r[26] <- ident(new_idx_left,                   o8$index_left)
    r[27] <- ident(new_idx_right,                  o8$index_right)
    r[23] <- ident(new_mean_score,                 o8$mean_score)
    r[24] <- ident(new_mean_sd,                    o8$sd_score)
    r[28] <- equal(new$peak$center[new$peak$high], o8$filtered_peaks + 1) # (1)
    r[29] <- ident(new$peak$right[new$peak$high],  o8$filtered_left_position + 1)
    r[30] <- ident(new$peak$left[new$peak$high],   o8$filtered_right_position + 1)
    r[31] <- ident(new$peak$score[new$peak$high],  o8$save_scores)

    # Lorent Curve init values
    o9    <- old$debuglist$parinit
    r[32] <- equal(new$lci$P$ic,   o9$filtered_peaks + 1) # (1)
    r[33] <- ident(new$lci$P$ir,   o9$filtered_left_position + 1)
    r[34] <- ident(new$lci$P$il,   o9$filtered_right_position + 1)
    r[35] <- equal(new$lci$A,      o9$A)
    r[36] <- equal(new$lci$lambda, o9$lambda)
    r[37] <- equal(new$lci$w,      o9$w)

    # Lorent Curve refined values
    o10 <- old$debuglist$parapprox
    r[38] <- ident(new$lcr$w,         o10$w_new         )
    r[39] <- equal(new$lcr$lambda,    o10$lambda_new    )
    r[40] <- equal(new$lcr$A,         o10$A_new         )
    r[41] <- equal(new$lcr$integrals, o10$integrals[1, ])

    # Return List
    old_sfr <- old$signal_free_region %||% new$ret$signal_free_region
    old_rws <- old$range_water_signal_ppm %||% new$ret$range_water_signal_ppm
    r[42] <- equal(new$ret$number_of_files,            old$number_of_files)
    r[43] <- ident(new$ret$filename,                   old$filename)
    r[44] <- ident(new$ret$x_values,                   old$x_values)
    r[45] <- ident(new$ret$x_values_ppm,               old$x_values_ppm)
    r[46] <- ident(new$ret$y_values,                   old$y_values)
    r[47] <- equal(new$ret$spectrum_superposition,     old$spectrum_superposition[1, ])
    r[48] <- ident(new$ret$mse_normed,                 old$mse_normed)
    r[49] <- equal(new$ret$index_peak_triplets_middle, old$index_peak_triplets_middle)
    r[50] <- equal(new$ret$index_peak_triplets_left,   old$index_peak_triplets_left)
    r[51] <- equal(new$ret$index_peak_triplets_right,  old$index_peak_triplets_right)
    r[52] <- ident(new$ret$peak_triplets_middle,       old$peak_triplets_middle)
    r[53] <- ident(new$ret$peak_triplets_left,         old$peak_triplets_left)
    r[54] <- ident(new$ret$peak_triplets_right,        old$peak_triplets_right)
    r[55] <- equal(new$ret$integrals,                  old$integrals[1, ])
    r[56] <- equal(new$ret$A,                          old$A)
    r[57] <- equal(new$ret$lambda,                     old$lambda)
    r[58] <- ident(new$ret$x_0,                        old$x_0)
    r[59] <- equal(new$ret$signal_free_region,         old_sfr)
    r[60] <- equal(new$ret$range_water_signal_ppm,     old_rws)

    # Print summary
    if (!silent) {
        msg <- "Identical: %s, Equal: %s, Different: %s, Error: %s, Empty: %s"
        logf(msg, sum(r == 0), sum(r == 1), sum(r == 2), sum(r == 3), sum(r == 4))
    }


    # (2) The original MetaboDecon1D implementation throws away NAs, so for
    #     comparsion we need to do the same

    # Return results
    r[is.na(r)] <- 4
    invisible(r)
}
# styler: on

#' @noRd
#' @description Helper of [compare_spectra()].
update_defaults <- function(func, ...) {
    kwargs <- list(...)
    defaults <- formals(func)
    for (name in names(kwargs)) {
        defaults[[name]] <- kwargs[[name]]
    }
    formals(func) <- defaults
    func
}

# Deconvolution #####

#' @noRd
#' @author Tobias Schmidt
MetaboDecon1D_silent <- function(# Passed on to [MetaboDecon1D()]
                                 filepath,
                                 filename = NA,
                                 file_format = "bruker",
                                 number_iterations = 10,
                                 range_water_signal_ppm = 0.1527692,
                                 signal_free_region = c(11.44494, -1.8828),
                                 smoothing_param = c(2, 5),
                                 delta = 6.4,
                                 scale_factor = c(1000, 1000000),
                                 debug = FALSE,
                                 store_results = NULL,
                                 # Passed on to [evalwith()]
                                 output = "captured",
                                 message = "captured",
                                 plot = "captured",
                                 # Passed on to [get_MetaboDecon1D_answers()]
                                 expno = 10,
                                 procno = 10) {
    answers <- get_MetaboDecon1D_answers(
        ns = if (is.na(filename)) length(list.dirs(filepath)) else 1,
        wshw = range_water_signal_ppm,
        sfr = signal_free_region,
        format = file_format,
        expno = expno,
        procno = procno
    )
    evalwith(
        answers = answers,
        output = output,
        message = output,
        plot = plot,
        decon0 <- MetaboDecon1D(
            filepath, filename, file_format, number_iterations,
            range_water_signal_ppm, signal_free_region, smoothing_param, delta,
            scale_factor, debug, store_results
        )
    )
    decon0
}

#' @noRd
#' @author Tobias Schmidt
MetaboDecon1D_silent_sim <- function(# Passed on to [MetaboDecon1D()]
                                     filepath,
                                     filename = NA,
                                     file_format = "bruker",
                                     number_iterations = 3,
                                     range_water_signal_ppm = 0,
                                     signal_free_region = c(3.55, 3.35),
                                     smoothing_param = c(2, 5),
                                     delta = 6.4,
                                     scale_factor = c(1000, 1000000),
                                     debug = FALSE,
                                     store_results = NULL,
                                     # Passed on to [evalwith()]
                                     output = "captured",
                                     message = "captured",
                                     plot = "captured",
                                     # Passed to [get_MetaboDecon1D_answers()]
                                     expno = 10,
                                     procno = 10) {
    MetaboDecon1D_silent(
        filepath, filename, file_format,
        number_iterations, range_water_signal_ppm, signal_free_region,
        smoothing_param, delta, scale_factor, debug, store_results,
        output, message, plot, expno, procno
    )
}

#' @noRd
#' @author Tobias Schmidt
#' @examples
#' sim <- metabodecon_file("bruker/sim_subset")
#' answers <- get_MetaboDecon1D_answers(ns = 1, wshw = 0, sfr = c(3.55, 3.35))
#' x <- evalwith(
#'     answers = answers,
#'     output = "captured",
#'     message = "captured",
#'     plot = "captured",
#'     expr = { decon_01 <- MetaboDecon1D(sim, "sim_01") }
#' )
#' str(decon_01, 1)
get_MetaboDecon1D_answers <- function(ns = 1, # Number of spectra
                                      wshw = 0.1527692,
                                      sfr = c(11.44494, -1.8828),
                                      format = "bruker",
                                      expno = 10,
                                      procno = 10) {
    answers <- c(
        ExpNo       = if (format == "bruker") expno else NULL,
        ProcNo      = if (format == "bruker") procno else NULL,
        SameParam   = if (ns > 1) "y" else NULL,
        AdjNo       = if (ns > 1) "1" else NULL,
        SFRok       = "n",
        Left        = max(sfr),
        Right       = min(sfr),
        SFRok       = "y",
        WSok        = "n",
        WSHW        = wshw,
        WSok        = "y",
        SaveResults = "n"
    )
}