# evalwith #####

#' @description Run expression with predefined global state
#' @param expr Expression to be evaluated.
#' @param testdir ID of the test directory. E.g. `"xyz/2"`. Will be created and populated with `inputs`. To clear, use `clear(testdir("xyz/2"))`.
#' @param answers Answers to be returned by readline().
#' @param output Path to the file where output stream should be redirected to. Use `"captured"` to capture the output.
#' @param message Path to the file where message stream be redirected to. Use `"captured"` to capture the messages.
#' @param plot Path to the pdf file where plots should be saved to.
#' @param datadir_temp State of the mocked temporary data directory. See details section.
#' @param datadir_persistent State of the mocked persistent data directory. See details section.
#' @param inputs Paths to be copied to the test directory before evaluting `expr`.
#' @param opts Named list of options to be set. See [options()].
#' @param pars Named list of parameters to be set. See [par()].
#' @param cache Logical indicating whether to cache the result of the expression.
#' @param overwrite Logical indicating whether to overwrite the cache file if it already exists.
#' @details The `datadir_temp` and `datadir_persistent` arguments accept values "missing", "filled" and "empty". Setting a value unequal NULL causes the functions [datadir_temp()] and/or [datadir_persistent()] to be replaced with mock functions pointing to fake directories. Functions depending on these functions will then use the fake directories instead of the real ones. When set to "missing" the returned mock directory does not exist. When set to "empty" it exists and is guaranteed to be empty. When set to "filled", it is populated with example datasets.
#' Attention: the mocked functions, i.e. [datadir_temp()] and [datadir_persistent()] cannot be used directly inside `expr` when called via `devtools::test()`. I'm not sure why, but it seems as if devtools and/or testthat have their own copies of the functions which are used when the expression is evaluated.
#' @noRd
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
    hash <- NULL
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
            src_dir <- download_example_datasets()
            bruker_dir <- file.path(src_dir, "bruker")
            jcampdx_dir <- file.path(src_dir, "jcampdx")
            inputpaths <- gsub("bruker", bruker_dir, inputs, fixed = TRUE)
            inputpaths <- gsub("jcampdx", jcampdx_dir, inputpaths, fixed = TRUE)
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
    if (!is.null(plot)) {
        pdf(plot)
        on.exit(dev.off(), add = TRUE)
    }
    if (!is.null(opts)) {
        oldopts <- options(opts)
        on.exit(options(oldopts), add = TRUE)
    }
    if (!is.null(pars)) {
        oldpars <- par(pars)
        on.exit(par(oldpars), add = TRUE)
    }
    withCallingHandlers(
        testthat::with_mocked_bindings(
            tryCatch(
                {
                    runtime <- system.time(rv <- expr)[["elapsed"]]
                },
                error = function(e) {
                    sink(NULL, type = "message")
                    stop(e)
                }
            ),
            datadir_temp = get_datadir_mock(type = "temp", state = datadir_temp),
            datadir_persistent = get_datadir_mock(type = "persistent", state = datadir_persistent),
            readline = get_readline_mock(answers)
        ),
        warning = function(w) {
            message("Warning: ", conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    retobj <- invisible(list(
        rv = rv, runtime = runtime,
        output = outvec, message = msgvec, plot = plot,
        testdir = testdir, inputs = inputs, hash = hash
    ))
    if (isTRUE(cache) && (!file.exists(cachefile) || isTRUE(overwrite))) saveRDS(retobj, cachefile)
    retobj
}

testdir <- function(p = NULL) {
    normPath(paste(tempdir(), "tests", p, sep = "/"))
    # use paste instead of file.path, because it can deal with NULL
}

mockdir <- function() {
    normPath(file.path(tempdir(), "mocks"))
}

#' @description Create and return cache dir. If existing, the persistent cache dir is returned, else the temp cache dir. To force creation of the persistent cache dir, call once with `persistent=TRUE`.
#' @noRd
cachedir <- function(persistent = NULL) {
    tcd <- file.path(tempdir(), "cache") # temporary cache dir
    pcd <- file.path(tools::R_user_dir("metabodecon", "cache")) # persistent cache dir
    cd <- if (isTRUE(persistent) || (is.null(persistent) && dir.exists(pcd))) pcd else tcd
    ncd <- normalizePath(cd, "/", mustWork = FALSE)
    mkdirs(ncd)
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
        message(prompt, appendLF = FALSE)
        if (env$readline_called > length(texts)) {
            msg <- "readline called %s times, but only %s answers were provided."
            stop(sprintf(msg, env$readline_called, length(texts)))
        }
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
#' popgs()
#' @noRd
get_datadir_mock <- function(type = "temp", state = "default") {
    type <- match.arg(type, c("temp", "persistent"))
    state <- match.arg(state, c("default", "missing", "empty", "filled"))
    if (state == "default" && type == "persistent") {
        return(datadir_persistent)
    }
    if (state == "default" && type == "temp") {
        return(datadir_temp)
    }
    p <- normPath(file.path(mockdir(), "datadir", type, state))
    if (state %in% c("missing", "empty")) unlink(p, recursive = TRUE, force = TRUE)
    if (state == "empty") mkdirs(p)
    if (state == "filled") download_example_datasets(dst_dir = p)
    function() p
}

# interactive testing #####

#' @title Compare two vectors
#' @description Checks if `x` and `y` are identical, all.equal or different vectors. If different, differences are printed.
#' @return 0==identical, 1==all.equal, 2==different, 3==error
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
    status <- if (identical(x, y)) 0 else if (isTRUE(all.equal(x, y))) 1 else if (isvec(x) && isvec(y)) 2 else 3
    if (!silent) {
        msg <- c("identical", "all.equal", "different", "error")[status + 1]
        if (status == 3) msg <- paste0(msg, ": class(x)==", class(x), ", class(y)==", class(y))
        col <- c(GREEN, YELLOW, RED, RED)[status + 1]
        cat2(callstr, " ", col, msg, RESET, sep = "")
        if (status %in% 1:2 && status != xpct) {
            i <- which(x != y)
            a <- x[i]
            b <- y[i]
            z <- a - b # nolint: object_usage_linter
            rows <- c("class", "length", "head", "tail")
            cols <- c("x", "y", "i", "a", "b", "z")
            objs <- list(x = x, y = y, i = i, a = a, b = b, z = z)
            df <- as.data.frame(matrix(".", length(rows), length(cols), dimnames = list(rows, cols)))
            df["class", c("x", "y")] <- c(class(x), class(y))
            df["length", ] <- sapply(cols, function(col) length(get(col)))
            df["head", ] <- sapply(objs, function(obj) dput2(head(obj)))
            df["tail", ] <- sapply(objs, function(obj) dput2(tail(obj)))
            print(df)
        }
    }
    invisible(status)
}

#' @description Compare a spectrum deconvoluted with [generate_lorentz_curves_v12()] with a spectrum deconvoluted with [MetaboDecon1D()].
#' @param x Result of [generate_lorentz_curves_v12()].
#' @param y Result of [MetaboDecon1D()].
#' @examples \donttest{
#' new <- glc_urine1_yy_ni3_dbg(overwrite = TRUE)$rv$urine_1
#' old <- MetaboDecon1D_urine1_1010yy_ni3_dbg()$rv
#' compare_spectra(new, old)
#' }
#' @noRd
compare_spectra <- function(new = glc_urine1_yy_ni3_dbg()$rv$urine_1,
                            old = MD1D()$rv,
                            silent = FALSE) {
    comp <- function(new, old, xpct = 1) vcomp(new, old, xpct = xpct, silent = silent)
    r <- numeric()
    # o1 <- old$debuglist$args
    o2 <- old$debuglist$data_read
    o3 <- old$debuglist$ws_rm
    o4 <- old$debuglist$neg_rm
    o5 <- old$debuglist$smoothed
    o6 <- old$debuglist$peaks_sel
    o7 <- old$debuglist$peaks_wob_rm
    o8 <- old$debuglist$peak_scores_calc
    o9 <- old$debuglist$params_init
    o10 <- old$debuglist$params_approx
    # o11 <- old$debuglist$params_saved

    # [1-5] spectra <- get_spectra(data_path, file_format, expno, procno, ask, sf)
    # [6-7] spectra <- add_sfrs(spectra, sfr, ask, adjno)
    # [8-9] spectra <- add_wsrs(spectra, wshw, ask, adjno)
    # [10] spec <- rm_water_signal_v12(spec)
    # [11] spec <- rm_negative_signals_v12(spec)
    # [12] spec <- smooth_signals_v12(spec, reps = smopts[1], k = smopts[2])
    r[1] <- comp(new$y_raw, o2$spectrum_y_raw)
    r[2] <- comp(new$y_scaled, o2$spectrum_y)
    r[3] <- comp(new$n, o3$spectrum_length)
    r[4] <- comp(new$sdp, o3$spectrum_x)
    r[5] <- comp(new$ppm, o3$spectrum_x_ppm)
    r[6] <- comp(new$sfr$left_sdp, o3$signal_free_region_left)
    r[7] <- comp(new$sfr$right_sdp, o3$signal_free_region_right)
    r[8] <- comp(new$wsr$left_dp, o3$water_signal_left)
    r[9] <- comp(new$wsr$right_dp, o3$water_signal_right)
    r[10] <- comp(new$y_nows, o3$spectrum_y)
    r[11] <- comp(new$y_pos, o4$spectrum_y)
    r[12] <- comp(new$y_smooth, o5$spectrum_y)

    # spec <- find_peaks_v12(spec)
    r[13] <- comp(new$n, o6$spectrum_length)
    r[14] <- comp(new$d, c(NA, o6$second_derivative[2, ], NA)) # (1)
    r[16] <- comp(new$peak$center, as.integer(o6$peaks_index + 1)) # (1)
    r[17] <- comp(new$peak$center, as.integer(o6$peaks_x + 1)) # (1)
    r[18] <- comp(new$peak$right, as.integer(o6$left_position[1, ]) + 1) # (1)
    r[19] <- comp(new$peak$left, as.integer(o6$right_position[1, ]) + 1) # (1)
    # (1) MetaboDecon1D did not store NAs at the border, which is bad, because you need to shift every index by one when you switch from `second_derivative` to any other vector like `x_ppm` or `y_au`.

    # spec <- rm_peaks_with_low_scores_v12(spec, delta)
    border_is_na <- which(is.na(new$peak$left) | is.na(new$peak$right)) # the original MetaboDecon1D implementation throws away NAs, so for comparsion we need to do the same
    new_peak2 <- if (length(border_is_na) > 0) new$peak[-border_is_na, ] else new$peak
    r[20] <- comp(new_peak2$center, as.integer(o7$peaks_index + 1)) # (1) see above
    r[21] <- comp(new_peak2$right, as.integer(o7$left_position) + 1) # (1)
    r[22] <- comp(new_peak2$left, as.integer(o7$right_position) + 1) # (1)
    r[23] <- comp(new$sdp[new_peak2$center], o7$peaks_x)
    r[24] <- comp(mean(new_peak2$score[new_peak2$region %in% c("sfrl", "sfrr")]), o8$mean_score)
    r[25] <- comp(sd(new_peak2$score[new_peak2$region %in% c("sfrl", "sfrr")]), o8$sd_score)
    r[26] <- comp(new_peak2$score, o8$scores[1, ])
    r[27] <- comp(which(new_peak2$region == "sfrl"), o8$index_left)
    r[28] <- comp(which(new_peak2$region == "sfrr"), o8$index_right)
    r[29] <- comp(new$peak$center[new$peak$high], as.integer(o8$filtered_peaks + 1)) # (1)
    r[30] <- comp(new$peak$right[new$peak$high], o8$filtered_left_position + 1) # (1)
    r[31] <- comp(new$peak$left[new$peak$high], o8$filtered_right_position + 1) # (1)
    r[32] <- comp(new$peak$score[new$peak$high], o8$save_scores)

    # spec <- init_lorentz_curves_v12(spec)
    r[33] <- comp(new$sdp, o9$spectrum_x, xpct = 0)
    r[34] <- comp(new$y_smooth, o9$spectrum_y, xpct = 0)
    r[35] <- comp(new$peak$center[new$peak$high], as.integer(o9$filtered_peaks + 1), xpct = 0) # (1)
    r[36] <- comp(new$peak$right[new$peak$high], o9$filtered_left_position + 1, xpct = 0) # (1)
    r[37] <- comp(new$peak$left[new$peak$high], o9$filtered_right_position + 1, xpct = 0) # (1)
    r[38] <- comp(new$lc$A, o9$A, xpct = 1)
    r[39] <- comp(new$lc$lambda, o9$lambda, xpct = 1)
    r[40] <- comp(new$lc$w, o9$w, xpct = 0)

    # spec <- refine_lorentz_curves_v12(spec, nfit)
    r[42] <- comp(new$lcr$w_new, o10$w_new, xpct = 0)
    r[43] <- comp(new$lcr$lambda_new, o10$lambda_new, xpct = 1)
    r[44] <- comp(new$lcr$A_new, o10$A_new, xpct = 1)
    r[45] <- comp(new$lcr$integrals[1, ], o10$integrals[1, ], xpct = 1)

    # spec <- add_return_list_v12(spec, n, nam, debug)
    r[46] <- comp(new$ret$number_of_files, as.integer(old$number_of_files)) # int 1
    r[47] <- comp(new$ret$filename, old$filename) # chr[1] "urine_1"
    r[48] <- comp(new$ret$x_values, old$x_values) # num[131072] 131.071 131.070 131.069 ...
    r[49] <- comp(new$ret$x_values_ppm, old$x_values_ppm) # num[131072] 14.80254 14.80239 14.80223 ...
    r[50] <- comp(new$ret$y_values, old$y_values) # num[131072] 0.000831 0.000783 0.000743 ...
    r[51] <- comp(new$ret$spectrum_superposition, old$spectrum_superposition[1, ], 1) # num [1:131072] 3.29e-05 3.29e-05 3.29e-05 3.29e-05 3.29e-05 ...
    r[52] <- comp(new$ret$mse_normed, old$mse_normed) # num 4.1e-11
    r[53] <- comp(new$ret$index_peak_triplets_middle, as.integer(old$index_peak_triplets_middle)) # int[1227] 36158 37148 37418 37434 38942 38959 39030 39046 39091 39109 ...
    r[54] <- comp(new$ret$index_peak_triplets_left, old$index_peak_triplets_left) # num[1227] 36160 37159 37422 37437 38948 ...
    r[55] <- comp(new$ret$index_peak_triplets_right, old$index_peak_triplets_right) # num[1227] 36155 37139 37414 37431 38937 ...
    r[59] <- comp(new$ret$integrals, old$integrals[1, ], 1) # num[1227] 0.000502 0.026498 0.000396 0.000379 0.00732 ...
    if (!is.null(old$signal_free_region)) r[60] <- comp(new$ret$signal_free_region, old$signal_free_region) # num[2] 109.1 21.9
    if (!is.null(old$range_water_signal_ppm)) r[61] <- comp(new$ret$range_water_signal_ppm, old$range_water_signal_ppm) # num 0.153
    r[62] <- comp(new$ret$A, old$A, 1) # num[1227] -0.00016 -0.008437 -0.000126 -0.000121 -0.00233 ...
    r[63] <- comp(new$ret$lambda, old$lambda, 1) # num[1227] -0.00775 -0.02188 -0.00672 -0.00566 -0.01252 ...
    r[64] <- comp(new$ret$x_0, old$x_0, 0) # num[1227] 94.9 93.9 93.7 93.6 92.1 ...

    r[56] <- comp(new$ret$peak_triplets_middle, old$peak_triplets_middle) # NULL
    r[57] <- comp(new$ret$peak_triplets_left, old$peak_triplets_left) # NULL
    r[58] <- comp(new$ret$peak_triplets_right, old$peak_triplets_right) # NULL

    r[is.na(r)] <- 4
    msg <- "Identical: %s, Equal: %s, Different: %s, Error: %s, Empty: %s"
    if (!silent) {
        catf(msg, sum(r == 0), sum(r == 1), sum(r == 2), sum(r == 3), sum(r == 4))
    }

    invisible(r)
}

#' @title Run tests with the option to skip slow tests
#' @description This function runs the tests in the current R package. If `all` is TRUE, it well set environment variable `RUN_SLOW_TESTS` to "TRUE" so that all tests are run. If `all` is FALSE, it will set `RUN_SLOW_TESTS` to "FALSE" so that slow tests are skipped.
#' @param all Logical. If TRUE, all tests are run. If FALSE, slow tests are skipped.
#' @return The result of devtools::test()
#' @noRd
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

#' @description Compare a spectrum deconvoluted with [generate_lorentz_curves_v12()] with a spectrum deconvoluted with [MetaboDecon1D()].
#' @param x Result of [generate_lorentz_curves_v12()].
#' @param y Result of [MetaboDecon1D()].
#' @examples \donttest{
#' obj <- glc_urine1_yy_ni3_dbg()$rv$urine_1
#' du(obj)
#' }
#' @noRd
du <- function(obj, pname = "", level = 0, max_level = 1, max_len = 50) {
    sizeMB <- function(x) paste(round(object.size(x) / 1e6, 2), "MB")
    for (name in names(obj)) {
        x <- obj[[name]]
        indent <- collapse(rep("    ", level), "")
        cat2(indent, name, " ", sizeMB(x), sep = "")
        if (is.list(x) && length(x) < max_len && level < max_level) {
            du(x, name, level + 1, max_level)
        }
    }
    if (level == 0) {
        cat2("Total:", sizeMB(obj))
    }
}

pairwise_identical <- function(x) {
    lapply(seq_len(length(x) - 1), function(i) identical(x[[i]], x[[i + 1]]))
}
