
# evalwith #####

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
#' @param opts Named list of options to be set. See [options()].
#' @param pars Named list of parameters to be set. See [par()].
#' @param cache Logical indicating whether to cache the result of the expression.
#' @param overwrite Logical indicating whether to overwrite the cache file if it already exists.
#' @details The `datadir_temp` and `datadir_persistent` arguments accept values "missing", "filled" and "empty". Setting a value unequal NULL causes the functions [datadir_temp()] and/or [datadir_persistent()] to be replaced with mocks functions pointing to fake directories. Functions depending on these functions will then use the fake directories instead of the real ones. When set to "missing" the returned mock directory does not exist. When set to "empty" it exists and is guaranteed to be empty. When set to "filled", it is populated with example datasets.
#' Attention: the mocked functions, i.e. [datadir_temp()] and [datadir_persistent()] cannot be used directly inside the expression when called via `devtools::test()`. I'm not sure why, but it seems as if devtools and/or testthat have their own copies of the functions which are used when the expression is evaluated.
#' @noRd
evalwith <- function(expr,
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
        if (file.exists(cachefile) && isFALSE(overwrite)) return(readRDS(cachefile))
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
    # readline <- get_readline_mock(answers) # (1)
    # datadir_temp <- get_datadir_mock(type = "temp", state = datadir_temp) # (1)
    # datadir_persistent <- get_datadir_mock(type = "persistent", state = datadir_persistent) # (1)
    # (1) mocks must be defined before calling with_mocked_bindings to prevent problems due to R's non standard evaluation mechanism
    withCallingHandlers(
        testthat::with_mocked_bindings(
            tryCatch(
                {
                    runtime <- system.time(rv <- {expr})[["elapsed"]]
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

test_evalwith <- function() {
    before_wd <- getwd()
    x <- evalwith(
        testdir = "with/1",
        answers = c("y", "n"),
        output = "captured", message = "captured", plot = "plots.pdf",
        inputs = "jcampdx/urine/urine_1.dx",
        expr = {
            cat2("Helloworld!") # output is captured
            readline("Continue?") # readline is mocked
            readline("Continue?") # readline is mocked
            message("Roar") # messages are captured
            warning("Blub") # warnings are transformed to messages and captured as well
            test_wd <- getwd() # working dir is set to '{testdir()}/with/1'
            y <- 2
            z <- 3 # vars can be assigned
            list(y=y, z=z) # return value of expression is captured in rv
        }
    )
    after_wd <- getwd()
    test_that("evalwith works", {
        expect_true(file.exists(file.path(x$testdir, "plots.pdf")))
        expect_true(file.exists(file.path(x$testdir, "urine_1.dx")))
        expect_equal(x$rv, list(y=2, z=3))
        expect_equal(z, 3)
        expect_true(x$runtime <= 1)
        expect_equal(x$output, "Helloworld!")
        expect_equal(x$message, c("Continue?y", "Continue?n", "Roar", "Warning: Blub"))
        expect_equal(test_wd, file.path(testdir(), "with/1"))
        expect_equal(after_wd, before_wd)
    })
}

testdir <- function() {
    p <- file.path(tempdir(), "tests")
    normalizePath(p, "/", mustWork = FALSE)
}

mockdir <- function() {
    p <- file.path(tempdir(), "mocks")
    normalizePath(p, "/", mustWork = FALSE)
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
    if (is.null(texts)) return(readline)
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
    if (state == "default" && type == "persistent") return(datadir_persistent)
    if (state == "default" && type == "temp") return(datadir_temp)
    p <- normPath(file.path(mockdir(), "datadir", type, state))
    if (state %in% c("missing", "empty")) unlink(p, recursive = TRUE, force = TRUE)
    if (state == "empty") mkdirs(p)
    if (state == "filled") download_example_datasets(dst_dir = p)
    function() p
}

# func results #####

#' @description The `metabodecon` package contains private copies of all functions from the precursor package `MetaboDecon1D`. The copies are functionally identical [^1] to the original versions expect for a few extra lines of code, which cause the functions to return additional information about intermediate calculations if argument `debug` is TRUE.
#'
#' This function uses that mechanism to call the original functions [MetaboDecon1D()] and [deconvolution] with different input parameters and stores the returned debug results incl. debug info as RDS files. This way, new or updated functions from `metabodecon` can be tested against the original function result to ensure correctness and/or backwards compatibility.
#'
#' [^1]: In fact, the functions are not 100% identical. All code parts modifying global state, such as writing to disk or changing the working directory have been altered to meet the CRAN requirements. But all calculations and the return value are completely identical.
#' @noRd
store_func_results <- function(overwrite = FALSE) {

    dstdir <- mkdirs(datadir("test_expects", warn = FALSE))
    xpcts <- list()
    store_as_rds2 <- function(name, expr) store_as_rds(name, dstdir, overwrite, expr)
    cat2("Storing test expects in:", dstdir)

    name <- "MetaboDecon1D__240319__urine1_1010yy_ni1"
    xpcts[[name]] <- store_as_rds2(name, {
        evalwith(
            testdir = name, inputs = "bruker/urine/urine_1", answers = c(10, 10, "y", "y"),
            output = "captured", message = "captured", plot = "plots.pdf",
            expr = MetaboDecon1D(filepath = ".", filename = "urine_1", file_format = "bruker", number_iterations = 1)
        )
    })

    name <- "MetaboDecon1D__240319__urine1_1010yy_ni1"
    xpcts[[name]] <- store_as_rds2(name, {
        evalwith(
            testdir = name, inputs = "bruker/urine/urine_1", answers = c(10, 10, "y", "y"),
            output = "captured", message = "captured", plot = "plots.pdf",
            expr = MetaboDecon1D(filepath = ".", filename = "urine_1", file_format = "bruker", number_iterations = 1)
        )
    })

    name <- "deconvolution_urine1_spF_ni10_cf1_nf2"
    xpcts[[name]] <- store_as_rds2(name, {
        evalwith(
            testdir = "deconvolution_v11/1", inputs = "bruker/urine/urine_1", answers = c("y", "y"),
            output = "captured", message = "captured", plot = "plots.pdf",
            expr = {
                set.seed(1234)
                deconvolution_v11(filepath = "urine_1/10", name = NULL, file_format = "bruker", same_parameter = FALSE, processing_value = 10, number_iterations = 10, range_water_signal_ppm = 0.1527692, signal_free_region = c(11.44494, -1.8828), smoothing_param = c(2, 5), delta = 6.4, scale_factor = c(1000, 1000000), current_filenumber = 1, number_of_files = 2, debug = TRUE)
            }
        )
    })

    name <- "deconvolution_urine1_spT_ni1_cf1_nf2"
    xpcts[[name]] <- store_as_rds2(name, {
        evalwith(
            testdir = "deconvolution_v11/2", inputs = c(urine = "bruker/urine/urine_1"), answers = c("y", "y"),
            output = "captured", message = "captured", plot = "plots.pdf",
            expr = {
                set.seed(1234)
                deconvolution_v11(filepath = "urine_1/10", name = NULL, file_format = "bruker", same_parameter = TRUE, processing_value = 10, number_iterations = 1, range_water_signal_ppm = 0.1527692, signal_free_region = c(11.44494, -1.8828), smoothing_param = c(2, 5), delta = 6.4, scale_factor = c(1000, 1000000), current_filenumber = 1, number_of_files = 2, debug = TRUE)
            }
        )
    })

    name <- "deconvolution_urine1_spT_ni1_cf2_nf2"
    xpcts[[name]] <- store_as_rds2(name, {
        evalwith(
            testdir = "deconvolution_v11/3", inputs = c(urine = "bruker/urine/urine_1"), answers = NULL,
            output = "captured", message = "captured", plot = "plots.pdf",
            expr = {
                set.seed(1234)
                deconvolution_v11(filepath = "urine_1/10", name = NULL, file_format = "bruker", same_parameter = TRUE, processing_value = 10, number_iterations = 1, range_water_signal_ppm = 0.1527692, signal_free_region = c(109.09458303373, 21.8529143006947), smoothing_param = c(2, 5), delta = 6.4, scale_factor = c(1000, 1000000), current_filenumber = 2, number_of_files = 2, debug = TRUE)
            }
        )
    })

    invisible(xpcts)
}

list_func_results <- function() {
    dstdir <- datadir("test_expects", warn = FALSE)
    if (!dir.exists(dstdir)) dir.create(dstdir, recursive = TRUE)
    dir(dstdir, full.names = TRUE)
}

get_func_result <- function(rds = "deconvolution_urine1_spF_ni1_cf1_nf2.rds") {
    dstdir <- datadir("test_expects", warn = FALSE)
    if (!dir.exists(dstdir)) dir.create(dstdir, recursive = TRUE)
    rds <- file.path(dstdir, rds)
    if (!file.exists(rds)) {
        text <- "File '%s' does not exist. Valid names are:\n%s\nIf you believe, the file should exist, run `store_func_results()` to create it."
        valid <- paste(list_func_results(), collapse = "\n")
        msg <- sprintf(text, rds, valid)
        stop(msg)
    }
    obj <- readRDS(rds)
    obj
}

#' @title Store the result of an expression as an RDS file
#' @description Evaluate an expression. Store the result as RDS file and return it. Overwrite existing files if `overwrite` is TRUE. Else skip execution of `expr` and instead read and return the RDS.
#' @param name The name of the expression.
#' @param dstdir The directory where the RDS file will be stored.
#' @param overwrite Logical indicating whether to overwrite the RDS file if it already exists.
#' @param expr The expression to evaluate and store the result of.
#' @return The result of the expression. The result is returned invisibly.
#' @noRd
store_as_rds <- function(name, dstdir, overwrite, expr) {
    rds <-  file.path(dstdir, paste0(name, ".rds"))
    exists <- file.exists(rds)
    status <- if (!exists) "generating" else if (overwrite) "overwriting" else "reading"
    cat3(paste0("\033[34m", name, "\033[0m"), status)
    if (!exists || overwrite) {
        x <- force(expr)
        saveRDS(x, file = rds)
    } else {
        x <- readRDS(rds)
    }
    invisible(x)
}

# interactive testing #####

loaded_via_devtools <- function() {
    nchar(pkg_file(".gitignore")) > 0 # .gitignore is part of .Rbuildignore, i.e. after installation it is no longer present and pkg_file will return ""
}

vcomp <- function(x, y) {
    callstr <- paste(deparse(sys.call()), collapse="")
    callstr <- gsub("\\s+", " ", callstr)
    o <- capture.output({
        r <- tryCatch({
            x <- vcomp_internal(x, y)
            m <- switch(as.character(x),
                "0" = paste0(GREEN, "identical", RESET),
                "1" = paste0(ORANGE, "all.equal", RESET),
                "2" = paste0(RED, "different", RESET)
            )
            list(x = x, m = m)
        }, error = function(e) {
            x <- 3
            m <- paste0(RED, e$message, RESET)
            list(x = x, m = m)
        })
    })
    cat2(callstr, ": ", r$m, sep = "")
    if (length(o) > 0) cat2(collapse(o, "\n"))
    r$x
}

vcomp_internal <- function(x, y) {
    if(!is.vector(x) || is.list(x)) stop("x must be a vector, but is: ", class(x))
    if(!is.vector(y) || is.list(y)) stop("y must be a vector, but is: ", class(y))
    if (identical(x, y)) {
        return(0)
    } else {
        capture.output(all_equal <- isTRUE(all.equal(x, y)))
        i <- which(x != y)
        a <- x[i]
        b <- y[i]
        z <- a - b
        cat2("Objects:")
        cat2("  length(x):", length(x))
        cat2("  length(y):", length(y))
        cat2("  head(x):", head(x))
        cat2("  head(y):", head(y))
        cat2("  tail(x):", tail(x))
        cat2("  tail(y):", tail(y))
        cat2("  class(x):", class(x))
        cat2("  class(y):", class(y))
        cat2("  attributes(x):", attributes(x))
        cat2("  attributes(y):", attributes(y))
        cat2("Differences:")
        cat2("  i <- which(x != y)")
        cat2("  a <- x[i]")
        cat2("  b <- y[i]")
        cat2("  z <- a - b")
        cat2("  length(i):", length(i))
        cat2("  head(i)", head(i))
        cat2("  tail(i)", tail(i))
        cat2("  head(a):", head(a))
        cat2("  head(b):", head(b))
        cat2("  tail(a):", head(a))
        cat2("  tail(b):", head(b))
        cat2("  head(z):", head(z))
        cat2("  tail(z):", tail(z))
        return(if (all_equal) 1 else 2)
    }
}

# list_func_results()
# # "deconvolution_urine1_spF_ni1_cf1_nf2.rds"
# # "deconvolution_urine1_spT_ni1_cf1_nf2.rds"
# # "deconvolution_urine1_spT_ni1_cf2_nf2.rds"
# spec <- deconvolution_v11()
# check_spec(spec, compare_against = "deconvolution_urine1_spF_ni1_cf1_nf2.rds")
check_spec <- function(spec, compare_against = "deconvolution_urine1_spF_ni10_cf1_nf2.rds") {
    func_result <- get_func_result(rds = compare_against)
    ref <- func_result$rv$debuglist
    blue <- function(x) cat2("\033[34m", x, "\033[0m", sep = "")
    x <- logical()

    blue("read_spectrum(path, type, sf, expno, procno)")
    x[1] <- vcomp(ref$data_read$spectrum_y_raw, spec$y_raw)
    x[2] <- vcomp(ref$data_read$spectrum_y, spec$y_scaled)

    blue("determine_signal_free_region(spec, sfr, ask)")
    blue("determine_water_signal(spec, hwidth_ppm = wshw, bwc, ask)")
    x[3] <- vcomp(ref$ws_rm$spectrum_length, spec$n)
    x[4] <- vcomp(ref$ws_rm$spectrum_x, spec$sdp)
    x[5] <- vcomp(ref$ws_rm$spectrum_x_ppm, spec$ppm)
    x[6] <- vcomp(ref$ws_rm$signal_free_region_left, spec$sfr$left_sdp)
    x[7] <- vcomp(ref$ws_rm$signal_free_region_right, spec$sfr$right_sdp)
    x[8] <- vcomp(ref$ws_rm$water_signal_left, spec$wsr$left_dp)
    x[9] <- vcomp(ref$ws_rm$water_signal_right, spec$wsr$right_dp)
    x[10] <- vcomp(ref$ws_rm$spectrum_y, spec$y_nows)

    blue("remove_negative_signals(spec)")
    x[11] <- vcomp(ref$neg_rm$spectrum_y, spec$y_pos)

    blue("smooth_signals(spec, reps = smopts[1], k = smopts[2], bwc)")
    x[12] <- vcomp(ref$smoothed$spectrum_y, spec$y_smooth)

    blue("find_peaks_v12(spec)")
    x[13] <- vcomp(ref$peaks_sel$spectrum_length, spec$n)
    x[14] <- vcomp(ref$peaks_sel$second_derivative[1, ], spec$sdp[2:(spec$n - 1)])
    x[15] <- vcomp(ref$peaks_sel$second_derivative[2, ], spec$d[2:(spec$n - 1)])
    x[16] <- vcomp(ref$peaks_sel$peaks_index, as.integer(spec$peak$center - 1))
    x[17] <- vcomp(ref$peaks_sel$peaks_x, as.integer(spec$peak$center - 1))
    x[18] <- vcomp(ref$peaks_sel$left_position[1, ], spec$peak$right - 1)
    x[19] <- vcomp(ref$peaks_sel$right_position[1, ], spec$peak$left - 1)
    x[20] <- vcomp(ref$peaks_wob_rm$peaks_x, spec$sdp[spec$peak$center])
    x[21] <- vcomp(ref$peaks_wob_rm$peaks_index, as.integer(spec$peak$center - 1))
    x[22] <- vcomp(ref$peaks_wob_rm$left_position[1, ], spec$peak$right - 1)
    x[23] <- vcomp(ref$peaks_wob_rm$right_position[1, ], spec$peak$left - 1)

    blue("rm_peaks_with_low_scores(spec)")
    peaks <- spec$peak$center
    scores <- spec$peak$score
    region <- spec$peak$region
    x[24] <- vcomp(ref$peak_scores_calc$mean_score, mean(scores[region %in% c("sfrl", "sfrr")]))
    x[25] <- vcomp(ref$peak_scores_calc$sd_score, sd(scores[region %in% c("sfrl", "sfrr")]))
    x[26] <- vcomp(ref$peak_scores_calc$scores[1,], spec$peak$score)
    x[27] <- vcomp(ref$peak_scores_calc$index_left, which(spec$peak$region == "sfrl"))
    x[28] <- vcomp(ref$peak_scores_calc$index_right, which(spec$peak$region == "sfrr"))
    x[29] <- vcomp(ref$peak_scores_calc$filtered_peaks, as.integer(spec$peak$center[spec$peak$high] - 1))
    x[30] <- vcomp(ref$peak_scores_calc$filtered_left_position, spec$peak$right[spec$peak$high] - 1)
    x[31] <- vcomp(ref$peak_scores_calc$filtered_right_position, spec$peak$left[spec$peak$high] - 1)
    x[32] <- vcomp(ref$peak_scores_calc$save_scores, spec$peak$score[spec$peak$high])

    blue("init_lorentz_curves(spec)")
    x[33] <- vcomp(ref$params_init$spectrum_x, spec$sdp)
    x[34] <- vcomp(ref$params_init$spectrum_y, spec$y_smooth)
    x[35] <- vcomp(ref$params_init$filtered_peaks, as.integer(spec$peak$center[spec$peak$high] - 1))
    x[36] <- vcomp(ref$params_init$filtered_left_position, spec$peak$right[spec$peak$high] - 1)
    x[37] <- vcomp(ref$params_init$filtered_right_position, spec$peak$left[spec$peak$high] - 1)
    x[38] <- vcomp(ref$params_init$A, spec$lc$A)
    x[39] <- vcomp(ref$params_init$lambda, spec$lc$lambda)
    x[40] <- vcomp(ref$params_init$w, spec$lc$w)
    x[41] <- vcomp(ref$params_init$number_iterations, nfit)

    blue("refine_lorentz_curves(spec, nfit)")
    x[42] <- vcomp(ref$params_approx$w_new, spec$lcr$w_new)
    x[43] <- vcomp(ref$params_approx$lambda_new, spec$lcr$lambda_new)
    x[44] <- vcomp(ref$params_approx$A_new, spec$lcr$A_new)
    x[45] <- vcomp(ref$params_approx$integrals[1, ], spec$lcr$integrals[1, ])

    blue("create_return_list(spec)")
    x[46] <- vcomp(ref$params_saved$index_peak_triplets_middle, spec$ret$index_peak_triplets_middle + 1)
    x[47] <- vcomp(ref$params_saved$index_peak_triplets_left, spec$ret$index_peak_triplets_left + 1)
    x[48] <- vcomp(ref$params_saved$index_peak_triplets_right, spec$ret$index_peak_triplets_right + 1)
    # x[49] <- vcomp(ref$params_saved$peak_triplets_middle, spec$ret$peak_triplets_middle)
    # x[50] <- vcomp(ref$params_saved$peak_triplets_left, spec$ret$peak_triplets_left)
    # x[51] <- vcomp(ref$params_saved$peak_triplets_right, spec$ret$peak_triplets_right)
    # x[52] <- vcomp(ref$params_saved$noise_threshold, spec$ret$noise_threshold)
    # x[54] <- vcomp(ref$params_saved$spectrum_info, spec$ret$spectrum_info)
    # x[55] <- vcomp(ref$params_saved$spectrum_output, spec$ret$spectrum_output)
    # x[56] <- vcomp(ref$params_saved$name_info_txt, spec$ret$name_info_txt)
    # x[57] <- vcomp(ref$params_saved$name_output_txt, spec$ret$name_output_txt)

    catf("Identical: %s, Numequal: %s, Different: %s, Error: %s", sum(x == 0), sum(x == 1), sum(x == 2), sum(x == 3))

    invisible(x)
}

list_tests <- function(testthat = FALSE) {
    ns <- asNamespace("metabodecon")
    grep("^test_", ls(ns), value = TRUE)
    # x <- capture.output(ls.str(tests))
    # x <- sapply(strsplit(x, ":"), function(x) x[1])
    # x <- gsub("^([^ ].*)", paste0(MAGENTA, "\\1", RESET), x)
    # cat2(x, sep = "\n")
}

make_test_files <- function() {
    if (!loaded_via_devtools()) stop("Only run during development while loading the package from source")
    tests <- list_tests()
    test_files <- dir(pkg_file("tests/testthat"), full.names = TRUE)
    to_rm <- setdiff(test_files, tests)
    to_add <- setdiff(tests, test_files)
    cat2("Removing:\n", collapse(to_rm, "\n"), sep = "")
    unlink(to_rm)
    lines <- paste0("# Generated from pkg source. Run `make_test_files()` after edits.\n", tests, "()")
    paths <- file.path(pkg_file("tests/testthat"), paste0(tests, ".R"))
    for (i in seq_along(paths)) {
        cat2("Writing:", paths[i])
        writeLines(lines[i], paths[i])
    }
    invisible(lines)
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
        skip("Slow tests are disabled. Set environment variable RUN_SLOW_TESTS to TRUE to enable them.")
    }
}
