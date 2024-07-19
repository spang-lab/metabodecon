# evalwith #####

#' @noRd
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
#' @return A list containing with following elements:
#' - `rv`: The return value of the expression.
#' - `runtime`: The "elapsed" runtime of the expression in seconds. Measured with [system.time()].
#' - `output`: The captured output.
#' - `message`: The captured messages.
#' - `plot`: The path to the saved plot.
#' - `testdir`: The path to the test directory.
#' - `inputs`: The paths to the copied input files.
#' @examples
#' x1 <- evalwith(output = "captured", cat("Helloworld\n"))
#' x2 <- evalwith(datadir_persistent = "missing", datadir())
#' x3 <- evalwith(testdir = "dummy", inputs = "bruker/urine/urine_1", dir())
#' x4 <- evalwith(Sys.sleep(0.02))
#' cat(sprintf("x1$output: '%s'", x1$output))
#' cat(sprintf("x2$rv: '%s'", x2$rv))
#' cat(sprintf("x3$rv: '%s'", x3$rv))
#' cat(sprintf("x4$runtime: '%s'", x4$runtime))
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
        testdir = testdir, inputs = inputs
    ))
    if (isTRUE(cache) && (!file.exists(cachefile) || isTRUE(overwrite))) saveRDS(retobj, cachefile)
    retobj
}

testdir <- function(p = NULL) {
    normPath(paste(tmpdir(), "tests", p, sep = "/"))
    # use paste instead of file.path, because it can deal with NULL
}

mockdir <- function() {
    normPath(file.path(tmpdir(), "mocks"))
}

#' @noRd
#' @description Create and return cache dir. If existing, the persistent cache dir is returned, else the temp cache dir. To force creation of the persistent cache dir, call once with `persistent=TRUE`.
cachedir <- function(persistent = NULL) {
    tcd <- file.path(tmpdir(), "cache") # temporary cache dir
    pcd <- file.path(tools::R_user_dir("metabodecon", "cache")) # persistent cache dir
    cd <- if (isTRUE(persistent) || (is.null(persistent) && dir.exists(pcd))) pcd else tcd
    ncd <- normalizePath(cd, "/", mustWork = FALSE)
    mkdirs(ncd)
}

#' @noRd
#' @title Creates a mock readline function for testing
#' @description Creates a mock readline function that returns the next element from a character vector each time it's called.
#' Used internally by [mock_readline()].
#' @param texts A character vector of responses to be returned by the readline function.
#' @return A function that mimics the readline function, returning the next element from `texts` each time it's called.
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
#' @description Returns a function that, when called, returns a path to a mock data directory. The type and state of the mock data directory can be specified. Used internally by [mock_datadir()].
#' @param type The type of data directory to mock. Can be "persistent" or "temp".
#' @param state The state of the data directory to mock. Can be "missing", "empty", or "filled".
#' @return A function that when called, returns a path to the mock data directory.
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
    p <- normPath(file.path(mockdir(), "datadir", type, state))
    if (state %in% c("missing", "empty")) unlink(p, recursive = TRUE, force = TRUE)
    if (state == "empty") mkdirs(p)
    if (state == "filled") download_example_datasets(dst_dir = p)
    function() p
}

# simulate #####

#' @noRd
#' @title Simulate multiple spectra based on real deconvolution results
#' @description Simulates multiple spectra based on the deconvolution results of real spectra. The simulated spectra are superpositions of Lorentzian functions from a small part of the original spectra.
#' @param pngdir Path to the directory where the PNG files of the simulated spectra should be saved.
#' @param pdfdir Path to the directory where the PDF files of the simulated spectra should be saved.
#' @param svgdir Path to the directory where the SVG files of the simulated spectra should be saved.
#' @param rdsdir Path to the directory where the RDS files of the simulated spectra should be saved.
#' @param brukerdir Path to the directory where the Bruker files of the simulated spectra should be saved.
#' @param verbose Logical indicating whether to print progress messages.
#' @param cache Logical indicating whether to cache the results.
#' @return A list containing the simulated spectra and the paths to the saved PNG, PDF, and RDS files.
#' @examples
#' simulate_spectra(cache = TRUE)
#' \dontrun{
#'      simulate_spectra(
#'          pngdir = "vignettes/Datasets/png",
#'          pdfdir = "vignettes/Datasets/pdf",
#'          svgdir = "vignettes/Datasets/svg",
#'          rdsdir = "inst/example_datasets/rds/sim",
#'          brukerdir = "inst/example_datasets/bruker/sim",
#'          verbose = TRUE
#'      )
#' }
simulate_spectra <- function(pngdir = NULL,
                             pdfdir = NULL,
                             svgdir = NULL,
                             rdsdir = NULL,
                             brukerdir = NULL,
                             verbose = TRUE) {
    deconvs <- glc("blood", debug = FALSE)$rv
    if (is.null(pngdir)) pngdir <- tmpdir("simulate_spectra/png", create = TRUE)
    if (is.null(pdfdir)) pdfdir <- tmpdir("simulate_spectra/pdf", create = TRUE)
    if (is.null(svgdir)) svgdir <- tmpdir("simulate_spectra/svg", create = TRUE)
    if (is.null(rdsdir)) rdsdir <- tmpdir("simulate_spectra/rds", create = TRUE)
    if (is.null(brukerdir)) brukerdir <- tmpdir("simulate_spectra/bruker", create = TRUE)
    logf("Saving pngs to %s", pngdir)
    logf("Saving pdfs to %s", pdfdir)
    logf("Saving svgs to %s", svgdir)
    logf("Saving rds to %s", rdsdir)
    logf("Saving bruker files to %s", brukerdir)
    logf("Simulating %s spectra, using following deconvolution results as base:", length(deconvs))
    sim <- lapply(deconvs, function(deconv) {
        logf(deconv$filename)
        old_name <- deconv$filename
        new_name <- gsub("blood", "sim", old_name)
        simulate_spectrum(
            deconv,
            show = FALSE,
            pngpath = file.path(pngdir, paste0(new_name, ".png")),
            pdfpath = file.path(pdfdir, paste0(new_name, ".pdf")),
            svgpath = file.path(svgdir, paste0(new_name, ".svg")),
            rdspath = file.path(rdsdir, paste0(new_name, ".rds")),
            brukerdir = file.path(brukerdir, new_name)
        )
    })
    invisible(named(sim, pngdir, pdfdir, rdsdir))
}

#' @noRd
#' @title Simulate a spectrum based on a real deconvolution result
#' @description Simulates a spectrum based on the deconvolution results of a real spectrum. The simulated spectrum is a superposition of Lorentzian functions from a small part of the original spectrum.
#' @param deconv A list containing the deconvolution result, as returned by [generate_lorentz_curves()].
#' @param show Logical indicating whether to display the simulated spectrum.
#' @param pngpath Path to save the simulated spectrum as a PNG file.
#' @param pdfpath Path to save the simulated spectrum as a PDF file.
#' @param svgpath Path to save the simulated spectrum as a SVG file.
#' @param rdspath Path to save the simulated spectrum as an RDS file.
#' @return A list containing the simulated spectrum and the parameters used for simulation.
#' @examples
#' simulate_spectrum()
#' \dontrun{
#'      deconv <- glc("blood_11")$rv$ret
#'      show = FALSE
#'      pngpath <- "vignettes/Datasets/png/sim_11.png"
#'      pdfpath <- "vignettes/Datasets/pdf/sim_11.pdf"
#'      svgpath <- "vignettes/Datasets/svg/sim_11.svg"
#'      rdspath <- "inst/example_datasets/rds/sim/sim_11.rds"
#'      brukerdir <- "inst/example_datasets/bruker/sim/sim_11"
#'      simulate_spectrum(deconv, show, pngpath, pdfpath, svgpath, rdspath, brukerdir)
#' }
simulate_spectrum <- function(deconv = glc("blood_01")$rv$blood_01$ret,
                              show = TRUE,
                              pngpath = NULL,
                              pdfpath = NULL,
                              svgpath = NULL,
                              rdspath = NULL,
                              brukerdir = NULL,
                              verbose = TRUE) {
    logv <- if (verbose) logf else function(...) NULL
    logv("Creating simulated spectrum based on deconvolution results of %s", deconv$filename)
    simspec <- create_sim_spec(deconv, verbose)
    if (show) plot_sim_spec(simspec)
    store_sim_spec(simspec, pngpath, pdfpath, svgpath, rdspath, brukerdir, verbose)
    logv("Finished creation of simulated spectrum")
    invisible(simspec)
}

#' @noRd
#' @title Create a simulated spectrum based on a real deconvolution result
#' @description Creates a simulated spectrum based on the deconvolution results of a real spectrum. The simulated spectrum is a superposition of Lorentzian functions from a small part of the original spectrum.
#' @param deconv A list containing the deconvolution result, as returned by [generate_lorentz_curves()].
#' @return A dataframe with following columns
#' `si_smooth`: Smoothed signal intensities.
#' `cs`: Chemical shifts in PPM.
#' `fq`: Frequencies in Hz.
#' `si_raw`: Raw signal intensities of base spectrum.
#' `si_smooth`: Smoothed simulated signal intensities of base spectrum.
#' `si_sim`: Simulated signal intensities of base spectrum (incl. noise).
#' `si_sim_raw`: Raw simulated signal intensities in raw format, i.e. `as.integer(si_sim * 1e6)`.
#' @examples
#' deconv = glc(dp = "blood_01", debug = FALSE)$rv$blood_01
#' simspec <- create_sim_spec(deconv)
#' str(simspec)
create_sim_spec <- function(deconv = glc(dp = "blood_01", debug = FALSE)$rv$blood_01,
                            verbose = TRUE) {

    logv <- if (verbose) logf else function(...) NULL
    cs <- deconv$x_values_ppm
    fq <- deconv$x_values_hz
    si_raw <- deconv$y_values_raw
    si_smooth <- deconv$y_values

    logv("Throwing away datapoints outside of 3.45 to 3.55 ppm range")
    ix <- which(cs >= 3.40 & cs <= 3.60)
    X <- data.frame(cs, fq, si_raw, si_smooth)[ix, ]

    logv("Throwing away lorentz curves outside of 3.45 to 3.55 ppm range")
    ix <- which(deconv$x_0_ppm >= 3.45 & deconv$x_0_ppm <= 3.55)
    P <- data.frame(
        A = -deconv$A_ppm[ix],
        w = +deconv$x_0_ppm[ix],
        l = -deconv$lambda_ppm[ix]
    )

    logv("Calculating simulated signal intensities (si_sim) as superposition of lorentz curves")
    rownames(X) <- NULL
    X$si_sim <- sapply(X$cs, function(csi) {
        sum(abs(P$A * (P$l / (P$l^2 + (csi - P$w)^2))))
    })
    colnames(P) <- c("A", "x_0", "lambda")

    logv("Adding noise to simulated data, with noise-SD being equal to SFR-SD.")
    # The true SFR covers approx. the first and last 20k datapoints. However, to
    # be on the safe side, only use the first and last 10k datapoints for the
    # calculation.
    sfr_sd <- sd(si_raw[c(1:10000, 121073:131072)] * 1e-6)
    X$si_sim <- X$si_sim + rnorm(n = length(X$si_sim), mean = 0, sd = sfr_sd)
    X$si_sim_raw <- as.integer(X$si_sim * 1e6) # This looses precision but will be the "ground truth" that we can also write to disk ==> we need to recalculate si_sim from si_sim_raw, so that si_sim and si_sim_raw are consistent.
    X$si_sim <- X$si_sim_raw / 1.e6
    named(X, P, filename = deconv$filename)
}

#' @noRd
#' @title Plots a spectrum simulated with [create_sim_spec()]
#' @param simspec A simulated spectrum as returned by [create_sim_spec()].
#' @examples
#' simspec <- create_sim_spec()
#' plot_sim_spec(simspec)
plot_sim_spec <- function(simspec = create_sim_spec()) {
    X <- simspec$X
    top_ticks <- seq(from = min(X$cs), to = max(X$cs), length.out = 5)
    top_labels <- round(seq(from = min(X$fq), to = max(X$fq), length.out = 5))
    line1_text <- sprintf("Name: %s", gsub("blood", "sim", simspec$filename))
    line2_text <- sprintf("Base: %s ", simspec$filename)
    line3_text <- sprintf("Range: 3.6-3.4 ppm")
    plot(
        x = X$cs, y = X$si_raw * 1e-6, type = "l",
        xlab = "Chemical Shift [PPM]", ylab = "Signal Intensity [AU]",
        xlim = c(3.6, 3.4), col = "black"
    )
    lines(x = X$cs, y = X$si_smooth, col = "blue")
    lines(x = X$cs, y = X$si_sim, col = "red")
    legend(
        x = "topleft", lty = 1,
        col = c("black", "blue", "red"),
        legend = c("Original Raw / 1e6", "Original Smoothed", "Simulated")
    )
    axis(3, at = top_ticks, labels = top_labels, cex.axis = 0.75)
    mtext(sprintf("Frequency [Hz]"), side = 3, line = 2)
    mtext(line1_text, side = 3, line = -1.1, col = "red", cex = 1, adj = 0.99)
    mtext(line2_text, side = 3, line = -2.1, col = "red", cex = 1, adj = 0.99)
    mtext(line3_text, side = 3, line = -3.1, col = "red", cex = 1, adj = 0.99)
}

#' @noRd
#' @title Stores a spectrum simulated with [create_sim_spec()]
#' @param simspec A simulated spectrum as returned by [create_sim_spec()].
#' @param pngpath Path to save the simulated spectrum as a PNG file.
#' @param pdfpath Path to save the simulated spectrum as a PDF file.
#' @param svgpath Path to save the simulated spectrum as a SVG file.
#' @param rdspath Path to save the simulated spectrum as an RDS file.
#' @param brukerdir Path to the directory where the Bruker files should be saved.
#' @examples
#' deconv = glc(dp = "blood_01")$rv$blood_01$ret
#' X <- create_sim_spec(deconv)
#' store_sim_spec(X)
store_sim_spec <- function(simspec = create_sim_spec(),
                           pngpath = NULL,
                           pdfpath = NULL,
                           svgpath = NULL,
                           rdspath = NULL,
                           brukerdir = NULL,
                           verbose = TRUE) {

    X <- simspec$X
    logv <- if (verbose) logf else function(...) NULL
    opts <- options(digits = 15); on.exit(options(opts), add = TRUE)
    pdfpath <- normPath(pdfpath)
    pngpath <- normPath(pngpath)
    rdspath <- normPath(rdspath)
    brukerdir <- normPath(brukerdir)

    logv("Storing simulated spectrum")
    if (is.null(pdfpath)) {
        logv("No PDF path provided. Not saving PDF.")
    } else {
        logv("Saving PDF to %s", pdfpath)
        if (!dir.exists(pdfdir <- dirname(normPath(pdfpath)))) mkdirs(pdfdir)
        pdf(file = pdfpath)
        tryCatch(plot_sim_spec(simspec), finally = dev.off())
    }

    if (is.null(pngpath)) {
        logv("No PNG path provided. Not saving PNG")
    } else {
        logv("Saving PNG to %s", pngpath)
        if (!dir.exists(pngdir <- dirname(normPath(pngpath)))) mkdirs(pngdir)
        png(filename = pngpath)
        tryCatch(plot_sim_spec(simspec), finally = dev.off())
    }

    if (is.null(svgpath)) {
        logv("No SVG path provided. Not saving SVG")
    } else {
        logv("Saving SVG to %s", svgpath)
        if (!dir.exists(svgdir <- dirname(normPath(svgpath)))) mkdirs(svgdir)
        svg(filename = svgpath)
        tryCatch(plot_sim_spec(simspec), finally = dev.off())
    }


    if (is.null(rdspath)) {
        logv("No RDS path provided. Not saving RDS.")
    } else {
        logv("Saving RDS to %s", rdspath)
        if (!dir.exists(rdsdir <- dirname(normPath(rdspath)))) mkdirs(rdsdir)
        saveRDS(simspec, file = rdspath)
    }

    if (is.null(brukerdir)) {
        logv("No bruker path provided. Not saving bruker.")
    } else {
        logv("Saving bruker files to %s", brukerdir)
        mkdirs(file.path(brukerdir, "10", "pdata", "10"))
        files <- list(
            acqus = file.path(brukerdir, "10", "acqus"),
            procs = file.path(brukerdir, "10", "pdata", "10", "procs"),
            one_r = file.path(brukerdir, "10", "pdata", "10", "1r")
        )
        strs <- list(
            procs = paste(
                sep = "\n",
                "##$BYTORDP=0", # Byte order (0 = Little endian)
                "##$NC_proc=0", # Exponent for intensity values (y_scaled = y_raw * 2^NC_proc)
                "##$DTYPP=0", # Data storage type (0=4-byte-integers, else=double)
                sprintf("##$SI=%d", length(X$si_sim)), # Number of data points
                sprintf("##$OFFSET=%.15f", max(X$cs)), # Offset in PPM, i.e. the maximum chemical shift
                ""
            ),
            acqus = paste(
                sep = "\n",
                sprintf("##$SW=%.15f", width(X$cs)), # Spectrum width in PPM
                sprintf("##$SFO1=%.15f", convert_pos(0, X$cs, X$fq) / 1e6), # Reference Frequency in MHz
                sprintf("##$SW_h=%.15f", width(X$fq)), # Spectrum width in Hz
                ""
            )
        )
        logf("Saving files to:\n%s", collapse(files, "\n"))
        cat(strs$acqus, file = files$acqus)
        cat(strs$procs, file = files$procs)
        writeBin(X$si_sim_raw, files$one_r)
        logf("Reading back files to validate them")
        cols <- c("cs", "fq", "si")
        spec_written <- X[, c("cs", "fq", "si_sim_raw")]
        colnames(spec_written) <- cols
        spec_read <- read_spectrum(brukerdir)[cols]
        if (!isTRUE(all.equal(spec_read, spec_written))) stop("The saved bruker files are not equal to the simulated spectrum")
        logf("Validation successful")
    }

    logv("Finished storing simulated spectrum")
    simspec
}

# testthat #####

skip_if_slow_tests_disabled <- function() {
    if (!Sys.getenv("RUN_SLOW_TESTS") == "TRUE") {
        testthat::skip("Slow tests are disabled. Use `Sys.setenv(RUN_SLOW_TESTS=TRUE)` or `run_tests(all=TRUE)` to enable.")
    }
}

#' @title Check if the size of each file in a directory is within a certain range
#' @description Check if the size of each file in a directory is within 90% to 110% of the expected size.
#' If a file size is not within this range, a message is printed and an error is thrown.
#' @param testdir A character string specifying the directory to check.
#' @param size_exp A named numeric vector where the names are filenames and the values are the expected file sizes.
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
#' @param expected_str The expected structure of the object as a string. Can be obtained by calling `dput(capture.output(str(obj)))`.
#' @return A logical value indicating whether the structure of the object matches the expected string
#' @examples
#' expect_str(list(a = 1, b = 2), c("List of 2", " $ a: num 1", " $ b: num 2"))
#' @noRd
expect_str <- function(obj, expected_str) {
    testthat::expect_identical(capture.output(str(obj)), expected_str)
}

#' @noRd
#' @title Run tests with the option to skip slow tests
#' @description Runs the tests in the current R package. If `all` is TRUE, it well set environment variable `RUN_SLOW_TESTS` to "TRUE" so that all tests are run. If `all` is FALSE, it will set `RUN_SLOW_TESTS` to "FALSE" so that slow tests are skipped.
#' @param all Logical. If TRUE, all tests are run. If FALSE, slow tests are skipped.
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

# compare #####

pairwise_identical <- function(x) {
    lapply(seq_len(length(x) - 1), function(i) identical(x[[i]], x[[i + 1]]))
}

#' @noRd
#' @title Compare two vectors
#' @description Checks if `x` and `y` are identical, all.equal or different vectors. If differences are greated than expected, the differing elements are printed.
#' @param x First vector.
#' @param y Second vector.
#' @param xpct Expected result. 0==identical, 1==all.equal, 2==different, 3==error.
#' @param silent Logical indicating whether to print the results.
#' @return 0==identical, 1==all.equal, 2==different, 3==error
#' @details The function compares the vectors `x` and `y` and prints the results. If the vectors are different, the differing elements are printed as follows:
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
    status <- if (identical(x, y)) 0 else if (isTRUE(all.equal(x, y))) 1 else if (isvec(x) && isvec(y)) 2 else 3
    if (!silent) {
        msg <- c("identical", "all.equal", "different", "error")[status + 1]
        if (status == 3) msg <- paste0(msg, ": class(x)==", class(x), ", class(y)==", class(y))
        col <- c(GREEN, YELLOW, RED, RED)[status + 1]
        cat2(callstr, " ", col, msg, RESET, sep = "")
        if (status %in% 1:2 && status != xpct) {
            line <- "----------------------------------------"
            cat2(col, line, RESET, sep = "")
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
            cat2(col, line, RESET, sep = "")
        }
    }
    invisible(status)
}

#' @noRd
#' @description Compares a spectrum deconvoluted with [generate_lorentz_curves_v12()] with a spectrum deconvoluted with [MetaboDecon1D()].
#' @param x Result of [generate_lorentz_curves_v12()].
#' @param y Result of [MetaboDecon1D()].
#' @examples \donttest{
#' new <- glc()$rv
#' old <- MetaboDecon1D_urine1_1010yy_ni3_dbg()$rv
#' compare_spectra(new, old)
#' }
compare_spectra <- function(new = glc()$rv,
                            old = md1d()$rv,
                            silent = FALSE) {
    # Define comparison function
    comp <- vcomp
    arglist <- formals(comp)
    arglist$silent <- silent
    formals(comp) <- arglist
    r <- numeric()

    # Extract old data
    o1 <- old$debuglist$args # nolint: object_usage_linter.
    o2 <- old$debuglist$data_read
    o3 <- old$debuglist$ws_rm
    o4 <- old$debuglist$neg_rm
    o5 <- old$debuglist$smoothed
    o6 <- old$debuglist$peaks_sel
    o7 <- old$debuglist$peaks_wob_rm
    o8 <- old$debuglist$peak_scores_calc
    o9 <- old$debuglist$params_init
    o10 <- old$debuglist$params_approx
    o11 <- old$debuglist$params_saved # nolint: object_usage_linter.

    # [1-5] spectra <- read_spectra(data_path, file_format, expno, procno, ask, sf, bwc = TRUE)
    # [6-7] spectra <- get_sfrs(spectra, sfr, ask, adjno)
    # [8-9] spectra <- get_wsrs(spectra, wshw, ask, adjno)
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

    # spec <- filter_peaks_v12(spec, delta)
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
        logf(msg, sum(r == 0), sum(r == 1), sum(r == 2), sum(r == 3), sum(r == 4))
    }

    invisible(r)
}

#' @noRd
#' @description Compares a spectrum deconvoluted with [generate_lorentz_curves_v12()] with a spectrum deconvoluted with [MetaboDecon1D()].
#' @param x Result of [generate_lorentz_curves_v12()].
#' @param y Result of [MetaboDecon1D()].
#' @examples
#' new <- glc(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE)$rv$urine_1
#' old <- md1d(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE)$rv
#' r <- compare_spectra_v13(new, old, silent = FALSE)
compare_spectra_v13 <- function(new = glc()$rv$urine_1,
                                old = md1d()$rv,
                                silent = FALSE) {
    # styler: off
    # Define comparison functions
    ident <- update_defaults(vcomp, xpct = 0, silent = silent)
    equal <- update_defaults(vcomp, xpct = 1, silent = silent)

    # Read and smooth spectrum
    o2 <- old$debuglist$data_read
    o3 <- old$debuglist$ws_rm
    o4 <- old$debuglist$neg_rm
    o5 <- old$debuglist$smoothed
    r    <-  ident(new$y_raw,         o2$spectrum_y_raw)
    r[2] <-  ident(new$y_scaled,      o2$spectrum_y)
    r[3] <-  ident(new$n,             o3$spectrum_length)
    r[4] <-  ident(new$sdp,           o3$spectrum_x)
    r[5] <-  ident(new$ppm,           o3$spectrum_x_ppm)
    r[6] <-  ident(new$sfr$left_sdp,  o3$signal_free_region_left)
    r[7] <-  ident(new$sfr$right_sdp, o3$signal_free_region_right)
    r[8] <-  ident(new$wsr$left_dp,   o3$water_signal_left)
    r[9] <-  ident(new$wsr$right_dp,  o3$water_signal_right)
    r[10] <- ident(new$y_nows,        o3$spectrum_y)
    r[11] <- ident(new$y_pos,         o4$spectrum_y)
    r[12] <- ident(new$y_smooth,      o5$spectrum_y)

    # Find peaks
    o6 <- old$debuglist$peaks_sel
    r[13] <- ident(new$n,           o6$spectrum_length)
    r[16] <- equal(new$peak$center, o6$peaks_x + 1) # (1)
    r[15] <- equal(new$peak$center, o6$peaks_index + 1)
    r[17] <- equal(new$peak$right,  o6$left_position[1, ] + 1)
    r[18] <- equal(new$peak$left,   o6$right_position[1, ] + 1)
    r[14] <- equal(new$d,           c(NA, o6$second_derivative[2, ], NA))

    # Filter peaks
    o7 <- old$debuglist$peaks_wob_rm
    o8 <- old$debuglist$peak_scores_calc
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
    o9    <- old$debuglist$params_init
    r[32] <- equal(new$lci$P$ic,   o9$filtered_peaks + 1) # (1)
    r[33] <- ident(new$lci$P$ir,   o9$filtered_left_position + 1)
    r[34] <- ident(new$lci$P$il,   o9$filtered_right_position + 1)
    r[35] <- equal(new$lci$A,      o9$A)
    r[36] <- equal(new$lci$lambda, o9$lambda)
    r[37] <- equal(new$lci$w,      o9$w)

    # Lorent Curve refined values
    o10 <- old$debuglist$params_approx
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

    # styler: on
    # (1) MetaboDecon1D did not store NAs at the border, which is bad, because you need to shift every index by one when you switch from `second_derivative` to any other vector like `x_ppm` or `y_au`.
    # (2) The original MetaboDecon1D implementation throws away NAs, so for comparsion we need to do the same

    # Return results
    r[is.na(r)] <- 4
    invisible(r)
}

#' @noRd
#' @description Helper of [compare_spectra_v13()].
update_defaults <- function(func, ...) {
    kwargs <- list(...)
    defaults <- formals(func)
    for (name in names(kwargs)) {
        defaults[[name]] <- kwargs[[name]]
    }
    formals(func) <- defaults
    func
}

# structures #####

str_urine_1_deconvoluted <- function(cf = 1, nf = 1, dx = FALSE, nested = TRUE, ni = 10) {
    fn <- if (dx) "urine_1.dx" else "urine_1"
    elemstr <- c(
        sprintf("number_of_files           : int %d", nf),
        sprintf('filename                  : chr "%s"', fn),
        sprintf("x_values                  : num [1:131072] 131 131 131 131 131 ..."),
        sprintf("x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ..."),
        sprintf("y_values                  : num [1:131072] 0.000831 0.000783 0.000743 0.000717 0.00065 ..."),
        sprintf("spectrum_superposition    : num [1, 1:131072] 3.51e-05 3.51e-05 3.51e-05 3.51e-05 3.52e-05 ..."),
        sprintf("mse_normed                : num 3.92e-11"),
        sprintf("index_peak_triplets_middle: num [1:1227] 36159 37149 37419 37435 38943 ..."),
        sprintf("index_peak_triplets_left  : num [1:1227] 36161 37160 37423 37438 38949 ..."),
        sprintf("index_peak_triplets_right : num [1:1227] 36156 37140 37415 37432 38938 ..."),
        sprintf("peak_triplets_middle      : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ..."),
        sprintf("peak_triplets_left        : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ..."),
        sprintf("peak_triplets_right       : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ..."),
        sprintf("integrals                 : num [1, 1:1227] 0.000501 0.026496 0.000402 0.000375 0.008274 ..."),
        # sprintf("signal_free_region        : num [1:2] 109.1 21.9"),
        # sprintf("range_water_signal_ppm    : num 0.153"),
        sprintf("A                         : num [1:1227] -0.00016 -0.008436 -0.000128 -0.000119 -0.002634 ..."),
        sprintf("lambda                    : num [1:1227] -0.00775 -0.02188 -0.00675 -0.00562 -0.01343 ..."),
        sprintf("x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ...")
    )
    ne <- length(elemstr)
    plainstr <- c(
        paste0("List of ", ne),
        paste0(" $ ", elemstr)
    )
    nestedstr <- c(
        paste0("List of ", nf),
        paste0(" $ ", fn, ":List of ", ne),
        paste0("  ..$ ", elemstr)
    )
    if (nested) nestedstr else plainstr
}

str_urine_2_deconvoluted <- function(cf = 2, nf = 2, dx = FALSE, nested = TRUE, ni = 10) {
    fn <- if (dx) "urine_2.dx" else "urine_2"
    elemstr <- c(
        sprintf("number_of_files           : int %d", nf),
        sprintf('filename                  : chr "%s"', fn),
        sprintf("x_values                  : num [1:131072] 131 131 131 131 131 ..."),
        sprintf("x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ..."),
        sprintf("y_values                  : num [1:131072] 0.00586 0.00578 0.00569 0.00557 0.00548 ..."),
        sprintf("spectrum_superposition    : num [1, 1:131072] 4.21e-05 4.21e-05 4.21e-05 4.21e-05 4.21e-05 ..."),
        sprintf("mse_normed                : num 2.86e-11"),
        sprintf("index_peak_triplets_middle: num [1:1393] 36290 37241 38346 38826 39025 ..."),
        sprintf("index_peak_triplets_left  : num [1:1393] 36297 37244 38349 38835 39028 ..."),
        sprintf("index_peak_triplets_right : num [1:1393] 36285 37234 38343 38823 39019 ..."),
        sprintf("peak_triplets_middle      : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ..."),
        sprintf("peak_triplets_left        : num [1:1393] 9.26 9.12 8.95 8.87 8.85 ..."),
        sprintf("peak_triplets_right       : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ..."),
        sprintf("integrals                 : num [1, 1:1393] 0.00683 0.00504 0.00322 0.0174 0.00274 ..."),
        # sprintf("signal_free_region        : num [1:2] 109.1 21.9"),
        # sprintf("range_water_signal_ppm    : num 0.153"),
        sprintf("A                         : num [1:1393] -0.002176 -0.001604 -0.001025 -0.005541 -0.000872 ..."),
        sprintf("lambda                    : num [1:1393] -0.0189 -0.0168 -0.014 -0.0291 -0.0146 ..."),
        sprintf("x_0                       : num [1:1393] 94.8 93.8 92.7 92.2 92 ...")
    )
    ne <- length(elemstr)
    plainstr <- c(
        paste0("List of ", ne),
        paste0(" $ ", elemstr)
    )
    nestedstr <- c(
        paste0("List of ", nf),
        paste0(" $ ", fn, ":List of ", ne),
        paste0("  ..$ ", elemstr)
    )
    if (nested) nestedstr else plainstr
}

str_urine_deconvoluted <- function(nf = 2, dx = FALSE, nested = TRUE, ni = 10) {
    u1 <- str_urine_1_deconvoluted(nf = 2, dx = dx, nested = TRUE, ni = ni)
    u2 <- str_urine_2_deconvoluted(nf = 2, dx = dx, nested = TRUE, ni = ni)
    c("List of 2", u1[2:length(u1)], u2[2:length(u2)])
}

# wrappers #####

# Wrappers for [generate_lorentz_curves()] and [MetaboDecon1D()] so we don't have to pass all arguments every time during development.

glc <- function(dp = "urine_1",
                ff = "bruker",
                nfit = 3,
                simple = TRUE,
                overwrite = FALSE,
                cout = TRUE,
                cplot = TRUE,
                cache = TRUE,
                debug = TRUE,
                nworkers = "auto",
                verbose = FALSE) {

    logv <- if (verbose) logf else function(...) NULL

    logv("Parsing GLC arguments")
    tid <- get_tid("glc", dp, ff, nfit, simple, debug)
    inputs <- if (dp %in% c("urine", "blood", "sim")) file.path(ff, dp) else file.path(ff, strsplit(dp, "_")[[1]][1], dp)
    answers <- get_glc_answers(dp, ff, simple, inputs)
    logv("Inputs: %s", collapse(inputs, "; "))
    logv("Answers: %s", collapse(answers, "; "))

    rds <- file.path(cachedir(), paste0(tid, ".rds"))
    if (file.exists(rds)) {
        logv("Reading %s", rds)
    } else {
        call <- substitute(generate_lorentz_curves(data_path = dp, file_format = ff, nfit = nfit, debug = debug, nworkers = nworkers))
        logv("Calling %s", collapse(format(call), ""))
    }

    invisible(evalwith(
        testdir = tid,
        inputs = inputs,
        answers = answers,
        cache = cache,
        overwrite = overwrite,
        plot = if (cplot) "plots.pdf" else NULL,
        output = if (cout) "captured" else NULL,
        message = if (cout) "captured" else NULL,
        expr = generate_lorentz_curves(data_path = dp, file_format = ff, nfit = nfit, debug = debug, nworkers = nworkers)
    ))
}

md1d <- function(dp = "urine_1",
                 ff = "bruker",
                 nfit = 3,
                 simple = TRUE,
                 overwrite = FALSE,
                 cout = TRUE,
                 cplot = TRUE,
                 cache = TRUE,
                 debug = TRUE,
                 verbose = FALSE) { # nolint: object_usage_linter.

    logv <- if (verbose) logf else function(...) NULL

    logv("Parsing GLC arguments")
    tid <- get_tid("md1d", dp, ff, nfit, simple, debug)
    if (dp %in% c("urine", "blood", "sim")) {
        inputs <- file.path(ff, dp) # e.g. 'bruker/urine'
        fp <- dp # e.g. 'urine', i.e. deconvolute all files in the 'urine' directory
        fn <- NA
    } else {
        pp <- strsplit(dp, "_")[[1]][1] # e.g. 'urine'
        inputs <- file.path(ff, pp, dp) # e.g. 'bruker/urine/urine_1', i.e. copy folder 'urine_1' from 'bruker/urine/urine_1' to testdir
        fp <- "." # i.e. deconvolute all files in the test directory, which is only 'urine_1'
        fn <- dp
    }
    answers <- get_md1d_answers(fn, ff, simple, inputs)
    logv("Inputs: %s", collapse(inputs, "; "))
    logv("Answers: %s", collapse(answers, "; "))

    rds <- file.path(cachedir(), paste0(tid, ".rds"))
    if (file.exists(rds)) {
        logv("Reading %s", rds)
    } else {
        call <- substitute(generate_lorentz_curves(data_path = dp, file_format = ff, nfit = nfit, debug = debug, nworkers = nworkers))
        logv("Calling %s", collapse(format(call), ""))
    }

    invisible(evalwith(
        testdir = tid,
        inputs = inputs,
        answers = answers,
        cache = cache,
        overwrite = overwrite,
        plot = if (cplot) "plots.pdf" else NULL,
        output = if (cout) "captured" else NULL,
        message = if (cout) "captured" else NULL,
        expr = MetaboDecon1D(filepath = fp, filename = fn, file_format = ff, number_iterations = nfit, debug = debug)
    ))
}

#' @noRd
#' @description Helper function for [md1d()].
get_md1d_answers <- function(fn, ff, simple, inputs) {
    if (simple) {
        if (any(grepl("sim", inputs))) {
           answers <- c(SFRok = "n", Left = "3.58", Right = "3.42", SFRok = "y", WSok = "n", WSHW = "0.0", WSok = "y", SaveResults = "n")
        } else {
            answers <- c(SFRok = "y", WSok = "y", SaveResults = "n")
        }
        if (is.na(fn)) answers <- c(SameParam = "y", AdjNo = "1", answers)
    } else {
        answers <- c(SFRok = "n", Left = "11", Right = "-1", SFRok = "y", WSok = "asdf", WSok = "n", WSHW = "0.13", WSok = "y", SaveResults = "n")
        if (is.na(fn)) answers <- c(SameParam = "n", answers, answers, SaveResults = "n")
    }
    if (ff == "bruker") {
        answers <- c(ExpNo = "10", ProcNo = "10", answers)
    }
    answers
}

#' @noRd
#' @description Helper function for [glc()].
get_glc_answers <- function(fn, ff, simple, inputs) {
    if (grepl("sim", inputs)) {
        answers <- c(SFRok = "n", Left = "3.58", Right = "3.42", SFRok = "y", WSok = "n", WSHW = "0.0", WSok = "y", SaveResults = "n")
    } else {
        answers <- c(SFRok = "y", WSok = "y", SaveResults = "n")
    }
    if (grepl("(urine|blood|sim)$", inputs)) {
        answers <- c(SameParam = "y", AdjNo = "1", answers)
    }
    answers
}

get_testmatrix <- function() {
    df <- expand.grid(dp = c("urine_1", "urine_2", "urine"), ff = c("bruker", "jcampdx"), nfit = c(1, 3, 10), simple = c(TRUE, FALSE), skip = TRUE, stringsAsFactors = FALSE)
    df$skip[df$ff == "bruker" | (df$ff == "jcampdx" & df$nfit == 3 & df$simple == TRUE)] <- FALSE
    x <- df$dp %in% c("urine_1", "urine_2") & df$ff == "jcampdx"
    df$dp[x] <- paste0(df$dp[x], ".dx")
    df
}

#' @description Generates a unique identifier for a test of `generate_lorentz_curves_v12` or `MetaboDecon1D`
#' @noRd
get_tid <- function(func, dp, ff, nfit, simple, debug) {
    paste(func, dp, ff, nfit, simple, debug, sep = "-")
}

#' @description Calls `func` for each row in `testmatrix` and caches the results
#' @param func Either "glc" or "md1d"
#' @param overwrite Logical indicating whether to overwrite cached results if they already exist
#' @noRd
cache_func_results <- function(func = "glc", overwrite = FALSE) {
    df <- get_testmatrix()
    cdir <- cachedir()
    tid <- get_tid(func, df$dp, df$ff, df$nfit, df$simple)
    rds <- file.path(cdir, paste0(tid, ".rds"))
    status <- ifelse(file.exists(rds), "cached", "todo")
    status[df$skip] <- "skip"
    callstr <- sprintf("%s(dp='%s', ff='%s', nfit=%d, simple=%s, overwrite=%s)", func, df$dp, df$ff, df$nfit, df$simple, overwrite)
    col <- ifelse(status == "cached", GREEN,  YELLOW)
    col[df$skip] <- BLUE
    df[, c("rds", "status", "col", "callstr")] <- list(rds, status, col, callstr)
    idxtodo <- which(status == "todo")
    idxdone <- which(status != "todo")
    process_row <- function(i) {
        row <- as.list(df[i, colnames(df)])
        cat2(row$callstr, " ", row$col, row$status, RESET, sep = "")
        x <- if (row$status == "todo") try(eval(parse(text = row$callstr))) else 0
        return(if (inherits(x, "try-error")) x else 0)
    }
    x <- lapply(idxdone, process_row) # dont spawn new processes for cached or skipped function calls
    y <- parallel::mclapply(idxtodo, process_row, mc.cores = ceiling(parallel::detectCores() / 2))
    return(unlist(c(x, y)))
}
