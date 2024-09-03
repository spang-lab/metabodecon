# Raw Spectrum Class #####

#' @export
#'
#' @title Metabodecon Classes
#'
#' @description
#' Objects of class `raw_spectrum` are used by 'metabodecon' to represent NMR spectra that have not been deconvoluted yet. The following functions can be used to create `raw_spectrum` objects:
#'
#' - `raw_spectrum()`: Create from scratch by specifying all elements directly. Usually not used directly, but through one of the following functions.
#' - `make_spectrum()`: Create by specifying signal intensities, chemical shifts and reference frequency.
#' - `read_spectrum()`: Create from local file or folder.
#' - `simulate_spectrum()`: Create from existing deconvolution result.
#' - `as_raw_spectrum()`: Create from existing compatible object.
#'
#' To test for `raw_spectrum` objects, `is_spectrum()` can be used.
#'
#' @details
#' A `raw_spectrum` object is a list with class attribute `"spectrum"` and the following elements in arbitrary order: `si`, `cs`, `fq`, `name`, `path`, `type`, `mfs`, `version`. Each element is described in 'Arguments' and must fulfill the following requirements:
#'
#' - `si`, `fq` and `cs` are numeric vectors of the same length.
#' - `name`, `type` and `path` are character vectors of length 1 or NULL.
#' - `version` is a character vector of length 1 of the form `x.y` where `x` and `y` are integers.
#' - `mfs` is a numeric vector of length 1 or NULL.
#'
#' Element access should always be done by name, e.g. using `spectrum_obj$si`, `spectrum_obj[["cs"]]` or `spectrum_obj[c("fq", "name")]`, because the order of elements is not guaranteed.
#'
#' @return
#' A `raw_spectrum` object as described in 'Details', except for `is_spectrum()`, which returns a locical.
#'
#' @param brukerdir Path to the directory where the Bruker files should be saved.
#' @param check_class If TRUE, class attribute "raw_spectrum" is required or the check will fail.
#' @param check_contents If TRUE, the contents of the provided object must comply to the expected structure of `raw_spectrum` objects, described in 'Details' or the check will fail.
#' @param cs Vector of "chemical shifts" in parts per pillion (ppm). Must be of the same length as `si` and `fq`.
#' @param cs_max The highest chemical shift value in ppm, usually shown as left end of the spectrum.
#' @param cs_width The width of the spectrum in ppm.
#' @param deconv A list containing the deconvolution result, as returned by [generate_lorentz_curves()].
#' @param expno The experiment number for the file. E.g. `"10"`. Only relevant if `file_format` equals `"bruker"`. For details about `expno` and `procno` see section [File Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure) in the metabodecon FAQ.
#' @param file_format The file format of the spectrum file. E.g. `"bruker"` or `"jcampdx"`.
#' @param force If `TRUE`, try to continue when encountering errors and print info messages instead. To hide these messages as well, set `silent = TRUE`.
#' @param fq Vector of "frequencies" in Hertz (Hz). Must be of the same length as `si` and `cs`.
#' @param fq_ref The frequency of the reference in Hz.
#' @param fq_width The width of the spectrum in Hz. If provided, the values calculated from `cs_max`, `cs_width` and `fq_ref` are compared against this value. If the values differs, an error is thrown. If set to `NULL`, this check will be skipped.
#' @param mfs Magnetic field strength in Tesla, e.g. `14.1`.
#' @param name The name of the spectrum, e.g. "Blood 1" or "Urine Mouse X23D".
#' @param noise_method Method to generate noise. Either `"RND"` or `"SFR"`.
#' @param path The path of the file/folder containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"` or `"example_datasets/bruker/urine/urine"`.
#' @param pdfpath Path to save the simulated spectrum as a PDF file.
#' @param pngpath Path to save the simulated spectrum as a PNG file.
#' @param procno The processing number for the file. E.g. `"10"`. Only relevant if `file_format` equals `"bruker"`.
#' @param raw If `FALSE`, scales the returned signal intensities based on information available in the spectrum metadata, in particular `NC_proc`. For details see `processing-reference.pdf`, available at <https://www.bruker.com/en.html> at section 'Services & Support > Documentation & Manuals > Magnetic Resonance > Acquisition & Processing > TopSpin Processing Commands and Parameters' (requires login).
#' @param rdspath Path to save the simulated spectrum as an RDS file.
#' @param show Logical indicating whether to display the simulated spectrum.
#' @param si Signal intensities in Arbitrary Units (au), ordered from highest to lowest corresponding chemical shift.
#' @param silent If `TRUE`, no output will be printed to the console.
#' @param svgpath Path to save the simulated spectrum as a SVG file.
#' @param type The type of experiment, e.g. "H1 CPMG" or "H1 NOESY".
#' @param verbose Logical indicating whether to print progress messages.
#' @param version The version of 'metabodecon' that was used to create the `spectrum` object.
#' @param x A [raw_spectrum] object.
#'
#' @examples
#'
#' ## Create from Scratch
#'
#' si <- c(1, 1, 3, 7, 8, 3, 8, 5, 2, 1)
#' cs_max <- 14.8
#' cs_width <- 20.0
#' fq_ref <- 600.25 * 1e6
#' fq_width <- 12005
#' spectrum1 <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width)
#' spectrum2 <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width = 12010, force = FALSE)
#' print(spectrum1)
#' print(spectrum2)
#' is_spectrum(spectrum1)
#'
#'
#' ## Read from Disk
#'
#' spectrum3 <- read_spectrum(metabodecon_file("urine_1"))
#' print(spectrum3)
#'
#' \dontrun{
#' # Reading files in JCAMP-DX format is very slow (about 30s on the development
#' # machine). So if possible, you should stick with the original Bruker format.
#' urine_1_dx <- system.file("example_datasets/jcampdx/urine/urine_1.dx", package = "metabodecon")
#' spectrum4 <- read_spectrum(urine_1_dx, file_format = "jcampdx")
#' stopifnot(all.equal(urine_1_spectrum, spectrum4))
#' }
#'
#'
#' ## Simulate
#'
#' spectrum5 <- simulate_spectrum()
#' print(spectrum5)
#'
#' \dontrun{
#'      deconv <- glc("blood_02")$rv$ret
#'      show = FALSE
#'      pngpath <- "vignettes/Datasets/png/sim_02.png"
#'      pdfpath <- "vignettes/Datasets/pdf/sim_02.pdf"
#'      svgpath <- "vignettes/Datasets/svg/sim_02.svg"
#'      rdspath <- "inst/example_datasets/rds/sim/sim_02.rds"
#'      brukerdir <- "inst/example_datasets/bruker/sim/sim_02"
#'      spectrum6 <- simulate_spectrum(deconv, show, pngpath, pdfpath, svgpath, rdspath, brukerdir)
#' }
NULL

#' @export
#' @rdname raw_spectrum
raw_spectrum <- function(si, cs, fq, name = NULL, path = NULL, type = NULL, mfs = NULL, version = packageVersion("metabodecon")) {
    x <- structure(named(si, cs, fq, name, path, type, mfs, version), class = "raw_spectrum")
    if (!is_spectrum(x, check_contents = TRUE)) stop("Invalid spectrum object")
    x
}

#' @export
#' @rdname raw_spectrum
make_spectrum <- function(si, cs_max, cs_width, fq_ref, fq_width = NULL, force = FALSE, silent = FALSE, name = NULL, path = NULL, type = NULL, mfs = NULL) {
    cs_min <- cs_max - cs_width # Lowest ppm value
    cs <- seq(cs_max, cs_max - cs_width, length.out = length(si)) # Chemical shift in parts per million
    fq_max <- fq_ref - (cs_min * 1e-6 * fq_ref)  # Highest frequency in Hz (corresponds to lowest ppm value)
    fq_min <- fq_ref - (cs_max * 1e-6 * fq_ref)  # Lowest frequency in Hz
    fq <- seq(fq_min, fq_max, length.out = length(si)) # Frequency in Hz
    fq_width_calc <- fq_max - fq_min
    if (!is.null(fq_width) && !isTRUE(all.equal(fq_width_calc, fq_width))) { # Check if calculated spectrum width in Hz matches the value provided by the user
        if (force) {
            stop(sprintf("Calculated spectrum width in Hz (%s) does not match the provided value (%s). Please read in the data manually or set `force = TRUE` to ignore this error. Please note that by doing so, all downstream calculations involving frequencies might be wrong, so be sure to double check the results.", round(fq_width_calc, 5), round(fq_width, 5)))
        } else if (!silent) {
            cat(sprintf("Calculated spectrum width in Hz (%s) does not match the provided value (%s). Continuing anyways, because `force` equals `TRUE`. Please note that all downstream calculations using frequencies might be wrong, so be sure to double check the results.", round(fq_width_calc, 5), round(fq_width, 5)))
        }
    }
    raw_spectrum(si, cs, fq, name, path, type, mfs, version = packageVersion("metabodecon"))
}

#' @export
#' @rdname raw_spectrum
read_spectrum <- function(path = metabodecon_file("bruker/sim/sim_01"),
                          file_format = "bruker",
                          expno = 10,
                          procno = 10,
                          raw = FALSE,
                          silent = TRUE,
                          force = FALSE) {
    file_format <- match.arg(file_format, c("bruker", "jcampdx"))
    if (file_format == "bruker") x <- read_bruker_spectrum(path, expno, procno, raw, silent, force)
    if (file_format == "jcampdx") x <- read_jcampdx_spectrum(path, silent, force)
    x
}

#' @export
#' @rdname raw_spectrum
as_raw_spectrum <- function(x, ...) {
    if (is.list(x)) {
        if (!"spectrum" %in% class(x)) class(x) <- c("spectrum", class(x))
        if (is_spectrum(x, check_contents = TRUE)) return(x)
    }
    stop("Cannot convert object to spectrum")
}

#' @export
#' @rdname raw_spectrum
is_spectrum <- function(x, check_class = TRUE, check_contents = FALSE) {
    if (!check_class && !check_contents) {
        stop("At least one of `check_class` and `check_contents` must be TRUE")
    }
    if (check_class) {
        if (!inherits(x, "raw_spectrum")) return(FALSE)
    }
    if (check_contents) {
        if (!is.list(x)) return(FALSE)
        mandatory <- c("si", "fq", "cs")
        if (!all(mandatory %in% names(x))) return(FALSE)
        lengths <- sapply(x[mandatory], length)
        if (length(unique(lengths)) != 1) return(FALSE)
        optional <- c("name", "type", "path", "mfs")
        ok <- sapply(x[optional], function(x) is.null(x) || (is.character(x) && length(x) == 1))
        if (!all(ok)) return(FALSE)
    }
    return(TRUE)
}

#' @export
#' @rdname raw_spectrum
simulate_spectrum <- function(deconv = glc("blood_01")$rv$ret,
                              show = TRUE,
                              pngpath = NULL,
                              pdfpath = NULL,
                              svgpath = NULL,
                              rdspath = NULL,
                              brukerdir = NULL,
                              verbose = TRUE,
                              noise_method = "SFR") {
    logv <- if (verbose) logf else function(...) NULL
    logv("Creating simulated spectrum from '%s' (noise method = '%s')", deconv$filename, noise_method)
    simspec <- create_sim_spec(deconv, verbose, noise_method)
    if (show) plot_sim_spec(simspec)
    store_sim_spec(simspec, pngpath, pdfpath, svgpath, rdspath, brukerdir, verbose)
    logv("Finished creation of simulated spectrum")
    invisible(simspec)
}

#' @export
#' @rdname raw_spectrum
print.spectrum <- function(x, ...) {
    catf("Spectrum with %d data points:\n", length(x$si))
    catf("- Signal Intensity Range: %.1f - %.1f\n", min(x$si), max(x$si))
    catf("- Frequency Range: %.1f - %.1f Hz\n", min(x$fq), max(x$fq))
    catf("- Chemical Shift Range: %.4f - %.4f ppm\n", max(x$cs), min(x$cs))
    catf("- Magnetic Field Strength: %s\n", paste(x$mfs, "T"))
    catf("- Name: %s\n", x$name %||% "NULL")
    catf("- Path: %s\n", x$path %||% "NULL")
    catf("- Type: %s\n", x$type %||% "NULL")
}

# Simulate Spectrum Helpers #####

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
#' deconv = glc(dp = "blood_01", debug = FALSE)$rv
#' simspec <- create_sim_spec(deconv)
#' str(simspec)
create_sim_spec <- function(deconv = glc(dp = "blood_01", debug = FALSE)$rv,
                            verbose = TRUE,
                            noise_method = "SFR") {

    noise_method <- match.arg(noise_method, c("RND", "SFR"))
    logv <- if (verbose) logf else function(...) NULL
    cs <- deconv$x_values_ppm
    fq <- deconv$x_values_hz
    si_raw <- deconv$y_values_raw
    si_smooth <- deconv$y_values

    logv("Throwing away datapoints outside of 3.45 to 3.55 ppm range")
    ix <- which(cs >= 3.40 & cs <= 3.60)
    X <- data.frame(cs, fq, si_raw, si_smooth)[ix, ]

    logv("Throwing away lorentz curves outside of 3.45 to 3.55 ppm range")
    if (deconv$filename == "blood_02") {
        ip <- which( deconv$x_0_ppm <= 3.54 & deconv$x_0_ppm >= 3.44) # Blood 02 is shifted approx. 0.01 ppm to the right, so we also need to shift the interval from which we pick our peaks so that we end up with signals from the same metabolites.
    } else {
        ip <- which(deconv$x_0_ppm <= 3.55 & deconv$x_0_ppm >= 3.45)
    }
    P <- data.frame(
        A      = -deconv$A_ppm[ip],
        x_0    = +deconv$x_0_ppm[ip],
        lambda = -deconv$lambda_ppm[ip]
    )

    logv("Calculating simulated signal intensities (si_sim) as superposition of lorentz curves")
    rownames(X) <- NULL
    X$si_sim <- sapply(X$cs, function(csi) {
        sum(abs(P$A * (P$lambda / (P$lambda^2 + (csi - P$x_0)^2))))
    })

    if (noise_method == "RND") { # ToSc, 2024-08-02: noise_method 'RND' turned out to be a bad idea, because the true noise is not normally distributed, but has long stretches of continuous increase or decrease that would be incredibly unlikely to occur with a normal distribution. This can be seen by running `analyze_noise_methods()`.
        logv("Adding noise to simulated data, with noise-SD being equal to SFR-SD.")
        sfr_sd <- sd(si_raw[c(1:10000, 121073:131072)] * 1e-6) # The true SFR covers approx. the first and last 20k datapoints. However, to be on the safe side, only use the first and last 10k datapoints for the calculation.
        noise <- rnorm(n = length(X$si_sim), mean = 0, sd = sfr_sd)
        X$si_sim <- X$si_sim + noise
    } else {
        logv("Adding noise to simulated data, with noise taken from SFR.")
        idx <- ix - min(ix) + 5000 # Use SI of datapoints 5000:6308 for noise
        noise <- si_raw[idx] * 1e-6
        X$si_sim <- X$si_sim + noise
    }
    X$si_sim_raw <- as.integer(X$si_sim * 1e6) # We need to convert to integers, because that's how the data gets written to disc. However, because this conversion loses precision, we need to recalculate `si_sim` from the less precise `si_sim_raw` as well, so that `si_sim` and `si_sim_raw` are consistent.
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

#' @noRd
#' @description Used during development of `simulate_spectra()` to find a realistic method for noise generation.
analyze_noise_methods <- function(ask = TRUE) {

    deconv = glc(dp = "blood_01", debug = FALSE)$rv
    si_raw <- deconv$y_values_raw
    sd_sfr <- sd(si_raw[c(1:10000, 121073:131072)] * 1e-6)
    siRND <- rnorm(n = 10000, mean = 0, sd = sd_sfr)
    siSFR <- si_raw[1:10000] * 1e-6

    logf("Visualizing raw SIs for noise methods RND and SFR")
    visualize <- function(siRND, siSFR, n = 300, start = 5000) {
        ymin <- min(c(min(siRND), min(siSFR)))
        ymax <- max(c(max(siRND), max(siSFR)))
        ylim <- c(ymin, ymax)
        opar <- par(mfrow = c(5, 1), mar = c(3, 4, 0, 1))
        on.exit(par(opar), add = TRUE)
        for (i in 1:5) {
            redT <- rgb(1, 0, 0, 0.1)
            bluT <- rgb(0, 0, 1, 0.1)
            idx <- ((i - 1) * n + 1):(i * n) + start
            ysiRND <- siRND[idx]
            ysiSFR <- siSFR[idx]
            plot(1:n, ysiRND, type = "n", ylim = ylim, ylab = "", xlab = "", xaxt = "n")
            axis(1, at = seq(1, n, by = 50), labels = idx[seq(1, n, by = 50)])
            points(1:n, ysiRND, col = "red", pch = 20)
            points(1:n, ysiSFR, col = "blue", pch = 20)
            lines(1:n, ysiRND, col = "red")
            lines(1:n, ysiSFR, col = "blue")
            lines(1:n, rep(0, n), col = "black", lty = 2)
            polygon(c(1:n, n:1), c(ysiRND, rep(0, n)), col = redT, border = NA)
            polygon(c(1:n, n:1), c(ysiSFR, rep(0, n)), col = bluT, border = NA)
            legend("topleft", NULL, c("RND", "SFR"), col = c("red", "blue"), lty = 1)
        }
    }
    visualize(siRND, siSFR)
    if (!get_yn_input("Continue?")) return()

    logf("Visualizing smoothed SIs for noise methods RND and SFR")
    siRND_sm <- smooth_signals(list(y_pos = siRND))$y_smooth
    siSFR_sm <- smooth_signals(list(y_pos = siSFR))$y_smooth
    visualize(siRND_sm, siSFR_sm)
    if (!get_yn_input("Continue?")) return()

    logf("Visualizing lengths of intervals of continuous increase and/or decrease")
    slRND <- count_stretches(siRND) # stretch lengths of siRND
    slSFR <- count_stretches(siSFR) # stretch lengths of SFR
    table(slRND)
    table(slSFR)
    opar <- par(mfrow = c(2, 1))
    on.exit(par(opar), add = TRUE)
    hist(slRND, breaks = 0:20 + 0.5, xlim = c(0, 20))
    hist(slSFR, breaks = 0:20 + 0.5, xlim = c(0, 20))

}

#' @title Count Stretches of Increases and Decreases
#' @description Counts the lengths of consecutive increases and decreases in a numeric vector.
#' @param vec A numeric vector.
#' @return A numeric vector containing the lengths of stretches of increases and decreases.
#' @examples
#' #
#' # Example Data (x)
#' # | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
#' # |-------------------------------|
#' # |   |   |   |###|   |   |   |###|
#' # |   |   |###|###|   |   |###|###|
#' # |   |###|###|###|###|   |###|###|
#' # |###|###|###|###|###|###|###|###|
#' # |-------------------------------|
#' # |   +   +   +   -   -   +   +   |
#' # |-------------------------------|
#' #
#' x <- c(1.0, 2.2, 3.0, 4.4, 2.0, 1.0, 3.0, 4.0)
#' count_stretches(x) # Returns c(3, 2, 2) because we have 3+, 2-, 2+
count_stretches <- function(x) {
  if (length(x) < 2) return(integer(0))
  ss <- numeric(length(x))
  s <- 1
  inc <- x[2] > x[1]
  for (i in 3:length(x)) {
    if ((inc && x[i] > x[i - 1]) || (!inc && x[i] < x[i - 1])) {
      s <- s + 1
    } else {
      ss[i - 1] <- s
      s <- 1
      inc <- x[i] > x[i - 1]
    }
  }
  ss[i] <- s
  return(ss[ss != 0])
}

# Read Spectrum Helpers #####

#' @noRd
#' @title Read single Bruker TopSpin 3 Spectrum
#' @description For params and return value see [read_spectrum()].
#' @examples
#' spldir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' X <- read_bruker_spectrum(spldir)
#' fq_ref <- X$fq[1] / (1 - (X$cs[1] / 1e6))
#' print(head(X))
#' cat("Frequency of reference in MHz:", fq_ref / 1e6)
read_bruker_spectrum <- function(spldir = file.path(download_example_datasets(), "bruker/urine/urine_1"),
                                 expno = 10,
                                 procno = 10,
                                 raw = FALSE,
                                 silent = TRUE,
                                 force = FALSE) {
    acqus <- read_acqus_file(spldir, expno)
    procs <- read_procs_file(spldir, expno, procno)
    one_r <- read_1r_file(spldir, expno, procno, procs, silent = TRUE)[c("raw", "scaled")]
    make_spectrum(
        si = if (raw) one_r$raw else one_r$scaled, # Signal intensities
        cs_max = procs$OFFSET, # Spectrum offset in PPM
        cs_width = acqus$SW, # Spectrum width in PPM
        fq_ref = acqus$SFO1 * 1e6, # Reference Frequency in Hz (better than procs$SF, because it fulfills `all.equal` check of `make_spectrum`)
        fq_width = acqus$SW_h, # Spectrum width in Hz
        force = force,
        silent = silent
    )
}

#' @noRd
#' @title Read single JCAMPDX Spectrum
#' @description For params and return value see [read_spectrum()].
#' \dontrun{
#' # Example Usage (took 30s on the development machine)
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1.dx")
#' spectrum_data <- read_jcampdx_spectrum(path)
#' str(spectrum_data, 1)
#' }
read_jcampdx_spectrum <- function(path, raw = FALSE, silent = TRUE, force = FALSE) {
    data <- readJDX::readJDX(file = path, SOFC = TRUE, debug = 0) # Example return: [dataGuide=df(3*3), metadata=chr(2180), commentLines=int(10053), real=df(131072*2), imaginary=df(131072*2)] with colnames(data$real) = c("x", "y"). Takes about 30s on machine r31 for urine_1.dx (1MB).
    meta <- parse_metadata_file(lines = data$metadata)
    si_raw <- data$real$y
    si_scaled <- if (meta$DTYPP == 0) data$real$y * 2 ^ meta$NC_proc else data$real$y
    make_spectrum(
        si = if (raw) si_raw else si_scaled, # Signal intensities
        cs_max = meta$OFFSET, # Spectrum offset in PPM
        cs_width = meta$SW, # Spectrum width in PPM
        fq_ref = meta$SFO1 * 1e6, # Reference Frequency in Hz (better than meta$SF, because it fulfills `all.equal` check of `make_spectrum`)
        fq_width = meta$SW_h, # Spectrum width in Hz
        force = force,
        silent = silent
    )
}

#' @noRd
#' @title Read Bruker TopSpin Acquistion Parameters
#' @param spldir The path of the directory holding the NMR measurements for a individual sample. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @return The signals acquisition parameters read from the file as named list.
#' @examples
#' blood1_dir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' acqus <- read_acqus_file(blood1_dir)
#' str(acqus, 0)
#' cat("spectrum width ppm:", as.numeric(acqus$SW))
#' cat("spectrum width Hz:", as.numeric(acqus$SW_h))
read_acqus_file <- function(spldir = pkg_file("example_datasets/bruker/urine/urine_1"),
                            expno = 10) {
    path <- file.path(spldir, expno, "acqus")
    acqus <- parse_metadata_file(path)
    acqus
}

#' @noRd
#' @title Read Bruker TopSpin Processing Parameters
#' @param spldir The path of the directory holding the NMR measurements for a individual sample. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @return The processing parameters read from the file as named list.
#' @examples
#' blood1_dir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' procs <- read_procs_file(blood1_dir)
read_procs_file <- function(spldir = pkg_file("example_datasets/bruker/urine/urine_1"),
                            expno = 10,
                            procno = 10) {
    path <- file.path(spldir, expno, "pdata", procno, "procs")
    procs <- parse_metadata_file(path)
    procs
}

#' @noRd
#' @title Read signal intensities from Bruker TopSpin 1r file
#' @param spldir spldir The path of the directory holding the NMR measurements for a individual sample. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @param procs The parsed content of the `procs` file as returned by `read_procs_file()`.
#' @param force If `TRUE`, the function will try to read the file as 64 bit floating point numbers if the processing parameter `DTYPP` is unequal 0. This behaviour is untested, so should be used with caution.
#' @param silent If `TRUE`, no output will be printed to the console.
#' @return A named list with following elements:
#' - `spldir` The path of the directory holding the NMR measurements for a individual sample.
#' - `expno` The experiment number for the file.
#' - `procno` The processing number for the file.
#' - `path_1r` The path of the 1r file.
#' - `path_procs` The path of the procs file.
#' - `procs` The parsed content of the `procs` file as returned by `read_procs_file()`.
#' - `byteordp` The byte ordering of the data. 0 = little endian, 1 = big endian.
#' - `dtypp` The data type of the data. 0 = integer, not 0 = double.
#' - `endian` The endianess of the data. Either "little" or "big".
#' - `nbytes` The number of bytes used to store a single value. Either 4 or 8.
#' - `ncproc` The exponent of the data. Only relevant if `dtypp` is 0, i.e. the data is stored as integer values.
#' - `type` The type of the data. Either "integer" or "double", derived from `dtypp`.
#' - `n` The number of data points.
#' - `raw` The raw signal intensity values.
#' - `scaled` The scaled signal intensity values.
#' The first 5 elements in the list (`spldir` - `path_procs`) are path related values. The next 7 elements (`procs` - `n`) are processing parameters and derived values. The last 2 elements (`raw` - `scaled`) are the raw and scaled signal intensity values.
#' @examples
#' spldir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' oneR <- read_1r_file(spldir, 10, 10)
#' str(oneR, 1)
read_1r_file <- function(spldir = pkg_file("example_datasets/bruker/urine/urine_1"),
                         expno = 10,
                         procno = 10,
                         procs = read_procs_file(spldir, expno, procno),
                         force = FALSE,
                         silent = FALSE) {
    # Bruker_NMR_Data_Formats.pdf:
    #
    # > The raw data files `fid` and `ser` contain one dimensional or multi-dimensional
    # > acquired data, respectively. They consist of a sequence of acquired data point
    # > values in binary format. The acquisition status parameter `DTYPA` defines, how
    # > the data values are stored. If the `DTYPA` is "int" the stored value represents
    # > a mantissa of the data point value, the acquisition parameter NC is the
    # > exponent. All data points share in this case the same exponent. If `DTYPA` is
    # > "double", the data points are stored as a double precision 64 bit floating
    # > number, parameter NC is not used.
    # >
    # > | FIGURE A (TopSpin 3)         | FIGURE B (TopSin 4)       |
    # > |------------------------------|---------------------------|
    # > | DTYPA/DTYPP = 0 (int) ==>    | DTYPA/DTYPP = ? (dbl) ==> |
    # > | Value = 32BitInt * 2^NC      | Value = 64BitDouble with  |
    # > |                              | bits 00 - 51 = fraction   |
    # > |                              | bits 52 - 62 = exponent   |
    # > |                              | bit  63      = sign       |
    # >
    # > The processing status parameter `DTYPP` defines how the data values are stored. If
    # > the `DTYPP` is 0 ("int"), the stored value represents a mantissa of the data
    # > point value, the processing status parameter `NC_proc` is the exponent. In this
    # > case all data points share the same exponent.
    # >
    # > Their format is given by the parameter `DTYPP`, the byte ordering is given by
    # > the parameter `BYTORDP`, both may be read from the processing status parameter
    # > file `procs`.
    #
    path_1r <- file.path(spldir, expno, "pdata", procno, "1r")
    path_procs <- file.path(spldir, expno, "pdata", procno, "procs")
    byteordp <- procs$BYTORDP
    endian <- if (byteordp == 0) "little" else "big"
    ncproc <- procs$NC_proc
    dtypp <- procs$DTYPP
    n <- procs$SI
    if (dtypp == 0) {
        msg <- sprintf("Reading '%s' as 32 bit integers", path_1r)
        type <- "integer"
        nbytes <- 4
    } else if (force) {
        msg <- sprintf("Reading '%s' as 64 bit floating point numbers, because processing parameter DTYPP equals '%s' and `force == TRUE`. This behaviour is untested, so please double check the returned values. For details see Bruker_NMR_Data_Formats.pdf.", path_1r,  dtypp)
        type <- "double"
        nbytes <- 8
    } else {
        msg <- sprintf("Processing parameter `DTYPP` has value '%s' but only '0' is supported. This indicates that the intensity values in file '%s' are stored doubles and not as integers. To read the file nonetheless, set `force = TRUE`, but note that this behaviour is completely untested, so please double check the returned values.", dtypp, path_1r)
        stop(msg)
    }
    if (!silent) logf(msg)
    con <- file(path_1r, "rb"); on.exit(close(con), add = TRUE)
    raw <- readBin(con, what = type, n = n, size = nbytes, signed = TRUE, endian = endian)
    scaled <- if (type == "integer") raw * 2 ^ ncproc else raw
    named(
        spldir, expno, procno, path_1r, path_procs, # path related variables
        procs, byteordp, dtypp, endian, nbytes, ncproc, type, n,# processing parameters and derived values
        raw, scaled # raw and scaled signal intensity values
    )
}

#' @noRd
#' @title Parse Metadata File
#' @description Parses a metadata file like Bruker's `acqu[s]` or `proc[s]` files and return the metadata as a named list.
#' @param path The path of the file containing the parameter data. E.g. `"example_datasets/bruker/urine/urine_1/10/acqus"` or `"example_datasets/bruker/urine/urine_1/10/pdata/10/procs"`.
#' @details For a detailed description of the format of burker parameter files, refer to 'Bruker_NMR_Data_Formats.pdf'.
#' @return A named list containing the metadata read from the file.
#' @examples
#' path <- pkg_file("example_datasets/bruker/urine/urine_1/10/acqus")
#' lines <- readLines(path)
#' acqus1 <- parse_metadata_file(path)
#' acqus2 <- parse_metadata_file(lines = lines)
#' stopifnot(all.equal(acqus1, acqus2))
parse_metadata_file <- function(path = NULL, lines = NULL) {
    if (is.null(path) && is.null(lines)) stop("Either `path` or `lines` must be provided.")
    if (is.null(lines)) lines <- readLines(path)
    lines <- lines[!startsWith(lines, "$$")] # strip comments
    content <- paste(lines, collapse = "") # Example: "##TITLE= Parameter file, TopSpin 3.6.2##JCAMPDX= 5.0"
    pattern <- "(##\\$?.+?= ?)([^#]*)"
    matches <- gregexpr(pattern, content, perl = TRUE)
    keyvals <- regmatches(content, matches)[[1]]
    tmp <- strsplit(keyvals, "= ?")
    keys <- sapply(tmp, "[", 1)
    keys <- gsub("^##\\$?", "", keys)
    vals <- sapply(tmp, "[", 2)
    ret <- structure(as.list(vals), names = keys)
    flt <- is_float_str(vals)
    int <- is_int_str(vals)
    ret[flt] <- as.numeric(vals[flt])
    ret[int] <- as.integer(vals[int])
    ret
}
