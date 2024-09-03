# Spectra #####

#' @export
#' @title Spectra Class
#' @description Each object of this class is a list of `spectrum` objects with class attribute `spectra`.
#' @examples
#' fq <- c(0.599771, 0.691471, 0.783171, 0.874871, 0.966571) + 600248506
#' cs <- c(7.164231, 7.164079, 7.163926, 7.163773, 7.163620)
#' si <- c(100861.5, 105161.5, 109285.3, 113199.5, 116719.7)
#' si1 <- si + rnorm(5, 0, 100)
#' si2 <- si + rnorm(5, 0, 100)
#' spectrum1 <- spectrum(si1, cs, fq)
#' spectrum2 <- spectrum(si2, cs, fq)
#' spectra <- spectra(spectrum1, spectrum2)
#' try(spectra(spectrum1, 2))
#' try(spectra(spectrum1, 2, spectrum2, "a", "b"))
spectra <- function(...) {
    spectra_obj <- list(...)
    class_list <- lapply(spectra_obj, class)
    is_ok <- sapply(class_list, function(x) "spectrum" %in% x)
    nok <- which(!is_ok)
    if (length(nok) == 1) stop("Argument ", nok, " is not of class spectrum")
    if (length(nok) >= 2) stop("Arguments ", collapse(nok, last = " and "), " are not of class spectrum")
    structure(spectra_obj, class = "spectra")
}

#' @export
#' @title Read one or more Spectra
#' @description Read one or more spectra files or folders from disk and return each parsed spectrum as dataframe.
#' `read_spectrum()` reads a single spectrum and returns its signal intensities, chemical shifts and frequencies as dataframe.
#' `read_spectra()` can be used to read multiple spectra at once and returns a list of dataframes in the above mentioned format.
#' @param path The path of the file/folder containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"` or `"example_datasets/bruker/urine/urine"`.
#' @param data_path The path of the directory holding the NMR measurements for a individual sample. E.g. `"example_datasets/bruker/urine"`.
#' @param file_format The file_format of the spectrum file. E.g. `"bruker"` or `"jcampdx"`.
#' @param expno The experiment number for the file. E.g. `"10"`. Only relevant if `file_format` equals `"bruker"`.
#' @param procno The processing number for the file. E.g. `"10"`. Only relevant if `file_format` equals `"bruker"`.
#' @param raw If `FALSE`, scales the returned signal intensities based on information available in the spectrum metadata, in particular `NC_proc`. For details see `processing-reference.pdf`, available at <https://www.bruker.com/en.html> at section 'Services & Support > Documentation & Manuals > Magnetic Resonance > Acquisition & Processing > TopSpin Processing Commands and Parameters' (requires login).
#' @param silent If `TRUE`, no output will be printed to the console.
#' @param force If `TRUE`, try to continue when encountering errors and print info messages instead. To hide these messages as well, set `silent = TRUE`.
#' @return For `read_spectrum` a 'spectrum' object as described in `make_spectrum()`. The object is a list with class `spectrum` and the following elements:
#'
#' - `si`: signal intensities in arbitrary units
#' - `cs`: chemical shifts in ppm
#' - `fq`: frequencies in Hz
#'
#' For `read_spectra` a named list of such dataframes, where the names are the file names of the spectra.
#' @details For details about `procno` and `expno` see section [File Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure) in the metabodecon FAQ.
#' @examples
#' urine <- system.file("example_datasets/bruker/urine", package = "metabodecon")
#' urine_1 <- file.path(urine, "urine_1")
#' urine_2 <- file.path(urine, "urine_2")
#' X1 <- read_spectrum(urine_1)
#' X2 <- read_spectrum(urine_2)
#' XX <- read_spectra(urine)
#' str(XX)
#' str(X1)
#' stopifnot(all.equal(X1, XX$urine_1))
#'
#' # Below code shows how a spectrum stored in JCAMP-DX format can be read.
#' # Reading files in this format is very slow (about 30s on the development
#' # machine). So if possible, you should stick with the original Bruker
#' # data storage format.
#' \dontrun{
#' urine_1_dx <- system.file("example_datasets/jcampdx/urine/urine_1.dx", package = "metabodecon")
#' X1_dx <- read_spectrum(urine_1_dx, file_format = "jcampdx")
#' stopifnot(all.equal(X1, X1_dx))
#' }
read_spectra <- function(data_path = pkg_file("example_datasets/bruker/urine"),
                         file_format = "bruker",
                         expno = 10,
                         procno = 10,
                         raw = FALSE,
                         silent = TRUE,
                         force = FALSE) {
    if (!file_format %in% c("bruker", "jcampdx")) {
        stop("Argument `file_format` should be either 'bruker' or 'jcampdx'")
    }
    dp <- normPath(data_path)
    jcampdx <- file_format == "jcampdx"
    bruker <- file_format == "bruker"
    r1_path <- file.path(dp, expno, "pdata", procno, "1r")
    r1_path_exists <- file.exists(r1_path)
    ends_with_dx <- grepl("\\.dx$", dp)
    if ((jcampdx && ends_with_dx) || (bruker && r1_path_exists)) {
        files <- basename(dp)
        paths <- dp
    } else if (jcampdx) {
        files <- dir(dp, pattern = "\\.dx$") # `.dx` files inside `path`
        paths <- dir(dp, pattern = "\\.dx$", full.names = TRUE)
    } else if (bruker) {
        files <- list.dirs(dp, recursive = FALSE, full.names = FALSE)
        paths <- list.dirs(dp, recursive = FALSE, full.names = TRUE) # folders inside `path`
        r1_paths <- file.path(paths, expno, "pdata", procno, "1r")
        r1_paths_exists <- file.exists(r1_paths)
        paths <- paths[r1_paths_exists]
        files <- files[r1_paths_exists]
    }
    if (length(files) == 0) {
        msg <- sprintf("No spectra found in directory '%s'.", data_path)
        if (file_format == "bruker") msg <- paste(msg, "Did you specify the correct `expno` and `procno`?")
        stop(msg)
    }
    spectra <- lapply(paths, function(path) {
        if (!silent) logf("Reading spectrum %s", path)
        read_spectrum(path, file_format, expno, procno, raw, silent, force)
    })
    names(spectra) <- files
    invisible(spectra)
}

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
#'          verbose = TRUE,
#'          noise_method = "SFR"
#'      )
#' }
simulate_spectra <- function(pngdir = NULL,
                             pdfdir = NULL,
                             svgdir = NULL,
                             rdsdir = NULL,
                             brukerdir = NULL,
                             verbose = TRUE,
                             noise_method = "SFR") {
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
            brukerdir = file.path(brukerdir, new_name),
            verbose = verbose,
            noise_method = noise_method
        )
    })
    invisible(named(sim, pngdir, pdfdir, rdsdir))
}
