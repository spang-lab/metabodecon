# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Public API #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @export
#' @family {spectra functions}
#' @inherit read_spectrum
#' @title Test, Convert or Read Spectra Files from Disk
#' @description
#' `read_spectra()` reads spectra files from disk.
#' `is_spectra()` tests if an object is a `spectra` object.
#' `as_spectra()` convert an object to a `spectra` object.
#' For details about `spectra` objects, see [metabodecon_classes].
#' @return
#' `read_spectra()` and `as_spectra` return a `spectra` object.
#' `is_spectra()` returns a logical.
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
    dp <- norm_path(data_path)
    jcampdx <- file_format == "jcampdx"
    bruker <- file_format == "bruker"
    r1_path <- file.path(dp, expno, "pdata", procno, "1r")
    r1_path_exists <- file.exists(r1_path)
    ends_with_dx <- grepl("\\.dx$", dp)
    if ((jcampdx && ends_with_dx) || (bruker && r1_path_exists)) {
        files <- basename(dp)
        paths <- dp
    } else if (jcampdx) {
        files <- dir(dp, pattern = "\\.dx$")
        paths <- dir(dp, pattern = "\\.dx$", full.names = TRUE)
    } else if (bruker) {
        files <- list.dirs(dp, recursive = FALSE, full.names = FALSE)
        paths <- list.dirs(dp, recursive = FALSE, full.names = TRUE)
        r1_paths <- file.path(paths, expno, "pdata", procno, "1r")
        r1_paths_exists <- file.exists(r1_paths)
        paths <- paths[r1_paths_exists]
        files <- files[r1_paths_exists]
    }
    if (length(files) == 0) {
        msg <- sprintf("No spectra found in directory '%s'.", data_path)
        msg2 <- "Did you specify the correct `expno` and `procno`?"
        if (file_format == "bruker") msg <- paste(msg, ms2)
        stop(msg)
    }
    spectra <- lapply(paths, function(path) {
        if (!silent) logf("Reading spectrum %s", path)
        read_spectrum(path, file_format, expno, procno, raw, silent, force)
    })
    names(spectra) <- files
    class(spectra) <- "spectra"
    invisible(spectra)
}

#' @export
#' @family {spectrum functions}
#'
#' @title Read a spectrum from Disk
#'
#' @description
#' Read a single spectrum file or folder from disk and return it as `spectrum`
#' object.
#'
#' @param data_path The path of the file/folder containing the spectrum data.
#' E.g. `"example_datasets/jcampdx/urine/urine_1.dx"` or
#' `"example_datasets/bruker/urine/urine"`.
#'
#' @param file_format The file_format of the spectrum file. E.g. `"bruker"` or
#' `"jcampdx"`.
#'
#' @param expno,procno The experiment/processing number for the file. E.g.
#' `"10"`. Only relevant if `file_format` equals `"bruker"`. For details see
#' section [File Structure](
#' https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure
#' ) in the metabodecon FAQ.
#'
#' @param raw If `FALSE`, scales the returned signal intensities based on
#' information available in the spectrum metadata, in particular `NC_proc`. For
#' details see `processing-reference.pdf`, available at
#' <https://www.bruker.com/en.html> at section 'Services & Support >
#' Documentation & Manuals > Magnetic Resonance > Acquisition & Processing >
#' TopSpin Processing Commands and Parameters' (requires login).
#'
#' @param silent If `TRUE`, no output will be printed to the console.
#'
#' @param force If `TRUE`, try to continue when encountering errors and print
#' info messages instead. To hide these messages as well, set `silent = TRUE`.
#'
#' @return
#' A `spectrum` object as described in [Metabodecon Classes](
#' https://spang-lab.github.io/metabodecon/articles/Metabodecon-Classes.html
#' ).
#'
#' @examples
#' relpath <- "example_datasets/bruker/urine"
#' urine <- system.file(relpath, package = "metabodecon")
#' urine_1 <- file.path(urine, "urine_1")
#' urine_2 <- file.path(urine, "urine_2")
#' x1 <- read_spectrum(urine_1)
#' x2 <- read_spectrum(urine_2)
#' xx <- read_spectra(urine)
#' str(xx)
#' str(x1)
#' stopifnot(all.equal(x1, xx$urine_1))
#'
#' ## Below code shows how a spectrum stored in JCAMP-DX format can be read.
#' ## Reading files in this format is very slow (about 30s on the development
#' ## machine). So if possible, you should stick with the original Bruker
#' ## data storage format.
#' \dontrun{
#' relpath <- "example_datasets/jcampdx/urine/urine_1.dx"
#' urine_1_dx <- system.file(relpath, package = "metabodecon")
#' x1_dx <- read_spectrum(urine_1_dx, file_format = "jcampdx")
#' stopifnot(all.equal(x1, x1_dx))
#' }
read_spectrum <- function(data_path = metabodecon_file("bruker/sim/sim_01"),
                          file_format = "bruker",
                          expno = 10,
                          procno = 10,
                          raw = FALSE,
                          silent = TRUE,
                          force = FALSE) {
    file_format <- match.arg(file_format, c("bruker", "jcampdx"))
    switch(file_format,
        "bruker" = read_bruker_spectrum(data_path, expno, procno, raw, silent, force),
        "jcampdx" = read_jcampdx_spectrum(data_path, silent, force)
    )
}

#'
#' @export
#'
#' @family {spectrum functions}
#'
#' @title Create a Spectrum Object
#'
#' @description
#' Creates a spectrum object from the provided signal intensities, frequencies
#' and chemical shifts.
#'
#' @param si Numeric vector of signal intensities, ordered from highest to
#' lowest corresponding chemical shift.
#'
#' @param cs_max The highest chemical shift value in ppm, usually shown as left
#' end of the spectrum.
#'
#' @param cs_width The width of the spectrum in ppm.
#'
#' @param fq_ref The reference frequency in Hz.
#'
#' @param fq_width The width of the spectrum in Hz. Only used to check whether
#' the values calculated from `cs_max`, `cs_width` and `fq_ref` match the
#' provided value. If `NULL`, this check will be skipped.
#'
#' @param force If `TRUE`, the function will not raise an error in case of
#' discrepancies between the calculated and the provided spectrum width in Hz,
#' but will print a info message instead. To hide this message as well, set
#' `silent = TRUE`.
#'
#' @param silent If `TRUE`, no output will be printed to the console.
#'
#' @param name The name of the spectrum, e.g. "Blood 1" or "Urine Mouse X23D".
#'
#' @param path The path to the spectrum file, e.g.
#' "/example_datasets/bruker/urine/urine_1".
#'
#' @param type The type of experiment, e.g. "H1 CPMG" or "H1 NOESY".
#'
#' @param mfs The magnetic field strength in Tesla.
#'
#' @return A `spectrum` object as described in [Metabodecon Classes](
#' https://spang-lab.github.io/metabodecon/articles/Metabodecon-Classes.html
#' ).
#'
#' @examples
#' si <- c(1, 1, 3, 7, 8, 3, 8, 5, 2, 1)
#' cs_max <- 14.8
#' cs_width <- 20.0
#' fq_ref <- 600.25 * 1e6
#' fq_width <- 12005
#' spectrum <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width)
#' spectrum2 <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width = 12010, force = FALSE)
make_spectrum <- function(si,
                          cs_max,
                          cs_width,
                          fq_ref,
                          fq_width = NULL,
                          force = FALSE,
                          silent = FALSE,
                          name = NULL,
                          path = NULL,
                          type = NULL,
                          simpar = NULL,
                          mfs = NULL) {
    cs_min <- cs_max - cs_width # Lowest ppm value
    cs <- seq(cs_max, cs_max - cs_width, length.out = length(si)) # Chemical shift in parts per million
    fq <- as_frequency(cs, fq_ref)
    fq_width_calc <- max(fq) - min(fq)
    if (!is.null(fq_width) && !is_equal(fq_width_calc, fq_width)) {
        if (force) {
            msg <- "Calculated spectrum width in Hz (%s) does not match the provided value (%s). Please read in the data manually or set `force = TRUE` to ignore this error. Please note that by doing so, all downstream calculations involving frequencies might be wrong, so be sure to double check the results."
            stopf(msg, round(fq_width_calc, 5), round(fq_width, 5))
        } else if (!silent) {
            msg <- "Calculated spectrum width in Hz (%s) does not match the provided value (%s). Continuing anyways, because `force` equals `TRUE`. Please note that all downstream calculations using frequencies might be wrong, so be sure to double check the results."
            logf(msg, round(fq_width_calc, 5), round(fq_width, 5))
        }
    }
    meta <- named(fq, name, path, type, simpar, mfs)
    structure(named(si, cs, meta), class = "spectrum")
}

#' @export
#'
#' @title Simulate a 1D NMR Spectrum
#'
#' @description
#' Simulates a 1D NMR spectrum based on the provided parameters.
#'
#' `r lifecycle::badge("experimental")`
#'
#' @param name The name of the spectrum.
#' @param seed The seed for the random number generator.
#' @param ndp The number of data points in the spectrum.
#' @param npk The number of peaks in the spectrum.
#' @param csres The chemical shift resolution in PPM.
#' @param cs The vector of chemical shifts in PPM.
#' @param pkr The start and stop of the peak region in PPM.
#' @param fqref The reference frequency in Hz.
#' @param x0 The peak center positions in PPM.
#' @param A The peak area parameter.
#' @param lambda The peak width parameter.
#' @param noise The noise to add to the spectrum.
#'
#' @return A `spectrum` object as described in [metabodecon_classes].
#'
#' @examples
#' simA <- simulate_spectrum("simA")
#' simA_copy <- simulate_spectrum("simA")
#' simB <- simulate_spectrum("simB")
#' simC <- simulate_spectrum("simC", npk = 20)
#' plot_spectrum(simC, foc_rgn = c(0.25, 0.75))
#' if (!identical(simA, simA_copy)) stop()
#' if (identical(simA, simB)) stop()
simulate_spectrum <- function(name = "sim_00",
                              seed = sum(utf8ToInt(name)), # (1)
                              ndp = 2048,
                              npk = 10,
                              csres = 0.00015,
                              cs = seq(from = 3.6, length.out = ndp, by = -csres),
                              pkr = quantile(cs, c(0.25, 0.75)),
                              fqref = 600252806.95,
                              x0 = sort(runif(npk, pkr[1], pkr[2])),
                              A = runif(npk, 2.5, 20) * 1e3,
                              lambda = runif(npk, 0.9, 1.3) / 1e3,
                              noise = rnorm(length(cs), sd = 1200)) {
    set.seed(seed) # (1)
    si <- round(lorentz_sup(cs, x0, A, lambda) + noise) # (2)
    fq <- as_frequency(cs, fqref)
    simpar <- named(name, cs, pkr, x0, A, lambda, noise)
    meta <- named(name, fq, simpar)
    spectrum <- structure(named(si, cs, meta), class = "spectrum")
    spectrum
    # (1) Setting the seed in the function also affects runif and rnorm calls in
    #     the function default values as they are only evaluated upon first use.
    # (2) We convert to integers to not lose precision later on when the data
    #     gets stored as integers on disk.
}



# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Private Helpers #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @noRd
#' @title Read single Bruker TopSpin 3 Spectrum
#' @inheritParams read_spectrum
#' @examples
#' spldir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' x <- read_bruker_spectrum(spldir)
#' fq_ref <- x$fq[1] / (1 - (x$cs[1] / 1e6))
#' print(head(x))
#' cat("Frequency of reference in MHz:", fq_ref / 1e6)
read_bruker_spectrum <- function(spldir,
                                 expno = 10,
                                 procno = 10,
                                 raw = FALSE,
                                 silent = TRUE,
                                 force = FALSE) {
    acqus <- read_acqus_file(spldir, expno)
    procs <- read_procs_file(spldir, expno, procno)
    simpar <- read_simpar_file(spldir, expno, procno)
    one_r <- read_one_r_file(spldir, expno, procno, procs, silent = silent)
    one_r <- one_r[c("raw", "scaled")]
    x <- make_spectrum(
        si = if (raw) one_r$raw else one_r$scaled, # Signal intensities
        cs_max = procs$OFFSET, # Spectrum offset in PPM
        cs_width = acqus$SW, # Spectrum width in PPM
        fq_ref = acqus$SFO1 * 1e6, # Reference Frequency in Hz (1)
        fq_width = acqus$SW_h, # Spectrum width in Hz
        force = force,
        silent = silent,
        name = basename(spldir),
        simpar = simpar,
        path = spldir
        # (1) Better than procs$SF, because it fulfills `all.equal` check of
        # `make_spectrum`.
    )
}

#' @noRd
#' @title Read single JCAMPDX Spectrum
#' @inheritParams read_spectrum
#' @examples
#' \dontrun{
#' # Example Usage (took 30s on the development machine)
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1.dx")
#' spectrum_data <- read_jcampdx_spectrum(path)
#' str(spectrum_data, 1)
#' }
read_jcampdx_spectrum <- function(path,
                                  raw = FALSE,
                                  silent = TRUE,
                                  force = FALSE) {
    data <- readJDX::readJDX(file = path, SOFC = TRUE, debug = 0)
    # Object `data` is a list with following elements:
    # - dataGuide = df (3*3)
    # - metadata = chr (2180)
    # - commentLines = int (10053),
    # - real = df (131072*2)
    # - imaginary = df (131072*2)]
    # Colnames(data$real) are c("x", "y")
    # Reading takes about 30s on machine r31 for urine_1.dx (1MB).
    meta <- parse_metadata_file(lines = data$metadata)
    si_raw <- data$real$y
    si_scaled <- if (meta$DTYPP == 0) data$real$y * 2^meta$NC_proc else data$real$y
    make_spectrum(
        si = if (raw) si_raw else si_scaled, # Signal intensities
        cs_max = meta$OFFSET, # Spectrum offset in PPM
        cs_width = meta$SW, # Spectrum width in PPM
        fq_ref = meta$SFO1 * 1e6, # Reference Frequency in Hz (1)
        fq_width = meta$SW_h, # Spectrum width in Hz
        force = force,
        silent = silent,
        path = path,
        name = basename(path)
        # (1) Better than meta$SF, because it fulfills `all.equal` check of
        #     [make_spectrum()]
    )
}

#' @noRd
#' @title Parse Metadata File
#' @description
#' Parses a metadata file like Bruker's `acqu[s]` or `proc[s]` files and return
#' the metadata as a named list.
#'
#' @param path The path of the file containing the parameter data. E.g.
#' `"example_datasets/bruker/urine/urine_1/10/acqus"` or
#' `"example_datasets/bruker/urine/urine_1/10/pdata/10/procs"`.
#'
#' @return A named list containing the metadata read from the file.
#'
#' @details
#' For a detailed description of the format of burker parameter files, refer to
#' 'Bruker_NMR_Data_Formats.pdf'.
#'
#' @examples
#' path <- pkg_file("example_datasets/bruker/urine/urine_1/10/acqus")
#' lines <- readLines(path)
#' acqus1 <- parse_metadata_file(path)
#' acqus2 <- parse_metadata_file(lines = lines)
#' stopifnot(all.equal(acqus1, acqus2))
parse_metadata_file <- function(path = NULL, lines = NULL) {
    msg <- "Either `path` or `lines` must be provided."
    if (is.null(path) && is.null(lines)) stop(msg)
    if (is.null(lines)) lines <- readLines(path)
    lines <- lines[!startsWith(lines, "$$")] # Strip comments
    content <- paste(lines, collapse = "") # (1)
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
    # (1) Example content of metadata file:
    # >>> ##TITLE= Param file, TopSpin 3.6
    # >>> ##JCAMPDX= 5.0
    # >>> ...
}

#' @noRd
#'
#' @title Read Bruker TopSpin Acquistion Parameters
#'
#' @inheritParams read_spectrum
#' @param spldir The path of the directory holding the NMR measurements for a
#' individual sample. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#'
#' @return
#' The signals acquisition parameters read from the file as named list.
#'
#' @examples
#' blood1_dir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' acqus <- read_acqus_file(blood1_dir)
#' str(acqus, 0)
#' cat("spectrum width ppm:", as.numeric(acqus$SW))
#' cat("spectrum width Hz:", as.numeric(acqus$SW_h))
read_acqus_file <- function(spldir, expno = 10) {
    path <- file.path(spldir, expno, "acqus")
    acqus <- parse_metadata_file(path)
    acqus
}

#' @noRd
#' @title Read Bruker TopSpin Processing Parameters
#' @inheritParams read_spectrum read_acqus_file
#' @return The processing parameters read from the file as named list.
#' @examples
#' blood1_dir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' procs <- read_procs_file(blood1_dir)
read_procs_file <- function(spldir, expno = 10, procno = 10) {
    path <- file.path(spldir, expno, "pdata", procno, "procs")
    procs <- parse_metadata_file(path)
    procs
}

#' @noRd
#' @title Read Simulation Parameters
#' @details
#' Usually, spectra stored in Bruker Format have three files that are relevant
#' for metabodecon: `acqus`, `procs` and `1r`. However, if a spectrum is
#' simulated using [simulate_spectrum()] and stored in Bruker format using
#' [save_spectrum()], an additional file called `simpar` is created, containing
#' the simulation parameters that were used to generate the spectrum. The file
#' contains R code, as generated by `dput()`, that can be used to recreate the
#' list object holding the simulation parameters.
#' @inheritParams read_spectrum read_acqus_file
#' @return The simulation parameters read from the file as named list.
#' @examples
#' sim_1_dir <- pkg_file("example_datasets/bruker/sim2/sim_01")
#' simpar <- read_simpar_file(blood1_dir)
read_simpar_file <- function(spldir, expno = 10, procno = 10) {
    path <- file.path(spldir, expno, "pdata", procno, "simpar")
    read <- switch(simpar_format, "rds" = readRDS, "ascii" = dget)
    if (file.exists(path)) read(path)
}

#' @noRd
#' @title Read signal intensities from Bruker TopSpin 1r file
#' @inheritParams read_spectrum read_acqus_file
#' @examples
#' spldir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' oneR <- read_one_r_file(spldir, 10, 10)
#' str(oneR, 1)
read_one_r_file <- function(spldir,
                            expno = 10,
                            procno = 10,
                            procs = read_procs_file(spldir, expno, procno),
                            force = FALSE,
                            silent = FALSE) {
    # Bruker_NMR_Data_Formats.pdf:
    #
    # > The raw data files `fid` and `ser` contain one dimensional or
    # > multi-dimensional acquired data, respectively. They consist of a
    # > sequence of acquired data point values in binary format. The acquisition
    # > status parameter `DTYPA` defines, how the data values are stored. If the
    # > `DTYPA` is "int" the stored value represents a mantissa of the data
    # > point value, the acquisition parameter NC is the exponent. All data
    # > points share in this case the same exponent. If `DTYPA` is "double", the
    # > data points are stored as a double precision 64 bit floating number,
    # > parameter NC is not used.
    # >
    # > | FIGURE A (TopSpin 3)         | FIGURE B (TopSin 4)       |
    # > |------------------------------|---------------------------|
    # > | DTYPA/DTYPP = 0 (int) ==>    | DTYPA/DTYPP = ? (dbl) ==> |
    # > | Value = 32BitInt * 2^NC      | Value = 64BitDouble with  |
    # > |                              | bits 00 - 51 = fraction   |
    # > |                              | bits 52 - 62 = exponent   |
    # > |                              | bit  63      = sign       |
    # >
    # > The processing status parameter `DTYPP` defines how the data values are
    # > stored. If the `DTYPP` is 0 ("int"), the stored value represents a
    # > mantissa of the data point value, the processing status parameter
    # > `NC_proc` is the exponent. In this case all data points share the same
    # > exponent.
    # >
    # > Their format is given by the parameter `DTYPP`, the byte ordering is
    # > given by the parameter `BYTORDP`, both may be read from the processing
    # > status parameter file `procs`.
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
        msg <- sprintf("Reading '%s' as 64 bit floating point numbers, because processing parameter DTYPP equals '%s' and `force == TRUE`. This behaviour is untested, so please double check the returned values. For details see Bruker_NMR_Data_Formats.pdf.", path_1r, dtypp)
        type <- "double"
        nbytes <- 8
    } else {
        msg <- sprintf("Processing parameter `DTYPP` has value '%s' but only '0' is supported. This indicates that the intensity values in file '%s' are stored doubles and not as integers. To read the file nonetheless, set `force = TRUE`, but note that this behaviour is completely untested, so please double check the returned values.", dtypp, path_1r)
        stop(msg)
    }
    if (!silent) logf(msg)
    con <- file(path_1r, "rb")
    on.exit(close(con), add = TRUE)
    raw <- readBin(con, what = type, n = n, size = nbytes, signed = TRUE, endian = endian)
    scaled <- if (type == "integer") raw * 2^ncproc else raw
    named(
        # -~-~-~-~-~-~-~-~-~-~-~
        # Path related variables
        # -~-~-~-~-~-~-~-~-~-~-~
        spldir, # Path of dir holding NMR measurements of a individual sample.
        expno, # Experiment number for the file.
        procno, # Processing number for the file.
        path_1r, # Path of the 1r file.
        path_procs, # Path of the procs file.
        # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
        # Processing parameters and derived values
        # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
        procs, # Parsed content of `procs` file from `read_procs_file()`.
        byteordp, # Byte ordering of the data. 0/1 = little/big endian.
        dtypp, # Type of the data. 0/not0 = integer/double.
        endian, # Endianess of the data. Either "little" or "big".
        nbytes, # Number of bytes used to store a single value. Either 4 or 8.
        ncproc, # Exponent of the data. Only relevant if `dtypp` is 0.
        type, # Type as string. Either "integer" or "double". Like `dtypp`.
        n, # Number of data points.
        # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
        # Raw and scaled signal intensity values
        # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
        raw, # The raw signal intensity values.
        scaled # The scaled signal intensity values.
    )
}

#' @noRd
#' @description
#' Converts a vector of chemical shifts given in ppm to a vector of frequencies
#' in Hz.
#' @param cs Vector of chemical shifts in ppm.
#' @param fqref Frequency of the reference molecule in Hz.
as_frequency <- function(cs, fqref) {
    fqmax <- fqref - (min(cs) * 1e-6 * fqref) # Highest frequency in Hz
    fqmin <- fqref - (max(cs) * 1e-6 * fqref) # Lowest frequency in Hz
    fq <- seq(fqmin, fqmax, length.out = length(cs))
    fq <- fq[rev(order(cs))] # (1)
    # (1) Sort frequencies in descending order of chemical shifts, as the
    #     highest chemical shift corresponds to the lowest frequency.
    fq
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Sim Dataset #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

update_sim_dataset <- function(dry_run = TRUE, verbose = TRUE) {
    parent <- if (dry_run) tmpdir() else pkg_file("example_datasets/bruker")
    path <- file.path(parent, "sim2")
    if (verbose) logf("Dry_run %s. Updating %s." , dry_run, path)
    x <- make_sim_dataset()
    save_spectra(x, path, force = FALSE, verbose = verbose)
    if (verbose) logf("Finished update of %s." , path)
}

make_sim_dataset <- function() {
    decons <- deconvolute_blood()
    simpars <- lapply(decons, get_sim_params, pkr = c(3.52, 3.37))
    simpars[[2]] <- get_sim_params(decons[[2]], pkr = c(3.51, 3.36)) # (1)
    # (1) Blood_02 is shifted approx. 0.01 ppm to the right, so we also need to
    # shift the interval from which we pick our peaks by 0.01 to end up with
    # signals from the same metabolites. This was determined by visual
    # inspection after deconvoluting the blood dataset.
    sim <- lapply(simpars, function(simpar) {
        simulate_spectrum(
            name = gsub("blood", "sim", simpar$name),
            ndp = 2048,
            x0 = simpar$x0,
            A = simpar$A,
            lambda = simpar$lambda,
            noise = rnorm(2048, sd = simpar$noiseSD)
        )
    })
    names(sim) <- gsub("blood", "sim", names(sim))
    class(sim) <- "spectra"
    sim
}

#' @noRd
#'
#' @title Deconvolute the Blood Dataset
#'
#' @description
#' Downloads and deconvolutes the Blood Dataset.
#'
#' @param read_cache If TRUE the function will try to read the results from a
#' cache file instead of downloading and deconvoluting the dataset again.
#'
#' @param write_cache If TRUE, the results will be written to the cache file
#' path. If the cache file already exists and the results from the deconvolution
#' and the stored object differ, a warning will be issued and the new results
#' will only be stored if force is TRUE.
#'
#' @param force If TRUE, the results will be written to the cache file path even
#' if the results from the deconvolution and the stored object differ.
#'
#' @return
#' A `decons2` object containing the deconvolution results. For details
#' see [deconvolute()].
#'
#' @examples
#' deconvolute_blood() # Use cache file if it exists
#' deconvolute_blood(read_cache = FALSE) # Deconvolute from scratch
#' deconvolute_blood(force = TRUE) # Force update of cache file
deconvolute_blood <- function(read_cache = TRUE,
                              write_cache = TRUE,
                              force = FALSE,
                              verbose = TRUE,
                              nworkers = 1) {
    rds <- file.path(cachedir(), "deconvolute_blood.rds")
    new <- old <- if (file.exists(rds)) readRDS(rds) else NULL # 62.6 MB ~= 0.6s
    if (!read_cache || is.null(new)) {
        download_example_datasets()
        path <- datadir("example_datasets/bruker/blood")
        new <- generate_lorentz_curves(
            data_path = path,
            nworkers = nworkers,
            verbose = verbose,
            ask = FALSE
        )
    }
    if (!(identical(new, old) || is.null(old))) {
        warning("Cache and deconvolution differ.", immediate. = TRUE)
    }
    if (write_cache && (is.null(old) || force)) {
        logf("Writing results to cache file %s", rds)
        saveRDS(new, rds)
    }
    return(new)
}

#' @noRd
#' @examples
#' xdir <- download_example_datasets()
#' path <- file.path(xdir, "bruker/blood/blood_01")
#' spec <- read_spectrum(path)
#' x <- deconvolute_gspec(spec, rtyp = "decon1")
#' sim_params <- get_sim_params(x)
get_sim_params <- function(x, pkr = c(3.4, 3.5)) {
    d <- as_decon1(x)
    p <- which(d$x_0_ppm <= max(pkr) & d$x_0_ppm >= min(pkr))
    A <- -d$A_ppm[p] * 1e6
    x0 <- d$x_0_ppm[p]
    lambda <- -d$lambda_ppm[p]
    noiseSD <- sd(d$y_values_raw[c(1:10000, 121073:131072)])
    name <- d$filename
    named(A, x0, lambda, noiseSD, name)
}

#' @noRd
#' @examples
#' x <- simulate_spectrum()
#' path <- save_spectrum(x)
#' y <- read_spectrum(path)
save_spectrum <- function(x,
                          path = tmpdir(),
                          force = FALSE,
                          verbose = TRUE) {

    # Check input args, init temp dir and log function
    stopifnot(is_spectrum(x))
    stopifnot(length(dir(path)) == 0 || isTRUE(force))
    stopifnot(is_bool(force))
    stopifnot(is_bool(verbose))
    temp <- tmpdir()
    logv <- get_logv(verbose)
    logv("Saving bruker files to %s", temp)

    # Prepare processing parameters to write to procs file
    BYTORDP <- "##$BYTORDP=0" # Byte order (0 = Little endian)
    NC_proc <- "##$NC_proc=0" # Exponent for intensity values (si = y * 2^NC_proc)
    DTYPP <- "##$DTYPP=0" # Data storage type (0=4-byte-integers, else=double)
    SI <- sprintf("##$SI=%d", length(x$si)) # Number of data points
    OFFSET <- sprintf("##$OFFSET=%.15f", max(x$cs)) # Maximum chemical shift in PPM
    procs_str <- paste(BYTORDP, NC_proc, DTYPP, SI, OFFSET, "", sep = "\n")
    procs_path <- file.path(temp, "10", "pdata", "10", "procs")

    # Prepare aquisition parameters to write to acqus file
    SW <- sprintf("##$SW=%.15f", width(x$cs)) # Spectrum width in PPM
    acqus_str <- paste(SW, "", sep = "\n")
    if (!is.null(x$meta$fq)) {
        fq_ref <- convert_pos(0, x$cs, x$meta$fq)
        SFO1 <- sprintf("##$SFO1=%.15f", fq_ref / 1e6) # Reference Frequency in MHz
        SW_h <- sprintf("##$SW_h=%.15f", width(x$meta$fq)) # Spectrum width in Hz
        acqus_str <- paste(SW, SFO1, SW_h, "", sep = "\n")
    }
    acqus_path <- file.path(temp, "10", "acqus")

    # Prepare signal intensities to write to one_r file
    if (!is_int(x$si)) stop("Signal intensities must be integers.")
    x$si <- as.integer(x$si)

    # Write files
    one_r_path <- file.path(temp, "10", "pdata", "10", "1r")
    mkdirs(file.path(temp, "10", "pdata", "10"))
    logv("Writing %s", procs_path); cat(procs_str, file = procs_path)
    logv("Writing %s", acqus_path); cat(acqus_str, file = acqus_path)
    logv("Writing %s", one_r_path); writeBin(x$si, one_r_path)
    if (!is.null(x$meta$simpar)) {
        simpar_path <- file.path(temp, "10", "pdata", "10", "simpar")
        logv("Writing %s", simpar_path)
        if (simpar_format == "rds") {
            saveRDS(x$meta$simpar, simpar_path)
        } else {
            simpar_str <- dput2(x$meta$simpar)
            cat(simpar_str, file = simpar_path)
        }
    }

    # Validate written files
    logv("Reading back files to validate them")
    y <- read_bruker_spectrum(spldir = temp)
    is_valid <- all(
        is_equal(x$si, y$si),
        is_equal(x$cs, y$cs),
        is_equal(x$meta$fq, y$meta$fq),
        is_equal(x$meta$simpar, y$meta$simpar)
    )
    if (!is_valid) stop("Stored spectrum is corrupted")

    # Move the files to the destination directory
    logv("Moving files to %s", path)
    mkdirs(path)
    from = file.path(temp, "10")
    file.copy(from, path, recursive = TRUE, overwrite = TRUE)

    return(path)
}

#' @noRd
#' @title Save Spectra to Disk in Bruker Format
save_spectra <- function(x, path, force = FALSE, verbose = TRUE) {
    stopifnot(is_spectra(x))
    stopifnot(length(dir(path)) == 0 || isTRUE(force))
    stopifnot(is_bool(force))
    stopifnot(is_bool(verbose))
    subpaths <- sapply(x, function(s) file.path(path, s$meta$name))
    mkdirs(path)
    for (i in seq_along(x)) {
        save_spectrum(
            x = x[[i]],
            path = subpaths[i],
            force = force,
            verbose = verbose
        )
    }
}

simpar_format <- "rds"

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Deprecated #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

simulate_from_decon <- function(x,
                                # Simulation Params
                                cs_min = 3.4,
                                cs_max = 3.6,
                                pk_min = 3.45,
                                pk_max = 3.55,
                                noise_method = "SFR",
                                # Output params
                                show = TRUE,
                                save = FALSE,
                                verbose = TRUE,
                                ...) {
    if (!isFALSE(verbose)) logf("Checking function inputs")
    stopifnot(is_num(cs_min, 1))
    stopifnot(is_num(cs_max, 1))
    stopifnot(is_num(pk_min, 1))
    stopifnot(is_num(pk_max, 1))
    stopifnot(is_char(noise_method, 1, "(RND|SFR)"))
    stopifnot(is_bool(show, 1))
    stopifnot(is_bool(save, 1))
    stopifnot(is_bool(verbose, 1))
    d <- as_decon1(x)
    logv <- if (verbose) logf else function(...) NULL

    logv("Simulating spectrum from '%s' (noise method = '%s')", d$filename, noise_method)
    s <- as_spectrum(d)

    logv("Throwing away datapoints outside of %.1f to %.1f", cs_min, cs_max)
    ix <- which(s$cs >= cs_min & s$cs <= cs_max)
    s$cs <- s$cs[ix]
    s$si <- s$si[ix]
    if (!is.null(s$meta$fq)) s$meta$fq <- s$meta$fq[ix]

    logv("Keeping only within %.1f to %.1f", pk_min, pk_max)
    ip <- which(d$x_0_ppm <= pk_max & d$x_0_ppm >= pk_min)
    s$lcpt <- list(
        A      = -d$A_ppm[ip],
        x_0    = d$x_0_ppm[ip],
        lambda = -d$lambda_ppm[ip]
    )

    logv("Calculating simulated signal intensities (si_sim) as superposition of lorentz curves")
    rownames(s) <- NULL
    si <- lorentz_sup(x = s$cs, lcp = s$lcpt)

    logv("Adding %s noise to simulated data", noise_method)
    if (noise_method == "RND") {
        # ToSc, 2024-08-02: noise_method 'RND' turned out to be a bad idea,
        # because the true noise is not normally distributed, but has long
        # stretches of continuous increase or decrease that would be incredibly
        # unlikely to occur with a normal distribution. This can be seen by
        # running `analyze_noise_methods()`.
        sfr_sd <- sd(d$y_values_raw[c(1:10000, 121073:131072)] * 1e-6) # (1)
        # (1) The true SFR covers approx. the first and last 20k datapoints.
        # However, to be on the safe side, only use the first and last 10k
        # datapoints for the calculation.
        noise <- rnorm(n = length(si), mean = 0, sd = sfr_sd)
        si <- si + noise
    } else {
        idx <- ix - min(ix) + 5000 # Use SI of datapoints 5000:6308 for noise
        noise <- d$y_values_raw[idx] * 1e-6
        si <- si + noise
    }

    logv("Discretize signal intensities to allow efficient storage as integers")
    s$si <- as.integer(si * 1e6) / 1.e6 # (2)
    # (2) Convert to integers and back. We do this to not lose precision later
    # on when the data gets written to disc as integers.

    if (show) plot_sim_spec(s)
    if (store) save_spectrum(s, verbose = verbose, ...)

    logv("Finished spectrum simulation")
    invisible(s)
}

#' @noRd
#' @title Count Stretches of Increases and Decreases
#' @description Counts the lengths of consecutive increases and decreases in a
#' numeric vector.
#' @param x A numeric vector.
#' @return A numeric vector containing the lengths of stretches of increases and
#' decreases.
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
#' # |   +   +   +   -   -   +   +   | + means "increase"
#' # |-------------------------------| - means "decrease"
#' #
#' x <- c(1.0, 2.2, 3.0, 4.4, 2.0, 1.0, 3.0, 4.0)
#' count_stretches(x) # Returns c(3, 2, 2) because we have + + + - - + +
count_stretches <- function(x) {
    if (length(x) < 2) {
        return(integer(0))
    }
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

#' @noRd
#' @description Used during development of `simulate_spectra()` to find a
#' realistic method for noise generation.
analyze_noise_methods <- function(ask = TRUE) {
    download_example_datasets()
    blood_1 <- datadir("example_datasets/bruker/blood/blood_01")
    deconv <- generate_lorentz_curves(blood_1, ask = FALSE)
    si_raw <- deconv$y_values_raw
    sd_sfr <- sd(si_raw[c(1:10000, 121073:131072)] * 1e-6)
    siRND <- rnorm(n = 10000, mean = 0, sd = sd_sfr)
    siSFR <- si_raw[1:10000] * 1e-6

    logf("Visualizing raw SIs for noise methods RND and SFR")
    plot_noise_methods(siRND, siSFR)
    if (!get_yn_input("Continue?")) {
        return()
    }

    logf("Visualizing smoothed SIs for noise methods RND and SFR")
    siRND_sm <- smooth_signals(list(y_pos = siRND))$y_smooth
    siSFR_sm <- smooth_signals(list(y_pos = siSFR))$y_smooth
    plot_noise_methods(siRND_sm, siSFR_sm)
    if (!get_yn_input("Continue?")) {
        return()
    }

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
