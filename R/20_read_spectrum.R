# Public API #####

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
#' @param raw If `TRUE`, scales the returned signal intensities based on information available in the spectrum metadata, in particular `NC_proc`. For details see `processing-reference.pdf`, available at <https://www.bruker.com/en.html> at section 'Services & Support > Documentation & Manuals > Magnetic Resonance > Acquisition & Processing > TopSpin Processing Commands and Parameters' (requires login).
#' @param silent If `TRUE`, no output will be printed to the console.
#' @param force If `TRUE`, try to continue when encountering errors and print info messages instead. To hide these messages as well, set `silent = TRUE`.
#' @return
#' For `read_spectrum`, a dataframe with following columns:
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
read_spectrum <- function(path = pkg_file("example_datasets/bruker/urine/urine_1"),
                          file_format = "bruker",
                          expno = 10,
                          procno = 10,
                          raw = FALSE,
                          silent = TRUE,
                          force = FALSE) {
    file_format <- match.arg(file_format, c("bruker", "jcampdx"))
    X <- switch(file_format,
        "bruker" = read_topspin3_spectrum(path, expno, procno, raw, silent, force),
        "jcampdx" = read_jcampdx_spectrum(path, silent, force)
    )
}

#' @export
#' @rdname read_spectrum
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

# Private Core #####

#' @noRd
#' @title Read single Bruker TopSpin 3 Spectrum
#' @description For params and return value see [read_spectrum()].
#' @examples
#' spldir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' X <- read_topspin3_spectrum(spldir)
#' fq_ref <- X$fq[1] / (1 - (X$cs[1] / 1e6))
#' print(head(X))
#' cat("Frequency of reference in MHz:", fq_ref / 1e6)
read_topspin3_spectrum <- function(spldir = file.path(download_example_datasets(), "bruker/urine/urine_1"),
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

# Make spectra #####

#' @noRd
#' @title Create a spectrum data frame from metadata
#' @description Creates a data frame from the signal intensities and metadata.
#' @param si The signal intensities.
#' @param cs_max The highest chemical shift value in ppm, i.e. the left border of the spectrum.
#' @param cs_width The width of the spectrum in ppm.
#' @param fq_ref The reference frequency in Hz.
#' @param fq_width The width of the spectrum in Hz. Only used to check whether the values calculated from `cs_max`, `cs_width` and `fq_ref` match the provided value. If `NULL`, this check will be skipped.
#' @param force If `TRUE`, the function will not raise an error in case of discrepancies between the calculated and the provided spectrum width in Hz, but will print a info message instead. To hide this message as well, set `silent = TRUE`.
#' @param silent If `TRUE`, no output will be printed to the console.
#' @return A data frame containing the following columns:
#' - `si`: signal intensities in arbitrary units
#' - `cs`: chemical shifts in ppm
#' - `fq`: frequencies in Hz
#' @examples
#' si <- c(1, 1, 3, 7, 8, 3, 8, 5, 2, 1)
#' cs_max <- 14.8
#' cs_width <- 20.0
#' fq_ref <- 600.25 * 1e6
#' fq_width <- 12005
#' df <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width)
#' df2 <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width = 12010, force = FALSE)
make_spectrum <- function(si,
                          cs_max,
                          cs_width,
                          fq_ref,
                          fq_width = NULL,
                          force = FALSE,
                          silent = FALSE) {
    cs_min <- cs_max - cs_width # lowest ppm value
    cs <- seq(cs_max, cs_max - cs_width, length.out = length(si)) # chemical shift in parts per million
    fq_max <- fq_ref - (cs_min * 1e-6 * fq_ref)  # highest frequency in Hz (corresponds to lowest ppm value)
    fq_min <- fq_ref - (cs_max * 1e-6 * fq_ref)  # lowest frequency in Hz
    fq <- seq(fq_min, fq_max, length.out = length(si)) # frequency in Hz
    fq_width_calc <- fq_max - fq_min
    if (!is.null(fq_width) && !isTRUE(all.equal(fq_width_calc, fq_width))) { # check if calculated spectrum width in Hz matches the value from the metadata
        if (!force) {
            stop(logf("Calculated spectrum width in Hz (%s) does not match the provided value (%s). Please read in the data manually or set `force = TRUE` to ignore this error. Please note that by doing so, all downstream calculations using frequencies might be wrong, so be sure to double check the results.", round(fq_width_calc, 5), round(fq_width, 5)))
        } else {
            if (!silent) logf(sprintf("Calculated spectrum width in Hz (%s) does not match the provided value (%s). Continuing anyways, because `force` equals `TRUE`. Please note that all downstream calculations using frequencies might be wrong, so be sure to double check the results.", round(fq_width_calc, 5), round(fq_width, 5)))
        }
    }
    data.frame(si, cs, fq)
}

#' @noRd
#' @title Convert spectrum data to format for generate_lorentz_curves
#' @description Takes a spectrum as returned by [read_topspin3_spectrum()] or [read_jcampdx_spectrum] and formats it so that it can be used as input for [generate_lorentz_curves()].
#' @param X The spectrum data as returned by [read_topspin3_spectrum()] or [read_jcampdx_spectrum].
#' @param sfx The scaling factor for the x-axis.
#' @param sfy The scaling factor for the y-axis.
#' @return A named list containing the following elements:
#' - `y_raw`: The raw signal intensities.
#' - `y_scaled`: The scaled signal intensities.
#' - `n`: The number of data points.
#' - `sfx`: The scaling factor for the x-axis.
#' - `sfy`: The scaling factor for the y-axis.
#' - `dp`: The data point numbers.
#' - `sdp`: The scaled data point numbers.
#' - `ppm`: The chemical shifts in ppm.
#' - `fq`: The frequencies in Hz.
#' - `ppm_min`: The minimum chemical shift in ppm.
#' - `ppm_max`: The maximum chemical shift in ppm.
#' - `ppm_range`: The range of the chemical shifts in ppm.
#' - `ppm_step`: The step size of the chemical shifts in ppm.
#' - `ppm_nstep`: The step size of the chemical shifts in ppm, calculated as `ppm_range / n`.
#' @examples
#' spldir <- pkg_file("example_datasets/bruker/urine/urine_1")
#' X <- read_topspin3_spectrum(spldir)
#' X_glc <- convert_spectrum(X, 1e3, 1e6)
convert_spectrum <- function(X, sfx, sfy) {
    y_raw <- X$si
    y_scaled <- y_raw / sfy
    n <- length(y_raw)
    ppm_range <- diff(range(X$cs))
    ppm_max <- max(X$cs)
    ppm_min <- min(X$cs)
    ppm_step <- ppm_range / (n - 1)
    ppm_nstep <- ppm_range / n
    # Example: data points in ppm = 1.1, 2.3, 3.5, 4.7
    # ==> ppm_step == 1.2
    # ==> ppm_nstep ~= 1.06 (not really useful, but we need it for backwards compatibility with MetaboDecon1D results)
    ppm <- X$cs # Parts per million
    hz <- X$fq # Frequency in Hz
    dp <- seq(n - 1, 0, -1) # Data point numbers
    sdp <- seq((n - 1) / sfx, 0, -1 / sfx) # Scaled data point numbers. (Same as `dp / sfx`, but with slight numeric differences, so we stick with the old calculation method for backwards compatibility)
    named(
        y_raw, y_scaled, # y-axis
        n, sfx, sfy, # misc
        dp, sdp, ppm, hz,# x-axis
        ppm_min, ppm_max, ppm_range, ppm_step, ppm_nstep # additional ppm info
    )
}

# File Reading #####

#' @noRd
#' @title Parse Metadata File
#' @description Parses a metadata file like Bruker's `acqu[s]` or `proc[s]` files and return the metadata as a named list.
#' @param path The path of the file containing the parameter data. E.g. `"example_datasets/bruker/urine/urine_1/10/acqus"` or `"example_datasets/bruker/urine/urine_1/10/pdata/10/procs"`.
#' @details For a detailed description of the format of burker parameter files, refer to 'Bruker_NMR_Data_Formats.pdf'.
#' @return A named list containing the metadata read from the file.
#' @examples
#' path <- pkg_file("example_datasets/bruker/blood/blood_01/10/acqus")
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
