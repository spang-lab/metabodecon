# Private Main #####

#' @noRd
#' @description Reads spectra from the user specified data path.
read_spectra <- function(data_path = file.path(download_example_datasets(), "bruker/urine"),
                         file_format = "bruker",
                         expno = 10,
                         procno = 10,
                         ask = FALSE,
                         sf = c(1e3, 1e6),
                         bwc = FALSE) {
    if (!file_format %in% c("bruker", "jcampdx")) {
        stop("Argument `file_format` should be either 'bruker' or 'jcampdx'")
    }
    dp <- normPath(data_path)
    jcampdx <- file_format == "jcampdx"
    bruker <- file_format == "bruker"
    r1_path <- file.path(dp, expno, "pdata", procno, "1r")
    r1_path_exists <- file.exists(r1_path)
    ends_with_dx <- grepl("\\.dx$", dp)
    if (jcampdx && ends_with_dx || bruker && r1_path_exists) {
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
        stop("No spectra found in directory", data_path)
    }
    spectra <- lapply(paths, function(path) {
        cat3("Reading spectrum", path)
        read_spectrum(path, file_format, sf, expno, procno, bwc)
    })
    names(spectra) <- files
    invisible(spectra)
}

#' @noRd
#' @title Read Spectrum
#' @description Reads a single spectrum file or folder and returns the spectrum data as list.
#' @param path The path of the file/folder containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"` or `"example_datasets/bruker/urine/urine"`.
#' @param type The type of the spectrum file. E.g. `"bruker"` or `"jcampdx"`.
#' @param sf A vector of two elements to scale the x and y axis values.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @param bwc If `TRUE`, the function will return the spectrum data in a backwards compatible format, suitable as input for [generate_lorentz_curves()]. If `FALSE`, the function will return spectrum data in the new format.
#' @return If `bwc` A list containing the spectrum data.
#' @details For details about `procno` and `expno` see section [File Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure) in the metabodecon FAQ.
#' @examples
#' \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1.dx")
#' spectrum_data <- read_spectrum(path, type = "jcampdx")
#' str(spectrum_data, 1)
#' }
read_spectrum <- function(path,
                          type = "bruker",
                          sf = c(1e3, 1e6),
                          expno = 10,
                          procno = 10,
                          bwc = FALSE) {
    type <- match.arg(type, c("bruker", "jcampdx"))
    switch(type,
        "bruker" = read_topspin3_spectrum_glc(path, sf[1], sf[2], expno, procno),
        "jcampdx" = read_jcampdx_spectrum_glc(path, sf[1], sf[2])
    )
}

# Private Helpers #####

#' @noRd
#' @title Read single JCAMPDX Spectrum
#' @description Reads a single JCAMPDX spectrum file and returns the spectrum data in ppm.
#' @param path The path of the file containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"`.
#' @param scale_factor A vector of two elements to scale the x and y axis values. Default is c(1, 1).
#' @return A list containing the spectrum data.
#' \dontrun{
#' # Example Usage (took 30s on the development machine)
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1.dx")
#' spectrum_data <- read_jcampdx_spectrum_glc(path)
#' str(spectrum_data, 1)
#' }
read_jcampdx_spectrum_glc <- function(path, sfx = 1e3, sfy = 1e6) {
    # List of 5
    # $ dataGuide   :'data.frame':   3 obs. of  3 variables:
    # ..$ Format   : chr [1:3] "metadata" "XRR" "XII"
    # ..$ FirstLine: num [1:3] 1 2181 10054
    # ..$ LastLine : num [1:3] 2180 10052 18572
    # $ metadata    : chr [1:2180] "##TITLE=GCKD_U_FA_1152_06_03_2020" ...
    # $ commentLines: int 10053
    # $ real        :'data.frame':   131072 obs. of  2 variables:
    # ..$ x: num [1:131072] 12019 12019 12019 12019 12019 ...
    # ..$ y: num [1:131072] 1265 1003 105 -937 -1062 ...
    # $ imaginary   :'data.frame':   131072 obs. of  2 variables:
    # ..$ x: num [1:131072] 12019 12019 12019 12019 12019 ...
    # ..$ y: num [1:131072] -35831 -36561 -36864 -36185 -34777 ...
    data <- readJDX::readJDX(file = path, SOFC = TRUE, debug = 0) # reading urine_1.dx (~1MB) takes ~30s on machine r31
    real <- data$real
    meta <- data$metadata
    n <- length(real$x) # number of data points (previously called `spectrum_length`)
    ppm_range <- as.numeric(sub("\\D+", "", meta[startsWith(meta, "##$SW=")]))
    ppm_max <- as.numeric(sub("\\D+", "", meta[startsWith(meta, "##$OFFSET=")]))
    ppm_min <- ppm_max - ppm_range
    ppm_step <- ppm_range / (n - 1) # Example: data points in ppm = 1.1, 2.3, 3.5, 4.7 --> ppm_step == 1.2
    ppm_nstep <- ppm_range / n # Don't really know what this is, but it's used in later calculations
    ppm <- seq(ppm_max, ppm_min, -ppm_step) # ppm (previously called `x_ppm`)
    dp <- seq(n - 1, 0, -1) # data points
    sdp <- seq((n - 1) / sfx, 0, -1 / sfx) # scaled data points (previously called `x`). Same as `dp / sfx`, but with slight numeric differences, so we stick with the old calculation method for backwards compatibility.
    return(list(
        y_raw = real$y, y_scaled = real$y / sfy,
        n = n, sfx = sfx, sfy = sfy, # misc
        dp = dp, sdp = sdp, ppm = ppm, # x-axis
        ppm_min = ppm_min, ppm_max = ppm_max, ppm_range = ppm_range, ppm_step = ppm_step, ppm_nstep = ppm_nstep # additional ppm info
        # , length = n, x_ppm = ppm, x = sdp, ppm_highest_value = ppm_max, ppm_lowest_value = ppm_min # backwards compatible names
    ))
}

#' @noRd
#' @title Read single Bruker TopSpin 3 Spectrum
#' @description Reads a single spectrum exported from bruker Topspin v3
#' @param spldir The path of the directory holding the NMR measurements for a individual sample. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @param raw If `TRUE`, the raw signal intensity values are returned. If `FALSE`, the scaled signal intensity values are returned.
#' @param silent Passed on to [read_1r_file()]. If `TRUE`, no output will be printed to the console.
#' @param force Passed on to [read_1r_file()]. If `TRUE`, the function will try to read the file as 64 bit floating point numbers if the processing parameter `DTYPP` is unequal 0. This behaviour is untested, so should be used with caution.
#' @examples
#' spldir <- pkg_file("example_datasets/bruker/blood/blood_01")
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
    # Read relevant files (0.01s)
    acqus <- read_acqus_file(spldir, expno)
    procs <- read_procs_file(spldir, expno, procno)
    one_r <- read_1r_file(spldir, expno, procno, procs, silent = TRUE)[c("raw", "scaled")]

    # Parse files (si: signal intensity, cs: chemical shift, fq: frequency)
    # https://www.bruker.com/protected/en/services/user-manuals/nmr/acquisition-processing.html
    si_raw <- one_r$raw
    si_scaled <- one_r$scaled
    cs_max <- procs$OFFSET # spectrum offset in ppm
    cs_width <- acqus$SW # SW = spectrum width in ppm
    fq_width <- acqus$SW_h # SW_h = spectrum width in Hz
    fq_ref <-  acqus$SFO1 * 1e6 # Reference Frequency in Hz (better than procs$SF, because it fulfills `all.equal` check below)

    # Calculate chemical shift and frequency verctors
    cs_min <- cs_max - cs_width # lowest ppm value
    cs <- seq(cs_max, cs_max - cs_width, length.out = length(si_raw)) # chemical shift in parts per million
    fq_max <- fq_ref - (cs_min * 1e-6 * fq_ref)  # highest frequency in Hz (corresponds to lowest ppm value)
    fq_min <- fq_ref - (cs_max * 1e-6 * fq_ref)  # lowest frequency in Hz
    fq <- seq(fq_min, fq_max, length.out = length(si_raw)) # frequency in Hz

    if (!all.equal(fq_max - fq_min, fq_width)) {
        if (force) {
            if (!silent) catf(sprintf("Calculated spectrum width in Hz (%s) does not match the value from the acqus file (%s). Continuing anyways, because `force` equals `TRUE`. Please note that all downstream calculations using frequencies might be wrong, so be sure to double check the results.", round(fq_width, 5), round(fq_max - fq_min, 5)))
        } else {
            stop(sprintf("Calculated spectrum width in Hz (%s) does not match the value from the acqus file (%s). Please read in the data manually or set `force = TRUE` to ignore this error. Please note that by doing so, all downstream calculations using frequencies might be wrong, so be sure to double check the results.", round(fq_width, 5), round(fq_max - fq_min, 5)))
        }
    }

    # Return spectrum data
    si <- if (raw) si_raw else si_scaled
    data.frame(si, cs, fq)
}

#' @noRd
#' @inherit read_topspin3_spectrum
#' @title Read single Bruker TopSpin 3 Spectrum
#' @description Reads a single spectrum exported from bruker Topspin v3 and returns the spectrum data in a format suitable for [generate_lorentz_curves()].
read_topspin3_spectrum_glc <- function(path = file.path(download_example_datasets(), "bruker/urine/urine_1"),
                                       sfx = 1e3,
                                       sfy = 1e6,
                                       expno = 10,
                                       procno = 10) {
    X <- read_topspin3_spectrum(path, expno, procno, raw = TRUE, silent = TRUE, force = FALSE)
    as_glc_spectrum(X, sfx, sfy)
}

#' @title Takes a spectrum as returned by [read_topspin3_spectrum()] and formats it so that it can be used as input for [generate_lorentz_curves()].
#' @noRd
as_glc_spectrum <- function(X, sfx, sfy) {
    y <- X$si
    n <- length(y)
    ppm_range <- diff(range(X$cs))
    ppm_max <- max(X$cs)
    ppm_min <- min(X$cs)
    ppm_step <- ppm_range / (n - 1)
    ppm_nstep <- ppm_range / n
    # Example: data points in ppm = 1.1, 2.3, 3.5, 4.7
    # ==> ppm_step == 1.2
    # ==> ppm_nstep ~= 1.06 (not really useful, but we need it for backwards compatibility with MetaboDecon1D results)
    ppm <- X$cs # Parts per million
    dp <- seq(n - 1, 0, -1) # Data point numbers
    sdp <- seq((n - 1) / sfx, 0, -1 / sfx) # Scaled data point numbers. (Same as `dp / sfx`, but with slight numeric differences, so we stick with the old calculation method for backwards compatibility)
    named(
        y_raw = y, y_scaled = y / sfy, # y-axis
        n, sfx, sfy, # misc
        dp, sdp, ppm, # x-axis
        ppm_min, ppm_max, ppm_range, ppm_step, ppm_nstep # additional ppm info
    )
}

# File Reading #####

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
#' spldir <- pkg_file("example_datasets/bruker/blood/blood_01")
#' oneR <- read_1r_file(spldir, 10, 10)
#' str(oneR, 1)
read_1r_file <- function(spldir = pkg_file("example_datasets/bruker/blood/blood_01"),
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
    if (!silent) catf(msg)
    con <- file(path_1r, "rb"); on.exit(close(con), add = TRUE)
    raw <- readBin(con, what = type, n = n, size = nbytes, signed = TRUE, endian = endian)
    scaled <- if (type == "integer") raw * 2 ^ procs$NC_proc else raw
    named(
        spldir, expno, procno, path_1r, path_procs, # path related variables
        procs, byteordp, dtypp, endian, nbytes, ncproc, type, n,# processing parameters and derived values
        raw, scaled # raw and scaled signal intensity values
    )
}

#' @noRd
#' @title Read Bruker TopSpin Acquistion Parameters
#' @param spldir The path of the directory holding the NMR measurements for a individual sample. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @return The signals acquisition parameters read from the file as named list.
#' @examples
#' blood1_dir <- pkg_file("example_datasets/bruker/blood/blood_01")
#' acqus <- read_acqus_file(blood1_dir)
#' str(acqus, 0)
#' cat("Spetrum width ppm:", as.numeric(acqus$SW))
#' cat("Spetrum width Hz:", as.numeric(acqus$SW_h))
read_acqus_file <- function(spldir = pkg_file("example_datasets/bruker/blood/blood_01"),
                            expno = 10) {
    path <- file.path(spldir, expno, "acqus")
    acqus <- parse_bruker_param_file(path)
    acqus
}

#' @noRd
#' @title Read Bruker TopSpin Processing Parameters
#' @param spldir The path of the directory holding the NMR measurements for a individual sample. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @return The processing parameters read from the file as named list.
#' @examples
#' blood1_dir <- pkg_file("example_datasets/bruker/blood/blood_01")
#' procs <- read_procs_file(blood1_dir)
read_procs_file <- function(spldir = pkg_file("example_datasets/bruker/blood/blood_01"),
                            expno = 10,
                            procno = 10) {
    path <- file.path(spldir, expno, "pdata", procno, "procs")
    procs <- parse_bruker_param_file(path)
    procs
}

#' @title Parse Bruker Parameter File
#' @description Parses a Bruker parameter file and returns the parameters as a named list.
#' @param path The path of the file containing the parameter data. E.g. `"example_datasets/bruker/urine/urine_1/10/acqus"` or `"example_datasets/bruker/urine/urine_1/10/pdata/10/procs"`.
#' @details For a detailed description of the format of burker parameter files, refer to 'Bruker_NMR_Data_Formats.pdf'.
parse_bruker_param_file <- function(path) {
    lines <- readLines(path)
    lines <- lines[!startsWith(lines, "$$")] # strip comments
    content <- paste(lines, collapse = "") # Example: "##TITLE= Parameter file, TopSpin 3.6.2##JCAMPDX= 5.0"
    keyvals <- strsplit(content, "(##\\$?|= )")[[1]] # c("", "TITLE", "Parameter file, TopSpin 3.6.2", "JCAMPDX", "5.0")
    keys <- keyvals[seq(2, length(keyvals), 2)] # c("TITLE", "JCAMPDX")
    vals <- keyvals[seq(3, length(keyvals), 2)] # c("Parameter file, TopSpin 3.6.2", "5.0")
    if (keys[length(keys)] == "END=" && length(keys) > length(vals)) keys <- keys[- length(keys)]
    ret <- structure(as.list(vals), names = keys)
    flt <- is_float_str(vals)
    int <- is_int_str(vals)
    ret[flt] <- as.numeric(vals[flt])
    ret[int] <- as.integer(vals[int])
    ret
}
