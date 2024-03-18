# Private API #####

#' @title Read Spectra
#' @description Reads multiple spectrum files or folders and returns the spectra data as a list.
#' @param path The path of the folder containing the spectrum data. E.g. `"example_datasets/jcampdx/urine"` or `"example_datasets/bruker/urine"`.
#' @param type The type of the spectrum files. E.g. `"bruker"` or `"jcampdx"`.
#' @param sf A vector of two elements to scale the x and y axis values.
#' @param expno The experiment number for the spectra files. E.g. `"10"`.
#' @param procno The processing number for the spectra. E.g. `"10"`.
#' @return A list containing the spectra data.
#' @details This function is experimental and has only been tested on outputs of Bruker's Topspin Software version 3.1.2. If it does not work for your spectra files, please create an issue at <https://github.com/spang-lab/metabodecon/issues>.
#' @examples
#' \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine")
#' spectra_data <- read_spectra(path, type = "jcampdx")
#' str(spectra_data, 1)
#' }
#' @noRd
read_spectra <- function(path = file.path(download_example_datasets(), "bruker/urine"),
                         type = "bruker",
                         expno = 10,
                         procno = 10,
                         ask = TRUE,
                         sf = c(1e3, 1e6)) {
    if (file_format == "jcampdx") {
        files <- dir(path, pattern = "\\.dx$") # `.dx` files inside `path`
        paths <- dir(path, pattern = "\\.dx$", full.names = TRUE)
    } else if (file_format == "bruker") {
        files <- list.dirs(path, recursive = FALSE, full.names = FALSE) # folders inside `path`
        paths <- list.dirs(path, recursive = FALSE, full.names = TRUE)
    }
    spectra <- lapply(paths, function(path) {
        msgf("Reading spectrum '%s'", path)
        read_spectrum(path, type, sf, expno, procno)
    })
    names(spectra) <- files
    return(spectra)
}

#' @title Read Spectrum
#' @description Reads a single spectrum file or folder and returns the spectrum data as list.
#' @param path The path of the file/folder containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"` or `"example_datasets/bruker/urine/urine"`.
#' @param type The type of the spectrum file. E.g. `"bruker"` or `"jcampdx"`.
#' @param sf A vector of two elements to scale the x and y axis values.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @return A list containing the spectrum data.
#' @details For details about `procno` and `expno` see section [File Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure) in the metabodecon FAQ.
#' @examples
#' \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1.dx")
#' spectrum_data <- load_spectrum(path, type = "jcampdx")
#' str(spectrum_data, 1)
#' }
#' @noRd
read_spectrum <- function(path, type = "bruker", sf = c(1e3, 1e6), expno = 10, procno = 10) {
    sfx <- sf[1]
    sfy <- sf[2]
    switch(type,
        "bruker" = read_bruker_spectrum_v11(path, sfx, sfy, expno, procno),
        "jcampdx" = read_jcampdx_spectrum_v11(path, sfx, sfy)
    )
}

# Helpers of `read_spectra` #####

#' @title Read single JCAMPDX Spectrum
#' @description Reads a single JCAMPDX spectrum file and returns the spectrum data in ppm.
#' @param path The path of the file containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"`.
#' @param scale_factor A vector of two elements to scale the x and y axis values. Default is c(1, 1).
#' @return A list containing the spectrum data.
#' \dontrun{
#' # Example Usage (took 30s on the development machine)
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1.dx")
#' spectrum_data <- load_jcampdx_spectrum_v10(path)
#' str(spectrum_data, 1)
#' }
#' @noRd
read_jcampdx_spectrum_v11 <- function(path, sfx = 1e3, sfy = 1e6) {
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

    # TODO: return as few values as possible, e.g.
    # >>> return(list(ppm = ppm, si = ss))
    # This should be sufficient to calculate everything else with super simple functions like
    # >>> max(ppm)
    # >>> seq(length(ppm) - 1, 0, -1)
    # etc.
}

#' @title Read single Bruker Spectrum
#' @description Reads a single Bruker spectrum and returns the spectrum data in ppm.
#' @param path The path of the directory holding the spectrum data. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param sfx Scaling factor for x axis in datapoints. E.g. `"1000"`. Only relevant for plotting.
#' @param sfy Scaling factor for y axis (signal strength). E.g. `"100000"`. Only relevant for plotting.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @return A list containing the spectrum data.
#' @details For details about `procno` and `expno` see section [File Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure) in the metabodecon FAQ.
#' @examples
#' \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "bruker/urine/urine_1")
#' spectrum_data <- read_bruker_spectrum_v11(path) # takes 3s on machine r31
#' str(spectrum_data, 1)
#' }
#' @noRd
read_bruker_spectrum_v11 <- function(path = file.path(download_example_datasets(), "bruker/urine/urine_1"),
                                     sfx = 1e3, sfy = 1e6,
                                     expno = 10, procno = 10) {
    acqus <- readLines(file.path(path, expno, "acqus"))
    procs <- readLines(file.path(path, expno, "pdata", procno, "procs"))
    ppm_range <- as.numeric(sub("\\D+", "", acqus[startsWith(acqus, "##$SW=")]))
    y <- read_1r_file(path = path, expno = expno, procno = procno, procs = procs) # signal strength (could also be called signal intensity)
    n <- length(y)
    ppm_max <- as.numeric(sub("\\D+", "", procs[startsWith(procs, "##$OFFSET=")]))
    ppm_min <- ppm_max - ppm_range
    ppm_step <- ppm_range / (n - 1) # Example: data points in ppm = 1.1, 2.3, 3.5, 4.7 --> ppm_step == 1.2
    ppm_nstep <- ppm_range / n # Not really useful, but we need it for backwards compatibility with MetaboDecon1D results
    ppm <- seq(ppm_max, ppm_max - ppm_range, by = -ppm_step) # parts per million
    dp <- seq(n - 1, 0, -1) # data points
    sdp <- seq((n - 1) / sfx, 0, -1 / sfx) # scaled data points (previously called `x`). Same as `dp / sfx`, but with slight numeric differences, so we stick with the old calculation method for backwards compatibility.
    return(list(
        y_raw = y, y_scaled = y / sfy, # y-axis
        n = n, sfx = sfx, sfy = sfy, # misc
        dp = dp, sdp = sdp, ppm = ppm, # x-axis
        ppm_min = ppm_min, ppm_max = ppm_max, ppm_range = ppm_range, ppm_step = ppm_step, ppm_nstep = ppm_nstep # additional ppm info
        # , length = n, x_ppm = ppm, x = sdp, ppm_highest_value = ppm_max, ppm_lowest_value = ppm_min # backwards compatible names
    ))
}

#' @title Read Bruker TopSpin spectrum from 1r file
#' @param path The path of the directory holding the spectrum data. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @param procs The content of the `procs` file. E.g. `readLines(file.path(path, expno, "pdata", procno, "procs"))`.
#' @examples \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "bruker/urine/urine_1")
#' procs <- readLines(file.path(path, 10, "pdata", 10, "procs"))
#' y <- read_1r_file(path, 10, 10, procs)
#' }
#' @return The signals read from the file as numeric vector
read_1r_file <- function(path, expno, procno, procs) {
    # TODO: see issue `Check: check use of DTYPP in load_spectrum` in TODOS.R
    n <- as.numeric(sub("\\D+", "", procs[startsWith(procs, "##$SI=")]))
    int_type <- as.numeric(sub("\\D+", "", procs[startsWith(procs, "##$DTYPP=")]))
    int_size <- if (int_type == 0) 4 else 8
    path_1r <- file.path(path, expno, "pdata", procno, "1r")
    spec_stream <- file(path_1r, "rb")
    on.exit(close(spec_stream), add = TRUE)
    y <- readBin(spec_stream, what = "int", size = int_size, n = n, signed = TRUE, endian = "little")
    return(y)
}
