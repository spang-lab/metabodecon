# Exported #####

#' @export
#'
#' @title Download metabodecon Example Datasets
#'
#' @description Downloads example datasets that can be used to test the
#' functionality of the metabodecon package. These datasets are not included in
#' the package by default due to size constraints. The datasets are downloaded
#' as zip file and extracted automatically, unless extraction is disabled by the
#' user.
#'
#' @param dst_dir The destination directory where the downloaded datasets will
#' be stored. If NULL, the function will return the path to the cached zip file.
#'
#' @param extract Logical. If TRUE, the downloaded zip file will be extracted.
#'
#' @param persistent Logical. If TRUE, the downloaded datasets will be cached at
#' [datadir_persistent()] to speed up future calls to
#' `download_example_datasets()`. If FALSE, the datasets will be cached at
#' [datadir_temp()]. If NULL, the function will check both paths for the cached
#' datasets but will return [datadir_temp()] if the cached file does not yet
#' exist.
#'
#' @param overwrite Logical. If TRUE, existing files with the same name in the
#' destination directory will be overwritten.
#'
#' @return The path to the downloaded (and possibly extracted) datasets.
#'
#' @seealso [datadir()]
#'
#' @examples
#' if (interactive()) {
#'      zip <- download_example_datasets(extract = FALSE, persistent = FALSE)
#'      dir <- download_example_datasets(extract = TRUE)
#' }
#'
download_example_datasets <- function(dst_dir = NULL,
                                      extract = TRUE,
                                      persistent = NULL,
                                      overwrite = FALSE,
                                      silent = FALSE) {

    # Example:
    # input   dst_dir    = C:/Users/max/Downloads
    # var     cached_zip = C:/Users/max/.local/share/R/metabodecon/example_datasets.zip
    # var     cached_xds = C:/Users/max/.local/share/R/metabodecon/example_datasets
    # return  dst_zip    = C:/Users/max/Downloads/example_datasets.zip
    # return  dst_xds    = C:/Users/max/Downloads/example_datasets
    cached_zip <- cache_example_datasets(persistent, extract, overwrite, silent)
    cached_xds <- file.path(dirname(cached_zip), "example_datasets")
    if (is.null(dst_dir)) {
        return(if (extract) cached_xds else cached_zip)
    } else {
        dst_zip <- norm_path(file.path(dst_dir, "example_datasets.zip"))
        dst_xds <- norm_path(file.path(dst_dir, "example_datasets"))
        if (overwrite || !file.exists(dst_zip) || isTRUE(file.size(dst_zip) != xds$zip_size)) {
            if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE)
            file.copy(from = cached_zip, to = dst_dir, overwrite = overwrite)
        }
        if (extract && (overwrite || !dir.exists(dst_xds))) extract_example_datasets(dst_zip)
        return(if (extract) dst_xds else dst_zip)
    }
}

#' @export
#'
#' @title Return Path to File or Directory in metabodecon Package
#'
#' @description Recursively searches for files or directories within the
#' 'metabodecon' package that match the given name.
#'
#' @param name The name to search for.
#'
#' @return The file or directory path.
#'
#' @examples
#' # Unambiguous paths
#' metabodecon_file("urine_1")
#' metabodecon_file("urine_1.dx")
#' metabodecon_file("sim/sim_01")
#'
#' # Ambiguous paths (i.e. multiple matches)
#' metabodecon_file("sim")
#' metabodecon_file("urine")
#'
#' # Non-existing paths (i.e. a character vector of length zero gets returned)
#' metabodecon_file("asdfasdf")
metabodecon_file <- function(name = "sim_01") {
    paths <- list.files(
        path = system.file(package = "metabodecon"),
        full.names = TRUE,
        recursive = TRUE,
        include.dirs = TRUE
    )
    paths[grepl(paste0(name, "$"), paths)]
}

#' @export
#'
#' @title Return path to metabodecon's data directory
#'
#' @description Returns the path to the directory where
#' [download_example_datasets()] stores metabodecon's example data sets or any
#' file within that directory. By default this directory is a subdirectory of
#' R's temporary session directory. If `persistent` is set to `TRUE`, the
#' directory equals the data directory returned by [tools::R_user_dir()]
#' instead.
#'
#' @param file Relative path to a file within the data directory.
#'
#' @param warn Print a warning message when the requested path does not yet
#' exist?
#'
#' @param persistent Return the path to the persistent data directory instead of
#' the temporary one?
#'
#' @return Path to the data directory or a file within it.
#'
#' @details
#' The decision to use a temporary data dir as default and a persistent
#' one only optionally was made to conform to CRAN package policies, which state
#' that:
#'
#'     Packages should not write in the user's home filespace (including
#'     clipboards), nor anywhere else on the file system apart from the R
#'     session's temporary directory \[...\] Limited exceptions may be allowed
#'     in interactive sessions if the package obtains confirmation from the
#'     user. For R version 4.0 or later \[...\] packages may store user-specific
#'     data, configuration and cache files in their respective user directories
#'     obtained from [tools::R_user_dir()] \[...\].
#'
#' Source:
#' [cran.r-project.org/web/packages/policies](
#' https://cran.r-project.org/web/packages/policies.html
#' ).
#'
#' @seealso
#' [download_example_datasets()],
#' [datadir_persistent()],
#' [datadir_temp()]
#'
#' @examples
#' # Get temporary datadir and persistent datadir
#' datadir(persistent = FALSE, warn = FALSE)
#' datadir(persistent = TRUE,  warn = FALSE)
#'
#' # Get persistent datadir if existing else temp datadir. Set `warn = TRUE`
#' # to raise a warning if none of the directories exist yet.
#' datadir(warn = FALSE)
#' if (interactive()) datadir()
#'
#' # Get PERSISTENT_DATADIR/bruker if existing else TEMP_DATADIR/bruker
#' datadir(file = "bruker/urine", warn = FALSE)
datadir <- function(file = NULL, warn = TRUE, persistent = NULL) {
    datadir <- datadir_temp <- datadir_temp()
    datadir_persistent <- datadir_persistent()
    zip_peristent <- zip_persistent()
    zip_persistent_has_correct_size <- isTRUE(file.size(zip_peristent) == xds$zip_size) # implies existence
    if (isTRUE(persistent) || (is.null(persistent) && zip_persistent_has_correct_size)) {
        datadir <- datadir_persistent
    }
    file_path <- if (is.null(file)) datadir else file.path(datadir, file)
    if (warn && !dir.exists(file_path)) {
        warning(file_path, " does not exist. Please call `download_example_datasets()` first.", call. = FALSE)
    }
    normalizePath(file_path, "/", mustWork = FALSE)
}

#' @export
#'
#' @title Persistent Data Directory
#'
#' @description Returns the path to the persistent data directory where
#' metabodecon's data sets are stored. This directory equals the data directory
#' returned by [tools::R_user_dir()] plus additional path normalization.
#'
#' @return Path to the persistent data directory.
#'
#' @seealso [datadir()], [datadir_temp()]
#'
#' @examples
#'
#' datadir_persistent()
datadir_persistent <- function() {
    p <- tools::R_user_dir("metabodecon", "data")
    normalizePath(p, "/", mustWork = FALSE)
}

#' @export
#'
#' @title Temporary Data Directory
#'
#' @description Returns the path to the temporary data directory where
#' metabodecon's data sets are stored. This directory equals subdirectory 'data'
#' of metabodecons temporary session directory [tmpdir()] plus additional path
#' normalization.
#'
#' @return Returns the path to the temporary data directory.
#'
#' @seealso [tmpdir()], [datadir()], [datadir_persistent()]
#'
#' @examples
#'
#' datadir_temp()
datadir_temp <- function() {
    p <- file.path(tmpdir(), "data")
    normalizePath(p, "/", mustWork = FALSE)
}

#' @export
#'
#' @title Temporary Session Directory
#'
#' @description
#' Returns the path to metabodecon's temporary session directory. This directory
#' equals subdirectory 'metabodecon' of R's temporary session directory
#' [base::tempdir()] plus additional path normalization.
#'
#' @param subdir Optional subdirectory within the temporary session directory.
#'
#' @param create Whether to create the directory if it does not yet exist.
#'
#' @return
#' Returns the path to the temporary session directory.
#'
#' @seealso
#' [datadir_temp()]
#' [datadir_persistent()]
#'
#' @examples
#' tmpdir()
#' tmpdir("simulate_spectra")
tmpdir <- function(subdir = NULL, create = FALSE) {
    p <- file.path(base::tempdir(), "metabodecon")
    if (isTRUE(subdir)) p <- tempfile("", p)
    if (is_char(subdir)) p <- file.path(p, subdir)
    if (create) base::dir.create(p, recursive = TRUE, showWarnings = FALSE)
    normalizePath(p, "/", mustWork = FALSE)
}

#' @export
#'
#' @title Retrieve directory path of an example dataset
#'
#' @description
#' Returns the path to the directory storing the example files shipped with
#' metabodecon.
#'
#' Deprecated since metabodecon v1.2.0. Please use [datadir()] instead. See
#' examples below for usage.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param dataset_name Either `""`, `"test"`, `"blood"` or `"urine"`.
#'
#' @param warn Whether to print a warning message when the example folders do
#' not yet exist, i.e. [download_example_datasets()] has not been called yet.
#'
#' @return Path to the directory storing the example files.
#'
#' @seealso [download_example_datasets()]
#'
#' @examples
#' x <- get_data_dir("urine")                     # Deprecated
#' y <- datadir("example_datasets/bruker/urine")  # Preferred
#' cat(x, y, sep = "\n")
get_data_dir <- function(dataset_name = c("", "blood", "test", "urine"), warn = TRUE) {
  dataset_name <- match.arg(dataset_name)
  base_data_dir <- datadir()
  data_dir <- file.path(base_data_dir, "example_datasets/bruker", dataset_name, fsep = "/")
  if (warn && !dir.exists(data_dir)) {
    warning(data_dir, " does not exist. Please call `download_example_datasets(extract = TRUE)` first.", call. = FALSE)
  }
  return(data_dir)
}

# Helpers (Private) #####

#' @noRd
#'
#' @title Example Datasets Information
#'
#' @description
#' This list contains information about the example datasets provided for users
#' to try the package.
#'
#' The datasets can be downloaded from the provided URL.
#'
#' @field url The URL from where the example datasets can be downloaded.
#'
#' @field zip_size The expected size of the zipped example datasets file in
#' bytes.
#'
#' @field dir_size The expected size of the directory containing the example
#' datasets in bytes.
#'
#' @field n_files The expected total number of files in the example_datasets
#' folder after extraction.
#'
#' @examples
#' str(xds)
#'
xds <- list(
    url = "https://github.com/spang-lab/metabodecon/releases/download/v1.1.0/example_datasets.zip",
    zip_size = 38425397,
    dir_size = 56378684,
    n_files = 1018
)

#' @noRd
#'
#' @title Cache Example Datasets
#'
#' @description
#' If file `"example_datasets.zip"` does not yet exist at [datadir()], it is
#' downloaded.
#'
#' If the zip exists but has a wrong size, it can be overwritten based on the
#' `overwrite` parameter. If `extract` is TRUE and either folder
#' `"example_datasets"` does not yet exist at [datadir()] or `overwrite` is
#' TRUE, the zip file is extracted. The path to the zip file is returned.
#'
#' @param persistent
#' If TRUE, the datasets are cached permanently.
#' If FALSE, the datasets are cached temporarily.
#' If NULL, the path to the permanent zip file is returned if it exists and has
#' the correct size else the path to the temporary zip file is returned.
#'
#' @param extract If TRUE, the datasets are extracted after being cached.
#'
#' @param overwrite
#' First effect: If TRUE and the cached file exists but has a wrong size, the
#' file is overwritten. If FALSE, an error is thrown in this case.
#'
#' Second effect (only relevant if `extract` is TRUE): If TRUE, the datasets are
#' extracted even if the corresponding folder already exists. If FALSE, files
#' are not extracted again.
#'
#' @return The path to the cached datasets.
#'
#' @examples
#' \donttest{
#' cache_example_datasets(persistent = FALSE, extract = FALSE, overwrite = FALSE)
#' }
cache_example_datasets <- function(persistent = NULL,
                                   extract = FALSE,
                                   overwrite = FALSE,
                                   silent = FALSE) {

    pdir <- datadir_persistent()
    pzip <- file.path(pdir, "example_datasets.zip")
    pzip_exists <- file.exists(pzip)
    pzip_is_missing <- !pzip_exists
    pzip_has_correct_size <- isTRUE(file.size(pzip) == xds$zip_size) # Implies existence
    pzip_has_wrong_size <- !pzip_has_correct_size

    tdir <- datadir_temp()
    tzip <- file.path(tdir, "example_datasets.zip")
    tzip_has_correct_size <- isTRUE(file.size(tzip) == xds$zip_size)
    tzip_has_wrong_size <- !tzip_has_correct_size

    persistent <- if (isTRUE(persistent) || (is.null(persistent) && pzip_has_correct_size)) TRUE else FALSE

    if (isTRUE(persistent)) {
        if (pzip_is_missing || (pzip_has_wrong_size && overwrite)) {
            download_example_datasets_zip(pzip, copyfrom = tzip, silent = silent)
        }
        if (pzip_exists && pzip_has_wrong_size && !overwrite) {
            msg <- "Path %s exists already, but has size %d instead of %d. Set `overwrite = TRUE` to overwrite."
            stop(sprintf(msg, pzip, file.size(pzip), xds$zip_size))
        }
    } else {
        if (tzip_has_wrong_size) {
            download_example_datasets_zip(tzip, copyfrom = pzip, silent = silent)
            # No need to ask for permission here because we only write to the temp dir
        }
    }

    zip <- if (persistent) pzip else tzip
    dstdir <- file.path(dirname(zip), "example_datasets")
    if (extract) {
        if (dir.exists(dstdir) && overwrite) unlink(dstdir, recursive = TRUE)
        if (!dir.exists(dstdir)) extract_example_datasets(zip)
    }
    zip
}

extract_example_datasets <- function(path = datadir("example_datasets.zip")) {
    utils::unzip(zipfile = path, exdir = dirname(path))
}

#' @noRd
#'
#' @title Download Example Datasets Zip
#'
#' @description This function downloads the example_datasets.zip file from the
#' specified URL and saves it to the destination directory provided. If the
#' directory does not exist, it will be created. If the file already exists in
#' the directory, it will be overwritten without asking for permission.
#'
#' @param path A string. The path where the downloaded file will be saved.
#'
#' @return The path where the downloaded file is saved.
#'
#' @examples
#'
#' \donttest{
#' download_example_datasets_zip("/path/to/your/directory/example_datasets.zip")
#' }
download_example_datasets_zip <- function(path, copyfrom = NULL, silent = FALSE) {
    mkdirs(dirname(path))
    if (!is.null(copyfrom) && file.exists(copyfrom) && file.size(copyfrom) == xds$zip_size) {
        if (!silent) message(sprintf("Copying cached archive %s to %s", copyfrom, path))
        file.copy(copyfrom, path, overwrite = TRUE)
    } else {
        if (!silent) message(sprintf("Downloading %s as %s", xds$url, path))
        utils::download.file(xds$url, path, quiet = TRUE)
    }
    path
}

zip_temp <- function() {
    p <- file.path(datadir_temp(), "example_datasets.zip")
    normalizePath(p, "/", mustWork = FALSE)
}

zip_persistent <- function() {
    p <- file.path(datadir_persistent(), "example_datasets.zip")
    normalizePath(p, "/", mustWork = FALSE)
}

tmpfile <- function(pattern = "file", fileext = "") {
    tempfile(pattern, tmpdir(create = TRUE), fileext)
}

testdir <- function(p = NULL) {
    norm_path(paste(tmpdir("tests"), p, sep = "/"))
    # use paste instead of file.path, because it can deal with NULL
}

mockdir <- function() {
    norm_path(file.path(tmpdir(), "mocks"))
}

#' @noRd
#' @description
#' Create and return cache dir. If existing, the persistent cache dir is
#' returned, else the temp cache dir. To force creation of the persistent cache
#' dir, call once with `persistent=TRUE`.
cachedir <- function(persistent = NULL) {
    tcd <- file.path(tmpdir(), "cache") # temporary cache dir
    pcd <- file.path(tools::R_user_dir("metabodecon", "cache")) # persistent cache dir
    cd <- if (isTRUE(persistent) || (is.null(persistent) && dir.exists(pcd))) pcd else tcd
    ncd <- normalizePath(cd, "/", mustWork = FALSE)
    mkdirs(ncd)
}

# Sap (Private) #####

`sap Documentation` <- NULL

#' @name sap
#'
#' @title Simple As Possible (SAP) Spectra
#'
#' @description
#' Create, update or get "Simple as Possible" (SAP) spectrum objects.
#'
#' Use
#' - `make_sap*()` to create the sap objects from scratch
#' - `update_sap*()` to save it to the package and
#' - `get_sap*()` to load it from the package.
#'
#' @details
#' SAP spectra should be so simple, that the expected deconvolution results can
#' be either:
#'
#' 1. Calculated by hand, i.e. without calling a deconvolution function, or
#' 2. Checked for correctness by inspecting the deconvolution object. Inspection
#'    can be done by looking at the object structure or by plotting the object.
#'
#' @section SAP1
#' The SAP1 spectrum contains 11 datapoints ranging from -0.5 to 0.5 ppm with one
#' peak at 0.2 ppm, i.e. at the 4th datapoint. At the far right a huge noise
#' peak is added to achieve a non-zero MSE. The spectrum has no water signal and
#' no signal free region. It's most important properties are shown below, with
#' the frequency (fq) is given as difference to the reference frequency (6e8).
#' I.e. fq* = +60 means fq = 600000060 Hz.
#'
#' ```
#'                  left  centr  rigt
#' si   1660  3448  9756  25000  9756  3448  1660   962   624   437   322
#' cs    0.5   0.4   0.3    0.2   0.1     0  -0.1  -0.2  -0.3  -0.4  -0.5
#' dpi     1     2     3      4     5     6     7     8     9    10    11
#' dpn    10     9     8      7     6     5     4     3     2     1     0
#' sdp  .010  .009  .008   .007  .006  .005  .004  .003  .002  .001  .000
#' fq*  -300  -240  -180   -120   -60     0   +60  +120  +180  +240  +300
#' ```
#'
#' On a Windows Machine with Samsung 980 Pro SSD and AMD Ryzen 7 5700U CPU the
#' following runtimes were observed:
#'
#' 1. Creation from scratch via [make_sap1()]: approx. 0.11 seconds.
#' 2. Reading from disk via [get_sap1()]: approx. 0.01 seconds
#'
#' @section SAP2
#' The SAP2 spectrum contains 128 datapoints ranging from -6.4 to 6.4 ppm with
#' four peaks. A rough sketch of the spectrum is shown below:
#'
#' ```
#' -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' |      SFR      |               w               |     SFR      |
#' |               |  x           www       p      |              |
#  |               | xxxa        wwwww     ppp     |              |
#  |               |xxxaaa      wwwwwww  ppppppp   |              |
#' |~-~-~-~-~-~-~-~|~-|-|-~-~-~-~-~|~-~-~-~-|-~-~-~-~-~-~-~-~-~-~-~
#' |               |  | |          |        |      |
#' 6.4             |  | 2.24       0.047    -2.22  -3.2
#'                 |  2.61
#'                 3.2
#' ```
NULL

make_sap1 <- function() {
    # Meta Data
    name   <- "sap1"
    sf     <- c(1e3, 1e6)
    fq_ref <- 6e8
    # Lorentz params
    x0_ppm     <- 0.2
    A_raw_ppm  <- 2000
    lambda_ppm <- 0.08
    # Axis values
    cs             <- (5:-5) / 10
    noise          <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10000)
    si_sup         <- lorentz_sup(cs, x0_ppm, A_raw_ppm, lambda_ppm)
    si_sup_norm    <- si_sup / sum(si_sup)
    si_sup_sc      <- si_sup / sf[2]
    si_sup_sc_norm <- si_sup_sc / sum(si_sup_sc)
    si             <- si_sup + noise
    si_norm        <- si / sum(si)
    si_sc          <- si / sf[2]
    si_sc_norm     <- si_sc / sum(si_sc)
    n              <- length(si)
    dpi            <- 1:n
    dpn            <- (n-1):0
    sdp            <- dpn / sf[1]
    fq             <- fq_ref - (fq_ref / 1e6) * cs
    cs_step        <- 0.1
    dpn_step       <- 1
    fq_step        <- cs_step * (fq_ref / 1e6)
    sdp_step       <- dpn_step / sf[1]
    # Peak values
    peak_idx <- data.frame(left = 3,     center = 4,     right = 5)
    peak_ppm <- data.frame(left = 0.3,   center = 0.2,   right = 0.1)
    peak_sdp <- data.frame(left = 0.008, center = 0.007, right = 0.006)
    # Lorentz params converted
    x0_dp  <- convert_pos(x0_ppm, cs, dpn) # x0_dp  = 7
    x0_sdp <- convert_pos(x0_ppm, cs, sdp) # x0_sdp = 0.007
    x0_hz  <- convert_pos(x0_ppm, cs, fq)  # x0_hz  = 599999880
    A_raw_dp  <- A_raw_ppm * (dpn_step  / cs_step) # A_raw_dp  = 20000
    A_raw_sdp <- A_raw_ppm * (sdp_step / cs_step)  # A_raw_sdp = 20
    A_raw_hz  <- A_raw_ppm * (fq_step  / cs_step)  # A_raw_hz  = 1200000
    A_sc_ppm  <- A_raw_ppm / sf[2] # A_sc_ppm  = 0.002
    A_sc_dp   <- A_raw_dp  / sf[2] # A_sc_dp   = 0.02
    A_sc_sdp  <- A_raw_sdp / sf[2] # A_sc_sdp  = 0.00002
    A_sc_hz   <- A_raw_hz  / sf[2] # A_sc_hz   = 1.2
    lambda_dp  <- convert_width(lambda_ppm, cs, dpn)     # lambda_dp  = 0.8
    lambda_sdp <- convert_width(lambda_ppm, cs, sdp)     # lambda_sdp = 0.0008
    lambda_hz  <- abs(convert_width(lambda_ppm, cs, fq)) # lambda_hz  = 48
    # Deconv params
    nfit <- 3
    smopts <- c(0, 3)
    delta <- 6.4
    wshw <- 0
    sfr_ppm <- c(1, -1) # outside of spectrum range
    sfr_dp  <- convert_pos(sfr_ppm, cs, dpn) # c(15, -5)
    sfr_sdp <- convert_pos(sfr_ppm, cs, sdp) # c(0.015, -0.005)
    sfr_sdp_0 <- sfr_in_sdp_bwc(sfr_ppm, cs, sf) # c(15, -5)
    sfr_hz  <- convert_pos(sfr_ppm, cs, fq)  # c(599999400, 600000600)
    args <- named(nfit, smopts, delta, sfr = sfr_ppm, wshw)
    # Integrals
    # (use n-step intervals instead of n-1-steps to reproduce (gy) results from v0.x)
    cs_range_0  <- c(min(cs ), max(cs ) + cs_step ) # cs_range_0  = c(-0.5, 0.6)
    dpn_range_0 <- c(min(dpn), max(dpn) + dpn_step) # dpn_range_0 = c(0, 11)
    sdp_range_0 <- c(min(sdp), max(sdp) + sdp_step) # sdp_range_0 = c(0, 0.011)
    fq_range_0  <- c(min(fq ), max(fq ) + fq_step ) # fq_range_0  = c(599999700, 600000360)
    integrals_rawsi_ppm   <- lorentz_int(x0_ppm, A_raw_ppm, lambda_ppm, limits = range(cs ))   # integrals_rawsi_ppm   = 5534.39650939749
    integrals_rawsi_dp    <- lorentz_int(x0_dp,  A_raw_dp,  lambda_dp,  limits = range(dpn))   # integrals_rawsi_dp    = 55343.9650939749
    integrals_rawsi_sdp   <- lorentz_int(x0_sdp, A_raw_sdp, lambda_sdp, limits = range(sdp))   # integrals_rawsi_sdp   = 55.3439650939749
    integrals_rawsi_hz    <- lorentz_int(x0_hz,  A_raw_hz,  lambda_hz,  limits = range(fq ))   # integrals_rawsi_hz    = 3320637.90563849
    integrals_rawsi_ppm_0 <- lorentz_int(x0_ppm, A_raw_ppm, lambda_ppm, limits = cs_range_0)   # integrals_rawsi_ppm_0 = 5660.81017319241
    integrals_rawsi_dp_0  <- lorentz_int(x0_dp,  A_raw_dp,  lambda_dp,  limits = dpn_range_0)  # integrals_rawsi_dp_0  = 56608.1017319241
    integrals_rawsi_sdp_0 <- lorentz_int(x0_sdp, A_raw_sdp, lambda_sdp, limits = sdp_range_0)  # integrals_rawsi_sdp_0 = 56.6081017319241
    integrals_rawsi_hz_0  <- lorentz_int(x0_hz,  A_raw_hz,  lambda_hz,  limits = fq_range_0)   # integrals_rawsi_hz_0  = 3337585.93122155
    integrals_scsi_ppm   <- integrals_rawsi_ppm   / sf[2] # integrals_sc_ppm   = 0.00553439650939
    integrals_scsi_dp    <- integrals_rawsi_dp    / sf[2] # integrals_sc_dp    = 0.05534396509397
    integrals_scsi_sdp   <- integrals_rawsi_sdp   / sf[2] # integrals_sc_sdp   = 0.00005534396509
    integrals_scsi_hz    <- integrals_rawsi_hz    / sf[2] # integrals_sc_hz    = 3.32063790563849
    integrals_scsi_ppm_0 <- integrals_rawsi_ppm_0 / sf[2] # integrals_sc_ppm_0 = 0.00566081017319
    integrals_scsi_dp_0  <- integrals_rawsi_dp_0  / sf[2] # integrals_sc_dp_0  = 0.05660810173192
    integrals_scsi_sdp_0 <- integrals_rawsi_sdp_0 / sf[2] # integrals_sc_sdp_0 = 0.00005660810173
    integrals_scsi_hz_0  <- integrals_rawsi_hz_0  / sf[2] # integrals_sc_hz_0  = 3.33758593122155
    # MSEs
    vae_raw_si  <- abs(si - si_sup)           # Vector of absolute errors for raw SIs
    vae_norm_si <- abs(si_norm - si_sup_norm) # Vector of absolute errors for normalized SIs
    vse_raw_si  <- vae_raw_si^2               # Vector of squared errors for raw SIs
    vse_norm_si <- vae_norm_si^2              # Vector of squared errors for normalized SIs
    mae_raw_si  <- mean(vae_raw_si)           # Mean absolute error for raw SIs
    mae_norm_si <- mean(vae_norm_si)          # Mean absolute error for normalized SIs
    mse_raw_si  <- mean(vse_raw_si)           # Mean squared  error for raw SIs
    mse_norm_si <- mean(vse_norm_si)          # Mean squared  error for normalized SIs
    # Compound Objects
    simpar <- named(name = "sap1", cs, pkr = NULL, x0 = x0_ppm, A = A_raw_ppm, lambda = lambda_ppm, noise)
    lcpar  <- named(A = A_sc_ppm, lambda = lambda_ppm, x0 = x0_ppm)
    args   <- named(nfit, smopts, delta, sfr = sfr_ppm, wshw)
    meta   <- named(name, fq, simpar)
    sit    <- named(wsrm = si, nvrm = si, sm = si, sup = si_sup, al = NULL)
    mse    <- named(raw = mse_raw_si, norm = mse_norm_si, sm = mse_raw_si, smnorm = mse_norm_si)
    peak   <- peak_idx
    # Return Objects
    spectrum <- structure(class = "spectrum", named(si, cs, meta))
    decon0 <- list(
        number_of_files            = 1L,
        filename                   = "sap1",
        x_values                   = sdp,
        x_values_ppm               = cs,
        y_values                   = si_sc,
        spectrum_superposition     = structure(si_sup_sc, dim = c(1L, 11L)),
        mse_normed                 = mse_norm_si,
        index_peak_triplets_middle = peak_idx[1, 2],
        index_peak_triplets_left   = peak_idx[1, 3],  # (1)
        index_peak_triplets_right  = peak_idx[1, 1],  # (1)
        peak_triplets_middle       = peak_ppm[1, 2],
        peak_triplets_left         = peak_ppm[1, 3],  # (1)
        peak_triplets_right        = peak_ppm[1, 1],  # (1)
        integrals                  = structure(integrals_scsi_sdp_0, dim = c(1, 1)), # (2)
        signal_free_region         = sfr_sdp_0,    # (3)
        range_water_signal_ppm     = 0,
        A                          = - A_sc_sdp,   # (4)
        lambda                     = - lambda_sdp, # (4)
        x_0                        = x0_sdp
    )
    decon1 <- structure(class = "decon1", c(decon0, list(
        y_values_raw           = si,
        x_values_hz            = fq,
        mse_normed_raw         = mse_norm_si,
        signal_free_region_ppm = sfr_ppm,
        x_0_hz                 = x0_hz,
        x_0_dp                 = x0_dp,
        x_0_ppm                = x0_ppm,
        A_hz                   = A_sc_hz,
        A_dp                   = - A_sc_dp,    # (4)
        A_ppm                  = - A_sc_ppm,   # (4)
        lambda_hz              = lambda_hz,
        lambda_dp              = - lambda_dp,  # (4)
        lambda_ppm             = - lambda_ppm  # (4)
    )))
    decon2 <- structure(class = "decon2", named(cs, si, meta, args, sit, peak, lcpar, mse))
    ispec <- structure(class = "ispec", list(
        y_raw     = si,
        y_scaled  = si_sc,
        n         = n,
        dpn       = dpn,
        sdp       = sdp,
        sf        = sf,
        ppm       = cs,
        hz        = fq,
        ppm_range = diff(range(cs)),
        ppm_max   = max(cs),
        ppm_min   = min(cs),
        ppm_step  = cs_step,
        ppm_nstep = diff(range(cs)) / length(cs),
        name      = "sap1",
        meta      = meta
    ))
    idecon <- structure(class = "idecon", c(ispec, list(
        TODO     = "TODO: Finish idecon implementation",
        sf       = sf,
        y_nows   = si,
        y_pos    = si,
        Z        = NULL,
        y_smooth = si,
        d        = calc_second_derivative(si),
        peak     = NULL,
        lci      = NULL,
        lca      = NULL,
        lcr      = NULL
    )))
    # (1): right and left are switched in decon0/decon1
    # (2): integrals are calculated over [0, n] instead of [0, n-1] in decon0/decon1
    # (3): conversion of ppm to sdp is slightly off in decon0/decon1
    # (4): x_0 and lambda are returned as negative values in decon0/decon1
    named(spectrum, decon0, decon1, decon2, ispec, idecon)
}

make_sap2 <- function() {
    sap2 <- simulate_spectrum(
        name   = "sap2",
        cs     = (64:-63) / 10,
        x0     = c(2.81, 2.10, 0.047, -2.22),
        A      = c(2000,  500,  5000,   1500),
        lambda = c( 0.4,  0.3,   0.6,   0.5),
        noise  = round(rnorm(128, sd = 50)),
        fqref  = 6e8,
        pkr    = c(3.2, -3.2)
    )
}

update_sap1 <- function(overwrite = FALSE, verbose = TRUE) {
    sap1 <- make_sap1()
    if (isTRUE(overwrite)) {
        xds_path  <- pkg_file("example_datasets")
        rds_path  <- file.path(xds_path, "rds") # might not exist
        sap1_path <- file.path(rds_path, "sap1.rds")
        if (!dir.exists(rds_path)) dir.create(rds_path)
        if (verbose) logf("Overwrite is TRUE. Updating %s." , sap1_path)
        saveRDS(sap1, sap1_path)
    } else {
        sap1_path <- file.path(tmpdir(subdir = TRUE, create = TRUE), "sap1.rds")
        if (verbose) logf("Overwrite is FALSE. Writing spectra to %s." , sap1_path)
        saveRDS(sap1, sap1_path)
    }
    sap1_path
}

update_sap2 <- function(overwrite = TRUE, verbose = TRUE) {
    sap2 <- make_sap2()
    if (isTRUE(overwrite)) {
        path <- file.path(pkg_file("example_datasets/bruker"), "sap/sap2")
        if (verbose) logf("Overwrite is TRUE. Updating %s." , path)
        save_spectrum(sap2, path, force = TRUE, verbose = verbose)
        usethis::use_data(sap2, overwrite = TRUE)
    } else {
        path <- file.path(tmpdir(subdir = TRUE, create = TRUE), "sap/sap2")
        if (verbose) logf("Overwrite is FALSE. Writing spectra to %s." , path)
        save_spectrum(sap2, path, force = FALSE, verbose = verbose)
        save(sap2, file = file.path(path, "sap2.rda"), compress = "bzip2")
    }
    path
}

get_sap1 <- function() {
    xds_path  <- pkg_file("example_datasets")
    rds_path  <- file.path(xds_path, "rds")
    sap1_path <- file.path(rds_path, "sap1.rds")
    sap1 <- readRDS(sap1_path)
    sap1
}

get_sap2_ispec <- function() {
    as_ispec(make_sap2())
}

get_sap2_idecon <- function(rmws = FALSE, bwc = 2) {
    if (rmws %!=% FALSE|| bwc %!=% 2) {
        call <- deparse(match.call())
        stopf("%s is not implemented yet. Use rmws=FALSE and bwc=2.", call, call. = FALSE)
    }
    sap_decon <- deconvolute_ispec(
        ispec = as_ispec(make_sap2()),
        sfr = c(3.2, -3.2),
        smopts = c(0, 3),
        delta = 1,
        verbose = FALSE
    )
}

get_optimal_values <- function(spectrum) {
    # Meta Data
    name   <- spectrum$sap1
    sf     <- c(1e3, 1e6)
    fq_ref <- convert_pos(0, spectrum$cs, spectrum$meta$fq)
    # Lorentz params
    x0_ppm     <- spectrum$meta$simpar$x0
    A_raw_ppm  <- spectrum$meta$simpar$A
    lambda_ppm <- spectrum$meta$simpar$lambda
    # Axis values
    cs             <- spectrum$cs
    noise          <- spectrum$meta$simpar$noise
    si_sup         <- lorentz_sup(cs, x0_ppm, A_raw_ppm, lambda_ppm)
    si_sup_norm    <- si_sup / sum(si_sup)
    si_sup_sc      <- si_sup / sf[2]
    si_sup_sc_norm <- si_sup_sc / sum(si_sup_sc)
    si             <- si_sup + noise
    si_norm        <- si / sum(si)
    si_sc          <- si / sf[2]
    si_sc_norm     <- si_sc / sum(si_sc)
    n              <- length(si)
    dpi            <- 1:n
    dpn            <- (n-1):0
    sdp            <- dpn / sf[1]
    fq             <- fq_ref - (fq_ref / 1e6) * cs
    cs_step        <- 0.1
    dpn_step       <- 1
    fq_step        <- cs_step * (fq_ref / 1e6)
    sdp_step       <- dpn_step / sf[1]
    # Peak values
    peak_ppm <- data.frame(
        left = x0_ppm + lambda_ppm,
        center = x0_ppm,
        right = x0_ppm - lambda_ppm
    )
    peak_idx <- data.frame(
        left = round(convert_pos(peak_ppm$left, cs, dpi)),
        center = round(convert_pos(peak_ppm$center, cs, dpi)),
        right = round(convert_pos(peak_ppm$right, cs, dpi))
    )
    peak_sdp <- data.frame(
        left = convert_pos(peak_ppm$left, cs, sdp),
        center = convert_pos(peak_ppm$center, cs, sdp),
        right = convert_pos(peak_ppm$right, cs, sdp)
    )
    # Lorentz params converted
    x0_dp  <- convert_pos(x0_ppm, cs, dpn) # x0_dp  = 7
    x0_sdp <- convert_pos(x0_ppm, cs, sdp) # x0_sdp = 0.007
    x0_hz  <- convert_pos(x0_ppm, cs, fq)  # x0_hz  = 599999880
    A_raw_dp  <- A_raw_ppm * (dpn_step  / cs_step) # A_raw_dp  = 20000
    A_raw_sdp <- A_raw_ppm * (sdp_step / cs_step)  # A_raw_sdp = 20
    A_raw_hz  <- A_raw_ppm * (fq_step  / cs_step)  # A_raw_hz  = 1200000
    A_sc_ppm  <- A_raw_ppm / sf[2] # A_sc_ppm  = 0.002
    A_sc_dp   <- A_raw_dp  / sf[2] # A_sc_dp   = 0.02
    A_sc_sdp  <- A_raw_sdp / sf[2] # A_sc_sdp  = 0.00002
    A_sc_hz   <- A_raw_hz  / sf[2] # A_sc_hz   = 1.2
    lambda_dp  <- convert_width(lambda_ppm, cs, dpn)     # lambda_dp  = 0.8
    lambda_sdp <- convert_width(lambda_ppm, cs, sdp)     # lambda_sdp = 0.0008
    lambda_hz  <- abs(convert_width(lambda_ppm, cs, fq)) # lambda_hz  = 48
    # Deconv params
    nfit <- 3
    smopts <- c(0, 3)
    delta <- 6.4
    wshw <- 0
    sfr_ppm <- c(1, -1) # outside of spectrum range
    sfr_dp  <- convert_pos(sfr_ppm, cs, dpn) # c(15, -5)
    sfr_sdp <- convert_pos(sfr_ppm, cs, sdp) # c(0.015, -0.005)
    sfr_sdp_0 <- sfr_in_sdp_bwc(sfr_ppm, cs, sf) # c(15, -5)
    sfr_hz  <- convert_pos(sfr_ppm, cs, fq)  # c(599999400, 600000600)
    args <- named(nfit, smopts, delta, sfr = sfr_ppm, wshw)
    # Integrals
    # (use n-step intervals instead of n-1-steps to reproduce (gy) results from v0.x)
    cs_range_0  <- c(min(cs ), max(cs ) + cs_step ) # cs_range_0  = c(-0.5, 0.6)
    dpn_range_0 <- c(min(dpn), max(dpn) + dpn_step) # dpn_range_0 = c(0, 11)
    sdp_range_0 <- c(min(sdp), max(sdp) + sdp_step) # sdp_range_0 = c(0, 0.011)
    fq_range_0  <- c(min(fq ), max(fq ) + fq_step ) # fq_range_0  = c(599999700, 600000360)
    integrals_rawsi_ppm   <- lorentz_int(x0_ppm, A_raw_ppm, lambda_ppm, limits = range(cs ))   # integrals_rawsi_ppm   = 5534.39650939749
    integrals_rawsi_dp    <- lorentz_int(x0_dp,  A_raw_dp,  lambda_dp,  limits = range(dpn))   # integrals_rawsi_dp    = 55343.9650939749
    integrals_rawsi_sdp   <- lorentz_int(x0_sdp, A_raw_sdp, lambda_sdp, limits = range(sdp))   # integrals_rawsi_sdp   = 55.3439650939749
    integrals_rawsi_hz    <- lorentz_int(x0_hz,  A_raw_hz,  lambda_hz,  limits = range(fq ))   # integrals_rawsi_hz    = 3320637.90563849
    integrals_rawsi_ppm_0 <- lorentz_int(x0_ppm, A_raw_ppm, lambda_ppm, limits = cs_range_0)   # integrals_rawsi_ppm_0 = 5660.81017319241
    integrals_rawsi_dp_0  <- lorentz_int(x0_dp,  A_raw_dp,  lambda_dp,  limits = dpn_range_0)  # integrals_rawsi_dp_0  = 56608.1017319241
    integrals_rawsi_sdp_0 <- lorentz_int(x0_sdp, A_raw_sdp, lambda_sdp, limits = sdp_range_0)  # integrals_rawsi_sdp_0 = 56.6081017319241
    integrals_rawsi_hz_0  <- lorentz_int(x0_hz,  A_raw_hz,  lambda_hz,  limits = fq_range_0)   # integrals_rawsi_hz_0  = 3337585.93122155
    integrals_scsi_ppm   <- integrals_rawsi_ppm   / sf[2] # integrals_sc_ppm   = 0.00553439650939
    integrals_scsi_dp    <- integrals_rawsi_dp    / sf[2] # integrals_sc_dp    = 0.05534396509397
    integrals_scsi_sdp   <- integrals_rawsi_sdp   / sf[2] # integrals_sc_sdp   = 0.00005534396509
    integrals_scsi_hz    <- integrals_rawsi_hz    / sf[2] # integrals_sc_hz    = 3.32063790563849
    integrals_scsi_ppm_0 <- integrals_rawsi_ppm_0 / sf[2] # integrals_sc_ppm_0 = 0.00566081017319
    integrals_scsi_dp_0  <- integrals_rawsi_dp_0  / sf[2] # integrals_sc_dp_0  = 0.05660810173192
    integrals_scsi_sdp_0 <- integrals_rawsi_sdp_0 / sf[2] # integrals_sc_sdp_0 = 0.00005660810173
    integrals_scsi_hz_0  <- integrals_rawsi_hz_0  / sf[2] # integrals_sc_hz_0  = 3.33758593122155
    # MSEs
    vae_raw_si  <- abs(si - si_sup)           # Vector of absolute errors for raw SIs
    vae_norm_si <- abs(si_norm - si_sup_norm) # Vector of absolute errors for normalized SIs
    vse_raw_si  <- vae_raw_si^2               # Vector of squared errors for raw SIs
    vse_norm_si <- vae_norm_si^2              # Vector of squared errors for normalized SIs
    mae_raw_si  <- mean(vae_raw_si)           # Mean absolute error for raw SIs
    mae_norm_si <- mean(vae_norm_si)          # Mean absolute error for normalized SIs
    mse_raw_si  <- mean(vse_raw_si)           # Mean squared  error for raw SIs
    mse_norm_si <- mean(vse_norm_si)          # Mean squared  error for normalized SIs
    # Compound Objects
    simpar <- named(name = "sap1", cs, pkr = NULL, x0 = x0_ppm, A = A_raw_ppm, lambda = lambda_ppm, noise)
    lcpar  <- named(A = A_sc_ppm, lambda = lambda_ppm, x0 = x0_ppm)
    args   <- named(nfit, smopts, delta, sfr = sfr_ppm, wshw)
    meta   <- named(name, fq, simpar)
    sit    <- named(wsrm = si, nvrm = si, sm = si, sup = si_sup, al = NULL)
    mse    <- named(raw = mse_raw_si, norm = mse_norm_si, sm = mse_raw_si, smnorm = mse_norm_si)
    peak   <- peak_idx
    # Return Objects
    locals()
}

# Sim #####

`sim Documentation` <- NULL

#' @title The Sim Dataset
#'
#' @description
#' A simulated dataset generated from the [Blood](
#' https://spang-lab.github.io/metabodecon/articles/Datasets.html#Blood)
#' dataset.
#'
#' @format A [spectra](metabodecon_classes) object consisting of 16 [spectrum]
#' objects, where each spectrum object contains 2048 datapoints ranging from
#' 3.60 to 3.29 ppm
"sim"

update_sim <- function(overwrite = FALSE,
                       verbose = TRUE,
                       nworkers = 1) {
    if (isTRUE(overwrite)) {
        path <- file.path(pkg_file("example_datasets/bruker"), "sim")
        if (verbose) {
            logf("Overwrite is TRUE. Updating %s." , path)
            if (!get_yn_input("Continue?")) return(invisible())
        }
        sim <- make_sim(nworkers = nworkers)
        save_spectra(sim, path, force = TRUE, verbose = verbose)
        usethis::use_data(sim, overwrite = TRUE)
    } else {
        path <- file.path(tmpdir(subdir = TRUE, create = TRUE), "sim")
        if (verbose) {
            logf("Overwrite is FALSE. Writing spectra to %s." , path)
            if (!get_yn_input("Continue?")) return(invisible())
        }
        sim <- make_sim(nworkers = nworkers)
        save_spectra(sim, path, force = FALSE, verbose = verbose)
        save(sim, file = file.path(path, "sim.rda"), compress = "bzip2")
    }
    path
}

# "C:/Users/tobi/AppData/Local/Temp/RtmpE3Ma4o/metabodecon/4b0857617065/sim/sim.rda"

make_sim <- function(nworkers = 1) {
    decons <- deconvolute_blood(nworkers = nworkers)
    simpars <- lapply(decons, get_sim_params, pkr = c(3.52, 3.37))
    simpars[[2]] <- get_sim_params(decons[[2]], pkr = c(3.51, 3.36)) # (1)
    # (1) Blood_02 is shifted approx. 0.01 ppm to the right, so we also need to
    # shift the interval from which we pick our peaks by 0.01 to end up with
    # signals from the same metabolites. This was determined by visual
    # inspection after deconvoluting the blood dataset.
    sim <- lapply(simpars, function(simpar) {
        simulate_spectrum(
            name = gsub("blood", "sim", simpar$name),
            cs = seq(from = 3.59, length.out = 2048, by = -0.00015),
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
#' x <- deconvolute_ispec(spec)
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

get_sim1_decon1 <- function() {
    generate_lorentz_curves(
        data_path = sim[[1]],
        sfr = c(3.55, 3.35),
        wshw = 0,
        ask = FALSE,
        verbose = FALSE
    )
}
