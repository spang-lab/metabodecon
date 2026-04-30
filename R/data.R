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
#' @param dst_dir
#' The destination directory where the downloaded datasets will be stored. If
#' NULL, the function will return the path to the cached zip file.
#'
#' @param extract
#' Logical. If TRUE, the downloaded zip file will be extracted.
#'
#' @param persistent
#' Logical. If TRUE, the downloaded datasets will be cached at
#' [metabodecon::datadir_persistent()] to speed up future calls to
#' `download_example_datasets()`. If FALSE, the datasets will be cached at
#' [metabodecon::datadir_temp()]. If NULL, the function will check both paths for the cached
#' datasets but will return [metabodecon::datadir_temp()] if the cached file does not yet
#' exist.
#'
#' @param overwrite
#' Logical. If TRUE, existing files with the same name in the destination
#' directory will be overwritten.
#'
#' @param silent
#' Logical. If TRUE, no output will be printed to the console.
#'
#' @return
#' The path to the downloaded (and possibly extracted) datasets.
#'
#' @seealso
#' [metabodecon::datadir()]
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' if (interactive()) {
#'      zip <- download_example_datasets(extract = FALSE, persistent = FALSE)
#'      dir <- download_example_datasets(extract = TRUE)
#' }
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
    cached_zip <- cache_example_datasets(persistent, extract, silent)
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
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' [metabodecon::download_example_datasets()] stores metabodecon's example data sets or any
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
#' *Packages should not write in the user's home filespace (including*
#' *clipboards), nor anywhere else on the file system apart from the R*
#' *session's temporary directory \[...\] Limited exceptions may be allowed*
#' *in interactive sessions if the package obtains confirmation from the*
#' *user. For R version 4.0 or later \[...\] packages may store user-specific*
#' *data, configuration and cache files in their respective user directories*
#' *obtained from `tools::R_user_dir()` \[...\].*
#'
#' Source:
#' [cran.r-project.org/web/packages/policies](
#' https://cran.r-project.org/web/packages/policies.html
#' ).
#'
#' @seealso
#' [metabodecon::download_example_datasets()],
#' [metabodecon::datadir_persistent()],
#' [metabodecon::datadir_temp()]
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' @seealso [metabodecon::datadir()], [metabodecon::datadir_temp()]
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
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
#' of metabodecons temporary session directory [metabodecon::tmpdir()] plus additional path
#' normalization.
#'
#' @return Returns the path to the temporary data directory.
#'
#' @seealso [metabodecon::tmpdir()], [metabodecon::datadir()], [metabodecon::datadir_persistent()]
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
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
#' [metabodecon::datadir_temp()]
#' [metabodecon::datadir_persistent()]
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' Deprecated since metabodecon v1.2.0. Please use [metabodecon::datadir()] instead. See
#' examples below for usage.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param dataset_name Either `""`, `"test"`, `"blood"`, `"urine"` or `"aki"`.
#'
#' @param warn Whether to print a warning message when the example folders do
#' not yet exist, i.e. [metabodecon::download_example_datasets()] has not been called yet.
#'
#' @return Path to the directory storing the example files.
#'
#' @seealso [metabodecon::download_example_datasets()]
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' x <- get_data_dir("urine")                     # Deprecated
#' y <- datadir("example_datasets/bruker/urine")  # Preferred
#' cat(x, y, sep = "\n")
get_data_dir <- function(dataset_name = c("", "blood", "test", "urine", "aki"), warn = TRUE) {
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
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' str(xds)
xds <- list(
    url = "https://github.com/spang-lab/metabodecon/releases/download/v1.6.3/example_datasets.zip",
    zip_size = 75391701,
    dir_size = 113207749,
    n_files = 1336
)

#' @noRd
#'
#' @title Cache Example Datasets
#'
#' @description
#' If file `"example_datasets.zip"` does not yet exist at [metabodecon::datadir()], it is
#' downloaded. If the zip exists but is outdated or corrupt, it gets overwritten
#' automatically. If `extract` is TRUE, the zip file is extracted, replacing an
#' existing `"example_datasets"` directory. The path to the zip file is
#' returned.
#'
#' @param persistent
#' If TRUE, the datasets are cached permanently.
#' If FALSE, the datasets are cached temporarily.
#' If NULL, the path to the permanent zip file is returned if it exists and has
#' the correct size else the path to the temporary zip file is returned.
#'
#' @param extract If TRUE, the datasets are extracted after being cached.
#'
#' @return The path to the cached datasets.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' \donttest{
#' cache_example_datasets(persistent = FALSE, extract = FALSE)
#' }
cache_example_datasets <- function(persistent = NULL,
                                   extract = FALSE,
                                   silent = FALSE) {

    pzip <- zip_persistent()
    pzip_exists <- file.exists(pzip)
    pzip_ok <- isTRUE(file.size(pzip) == xds$zip_size) # Implies existence
    pzip_nok <- !pzip_ok

    tzip <- zip_temp()
    tzip_ok <- isTRUE(file.size(tzip) == xds$zip_size)
    tzip_nok <- !tzip_ok

    use_pzip <- isTRUE(persistent) || (is.null(persistent) && pzip_exists)
    use_tzip <- !use_pzip

    dst <- if (use_pzip) pzip else tzip
    cf <- if (use_pzip && tzip_ok) tzip else if (use_tzip && pzip_ok) pzip else NULL

    if (use_tzip && tzip_nok && pzip_exists && pzip_nok) {
        fmt <- paste0("\n",
            "Found outdated or corrupt persistent cache file at:\n%s\n",
            "To update/repair, use 'download_example_datasets(persistent = TRUE)'\n",
            "To remove, use 'unlink(\"%s\", recursive = TRUE)'\n"
        )
        msg <- sprintf(fmt, pzip, pzip)
        warning(msg, immediate. = TRUE, call. = FALSE)
    }

    if ((use_tzip && tzip_nok) || (use_pzip && pzip_nok)) {
        download_example_datasets_zip(dst, copyfrom = cf, silent = silent)
    }

    dstdir <- file.path(dirname(dst), "example_datasets")
    if (extract) {
        files <- dir(dstdir, recursive=TRUE, full.names=TRUE) # 0.6s on r4
        size <- sum(file.info(files)$size) # 0.3s on r4
        if (size != xds$dir_size) { # captures missing and outdated dstdir
            unlink(dstdir, recursive = TRUE) # 1s on r4
            system.time(extract_example_datasets(dst)) # 6.9s on r4
        }
    }
    dst
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' \donttest{
#' download_example_datasets_zip("/path/to/your/directory/example_datasets.zip")
#' }
download_example_datasets_zip <- function(path,
                                          copyfrom = NULL,
                                          silent = FALSE,
                                          url = xds$url,
                                          zip_size = xds$zip_size) {
    mkdirs(dirname(path))
    path <- normalizePath(path, "/", mustWork = FALSE)
    if (!is.null(copyfrom) && file.exists(copyfrom) && file.size(copyfrom) == zip_size) {
        if (!silent) message(sprintf("Copying cached archive %s to %s", copyfrom, path))
        file.copy(copyfrom, path, overwrite = TRUE)
    } else {
        if (!silent) message(sprintf("Downloading %s as %s", url, path))
        utils::download.file(url, path, quiet = TRUE)
    }
    if (!is.null(zip_size) && isTRUE(file.size(path) != zip_size)) {
        msg <- "Downloaded zip at %s has size %d instead of expected %d"
        stop(sprintf(msg, path, file.size(path), zip_size))
    }
    path
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
zip_temp <- function() {
    p <- file.path(datadir_temp(), "example_datasets.zip")
    normalizePath(p, "/", mustWork = FALSE)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
zip_persistent <- function() {
    p <- file.path(datadir_persistent(), "example_datasets.zip")
    normalizePath(p, "/", mustWork = FALSE)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
tmpfile <- function(pattern = "file", fileext = "") {
    tempfile(pattern, tmpdir(create = TRUE), fileext)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
testdir <- function(p = NULL) {
    norm_path(paste(tmpdir("tests"), p, sep = "/"))
    # use paste instead of file.path, because it can deal with NULL
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
mockdir <- function() {
    norm_path(file.path(tmpdir(), "mocks"))
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
#' @description
#' Create and return cache dir. If existing, the persistent cache dir is
#' returned, else the temp cache dir. To force creation of the persistent cache
#' dir, call once with `persistent=TRUE`.
cachedir <- function(subdir = NULL, persistent = NULL) {
    tcd <- file.path(tmpdir(), "cache") # temporary cache dir
    pcd <- file.path(tools::R_user_dir("metabodecon", "cache")) # persistent cache dir
    cd <- if (isTRUE(persistent) || (is.null(persistent) && dir.exists(pcd))) pcd else tcd
    if (!is.null(subdir)) cd <- file.path(cd, subdir)
    ncd <- normalizePath(cd, "/", mustWork = FALSE)
    mkdirs(ncd)
}

#' @export
#'
#' @title Get Default Cache Directory for Deconvolution
#'
#' @description
#' Returns the temporary session-scoped cache directory used by deconvolution
#' functions such as [metabodecon::deconvolute()] and the internal
#' `grid_deconvolute_spectrum()` helper.
#'
#' @return Path to the shared temporary cache directory for deconvolution.
#'
#' @author 2026 Tobias Schmidt: initial version.
#'
#' @examples
#' decon_cachedir()
decon_cachedir <- function() {
    x <- cachedir("deconvs", persistent = FALSE)
    if (!dir.exists(x)) dir.create(x, recursive = TRUE)
    x
}

# Sap (Public) #####

sap_docs <- NULL # To get a symbol in the outline

#' @name sap
#'
#' @title The SAP Dataset
#'
#' @description
#' The SAP Dataset consists of a single 'Simple-As-Possible' (SAP) spectrum. The
#' purpose of the SAP spectrum is to provide a straightforward example that can
#' be used to test and understand the deconvolution algorithm in detail.
#'
#' @details
#' The first (and only) spectrum within the SAP dataset contains 128 datapoints
#' ranging from -6.4 to 6.4 ppm with four peaks. A rough sketch of the spectrum
#' is shown below:
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
#'
"sap" # To regenerate this dataset, see `data-raw/data.R`.

# Sim (Public) #####

sim_docs <- NULL # To get a symbol in the outline

#' @title The Sim Dataset
#'
#' @description
#' A simulated dataset generated from the [Blood](
#' https://spang-lab.github.io/metabodecon/articles/Datasets.html#Blood)
#' dataset.
#'
#' @format
#' A `spectra` object consisting of 16 `spectrum` objects, where each spectrum
#' contains 2048 datapoints ranging from 3.60 to 3.29 ppm. For details about
#' `spectrum` and `spectra` objects see [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
"sim" # To regenerate this dataset, see `data-raw/data.R`.

# AKI #####

read_aki_metadata <- function(aki_path) {
    meta <- read.csv(file.path(aki_path, "s_MTBLS24.txt"), sep = "\t")
    meta <- meta[, c("Sample.Name", "Factor.Value.Acute.Kidney.Injury.")]
    meta <- setNames(meta, c("sid", "type"))
    meta$type <- ifelse(meta$type == "Acute Kidney Injury", "AKI", "Control")
    # Fix samples with wrong date in their sample
    meta$sid[meta$sid == "AKI_8_24_105_110812"] <- "AKI_8_24_105_110816"
    meta$sid[meta$sid == "AKI_8_24_106_110812"] <- "AKI_8_24_106_110816"
    meta$sid[meta$sid == "AKI_8_24_107_110812"] <- "AKI_8_24_107_110816"
    meta$sid[meta$sid == "AKI_8_24_108_110812"] <- "AKI_8_24_108_110816"
    meta$sid[meta$sid == "AKI_8_24_109_110812"] <- "AKI_8_24_109_110816"
    meta$sid[meta$sid == "AKI_8_24_110_110812"] <- "AKI_8_24_110_110816"
    # Sort in the same order as the spectra on disk
    rownames(meta) <- meta$sid
    spectra_dirs <- grep("AKI_8_24", dir(aki_path), value = TRUE)
    meta <- meta[spectra_dirs, , drop = FALSE]
    meta
}

read_aki_data <- function() {
    # Read in data
    aki_path <- datadir("example_datasets/bruker/aki")
    if (!dir.exists(aki_path)) metabodecon::download_example_datasets()
    meta <- read_aki_metadata(aki_path)
    spectra_raw <- metabodecon::read_spectra(aki_path)
    stopifnot(all.equal(names(spectra_raw), meta$sid))

    # Normalize spectra
    spectra <- creatinine_normalize(spectra_raw)

    list(spectra = spectra, meta = meta)
}

creatinine_normalize <- function(spectra, cr = c(3.053, 3.011)) {
    stopifnot(length(cr) == 2, cr[1] > cr[2])
    spectra_normed <- spectra # copy spectra to new object for normalization
    for (i in seq_along(spectra)) {
        s <- spectra[[i]]
        idx <- which(s$cs >= cr[2] & s$cs <= cr[1])
        ci <- sum(s$si[idx]) # creatinine intensity
        spectra_normed[[i]]$si <- s$si / ci
    }
    spectra_normed
}