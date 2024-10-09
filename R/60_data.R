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
                                      overwrite = FALSE) {

    # Example:
    # input   dst_dir    = C:/Users/max/Downloads
    # var     cached_zip = C:/Users/max/.local/share/R/metabodecon/example_datasets.zip
    # var     cached_xds = C:/Users/max/.local/share/R/metabodecon/example_datasets
    # return  dst_zip    = C:/Users/max/Downloads/example_datasets.zip
    # return  dst_xds    = C:/Users/max/Downloads/example_datasets
    cached_zip <- cache_example_datasets(persistent = persistent, extract = extract, overwrite = overwrite)
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

# Helpers #####

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
cache_example_datasets <- function(persistent = NULL, extract = FALSE, overwrite = FALSE) {

    pdir <- datadir_persistent()
    pzip <- file.path(pdir, "example_datasets.zip")
    pzip_exists <- file.exists(pzip)
    pzip_is_missing <- !pzip_exists
    pzip_has_correct_size <- isTRUE(file.size(pzip) == xds$zip_size) # implies existence
    pzip_has_wrong_size <- !pzip_has_correct_size

    tdir <- datadir_temp()
    tzip <- file.path(tdir, "example_datasets.zip")
    tzip_has_correct_size <- isTRUE(file.size(tzip) == xds$zip_size)
    tzip_has_wrong_size <- !tzip_has_correct_size

    persistent <- if (isTRUE(persistent) || (is.null(persistent) && pzip_has_correct_size)) TRUE else FALSE

    if (isTRUE(persistent)) {
        if (pzip_is_missing || (pzip_has_wrong_size && overwrite)) download_example_datasets_zip(pzip, copyfrom = tzip)
        if (pzip_exists && pzip_has_wrong_size && !overwrite) {
            msg <- "Path %s exists already, but has size %d instead of %d. Set `overwrite = TRUE` to overwrite."
            stop(sprintf(msg, pzip, file.size(pzip), xds$zip_size))
        }
    } else {
        if (tzip_has_wrong_size) download_example_datasets_zip(tzip, copyfrom = pzip) # no need to ask for permission here because we only write to the temp dir
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
download_example_datasets_zip <- function(path, copyfrom = NULL) {
    mkdirs(dirname(path))
    if (!is.null(copyfrom) && file.exists(copyfrom) && file.size(copyfrom) == xds$zip_size) {
        message(sprintf("Copying cached archive %s to %s", copyfrom, path))
        file.copy(copyfrom, path, overwrite = TRUE)
    } else {
        message(sprintf("Downloading %s as %s", xds$url, path))
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
