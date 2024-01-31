# Exports #####

#' @title Return path to metabodecon's data directory
#' @description Returns the path to the directory where [download_example_datasets()] stores metabodecon's example data sets or any file within that directory. By default this directory is a subdirectory of R's temporary session directory. If `persistent` is set to `TRUE`, the directory equals the data directory returned by [tools::R_user_dir()] instead.
#' @param file Relative path to a file within the data directory.
#' @param warn Print a warning message when the requested path does not yet exist?
#' @param persistent Return the path to the persistent data directory instead of the temporary one?
#' @examples
#' # Return path to persistent data directory if it exists, otherwise return path to temporary data directory
#' datadir()
#'
#' # Return path to temporary data directory
#' datadir(persistent = FALSE)
#'
#' # Return path to persistent data directory
#' datadir(persistent = TRUE)
#'
#' # Return path to "<persistent-data-dir>/bruker/urine" if it exists.
#' # Else return path to "<temp-data-dir>/bruker/urine"
#' datadir(path = "bruker/urine")
#' @seealso [download_example_datasets()]
#' @details The decision to use a temporary data dir as default and a persistent one only optionally was made to conform to CRAN package policies, which state that:
#' *Packages should not write in the user's home filespace (including clipboards), nor anywhere else on the file system apart from the R session's temporary directory (or during installation in the location pointed to by TMPDIR: and such usage should be cleaned up). Installing into the system's R installation (e.g., scripts to its bin directory) is not allowed.*
#'
#' *Limited exceptions may be allowed in interactive sessions if the package obtains confirmation from the user.*
#'
#' *For R version 4.0 or later (hence a version dependency is required or only conditional use is possible), packages may store user-specific data, configuration and cache files in their respective user directories obtained from tools::R_user_dir(), provided that by default sizes are kept as small as possible and the contents are actively managed (including removing outdated material).*
#'
#' Source: [cran.r-project.org/web/packages/policies](https://cran.r-project.org/web/packages/policies.html)
#' @export
datadir <- function(file = "", warn = TRUE, persistent = NULL) {
    datadir_persistent <- datadir_persistent()
    datadir <- datadir_temp <- datadir_temp()
    if (isTRUE(persistent) || is.null(persistent) && dir.exists(datadir_persistent)) {
        datadir <- datadir_persistent
    }
    file_path <- file.path(datadir, file, fsep = "/")
    if (warn && !dir.exists(file_path)) {
        warning(file_path, " does not exist. Please call `download_example_datasets()` first.", call. = FALSE)
    }
    normalizePath(file_path, "/", mustWork = FALSE)
}

#' @name download_example_datasets
#' @title Download metabodecon example datasets
#' @description Downloads example datasets for testing metabodecon's functionality. Due to the size constraints for R packages, these datasets are not included by default when the package is installed.
#' Datasets are download and extracted to the directory returned by [datadir()].
#' By default this directory is a subdirectory of R's temporary session directory.
#' If `persistent` is set to `FALSE`, the directory equals the one returned by [datadir_temp()]
#' If `persistent` is set to `TRUE`, the directory equals the one returned by [datadir_persistent()].
#' If `persistent` is set to `NULL` and `ask` is TRUE, the user is asked which one he wants to use.
#' @examples \dontrun{
#' download_example_datasets()
#' }
#' @seealso [datadir()]
#' @export
download_example_datasets <- function(dst_dir = NULL,
                                      extract = FALSE,
                                      cache_persistent = FALSE) {
    cached_zip <- cache_example_datasets(persistent = cache_persistent, extract = extract)
    if (is.null(dst_dir)) {
        return(cached_zip)
    } else {
        dst_zip <- file.path(dst_dir, basename(cached_zip))
        message(sprintf("Copying cached archive %s to %s", cached_zip, dst_dir))
        file.copy(cached_zip, dst_zip)
        if (extract) extract_example_datasets(dst_zip)
        return(dst_zip)
    }
}

#' @title Cache Example Files
#' @description Simpler Version of [cache_example_datasets()] that does not prompt the user for input.
#' @param persistent Logical. If TRUE, the datasets are cached permanently. If NULL, the function prompts the user. Default is NULL.
#' @return The path to the cached datasets.
#' @examples
#' cache_example_datasets(persistent = FALSE)
#' @export
cache_example_datasets <- function(persistent = FALSE, extract = TRUE, overwrite = FALSE) {
    dir <- if (persistent) datadir_persistent() else datadir_temp()
    zip <- file.path(dir, "example_datasets.zip")
    zip_exists <- file.exists(zip)
    zip_has_wrong_size <- isTRUE(file.size(zip) != 38253468)
    if ((!zip_exists) || (zip_exists && zip_has_wrong_size && overwrite)) {
        dir2 <- if (persistent) datadir_temp() else datadir_persistent()
        zip2 <- file.path(dir2, "example_datasets.zip")
        zip <- download_example_datasets_zip(zip, copyfrom = zip2)
    } else if (zip_exists && zip_has_wrong_size && !overwrite) {
        stop(paste("Path", zip, "exists already, but has size", file.size(zip), "instead of 38253468. Set `overwrite = TRUE` to overwrite."))
    }
    zip
}

extract_example_datasets <- function(path = datadir("example_datasets.zip")) {
    utils::unzip(zipfile = path, exdir = dirname(path))
}

# Helpers #####

#' @title Download Example Datasets Zip
#' @description This function downloads the example_datasets.zip file from the specified URL and saves it to the destination directory provided. If the directory does not exist, it will be created. If the file already exists in the directory, it will be overwritten without asking for permission.
#' @param path A string. The path where the downloaded file will be saved.
#' @return The path where the downloaded file is saved.
#' @examples \dontrun{
#' download_example_datasets_zip("/path/to/your/directory/example_datasets.zip")
#' }
#' @noRd
download_example_datasets_zip <- function(path, copyfrom = NULL) {
    mkdirs(dirname(path))
    url <- "https://github.com/spang-lab/metabodecon/releases/download/v1.0.2/example_datasets.zip"
    if (!is.null(copyfrom) && file.exists(copyfrom) && file.size(copyfrom) == 38253468) {
        message(sprintf("Copying cached archive %s to %s", copyfrom, path))
        file.copy(copyfrom, path, overwrite = TRUE)
    } else {
        message(sprintf("Downloading %s as %s", url, path))
        utils::download.file(url, path, quiet = TRUE)
    }
    path
}

readline <- function(...) {
    base::readline(...) # we must have our own copy of readline in the package namespace so we can mock it in tests
}

# Temporary Directory
tempdir <- function() {
    p <- file.path(base::tempdir(), "metabodecon")
    normalizePath(p, "/", mustWork = FALSE)
}

# Persistent Data Directory
datadir_persistent <- function() {
    p <- tools::R_user_dir("metabodecon", "data")
    normalizePath(p, "/", mustWork = FALSE)
}

# Temporary Data Directory
datadir_temp <- function() {
    p <- file.path(tempdir(), "data")
    normalizePath(p, "/", mustWork = FALSE)
}

testdir <- function() {
    p <- file.path(tempdir(), "tests")
    normalizePath(p, "/", mustWork = FALSE)
}

mockdir <- function() {
    p <- file.path(tempdir(), "mocks")
    normalizePath(p, "/", mustWork = FALSE)
}

# Deprecated #####

#' @title Cache Example Files
#' @description This function caches example datasets to speed up future calls. It checks if the file exists and has the correct size. If not, it prompts the user for action.
#' @param persistent Logical. If TRUE, the datasets are cached permanently. If NULL, the function prompts the user. Default is NULL.
#' @param ask Logical. If TRUE, the function operates in interactive mode and prompts the user for input. Default is interactive().
#' @return The path to the cached datasets.
#' @examples
#' cache_example_datasets(persistent = FALSE)
#' @export
cache_example_datasets_bak <- function(persistent = NULL, ask = interactive()) {
    persistent_path <- file.path(datadir_persistent(), "example_datasets.zip")
    persistent_path_exists <- file.exists(persistent_path)
    if (is.null(persistent)) {
        persistent <- if (persistent_path_exists) TRUE else if (!ask) FALSE else {
            prompt <- paste("Cache example datasets (38.3 MB) permanently on disk at", persistent_path, "to speed up future calls to this function?")
            persistent <- get_yn_input(prompt)
        }
    }
    path <- if (persistent) persistent_path else file.path(datadir_temp(), "example_datasets.zip")
    exists <- file.exists(path)
    size <- file.size(path)
    has_correct_size <- isTRUE(size == 38253468)
    if (exists && has_correct_size) {
        return(path)
    }
    if (persistent && exists && !has_correct_size) {
        prompt <- paste("Path", path, "exists already, but has size", size, "instead of 38253468. Overwrite?")
        if (!(get_yn_input(prompt))) {
            return(NULL)
        }
    }
    path <- download_example_datasets_zip(path)
    path
}