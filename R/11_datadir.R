# Exported Main #####

#' @export
#' @title Return path to metabodecon's data directory
#' @description Returns the path to the directory where [download_example_datasets()] stores metabodecon's example data sets or any file within that directory. By default this directory is a subdirectory of R's temporary session directory. If `persistent` is set to `TRUE`, the directory equals the data directory returned by [tools::R_user_dir()] instead.
#' @param file Relative path to a file within the data directory.
#' @param warn Print a warning message when the requested path does not yet exist?
#' @param persistent Return the path to the persistent data directory instead of the temporary one?
#' @return Path to the data directory or a file within it.
#' @details The decision to use a temporary data dir as default and a persistent one only optionally was made to conform to CRAN package policies, which state that: *Packages should not write in the user's home filespace (including clipboards), nor anywhere else on the file system apart from the R session's temporary directory \[...\] Limited exceptions may be allowed in interactive sessions if the package obtains confirmation from the user. For R version 4.0 or later \[...\] packages may store user-specific data, configuration and cache files in their respective user directories obtained from [tools::R_user_dir()] \[...\]*. Source: [cran.r-project.org/web/packages/policies](https://cran.r-project.org/web/packages/policies.html).
#' @seealso [download_example_datasets()], [datadir_persistent()], [datadir_temp()]
#' @examples
#' # Get temporary datadir and persistent datadir
#' datadir(persistent = FALSE, warn = FALSE)
#' datadir(persistent = TRUE, warn = FALSE)
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
#' @title Return Path to File or Directory in metabodecon Package
#' @description Recursively searches for files or directories within the 'metabodecon' package that match the given name.
#' @param name The name to search for.
#' @return The file or directory path.
#' @examples
#' # Unambiguous paths
#' metabodecon_file("sim_01")
#' metabodecon_file("urine_1")
#' metabodecon_file("urine_1.dx")
#'
#' # Ambiguous paths (i.e. multiple matches)
#' metabodecon_file("sim")
#' metabodecon_file("urine")
#'
#' # Non-existing paths (i.e. a charactr vector of length zero gets returned)
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

# Exported Helpers #####

#' @export
#' @title Persistent Data Directory
#' @description Returns the path to the persistent data directory where metabodecon's data sets are stored. This directory equals the data directory returned by [tools::R_user_dir()] plus additional path normalization.
#' @return Path to the persistent data directory.
#' @seealso [datadir()], [datadir_temp()]
#' @examples
#' datadir_persistent()
datadir_persistent <- function() {
    p <- tools::R_user_dir("metabodecon", "data")
    normalizePath(p, "/", mustWork = FALSE)
}

#' @export
#' @title Temporary Data Directory
#' @description Returns the path to the temporary data directory where metabodecon's data sets are stored. This directory equals subdirectory `"data"` of metabodecons temporary session directory `[tmpdir()]` plus additional path normalization.
#' @return Returns the path to the temporary data directory.
#' @seealso [tmpdir()], [datadir()], [datadir_persistent()]
#' @examples
#' datadir_temp()
datadir_temp <- function() {
    p <- file.path(tmpdir(), "data")
    normalizePath(p, "/", mustWork = FALSE)
}

#' @export
#' @title Temporary Session Directory
#' @description Returns the path to metabodecon's temporary session directory.
#' This directory equals subdirectory `"metabodecon"` of R's temporary session directory `[base::tempdir()]` plus additional path normalization.
#' @param subdir Optional subdirectory within the temporary session directory.
#' @param create Whether to create the directory if it does not yet exist.
#' @return Returns the path to the temporary session directory.
#' @seealso [datadir_temp()], [datadir_persistent()]
#' @examples
#' tmpdir()
#' tmpdir("simulate_spectra")
tmpdir <- function(subdir = NULL, create = FALSE) {
    p <- file.path(base::tempdir(), "metabodecon")
    if (!is.null(subdir)) p <- file.path(p, subdir)
    if (create) base::dir.create(p, recursive = TRUE, showWarnings = FALSE)
    normalizePath(p, "/", mustWork = FALSE)

}

# Private Helpers #####

zip_temp <- function() {
    p <- file.path(datadir_temp(), "example_datasets.zip")
    normalizePath(p, "/", mustWork = FALSE)
}

zip_persistent <- function() {
    p <- file.path(datadir_persistent(), "example_datasets.zip")
    normalizePath(p, "/", mustWork = FALSE)
}

# Exported Deprecated #####

#' @export
#' @title Retrieve directory path of an example dataset
#' @description
#' Returns the path to the directory storing the example files shipped with metabodecon.
#'
#' Deprecated since metabodecon v1.2.0. Please use [datadir()] instead. See examples below for usage.
#'
#' `r lifecycle::badge("deprecated")`
#' @param dataset_name Either `""`, `"test"`, `"blood"` or `"urine"`.
#' @param warn Whether to print a warning message when the example folders do not yet exist, i.e. [download_example_datasets()] has not been called yet.
#' @return Path to the directory storing the example files.
#' @seealso [download_example_datasets()]
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
