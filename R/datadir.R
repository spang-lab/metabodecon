# Exported #####

#' @title Return path to metabodecon's data directory
#' @description Returns the path to the directory where [download_example_datasets()] stores metabodecon's example data sets or any file within that directory.
#' By default this directory is a subdirectory of R's temporary session directory.
#' If `persistent` is set to `TRUE`, the directory equals the data directory returned by [tools::R_user_dir()] instead.
#' @param file Relative path to a file within the data directory.
#' @param warn Print a warning message when the requested path does not yet exist?
#' @param persistent Return the path to the persistent data directory instead of the temporary one?
#' @examples
#' # Return path to persistent data dir if it exists, else path to temp data dir
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
#' datadir(file = "bruker/urine")
#' @seealso [download_example_datasets()]
#' @details The decision to use a temporary data dir as default and a persistent one only optionally was made to conform to CRAN package policies, which state that:
#' *Packages should not write in the user's home filespace (including clipboards), nor anywhere else on the file system apart from the R session's temporary directory (or during installation in the location pointed to by TMPDIR: and such usage should be cleaned up). Installing into the system's R installation (e.g., scripts to its bin directory) is not allowed.*
#'
#' *Limited exceptions may be allowed in interactive sessions if the package obtains confirmation from the user.*
#'
#' *For R version 4.0 or later (hence a version dependency is required or only conditional use is possible), packages may store user-specific data, configuration and cache files in their respective user directories obtained from tools::R_user_dir(), provided that by default sizes are kept as small as possible and the contents are actively managed (including removing outdated material).*
#'
#' Source: [cran.r-project.org/web/packages/policies](https://cran.r-project.org/web/packages/policies.html)
#' @seealso [datadir_persistent()], [datadir_temp()]
#' @export
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

#' @title Persistent Data Directory
#' @description Returns the path to the persistent data directory where metabodecon's data sets are stored.
#' This directory equals the data directory returned by [tools::R_user_dir()] plus additional path normalization.
#' @return Returns the path to the persistent data directory.
#' @seealso [datadir()], [datadir_temp()]
#' @export
datadir_persistent <- function() {
    p <- tools::R_user_dir("metabodecon", "data")
    normalizePath(p, "/", mustWork = FALSE)
}

#' @title Temporary Data Directory
#' @description Returns the path to the temporary data directory where metabodecon's data sets are stored.
#' This directory equals subdirectory `"data"` of metabodecons temporary session directory `[tempdir()]` plus additional path normalization.
#' @return Returns the path to the temporary data directory.
#' @export
#' @seealso [tempdir()], [datadir()], [datadir_persistent()]
datadir_temp <- function() {
    p <- file.path(tempdir(), "data")
    normalizePath(p, "/", mustWork = FALSE)
}

#' @title Temporary Session Directory
#' @description Returns the path to metabodecon's temporary session directory.
#' This directory equals subdirectory `"metabodecon"` of R's temporary session directory `[tempdir()]` plus additional path normalization.
#' @return Returns the path to the temporary session directory.
#' @export
tempdir <- function() {
    p <- file.path(base::tempdir(), "metabodecon")
    normalizePath(p, "/", mustWork = FALSE)
}

# Helpers #####

zip_temp <- function() {
    p <- file.path(datadir_temp(), "example_datasets.zip")
    normalizePath(p, "/", mustWork = FALSE)
}

zip_persistent <- function() {
    p <- file.path(datadir_persistent(), "example_datasets.zip")
    normalizePath(p, "/", mustWork = FALSE)
}

# Deprecated #####

#' @export
#' @name get_data_dir
#' @title Retrieve directory path of an example dataset
#' @description Returns the path to the directory storing the example files
#' shipped with metabodecon.
#'
#' Deprecated. Please use [datadir()] instead. See examples below for usage.
#' ``
#' @param dataset_name Either `""`, `"test"`, `"blood"` or `"urine"`.
#' @param warn Whether to print a warning message when the example folders do
#' not yet exist, i.e. [download_example_datasets()] has not been called yet.
#' @examples
#' x <- get_data_dir("urine")                     # Deprecated
#' y <- datadir("example_datasets/bruker/urine")  # Preferred
#' @seealso [download_example_datasets()]
get_data_dir <- function(dataset_name = c("", "blood", "test", "urine"), warn = TRUE) {
  dataset_name <- match.arg(dataset_name)
  base_data_dir <- datadir()
  data_dir <- file.path(base_data_dir, "example_datasets/bruker", dataset_name, fsep = "/")
  if (warn && !dir.exists(data_dir)) {
    warning(data_dir, " does not exist. Please call `download_example_datasets(extract = TRUE)` first.", call. = FALSE)
  }
  return(data_dir)
}
