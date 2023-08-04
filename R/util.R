#' @export
#' @name get_data_dir
#' @title Retrieve directory path of an example dataset
#' @description This function return the path to the directory storing the
#' example files shipped with metabodecon.
#' @param dataset_name Either an string, "test" or "urine".
#' @param warn Whether to print a warning message when the example folders do
#' not yet exist, i.e. [download_example_datasets()] has not been called yet.
#' @examples get_data_dir("urine")
#' @seealso [download_example_datasets()]
get_data_dir <- function(dataset_name = c("", "test", "urine"), warn = TRUE) {
  dataset_name = match.arg(dataset_name)
  base_data_dir <- normalizePath(
    path = tools::R_user_dir(package = "metabodecon", which = "data"),
    winslash = "/",
    mustWork = FALSE
  )
  data_dir <- file.path(base_data_dir, dataset_name, fsep = "/")
  if (warn && !dir.exists(data_dir)) {
    warning(
      data_dir, " does not exist. Please call `download_example_datasets()` first.",
      call. = FALSE
    )
  }
  return(data_dir)
}

#' @export
#' @name download_example_datasets
#' @title Download metabodecon example datasets
#' @description Downloads the example files shipped wit metabodecon.  Due to the size
#' constraints for R packages, these datasets are not included by default when
#' the package is installed, but must be downloaded explicitly afterwards.
#' @examples \dontrun{download_example_datasets()}
#' @seealso [get_data_dir()]
download_example_datasets <- function() {
  data_dir = get_data_dir(warn = FALSE)
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  url <- "https://gitlab.spang-lab.de/api/v4/projects/690/repository/archive.zip"
  repo_zip <- file.path(data_dir, "repo.zip")
  repo_dir <- gsub(".zip", "", repo_zip)
  headers <- c(`PRIVATE-TOKEN` = "glpat-ndxyfy5Ty7yksgy9MAFs")
  if (!file.exists(repo_zip)) {
    message(sprintf("Downloading %s as %s", url, repo_zip))
    utils::download.file(url, repo_zip, headers = headers, quiet = TRUE)
  }
  if (!dir.exists(repo_dir)) {
    message("Extracting ", repo_zip)
    utils::unzip(zipfile = repo_zip, exdir = repo_dir)
  }
  for (dataset_name in c("urine", "test")) {
    dataset_dir = get_data_dir(dataset_name, warn = FALSE)
    src_dir <- Sys.glob(file.path(repo_dir, "*/misc/datasets", dataset_name))
    if (!dir.exists(dataset_dir)) {
      message("Copying\n", src_dir, "to\n", data_dir)
      file.copy(src_dir, data_dir, recursive = TRUE)
    }
  }
}
