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
