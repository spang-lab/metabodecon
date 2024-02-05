# Exported #####

#' @name download_example_datasets
#' @title Download metabodecon Example Datasets
#' @description This function downloads example datasets that can be used to test the functionality of the metabodecon package.
#' These datasets are not included in the package by default due to size constraints.
#' The datasets are downloaded as a zip file.
#' @param dst_dir The destination directory where the downloaded datasets will be stored. If NULL, the function will return the path to the cached zip file.
#' @param extract Logical. If TRUE, the downloaded zip file will be extracted.
#' @param persistent Logical. If TRUE, the downloaded datasets will be cached at [datadir_persistent()] to speed up future calls to download_example_datasets().
#' If FALSE, the datasets will be cached at [datadir_temp()].
#' If NULL, the function will check both paths for the cached datasets but will use [datadir_temp()] if the cached file does not yet exist.
#' @param overwrite Logical. If TRUE, existing files with the same name in the destination directory will be overwritten.
#' @return The path to the downloaded (and possibly extracted) datasets.
#' @examples \dontrun{
#' download_example_datasets(dst_dir = ".", extract = TRUE, persistent = TRUE, overwrite = TRUE)
#' }
#' @seealso [datadir()]
#' @export
download_example_datasets <- function(dst_dir = NULL,
                                      extract = FALSE,
                                      persistent = NULL,
                                      overwrite = FALSE) {
    cached_zip <- cache_example_datasets(persistent = persistent, extract = extract)
    if (is.null(dst_dir)) {
        return(cached_zip)
    } else {
        zip <- file.path(dst_dir, basename(cached_zip))
        zip_is_missing <- !file.exists(zip)
        zip_has_wrong_size <- isTRUE(file.size(zip) != xds$zip_size)
        if (zip_is_missing || zip_has_wrong_size) file.copy(from = cached_zip, to = zip, overwrite = overwrite)
        if (extract) extract_example_datasets(zip)
        return(zip)
    }
}



# Helpers #####

#' @title Example Datasets Information
#' @description This list contains information about the example datasets provided for users to try the package.
#' The datasets can be downloaded from the provided URL.
#' @field url The URL from where the example datasets can be downloaded.
#' @field zip_size The expected size of the zipped example datasets file in bytes.
#' @field dir_size The expected size of the directory containing the example datasets in bytes.
#' @field n_files The expected total number of files in the example_datasets folder after extraction.
#' @examples
#' \dontrun{
#' print(xds$url)
#' print(xds$zip_size)
#' print(xds$dir_size)
#' print(xds$n_files)
#' }
#' @noRd
xds <- list(
    url = "https://github.com/spang-lab/metabodecon/releases/download/v1.0.2/example_datasets.zip",
    zip_size = 38253468,
    dir_size = 56378684,
    n_files = 1018
)

#' @title Cache Example Datasets
#' @description This function caches example datasets.
#' It checks if the datasets are already cached and if not, it downloads and caches them.
#' If the cached file has a wrong size, it can be overwritten based on the `overwrite` parameter.
#' The path to the cached datasets is returned.
#' @param persistent If TRUE, the datasets are cached permanently.
#' If FALSE, the datasets are cached temporarily.
#' If NULL, the path to the permanent zip file is returned if it exists and has the correct size else the path to the temporary zip file is returned.
#' @param extract Logical. If TRUE, the datasets are extracted after being cached. Default is TRUE.
#' @param overwrite Logical. If TRUE and the cached file exists but has a wrong size, the file is overwritten. If FALSE, an error is thrown in this case. Default is FALSE.
#' @return The path to the cached datasets.
#' @examples
#' \dontrun{
#' cache_example_datasets(persistent = FALSE, extract = FALSE, overwrite = FALSE)
#' }
#' @noRd
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
    if (extract) extract_example_datasets(zip)
    zip
}

extract_example_datasets <- function(path = datadir("example_datasets.zip")) {
    utils::unzip(zipfile = path, exdir = dirname(path))
}

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
    xds$url <- "https://github.com/spang-lab/metabodecon/releases/download/v1.0.2/example_datasets.zip"
    if (!is.null(copyfrom) && file.exists(copyfrom) && file.size(copyfrom) == xds$zip_size) {
        message(sprintf("Copying cached archive %s to %s", copyfrom, path))
        file.copy(copyfrom, path, overwrite = TRUE)
    } else {
        message(sprintf("Downloading %s as %s", xds$url, path))
        utils::download.file(xds$url, path, quiet = TRUE)
    }
    path
}
