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
#' [datadir_persistent()] to speed up future calls to
#' `download_example_datasets()`. If FALSE, the datasets will be cached at
#' [datadir_temp()]. If NULL, the function will check both paths for the cached
#' datasets but will return [datadir_temp()] if the cached file does not yet
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
#' [datadir()]
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
#' [download_example_datasets()],
#' [datadir_persistent()],
#' [datadir_temp()]
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
#' @seealso [datadir()], [datadir_temp()]
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
#' of metabodecons temporary session directory [tmpdir()] plus additional path
#' normalization.
#'
#' @return Returns the path to the temporary data directory.
#'
#' @seealso [tmpdir()], [datadir()], [datadir_persistent()]
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
#' [datadir_temp()]
#' [datadir_persistent()]
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
#' Deprecated since metabodecon v1.2.0. Please use [datadir()] instead. See
#' examples below for usage.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param dataset_name Either `""`, `"test"`, `"blood"`, `"urine"` or `"aki"`.
#'
#' @param warn Whether to print a warning message when the example folders do
#' not yet exist, i.e. [download_example_datasets()] has not been called yet.
#'
#' @return Path to the directory storing the example files.
#'
#' @seealso [download_example_datasets()]
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
#' file is overwritten. If FALSE, an error is thrown in this case. Important: if
#' the zip exists and has the correct size, overwrite has no effect, i.e. the
#' file will not be downloaded again. If you really want to force a re-download,
#' delete the zip file manually before calling this function.
#'
#' Second effect (only relevant if `extract` is TRUE): If TRUE, the datasets are
#' extracted even if the corresponding folder already exists. If FALSE, files
#' are not extracted again.
#'
#' @return The path to the cached datasets.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' \donttest{
#' cache_example_datasets(persistent = FALSE, extract = FALSE, overwrite = FALSE)
#' }
cache_example_datasets <- function(persistent = NULL,
                                   extract = FALSE,
                                   overwrite = FALSE,
                                   silent = FALSE) {

    pzip <- zip_persistent()
    pzip_exists <- file.exists(pzip)
    pzip_is_missing <- !pzip_exists
    pzip_has_correct_size <- isTRUE(file.size(pzip) == xds$zip_size) # Implies existence
    pzip_has_wrong_size <- !pzip_has_correct_size

    tzip <- zip_temp()
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

    zip <- normalizePath(if (persistent) pzip else tzip, "/", mustWork = FALSE)
    dstdir <- file.path(dirname(zip), "example_datasets")
    if (extract) {
        if (dir.exists(dstdir) && overwrite) unlink(dstdir, recursive = TRUE)
        if (!dir.exists(dstdir)) extract_example_datasets(zip)
    }
    zip
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
cachedir <- function(persistent = NULL) {
    tcd <- file.path(tmpdir(), "cache") # temporary cache dir
    pcd <- file.path(tools::R_user_dir("metabodecon", "cache")) # persistent cache dir
    cd <- if (isTRUE(persistent) || (is.null(persistent) && dir.exists(pcd))) pcd else tcd
    ncd <- normalizePath(cd, "/", mustWork = FALSE)
    mkdirs(ncd)
}

# Sap (Public) #####

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
"sap"

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
make_sap <- function() {
    sap_01 <- simulate_spectrum(
        name   = "sap_01",
        cs     = (64:-63) / 10,
        x0     = c(2.81, 2.10, 0.047, -2.22),
        A      = c(2000,  500,  5000,   1500),
        lambda = c( 0.4,  0.3,   0.6,   0.5),
        noise  = round(rnorm(128, sd = 50)),
        fqref  = 6e8,
        pkr    = c(3.2, -3.2)
    )
    as_spectra(sap_01)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
update_sap <- function() {
    path <- file.path(pkg_file("example_datasets/bruker"), "sap")
    logf("Updating %s." , path)
    if (!get_yn_input("Continue?")) return(invisible())
    sap <- make_sap()
    save_spectra(sap, path, force = TRUE)
    usethis::use_data(sap, overwrite = TRUE)
    path
}

# Sim (Public) #####

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
"sim"

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' @author 2024-2025 Tobias Schmidt: initial version.
update_sim <- function(nworkers = 1) {
    path <- file.path(pkg_file("example_datasets/bruker"), "sim")
    logf("Overwrite is TRUE. Updating %s." , path)
    if (!get_yn_input("Continue?")) return(invisible())
    sim <- make_sim(nworkers = nworkers)
    save_spectra(sim, path, force = TRUE)
    usethis::use_data(sim, overwrite = TRUE)
    path
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
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' @author 2024-2025 Tobias Schmidt: initial version.
#' @examples
#' xdir <- download_example_datasets()
#' path <- file.path(xdir, "bruker/blood/blood_01")
#' spec <- read_spectrum(path)
#' decon <- deconvolute_spectrum(spec)
#' sim_params <- get_sim_params(decon)
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

# AKI (Public) #####

#' @noRd
#' @title Downloads the AKI dataset from the Metabolights Database
#' @description
#' Recursively downloads and extracts all files of the AKI dataset (MTBLS24)
#' from the MetaboLights database, if not already available locally.
#' @return Path to the local directory containing the AKI dataset.
#' @examples
#' \dontrun{
#' aki_dir <- download_aki_data()
#' }
download_aki_data_from_metabolights <- function() {
    url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS24/"
    dest <- datadir("aki_raw", warn = FALSE, persistent = TRUE)
    n_files <- length(dir(dest, recursive = TRUE))
    if (n_files >= 2046) {
        # Return early if all files have already been downloaded and extracted
        return(dest)
    }
    ftp_download(url, dest, recursive = TRUE)
    spectra_dir <- file.path(dest, "FILES")
    zips <- dir(spectra_dir, ".zip$", full.names = TRUE)
    lapply(zips, unzip, exdir = spectra_dir)
    dest
}

#' @noRd
#' @title Update AKI Example Dataset (Development Only)
#' @description
#' Downloads the AKI urine NMR dataset (`MTBLS24`) from MetaboLights and copies
#' the Bruker files required by metabodecon (`1r`, `procs`, `acqus`) to
#' `misc/example_datasets/bruker/aki`.
#'
#' For details about the AKI dataset see `Dataset.Rmd`.
#'
#' The full source dataset can be downloaded directly from MetaboLights:
#' <https://www.ebi.ac.uk/metabolights/MTBLS24>. The copy prepared here is a
#' convenience subset and includes `s_MTBLS24.txt` plus only the Bruker files
#' required to read spectra (`1r`, `procs`, `acqus`). Other spectra-related
#' outputs and early analysis/classifier result files from Zacharias et al.
#' (2012) are excluded to keep archive size and download time low.
#'
#' To make the dataset available to package users, the `example_datasets.zip` in
#' GitHub must be updated as well. See `update_example_datasets()` for details.
#'
#' @return Path to the updated AKI example dataset directory.
#' @author 2026 Tobias Schmidt: initial version.
update_aki <- function() {

    stopifnot(loaded_via_devtools())

    # Download and extract AKI raw data if not already available locally
    raw <- download_aki_data_from_metabolights()
    src <- file.path(raw, "FILES")
    dst <- file.path(pkg_file("misc/example_datasets/bruker"), "aki")
    if (dir.exists(dst)) {
        cont <- get_yn_input(sprintf("Overwrite %s ?", dst))
        if (!cont) return(invisible())
        unlink(dst, recursive = TRUE, force = TRUE)
    }

    # Copy required Bruker files from all AKI sample folders
    spectra <- list.files(src, pattern = "^AKI_8", full.names = TRUE)
    spectra <- spectra[file.info(spectra)$isdir]
    one_r_files <- file.path(spectra, "10/pdata/10/1r")
    procs_files <- file.path(spectra, "10/pdata/10/procs")
    acqus_files <- file.path(spectra, "10/acqus")
    src_files <- c(one_r_files, procs_files, acqus_files)
    dst_files <- gsub(src, dst, src_files)
    dst_dirs <- unique(dirname(dst_files))
    lapply(dst_dirs, mkdirs)
    copied_successful <- file.copy(src_files, dst_files, overwrite = TRUE)
    stopifnot(all(copied_successful))

    # Copy study metadata file
    src_meta <- file.path(raw, "s_MTBLS24.txt")
    dst_meta <- file.path(dst, "s_MTBLS24.txt")
    copied_successful <- file.copy(src_meta, dst_meta, overwrite = TRUE)
    stopifnot(copied_successful)

}

# Example Datasets (Development) #####

#' @noRd
#' @title Update All Example Datasets (Development Only)
#' @description
#' Rebuilds `example_datasets.zip` from `misc/example_datasets/` and prints the
#' values needed to update `xds` (Example Datasets Summary). The actual content
#' of `misc/example_datasets/` must be updated manually before calling this
#' function.
#'
#' The subdirectories `blood`, `test` and `urine` contain raw spectra files as
#' created by the NMR spectrometer in the Oefner Lab and are not generated by
#' any function in the package.
#'
#' The `aki` subdirectory contains raw spectra files downloaded from the
#' MetaboLights database and can be recreated by calling [update_aki()].
#'
#' @return The path to the updated `example_datasets.zip` file.
#'
#' @author 2026 Tobias Schmidt: initial version.
update_example_datasets <- function() {
    stopifnot(loaded_via_devtools())
    misc <- pkg_file("misc")
    src <- file.path(misc, "example_datasets")
    zip <- file.path(misc, "example_datasets.zip")
    if (file.exists(zip)) unlink(zip)
    logf("Creating %s", zip)
    withr::with_dir(misc, utils::zip(zip, "example_datasets", "-rq9X"))
    logf("Done")
    files <- dir(src, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
    sizes <- file.info(files)$size
    zip_size <- file.info(zip)$size
    dir_size <- sum(sizes, na.rm = TRUE)
    n_files <- length(files)
    logf("")
    logf("Next steps:")
    logf("")
    logf("1. Upload example_datasets.zip to GitHub")
    logf("2. Update the xds object in data.R as follows:")
    logf("")
    logf("    xds <- list(")
    logf("        url = NEW_URL [^1],")
    logf("        zip_size = %d,", zip_size)
    logf("        dir_size = %d,", dir_size)
    logf("        n_files = %d", n_files)
    logf("    )")
    logf("")
    logf("[^1] OLD_URL was %s", xds$url)
    logf("")
    invisible(zip)
}
