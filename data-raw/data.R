# Data-generation scripts for metabodecon package datasets.
#
# These scripts are for development use only and are NOT part of the package.
# Run `devtools::load_all()` before sourcing this file.
#
# To regenerate the sap and sim datasets:
#   devtools::load_all()
#   source("data-raw/ftp.R")
#   source("data-raw/data.R")
#   update_sap()     # regenerates data/sap.rda and inst/example_datasets/bruker/sap
#   update_sim()     # regenerates data/sim.rda and inst/example_datasets/bruker/sim
#   update_aki()     # downloads MTBLS24, copies to inst/example_datasets/bruker/aki
#   update_example_datasets()  # rebuilds misc/example_datasets.zip

# Sap #####

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
    logf("Updating %s.", path)
    if (!get_yn_input("Continue?")) return(invisible())
    sap <- make_sap()
    save_spectra(sap, path, force = TRUE)
    usethis::use_data(sap, overwrite = TRUE)
    path
}

# Sim #####

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
    logf("Overwrite is TRUE. Updating %s.", path)
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
#' see [metabodecon::deconvolute()].
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

# AKI #####

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

# Example Datasets #####

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
#' MetaboLights database and can be recreated by calling `update_aki()`.
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
    logf("2. Update the xds object in R/data.R as follows:")
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
