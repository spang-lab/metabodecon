#!/usr/bin/env Rscript

# PURPOSE: test installation of metabodecon using different methods
# USAGE: Rscript test-install.R <method>
# PARAMS: <method>: Either 'CRAN-Modern', 'CRAN-Old' or 'Github'

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript test-install.R {CRAN-Modern|CRAN-Old|Github}")
method <- args[1]
if (!method %in% c("CRAN-Modern", "CRAN-Old", "Github")) stop(sprintf("Invalid method: %s\nMust be 'CRAN-Modern', 'CRAN-Old' or 'Github'"), method)
inst_pkgs <- rownames(installed.packages())

remove_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("Removing package:", pkg, "\n")
    remove.packages(pkg)
  } else {
    cat("Already uninstalled:", pkg, "\n")
  }
}

install_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("Package already installed:", pkg, "\n")
  } else {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg)
  }
}

remove_rtools <- function() {
  if (.Platform$OS.type == "windows") {
    cat("Removing Rtools from PATH on Windows...\n")
    path <- Sys.getenv("PATH")
    rtools_paths <- grep("rtools|Rtools", strsplit(path, ";")[[1]], value = TRUE, ignore.case = TRUE)
    if (length(rtools_paths) > 0) {
      new_path <- path
      for (rpath in rtools_paths) {
        new_path <- gsub(paste0(rpath, ";"), "", new_path, fixed = TRUE)
        new_path <- gsub(paste0(";", rpath), "", new_path, fixed = TRUE)
        new_path <- gsub(rpath, "", new_path, fixed = TRUE)
      }
      Sys.setenv(PATH = new_path)
      cat("Rtools removed from PATH\n")
    } else {
      cat("No Rtools found in PATH\n")
    }
  } else {
    cat("Not on Windows, skipping Rtools removal\n")
  }
}

install_metabodecon <- function(method) {
  if (method == "CRAN-Modern") {
    pak::pkg_install("metabodecon")
  } else if (method == "CRAN-Old") {
    install.packages("BiocManager")
    BiocManager::install(c("MassSpecWavelet", "impute"))
    install.packages("metabodecon")
  } else if (method == "Github") {
    pak::pkg_install("spang-lab/metabodecon")
  }
}

verify_install <- function() {
  if (require("metabodecon", quietly = TRUE)) {
    cat("metabodecon installed successfully\n")
    return(0)
  } else {
    cat("metabodecon installation failed\n")
    return(1)
  }
}

main <- function() {

  cat("Removing previous installation of metabodecon\n")
  remove_pkg("metabodecon")

  cat("Removing previous installations of Bioconductor dependencies\n")
  remove_pkg("BiocManager") # Bioconductor dependency
  remove_pkg("MassSpecWavelet") # Bioconductor dependency
  remove_pkg("impute") # Bioconductor dependency

  cat("Removing previous installations of R-Universe dependencies\n")
  remove_pkg("mdrb") # R-Universe dependency

  cat("Removing RTools\n")
  remove_rtools()

  cat("Installing dev dependencies\n")
  install_pkg("pak")
  install_pkg("pkgbuild")

  cat("Installing metabodecon using method: ", method, "\n")
  install_metabodecon(method)

  cat("Verifying installation\n")
  verify_install()
}

quit(status = main())
