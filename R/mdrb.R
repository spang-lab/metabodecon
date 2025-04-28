# Public #####

#' @export
#'
#' @title Check Rust Backend Requirements
#'
#' @description
#' `check_mdrb()` returns a boolean indicating whether a suitable version of the
#' metabodecon Rust backend [mdrb](https://github.com/spang-lab/mdrb) is
#' currently installed.
#'
#' `check_mdrb_deps()` returns a list with information about the
#' installation status of mdrb system dependencies.
#'
#' @param stop_on_fail
#' If TRUE, an error is thrown if the check fails, providing instructions on how
#' to install or upgrade mdrb.
#'
#' @param type
#' The planned installation type. Either "binary" or "source". Binary
#' installation of `mdrb` requires R >= 4.2. Source installation requires
#' R >= 4.2, R Buildtools, cargo >= 1.80 and rustc >= 1.80. Parameter added
#' with metabodecon version 1.5.0.
#'
#' @details
#' `check_mdrb_deps()` uses [pkgbuild::has_build_tools] to determine whether
#' buildtools are available on the current system. On Windows,
#' [pkgbuild::has_build_tools] searches for RTools. On other operating systems
#' it looks for a suitable compiler. If this fails, the function tries to build
#' a package from source, to check whether buildtools are available. Depending
#' on the environment this can take several seconds. This step can be disabled
#' by setting option `"buildtools.check"` to `FALSE`, e.g. via
#' `options(buildtools.check = FALSE)`.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @return
#' `check_mdrb()` returns TRUE if a suitable version of mdrb is installed,
#' else FALSE.
#'
#' `check_mdrb_deps()` returns a data.frame with the following columns:
#'
#' - `check`: description of the checked dependency as a string
#' - `passed`: boolean indicating whether the check passed
#' - `comment`: free text string describing the check result
#'
#' The rownames of the dataframe are: `r`, `rtools`, `cargo` and `rustc` and
#' correspond to the checked dependencies.
#'
#' Example:
#'
#' ```R
#' check_mdrb_deps()
#'
#' ##                check passed                                      comment
#' ## r           R >= 4.2   TRUE                             Current: R 4.4.2
#' ## rtools  Rtools exist   TRUE     Tested with: pkgbuild::has_build_tools()
#' ## cargo  cargo >= 1.80   TRUE Current: cargo 1.84.1 (66221abde 2024-11-19)
#' ## rustc  rustc >= 1.80   TRUE Current: rustc 1.84.1 (e71f9a9a9 2025-01-27)
#' ```
#'
#' @examples
#' check_mdrb()
#' withr::with_options(
#'     list(buildtools.check = FALSE), # (1)
#'     check_mdrb_deps()
#'     # (1) # Use faster, but less reliable check. See 'Details'
#' )
check_mdrb <- function(stop_on_fail = FALSE) {
    stopifnot(is_bool(stop_on_fail, 1))
    mdrb_version <- get_mdrb_version()
    req_version <- package_version("0.0.1")
    mdrb_is_ok <- mdrb_version >= req_version
    if (mdrb_is_ok || !stop_on_fail) return(mdrb_is_ok)
    err_msg <- paste(sep = "\n",
        "Using the Rust backend requires mdrb >= 0.0.1.\n",
        "To install or upgrade mdrb run: install_mdrb()",
        "To check system requirements run: check_mdrb_deps()\n",
        "For more information see: https://github.com/spang-lab/mdrb"
    )
    stop(err_msg, call. = FALSE)
}

#' @export
#' @rdname check_mdrb
check_mdrb_deps <- function(type = "binary") {

    stopifnot(is_char(type, 1, "(binary|source)"))

    check <- c("R >= 4.2", "Rtools exist", "cargo >= 1.80", "rustc >= 1.80")
    df <- data.frame(check, passed = NA, comment = NA)
    rownames(df) <- c("r", "rtools", "cargo", "rustc")

    # Check R version
    logf("Check R version...")
    r_cur <- package_version(getRversion())
    r_req <- package_version("4.2.0")
    df["r", "passed"] <- r_cur >= r_req
    df["r", "comment"] <- paste("Current: R", r_cur)

    # Skip most checks if type is binary (added in version 1.5.0)
    if (type == "binary") {
        skip <- c("rtools", "cargo", "rustc")
        df[skip, "passed"]  <- TRUE
        df[skip, "comment"] <- "Not required for binary installation"
        return(df)
    }

    # Check buildtools exist
    logf("Check buildtools exist...")
    df["rtools", "passed"] <- pkgbuild::has_build_tools()
    df["rtools", "comment"] <- "Testcall: pkgbuild::has_build_tools()"

    # Add $HOME/.cargo/bin to PATH before checking cargo and rustc
    logf("Add $HOME/.cargo/bin to PATH before checking cargo and rustc...")
    win <- .Platform$OS.type == "windows"
    home_cargo_bin <- paste(home(), ".cargo", "bin", sep = .Platform$file.sep)
    PATH <- Sys.getenv("PATH")
    PATH <- paste(PATH, home_cargo_bin, sep = if (win) ";" else ":")
    Sys.setenv(PATH = PATH)

    # Check cargo version
    logf("Check cargo version...")
    cargo_out <- "not found in PATH"
    cargo_str <- "0.0.0"
    try(silent = TRUE, {
        local_options(warn = -1)
        cargo_out <- system("cargo --version", intern = TRUE, ignore.stderr = TRUE)
        cargo_str <- strsplit(cargo_out, " ")[[1]][2] # cargo_out == 'cargo x.y.z (commitsha yyyy-mm-dd)'
    })
    cargo_cur <- package_version(cargo_str)
    cargo_req <- package_version("1.80.0")
    df["cargo", "passed"] <- cargo_cur >= cargo_req
    df["cargo", "comment"] <- paste("Current:", cargo_out)

    # Check rustc version
    logf("Check rustc version...")
    rustc_out <- "not found in PATH"
    rustc_str <- "0.0.0"
    try(silent = TRUE, {
        local_options(warn = -1)
        rustc_out <- system("rustc --version", intern = TRUE, ignore.stderr = TRUE)
        rustc_str <- strsplit(rustc_out, " ")[[1]][2] # rustc_out == 'rustc x.y.z (commitsha yyyy-mm-dd)'
    })
    rustc_cur <- package_version(rustc_str)
    rustc_req <- package_version("1.80.0")
    df["rustc", "passed"] <- rustc_cur >= rustc_req
    df["rustc", "comment"] <- paste("Current:", rustc_out)

    df
}

#' @export
#'
#' @title Install Rust Backend
#'
#' @description
#' Installs metabodecon's Rust backend [mdrb](https://github.com/spang-lab/mdrb)
#' from [R-Universe](https://spang-lab.r-universe.dev/mdrb).
#'
#' lifecycle::badge("experimental")
#'
#' @param ask
#' Whether to ask for confirmation before attempting installation. Default is
#' TRUE.
#'
#' @param type
#' Preffered installation type. Either "binary", "source" or "both". Passed to
#' [utils::install.packages()].
#'
#' @param ...
#' Additional arguments to pass to [utils::install.packages()] when attempting
#' installation of mdrb.
#'
#' @return NULL. Called for side effect of installing the Rust backend.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' if (interactive()) try(install_mdrb())
install_mdrb <- function(ask = TRUE, type = "both", ...) {
    msg <- "Proceeding will install package 'mdrb'. Continue?"
    cont <- if (isTRUE(ask) && isFALSE(get_yn_input(msg))) return()
    invisible(install.packages(
        'mdrb',
        repos = 'https://spang-lab.r-universe.dev',
        type = type,
        ...
    ))
}

# Internal #####

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
get_mdrb_version <- function() {
    tryCatch(
        packageVersion("mdrb"),
        error = function(e) package_version("0.0.0")
    )
}