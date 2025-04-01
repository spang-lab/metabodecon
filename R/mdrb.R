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
#' `check_mdrb_deps()` returns a list with detailed information about the
#' installation status of mdrb and its dependencies.
#'
#' @param stop_on_fail
#' If TRUE, an error is thrown if the check fails, providing instructions on how
#' to install or upgrade mdrb.
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
#' ## cargo  cargo >= 1.78   TRUE Current: cargo 1.84.1 (66221abde 2024-11-19)
#' ## rustc  rustc >= 1.78   TRUE Current: rustc 1.84.1 (e71f9a9a9 2025-01-27)
#' ```
#'
#' @examples
#' check_mdrb()
#' check_mdrb_deps()
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
check_mdrb_deps <- function() {
    check <- c("R >= 4.2", "Rtools exist", "cargo >= 1.78", "rustc >= 1.78")
    df <- data.frame(check, passed = NA, comment = NA)
    rownames(df) <- c("r", "rtools", "cargo", "rustc")

    # Check R version
    r_cur <- package_version(getRversion())
    r_req <- package_version("4.2.0")
    df["r", "passed"] <- r_cur >= r_req
    df["r", "comment"] <- paste("Current: R", r_cur)

    # Check buildtools exist
    df["rtools", "passed"] <- pkgbuild::has_build_tools()
    df["rtools", "comment"] <- "Testcall: pkgbuild::has_build_tools()"

    # Add $HOME/.cargo/bin to PATH before checking cargo and rustc
    win <- .Platform$OS.type == "windows"
    home_cargo_bin <- paste(home(), ".cargo", "bin", sep = .Platform$file.sep)
    PATH <- Sys.getenv("PATH")
    PATH <- paste(PATH, home_cargo_bin, sep = if (win) ";" else ":")
    Sys.setenv(PATH = PATH)

    # Check cargo version
    cargo_out <- "not found in PATH"
    cargo_str <- "0.0.0"
    try(silent = TRUE, {
        local_options(warn = -1)
        cargo_out <- system("cargo --version", intern = TRUE, ignore.stderr = TRUE)
        cargo_str <- strsplit(cargo_out, " ")[[1]][2] # cargo_out == 'cargo x.y.z (commitsha yyyy-mm-dd)'
    })
    cargo_cur <- package_version(cargo_str)
    cargo_req <- package_version("1.78.0")
    df["cargo", "passed"] <- cargo_cur >= cargo_req
    df["cargo", "comment"] <- paste("Current:", cargo_out)

    # Check rustc version
    rustc_out <- "not found in PATH"
    rustc_str <- "0.0.0"
    try(silent = TRUE, {
        local_options(warn = -1)
        rustc_out <- system("rustc --version", intern = TRUE, ignore.stderr = TRUE)
        rustc_str <- strsplit(rustc_out, " ")[[1]][2] # rustc_out == 'rustc x.y.z (commitsha yyyy-mm-dd)'
    })
    rustc_cur <- package_version(rustc_str)
    rustc_req <- package_version("1.78.0")
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
#' via [remotes::install_github()]. If remotes is not available, it is installed
#' first via [utils::install.packages()]
#'
#' @param ask
#' Whether to ask for confirmation before attempting installation. Default is
#' TRUE.
#'
#' @param args_remotes
#' Additional arguments to pass to [utils::install.packages()] when attempting
#' installation of [remotes](https://remotes.r-lib.org/index.html).
#'
#' @param args_mdrb
#' Additional arguments to pass to [remotes::install_github()] when attempting
#' installation of mdrb.
#'
#' @param verbose Whether to print messages to the console. Default is TRUE.
#'
#' @return NULL. Called for side effect of installing the Rust backend.
#'
#' @examples
#' if (interactive()) try(install_mdrb())
install_mdrb <- function(ask = TRUE,
                         args_remotes = list(),
                         args_mdrb = list(),
                         verbose = TRUE) {
    if (check_mdrb()) {
        if (verbose) logf("mdrb is already installed.")
        return(invisible())
    }
    checks <- check_mdrb_deps()
    if (!all(checks$passed)) {
        tbl <- capture.output2(print(checks, row.names = FALSE, right = FALSE))
        msg <- paste(
            "cannot install rust backend due to missing system dependencies:\n",
            tbl, "",
            "Please install the missing dependencies and try again.",
            "For installation instructions, see the following links:\n",
            "R & Rtools:    https://cran.r-project.org/",
            "cargo & rustc: https://www.rust-lang.org/tools/install",
            sep = "\n"
        )
        stop(msg, call. = FALSE)
    }
    remotes_available <- requireNamespace("remotes", quietly = TRUE)
    pkg_vec <- if (remotes_available) "mdrb" else c("remotes", "mdrb")
    pkg_str <- paste(pkg_vec, collapse = " and ")
    pkg_word <- if (length(pkg_vec) == 1) "package" else "packages"
    msg <- "Proceeding will install the following %s: %s. Continue?"
    msg <- sprintf(msg, pkg_word, pkg_str)
    cont <- if (isTRUE(ask) && isFALSE(get_yn_input(msg))) return()
    args_remotes$pkgs <- "remotes"
    if (!remotes_available) do.call(utils::install.packages, args_remotes)
    args_mdrb$repo <- "spang-lab/mdrb"
    do.call(remotes::install_github, args_mdrb)
    invisible()
}

# Internal #####

get_mdrb_version <- function() {
    tryCatch(
        packageVersion("mdrb"),
        error = function(e) package_version("0.0.0")
    )
}