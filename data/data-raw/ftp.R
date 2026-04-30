# FTP helpers used by data-raw/data.R to download datasets from MetaboLights.
# Run `devtools::load_all()` before sourcing this file.

# Internal API #####

#' @noRd
#' @title List Entries in a Remote FTP Directory
#' @description
#' Lists entries in a remote FTP directory. Supports both HTML-style directory
#' listings containing `href` links and plain-text FTP listings.
#' @param url URL of the remote FTP directory.
#' @return Character vector of entry names (files and/or directories).
#' @examples
#' \dontrun{
#' aki_url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS24/"
#' head(ftp_list(aki_url), 20)
#' }
ftp_list <- function(url) {
    ftp_list_info(url)$entry
}

#' @noRd
#' @title Download from FTP
#' @description
#' Downloads data from FTP into a local directory. By default this downloads a
#' single file. If `recursive = TRUE`, the function walks the remote directory
#' tree and downloads all files while preserving the remote subdirectory layout.
#' @param url URL of the remote FTP source (file or directory).
#' @param dest Path to the local destination directory.
#' @param recursive Logical indicating whether to recurse into subdirectories.
#' @param overwrite Logical. If `FALSE`, existing non-empty files are skipped.
#' @param quiet Logical passed to [utils::download.file()].
#' @param retries Integer number of retries for each file transfer.
#' @param verbose Logical. If `TRUE`, prints progress messages while traversing
#' and downloading.
#' @return Named list with elements:
#' \describe{
#'   \item{dest}{Normalized destination directory path.}
#'   \item{n_downloaded}{Number of files downloaded in this call.}
#' }
#' @examples
#' \dontrun{
#' # Download the complete AKI dataset
#' url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS24/"
#' dest <- datadir("aki", persistent = TRUE)
#' ret <- ftp_download(url, dest, recursive = TRUE)
#' ret
#' }
ftp_download <- function(url,
                         dest,
                         recursive = FALSE,
                         overwrite = FALSE,
                         quiet = TRUE,
                         retries = 3L,
                         verbose = TRUE) {

    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    dest <- normalizePath(dest, winslash = "/", mustWork = FALSE)
    if (verbose) message(sprintf("FTP Source: %s", url))
    if (verbose) message(sprintf("Local Destination: %s", dest))

    if (!recursive) {
        dest <- file.path(dest, basename(sub("/$", "", url)))
        ok <- ftp_download_one(url, dest, overwrite, quiet, retries, verbose)
        if (verbose && ok) message(sprintf("Downloaded: %s", dest))
        return(list(dest = dest, n_downloaded = as.integer(ok)))
    }

    queue <- new.env(parent = emptyenv())
    head <- 1L
    tail <- 2L
    queue[["1"]] <- ""
    n <- 0L
    n_seen <- 0L

    while (head < tail) {
        rel_dir <- queue[[as.character(head)]]
        rm(list = as.character(head), envir = queue)
        head <- head + 1L

        dir_url <- paste0(url, rel_dir)
        info <- ftp_list_info(dir_url)
        if (!nrow(info)) next

        if (verbose) {
            pretty_dir <- if (nzchar(rel_dir)) rel_dir else "/"
            message(sprintf("Scanning %s (%d entries)", pretty_dir, nrow(info)))
        }

        for (k in seq_len(nrow(info))) {
            entry <- info$entry[k]
            is_dir <- isTRUE(info$is_dir[k])
            n_seen <- n_seen + 1L
            entry_dir <- sub("/$", "", entry)
            rel_entry <- paste0(rel_dir, entry_dir)
            if (is_dir) {
                queue[[as.character(tail)]] <- paste0(rel_entry, "/")
                tail <- tail + 1L
            } else {
                file_url <- paste0(url, rel_entry)
                out <- file.path(dest, rel_entry)
                if (ftp_download_one(file_url, out, overwrite, quiet, retries, verbose)) {
                    n <- n + 1L
                    if (verbose) {
                        message(sprintf("Downloaded %d files (last: %s)", n, rel_entry))
                    }
                }
            }
        }
    }

    if (verbose) message(sprintf("Done. Seen %d entries, downloaded %d files.", n_seen, n))
    list(dest = dest, n_downloaded = n)
}

# Helpers #####

#' @noRd
#' @title Download a Single Remote FTP File
#' @description
#' Downloads a single remote file via FTP to a local destination path with
#' retry logic.
#' @param url URL of the remote file.
#' @param dest Path to the local destination file.
#' @param overwrite Logical. If `FALSE`, existing non-empty destination files
#' are kept and no new download is attempted.
#' @param quiet Logical passed to [utils::download.file()].
#' @param retries Integer number of retries if the transfer fails.
#' @param verbose Logical. If `TRUE`, prints skip messages for pre-existing
#' destination files.
#' @return Logical. `TRUE` if a new file was downloaded successfully, otherwise
#' `FALSE`.
#' @examples
#' \dontrun{
#' url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS24/i_Investigation.txt"
#' dest <- tempfile(fileext = ".txt")
#' ftp_download_one(url, dest, overwrite = TRUE, quiet = TRUE, retries = 1)
#' }
ftp_download_one <- function(url,
                             dest,
                             overwrite = FALSE,
                             quiet = TRUE,
                             retries = 1L,
                             verbose = FALSE) {
    if (!overwrite && file.exists(dest) && isTRUE(file.info(dest)$size > 0)) {
        if (verbose) message(sprintf("Skipped existing: %s", dest))
        return(FALSE)
    }
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    for (j in seq_len(retries)) {
        status <- tryCatch(
            suppressWarnings(utils::download.file(
                url = url,
                destfile = dest,
                quiet = quiet,
                mode = "wb"
            )),
            error = function(e) 1
        )
        if (status == 0 && file.exists(dest) && isTRUE(file.info(dest)$size > 0)) {
            return(TRUE)
        }
        Sys.sleep(j)
    }
    FALSE
}

#' @noRd
#' @title Parse FTP Directory Listing
#' @description
#' Returns entry names and a directory flag for a remote FTP directory.
#' Supports HTML-style listings (via `href`) and Unix-like FTP listings
#' (directory rows start with `d`).
#' @param url URL of the remote FTP directory.
#' @return A data.frame with columns `entry` and `is_dir`.
#' @examples
#' \dontrun{
#' url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS24/"
#' head(ftp_list_info(url), 20)
#' }
ftp_list_info <- function(url) {
    txt <- tryCatch({
        con <- base::url(url, open = "rb")
        on.exit(close(con), add = TRUE)
        readLines(con, warn = FALSE)
    }, warning = function(w) character(0), error = function(e) character(0))

    txt <- trimws(txt)
    txt <- txt[nzchar(txt)]
    if (!length(txt)) {
        return(data.frame(entry = character(0), is_dir = logical(0), stringsAsFactors = FALSE))
    }

    has_href <- grepl("href=\"", txt)
    if (any(has_href)) {
        href <- sub(".*href=\"([^\"]+)\".*", "\\1", txt[has_href])
        href <- unique(href)
        href <- href[!href %in% c("../", "./", "/")]
        is_dir <- grepl("/$", href)
        entry <- sub("/$", "", href)
        return(data.frame(entry = entry, is_dir = is_dir, stringsAsFactors = FALSE))
    }

    entry <- sub("^.*\\s", "", txt)
    entry <- entry[!entry %in% c(".", "..")]
    dir_flag <- substr(txt, 1, 1) == "d"
    keep <- !sub("^.*\\s", "", txt) %in% c(".", "..")
    data.frame(entry = entry, is_dir = dir_flag[keep], stringsAsFactors = FALSE)
}

#' @noRd
#' @title Check Whether a Remote FTP Path is a Directory
#' @description
#' Checks whether a remote FTP URL points to a directory by attempting to list
#' entries below that URL.
#' @param url URL of the remote directory candidate.
#' @return Logical. `TRUE` if the URL appears to be a directory, else `FALSE`.
#' @examples
#' \dontrun{
#' ftp_is_remote_dir("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS24/FILES/")
#' ftp_is_remote_dir("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS24/i_Investigation.txt/")
#' }
ftp_is_remote_dir <- function(url) {
    nrow(ftp_list_info(url)) > 0
}
