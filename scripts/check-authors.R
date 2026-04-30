#!/usr/bin/env Rscript

# PURPOSE: Check that all functions in the package have an author
# USAGE: Rscript check-authors.R [--missing-only] [--strict]

ignore <- c(
    # Files to ignore
)
no_author_needed <- c(
    # Functions to ignore
    "check_mdrb_deps", # mdrb.R
    # class.R
    "print.spectrum",
    "print.spectra",
    "format.spectrum",
    "format.spectra",
    "summary.spectrum",
    "summary.spectra",
    "c.spectrum",
    "is_spectrum",
    "is_spectra",
    "as_spectra",
    "as_decon2",
    "as_decons2",
    # mdm.R
    "print.mdm",
    "print.summary.mdm",
    "predict.mdm",
    "coef.mdm",
    "plot.mdm",
    "summary.mdm",
    # spec.R
    "read_spectra",
    # util.R
    "convert_width",
    NULL
)
args <- commandArgs(trailingOnly = TRUE)
missing_only <- "--missing-only" %in% args
strict <- "--strict" %in% args
total_missing <- 0
fg <- list(reset="\033[0m", warn="\033[91m", ok="\033[92m", file="\033[94m")
if (!isatty(stdout())) fg <- list(reset="", warn="", ok="", file="")
r_dir <- if (dir.exists("./R/")) "./R/" else "..R/"
r_files <- dir(r_dir, full.names=TRUE)
r_basenames <- basename(r_files)
r_files <- r_files[!r_basenames %in% ignore]
for (r_file in r_files) {
    r_file_colored <- sprintf("%s%s%s", fg$file, basename(r_file), fg$reset)
    lines_all <- readLines(r_file, warn=FALSE)
    pattern <- "^(#' @author|[a-zA-Z0-9._]+ *<- *function\\()"
    lines <- grep(pattern, lines_all, value=TRUE)
    n_ok <- n_author_missing <- 0
    for (i in seq_along(lines)) {
        line <- lines[i]
        if (grepl("^#' @author", line)) next
        line_before <- if (i == 1) "" else lines[i-1]
        fn <- trimws(sub("\\s*<- function\\(.*", "", line))
        fn_has_author <- grepl("^#' @author", line_before) || fn %in% no_author_needed
        if (fn_has_author) {
            n_ok <- n_ok + 1
            if (missing_only) next
        } else {
            n_author_missing <- n_author_missing + 1
        }
        state <- if (fn_has_author) "[ok]" else "[missing]"
        color <- if (fn_has_author) fg$ok else fg$warn
        state_colored <- sprintf("%s%s%s", color, state, fg$reset)
        cat(sprintf("%s %s %s\n", r_file_colored, state_colored, fn))
    }
    color <- if (n_author_missing == 0) fg$ok else fg$warn
    n_total <- n_ok + n_author_missing
    cat(sprintf("%s %s%d/%d%s\n", r_file_colored, color, n_ok, n_total, fg$reset))
    total_missing <- total_missing + n_author_missing
}
if (strict && total_missing > 0) {
    cat(sprintf("\n%d function(s) missing @author tag\n", total_missing))
    quit(status = 1)
}
