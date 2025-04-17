#!/usr/bin/env Rscript

# PURPOSE: Check that all functions in the package have an author
# USAGE: Rscript check-authors.R [--missing-only]

ignore <- c(
    # Files to ignore
)
no_author_needed <- c(
    # Functions to ignore
    "generate_lorentz_curves", # decon.R
    "generate_lorentz_curves_sim", # decon.R
    "check_mdrb_deps", # mdrb.R
    "print.decon1", # class.R
    "print.decon2", # class.R
    "print.align", # class.R
    "print.spectra", # class.R
    "print.decons1", # class.R
    "print.decons2", # class.R
    "print.aligns", # class.R
    "print.ispec", # class.R
    "print.idecon", # class.R
    "print.rdecon", # class.R
    "print.ispecs", # class.R
    "print.idecons", # class.R
    "print.rdecons", # class.R
    "is_decon0", # class.R
    "is_decon1", # class.R
    "is_decon2", # class.R
    "is_align", # class.R
    "is_spectra", # class.R
    "is_decons0", # class.R
    "is_decons1", # class.R
    "is_decons2", # class.R
    "is_aligns", # class.R
    "is_spectrum_or_spectra", # class.R
    "is_ispec", # class.R
    "is_idecon", # class.R
    "is_rdecon", # class.R
    "is_ispecs", # class.R
    "is_idecons", # class.R
    "is_rdecons", # class.R
    "as_decon0", # class.R
    "as_decon1", # class.R
    "as_decon2", # class.R
    "as_spectra", # class.R
    "as_decons0", # class.R
    "as_decons1", # class.R
    "as_decons2", # class.R
    "read_spectra", # spectrum.R
    "convert_width", # utils.R
    NULL
)
args <- commandArgs(trailingOnly = TRUE)
missing_only <- "--missing-only" %in% args
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
}
