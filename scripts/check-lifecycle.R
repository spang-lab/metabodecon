#!/usr/bin/env Rscript

# PURPOSE: Check that all exported functions in the package have a defined lifecycle.
# USAGE: Rscript check-lifecycle.R [--missing-only]

file_ignore <- c()
fn_ignore <- c(
    # align.R
    "align", # stable
    "get_si_mat", # stable
    # class.R
    "print.spectrum", # s3
    "print.decon1", # s3
    "print.decon2", # s3
    "print.align", # s3
    "print.spectra", # s3
    "print.decons1", # s3
    "print.decons2", # s3
    "print.aligns", # s3
    "print.ispec", # s3
    "print.idecon", # s3
    "print.rdecon", # s3
    "print.ispecs", # s3
    "print.idecons", # s3
    "print.rdecons", # s3
    "is_spectrum", # stable
    "is_decon0", # stable
    "is_decon1", # stable
    "is_decon2", # stable
    "is_align", # stable
    "is_spectra", # stable
    "is_decons0", # stable
    "is_decons1", # stable
    "is_decons2", # stable
    "is_aligns", # stable
    "as_spectrum", # stable
    "as_decon0", # stable
    "as_decon1", # stable
    "as_decon2", # stable
    "as_spectra", # stable
    "as_decons0", # stable
    "as_decons1", # stable
    "as_decons2", # stable
    # data.R
    "download_example_datasets", # stable
    "metabodecon_file", # stable
    "datadir", # stable
    "datadir_persistent", # stable
    "datadir_temp", # stable
    "tmpdir", # stable
    # decon.R
    "deconvolute", # stable
    "generate_lorentz_curves_sim", # documented with "generate_lorentz_curves"
    # plot.R
    "plot_spectra", # stable
    # spectrum.R
    "read_spectrum", # stable
    "read_spectra", # stable
    "make_spectrum", # stable
    # test.R
    "evalwith", # stable
    # util.R
    "convert_pos", # stable
    "convert_width", # stable
    "width", # stable
    "tree", # stable
    "transp", # stable
    # zzz.R
    "aaa_Get_Started" # stable
)
args <- commandArgs(trailingOnly = TRUE)
missing_only <- "--missing-only" %in% args
fg <- list(reset="\033[0m", warn="\033[91m", ok="\033[92m", file="\033[94m")
if (!isatty(stdout())) fg <- list(reset="", warn="", ok="", file="")
r_dir <- if (dir.exists("./R/")) "./R/" else "..R/"
r_files <- dir(r_dir, full.names=TRUE)
r_basenames <- basename(r_files)
r_files <- r_files[!r_basenames %in% file_ignore]
for (r_file in r_files) {
    r_file_colored <- sprintf("%s%s%s", fg$file, basename(r_file), fg$reset)
    lines_all <- readLines(r_file, warn=FALSE)
    idx_fns <- grep("^[a-zA-Z0-9._]+ *<- *function\\(", lines_all)
    idx_exports <- grep("^#' @export", lines_all)
    idx_lifecycles <- grep("#'.*lifecycle::", lines_all)
    lines <- lines_all[sort(c(idx_fns, idx_exports, idx_lifecycles))]
    types <- rep(NA, length(lines))
    types[idx_fns] <- "function"
    types[idx_exports] <- "export"
    types[idx_lifecycles] <- "lifecycle"
    types <- types[!is.na(types)]
    n_ok <- n_miss <- 0
    for (i in seq_along(lines)) {
        type0 <- types[i]
        type1 <- if (i > 1) types[i-1] else ""
        type2 <- if (i > 2) types[i-2] else ""
        if (type0 != "function") next
        if (type1 != "export") next
        if (type1 == "lifecycle" && type2 != "export") next
        fn <- trimws(sub("\\s*<- function\\(.*", "", lines[i]))
        fn_ok <- (type1 == "lifecycle") || fn %in% fn_ignore
        if (fn_ok) n_ok <- n_ok + 1 else n_miss <- n_miss + 1
        if (fn_ok && missing_only) next
        state <- if (fn_ok) "[ok]" else "[missing]"
        color <- if (fn_ok) fg$ok else fg$warn
        state_colored <- sprintf("%s%s%s", color, state, fg$reset)
        cat(sprintf("%s %s %s\n", r_file_colored, state_colored, fn))
    }
    color <- if (n_miss == 0) fg$ok else fg$warn
    n_total <- n_ok + n_miss
    cat(sprintf("%s %s%d/%d%s\n", r_file_colored, color, n_ok, n_total, fg$reset))
}
