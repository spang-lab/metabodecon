check_args <- function(func,
                       envir = parent.frame(),
                       asenv = TRUE) {
    name <- deparse(substitute(func))
    check <- check[[name]]
    if (is.null(check)) stop(sprintf("No check defined for [%s()]", name))
    check(envir)
    if (asenv) return(invisible(envir))
    keys <- names(formals(func))
    sapply(keys, get, envir = envir, simplify = FALSE)
}

check_args_generate_lorentz_curves <- function(e = parent.frame()) {
    x <- e$data_path
    stopifnot(is_char(x, 1) || is_spectrum(x) || is_spectra(x) || is_gspec(x) || is_gspecs(x))
    stopifnot(is_char(e$file_format, 1, "(bruker|jcampdx)"))
    stopifnot(is_bool(e$make_rds, 1) || is_char(e$make_rds, 1))
    stopifnot(is_int(e$expno, 1))
    stopifnot(is_int(e$procno, 1))
    stopifnot(is_int(e$nfit, 1))
    stopifnot(is_int(e$smopts, 2))
    stopifnot(is_num(e$delta, 1))
    stopifnot(is_num(e$sfr, 2))
    stopifnot(is_num(e$wshw, 1))
    stopifnot(is_bool(e$ask, 1))
    stopifnot(is_bool(e$force, 1))
    stopifnot(is_bool(e$verbose, 1))
    stopifnot(is_int(e$nworkers, 1))
}

# checks #####
check <- list()

## rm_water_signal #####
check$rm_water_signal <- function(e) {
    stopifnot(is_gspec(e$x))
    stopifnot(is_num(e$bwc, 1))
}

## deconvolute_gspec #####
check$deconvolute_gspec <- function(e) {
    stopifnot(is_gspec(e$gspec))
    stopifnot(is_int(e$nfit, 1))
    stopifnot(is_int(e$smopts, 2))
    stopifnot(is_num(e$delta, 1))
    stopifnot(is_num(e$sfr, 2))
    stopifnot(is_num(e$wshw, 1))
    stopifnot(is_bool(e$force, 1))
    stopifnot(is_num(e$bwc, 1))
    stopifnot(is_char(e$rtyp, 1, "(gdecon|decon0|decon1|decon2)"))
}

## deconvolute_gspecs #####
check$deconvolute_gspecs <- function(e) {
    n <- length(e$gspecs)
    stopifnot(is_gspecs(e$gspecs) || is_gspec(e$gspecs))
    stopifnot(is_int(e$nfit, 1))
    stopifnot(is_int(e$smopts, 2))
    stopifnot(is_num(e$delta, 1))
    stopifnot(is_num(e$sfr, 2) || is_list_of_nums(sfr, n, 2))
    stopifnot(is_num(e$wshw, 1) || is_list_of_nums(wshw, n, 1))
    stopifnot(is_bool(e$ask, 1))
    stopifnot(is_bool(e$force, 1))
    stopifnot(is_bool(e$verbose, 1))
    stopifnot(is_num(e$bwc, 1))
    stopifnot(is_char(e$rtyp, 1, "(gdecons|decons0|decons1)"))
}

## filter_peaks #####
check$filter_peaks <- function(e) {
    stopifnot(is_gspec(e$gspec))
}
