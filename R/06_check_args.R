check_args <- function(fn = caller(), env = parent.frame()) {
    browser()
    cn <- paste0("check_args_", fn)
    cf <- get(cn, envir = environment(check_args))
    cf(env)
}

check_args_generate_lorentz_curves <- function(e = parent.frame()) {
    stopifnot(
        is_char(data_path, 1) ||
        is_spectrum(data_path) || is_spectra(data_path) ||
        is_gspec(data_path) || is_gspecs(data_path)
    )
    stopifnot(is_char(file_format, 1, "(bruker|jcampdx)"))
    stopifnot(is_bool(make_rds, 1) || is_char(make_rds, 1))
    stopifnot(is_int(expno, 1))
    stopifnot(is_int(procno, 1))
    stopifnot(is_int(e$nfit, 1))
    stopifnot(is_int(e$smopts, 2))
    stopifnot(is_num(e$delta, 1))
    stopifnot(is_num(e$sfr, 2))
    stopifnot(is_num(e$wshw, 1))
    stopifnot(is_bool(e$bwc, 1))
    stopifnot(is_int(e$rm_ws_version, 1))
    stopifnot(is_bool(e$force, 1))
}

check_args_deconvolute_gspec <- function(e = parent.frame()) {
    stopifnot(is_gspec(e$gspec))
    stopifnot(is_int(e$nfit, 1))
    stopifnot(is_int(e$smopts, 2))
    stopifnot(is_num(e$delta, 1))
    stopifnot(is_num(e$sfr, 2))
    stopifnot(is_num(e$wshw, 1))
    stopifnot(is_bool(e$bwc, 1))
    stopifnot(is_int(e$rm_ws_version, 1))
    stopifnot(is_bool(e$force, 1))
}

check_args_deconvolute_gspecs <- function(e = parent.frame()) {
    stopifnot(is_gspecs(e$gspecs))
    stopifnot(is_int(e$nfit, 1))
    stopifnot(is_int(e$smopts, 2))
    stopifnot(is_num(e$delta, 1))
    stopifnot(is_num(e$sfr, 2) || is_list_of_nums(sfr, length(e$x), 2))
    stopifnot(is_num(e$wsr, 2) || is_list_of_nums(wsr, length(e$x), 2))
    stopifnot(is_bool(e$bwc, 1))
    stopifnot(is_int(e$rm_ws_version, 1))
    stopifnot(is_bool(e$force, 1))
}

check_args_filter_peaks <- function(e = parent.frame()) {
    stopifnot(is_gspec(e$gspec))
}
