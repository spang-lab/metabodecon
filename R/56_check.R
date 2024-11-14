check_args_deconvolute <- function(e = parent.frame()) {
    stopifnot(is_char(e$x, 1) || is_spectrum(e$x) || is_spectra(e$x) || is_ispec(e$x) || is_ispecs(e$x))
    stopifnot(is_int(e$nfit, 1) || is.null(e$nfit))
    stopifnot(is_int(e$smopts, 2) || is.null(e$smopts))
    stopifnot(is_num(e$delta, 1) || is.null(e$delta))
    stopifnot(is_num(e$sfr, 2) || is.null(e$sfr))
    stopifnot(is_num(e$wshw, 1) || is.null(e$wshw))
    stopifnot(is_bool(e$ask, 1))
    stopifnot(is_bool(e$force, 1))
    stopifnot(is_bool(e$verbose, 1))
    stopifnot(is_int(e$nworkers, 1))
    return(e)
}

check_args_generate_lorentz_curves <- function(e = parent.frame()) {
    with(e, {
        stopifnot(
            is_existing_path(data_path) ||
            is_spectrum(data_path) ||
            is_spectra(data_path) ||
            is_ispec(data_path) ||
            is_ispecs(data_path)
        )
        stopifnot(is_char(file_format, 1, "(bruker|jcampdx)"))
        stopifnot(is_bool(make_rds, 1) || is_char(make_rds, 1))
        stopifnot(is_int(expno, 1))
        stopifnot(is_int(procno, 1))
        stopifnot(is_int(nfit, 1))
        stopifnot(is_int(smopts, 2))
        stopifnot(is_num(delta, 1))
        stopifnot(is_num(sfr, 2))
        stopifnot(is_num(wshw, 1))
        stopifnot(is_bool(ask, 1))
        stopifnot(is_bool(force, 1))
        stopifnot(is_bool(verbose, 1))
        stopifnot(is_int(nworkers, 1))
    })
    return(e)
}

check_args_rm_water_signal <- function(e = parent.frame()) {
    stopifnot(is_ispec(e$x))
    stopifnot(is_num(e$bwc, 1))
    return(e)
}

check_args_deconvolute_ispec <- function(e = parent.frame()) {
    stopifnot(is_ispec(e$ispec))
    stopifnot(is_int(e$nfit, 1))
    stopifnot(is_int(e$smopts, 2))
    stopifnot(is_num(e$delta, 1))
    stopifnot(is_num(e$sfr, 2))
    stopifnot(is_num(e$wshw, 1))
    stopifnot(is_bool(e$force, 1))
    stopifnot(is_num(e$bwc, 1))
    return(as.list(e))
}

check_args_deconvolute_ispecs <- function(e = parent.frame()) {
    stopifnot(is_ispecs(e$ispecs) || is_ispec(e$ispecs))
    stopifnot(is_int(e$nfit, 1))
    stopifnot(is_int(e$smopts, 2))
    stopifnot(is_num(e$delta, 1))
    stopifnot(is_num(e$sfr, 2) || is_list_of_nums(e$sfr, length(e$ispecs), 2))
    stopifnot(is_num(e$wshw, 1) || is_list_of_nums(e$wshw, length(e$ispecs), 1))
    stopifnot(is_bool(e$ask, 1))
    stopifnot(is_bool(e$force, 1))
    stopifnot(is_bool(e$verbose, 1))
    stopifnot(is_num(e$bwc, 1))
    return(e)
}

check_args_filter_peaks <- function(e = parent.frame()) {
    stopifnot(is_ispec(e$ispec))
    return(e)
}
