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

draw_spectrum_standardize_arguments <- quote({
    if (isFALSE(si_line)  || is.null(si_line))  si_line  <- list(show = FALSE)
    if (isFALSE(sm_line)  || is.null(sm_line))  sm_line  <- list(show = FALSE)
    if (isFALSE(lc_lines) || is.null(lc_lines)) lc_lines <- list(show = FALSE)
    if (isFALSE(sp_line)  || is.null(sp_line))  sp_line  <- list(show = FALSE)
    if (isFALSE(d2_line)  || is.null(d2_line))  d2_line  <- list(show = FALSE)
    if (isFALSE(cent_pts) || is.null(cent_pts)) cent_pts <- list(show = FALSE)
    if (isFALSE(bord_pts) || is.null(bord_pts)) bord_pts <- list(show = FALSE)
    if (isFALSE(norm_pts) || is.null(norm_pts)) norm_pts <- list(show = FALSE)
    if (isFALSE(bg_rect)  || is.null(bg_rect))  bg_rect  <- list(show = FALSE)
    if (isFALSE(foc_rect) || is.null(foc_rect)) foc_rect <- list(show = FALSE)
    if (isFALSE(lc_rects) || is.null(lc_rects)) lc_rects <- list(show = FALSE)
    if (isFALSE(bt_axis)  || is.null(bt_axis))  bt_axis  <- list(show = FALSE)
    if (isFALSE(lt_axis)  || is.null(lt_axis))  lt_axis  <- list(show = FALSE)
    if (isFALSE(tp_axis)  || is.null(tp_axis))  tp_axis  <- list(show = FALSE)
    if (isFALSE(rt_axis)  || is.null(rt_axis))  rt_axis  <- list(show = FALSE)
})

generate_lorentz_curves_type_checks <- quote(
    stopifnot(
    is_existing_path(data_path) ||
        is_spectrum(data_path) ||
        is_spectra(data_path) ||
        is_ispec(data_path) ||
        is_ispecs(data_path),
    is_char(file_format, 1, "(bruker|jcampdx)"),
    is_bool(make_rds, 1) || is_char(make_rds, 1),
    is_int(expno, 1),
    is_int(procno, 1),
    is_int(nfit, 1),
    is_int(smopts, 2),
    is_num(delta, 1),
    is_num(sfr, 2),
    is_num(wshw, 1),
    is_bool(ask, 1),
    is_bool(force, 1),
    is_bool(verbose, 1),
    is_int(nworkers, 1)
))
