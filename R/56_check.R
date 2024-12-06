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

draw_spectrum_check_arguments <- quote(
    stopifnot(
        is_num(foc_rgn, 2)  || is.null(foc_rgn),
        is_num(foc_frac, 2) || is.null(foc_frac),
        is_bool(foc_only),
        is_bool(add),
        is_num(fig, 4) || is.null(fig),
        is_num(mar, 4) || is.null(mar),
        is_str(main)   || is.null(main),
        is.list(si_line),  is.list(sm_line),  is.list(lc_lines),
        is.list(sp_line),  is.list(d2_line),
        is.list(cent_pts), is.list(bord_pts), is.list(norm_pts),
        is.list(bg_rect),  is.list(foc_rect), is.list(lc_rects),
        is.list(bt_axis),  is.list(lt_axis),  is.list(tp_axis) ,
        is.list(rt_axis)
    )
)

generate_lorentz_curves_type_checks <- quote(stopifnot(
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
