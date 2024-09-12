
#' @description Return number of spectrum to adjust all others or 0 if every spectrum should be adjusted individually.
#' @noRd
get_adjno <- function(gspecs, sfr, wshw, ask) {
    if (isFALSE(ask) || length(gspecs) == 1) {
        return(0)
    }
    same_param <- get_yn_input("Use same parameters for all spectra?")
    if (!same_param) {
        return(0)
    }
    namestr <- paste(seq_along(gspecs), names(gspecs), sep = ": ", collapse = ", ")
    prompt <- sprintf("Number of spectrum for adjusting parameters? (%s)", namestr)
    get_num_input(prompt, min = 1, max = length(gspecs), int = TRUE)
}

#' @description Convert SFR to list of correct length and let user confirm.
#' @noRd
get_sfr <- function(gspecs, sfr, ask, adjno) {
    n <- length(gspecs)
    if (is_num(sfr, 2)) sfr <- rep(list(sfr), n)
    if (!is_list_of_nums(sfr, n, 2)) stop("sfr should be a [list of] num(2)")
    if (ask && adjno == 0) { # Different SFR for each spectrum.
        sfr <- lapply(seq_len(n), function(i) confirm_sfr(gspecs[[i]], sfr[[i]]))
    }
    if (ask && adjno >= 1) { # Same SFR for each spectrum.
        sfr_adjno <- confirm_sfr(gspecs[[adjno]], sfr[[adjno]])
        sfr <- rep(list(sfr_adjno), n)
    }
    names(sfr) <- names(gspecs)
    sfr
}

#' @description Convert WSHW to list of correct length and let user confirm.
#' @noRd
get_wshw <- function(gspecs, wshw, ask, adjno) {
    n <- length(gspecs)
    if (is_num(wshw, 1)) wshw <- rep(list(wshw), n)
    if (is_num(wshw, n)) wshw <- as.list(wshw)
    if (!is_list_of_nums(wshw, n, 1)) stop("wshw should be a [list of] num(1)")
    if (ask && adjno == 0) wshw <- mapply(confirm_wshw, gspecs, wshw)
    if (ask && adjno >= 1) wshw <- rep(list(confirm_wshw(gspecs[[adjno]], wshw[[adjno]])), n)
    names(wshw) <- names(gspecs)
    wshw
}

get_smopts <- function(gspecs, smopts) {
    n <- length(gspecs)
    if (is_int(smopts, 2)) smopts <- rep(list(smopts), n)
    if (!is_list_of_nums(smopts, n, 2)) stop("smopts should be a [list of] int(2)")
    names(smopts) <- names(gspecs)
    smopts
}

#' @description Repeatedly ask the user to confirm/refine SFR borders.
#' @noRd
confirm_sfr <- function(gspec, sfr = c(11.44494, -1.8828)) {
    plot_sfr(gspec, sfr[1], sfr[2])
    sfr_ok <- get_yn_input("Signal free region correctly selected?")
    while (!sfr_ok) {
        sfr[1] <- get_num_input("Choose another left border: [e.g. 12]", min = gspec$ppm_min, max = gspec$ppm_max)
        sfr[2] <- get_num_input("Choose another right border: [e.g. -2]", min = gspec$ppm_min, max = gspec$ppm_max)
        plot_sfr(gspec, sfr[1], sfr[2])
        sfr_ok <- get_yn_input("Signal free region correctly selected?")
    }
    sfr
}

#' @description Repeatly ask the user to confirm/refine the WSHW.
#' @noRd
confirm_wshw <- function(gspec, wshw) {
    plot_ws(gspec, wshw)
    ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    while (!ws_ok) {
        wshw <- get_num_input("Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154]")
        plot_ws(gspec, wshw)
        ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    }
    wshw
}

#' @description Takes the SFR in PPM and returns the SFR in PPM, DP and SDP.
#' @note Because the conversion from PPM to DP/SDP is slightly off (by 1-2 data points), the returned borders in DP and SDP are also incorrect. However, to maintain backwards compatibility with the old MetaboDecon1D function, these values are still returned. To only use the correct ppm values, set `version` to 2 in [filter_peaks()]. For details see `CHECK-2: signal free region calculation` in `TODOS.md`.
#' @noRd
enrich_sfr <- function(gspec, sfr) {
    left_ppm <- sfr[1]
    right_ppm <- sfr[2]
    left_dp <- (gspec$n + 1) - (gspec$ppm_max - left_ppm) / gspec$ppm_nstep
    left_sdp <- left_dp / gspec$sf[1]
    right_dp <- (gspec$n + 1) - (gspec$ppm_max - right_ppm) / gspec$ppm_nstep
    right_sdp <- right_dp / gspec$sf[1]
    named(left_ppm, right_ppm, left_dp, right_dp, left_sdp, right_sdp)
}

#' @description Calculates the WSR in dp and ppm from the WSHW in ppm.
#' @note Because the conversion from PPM to DP is slightly off (by 1-2 data points), the returned WSR in dp is also incorrect. However, to maintain backwards compatibility with the old MetaboDecon1D function, these values are still returned. To only use the correct ppm values, set `version` to 2 in [rm_water_signal()]. For details see `CHECK-3: water signal calculation` in `TODOS.md`.
#' @noRd
enrich_wshw <- function(gspec, wshw) {
    if (!is_gspec(gspec)) stop("Input must be a `gspec`, not ", class(gspec))
    gspec <- as_gspec(gspec)
    hwidth_ppm <- wshw
    hwidth_dp <- hwidth_ppm / gspec$ppm_nstep
    center_dp <- gspec$n / 2
    right_dp <- center_dp + hwidth_dp
    left_dp <- center_dp - hwidth_dp
    center_ppm <- gspec$ppm[center_dp]
    right_ppm <- gspec$ppm[right_dp]
    left_ppm <- gspec$ppm[left_dp]
    if (left_dp <= 1 || right_dp >= gspec$n) stop("WSR is out of range")
    named(left_ppm, right_ppm, center_ppm, hwidth_ppm, left_dp, right_dp, center_dp, hwidth_dp)
}
