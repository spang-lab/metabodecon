# Helpers of generate_lorentz_curves #####

#' @title Get the spectrum number for adjusting parameters
#' @description Asks the user which spectrum should be used for adjusting signal free region (SFR) and water signal half width (WSHW) and returns the corresponding spectrum number. If `ask` is `FALSE` or the user chooses to use different parameters for each spectrum, 0 is returned.
#' @return The spectrum number chosen by the user for adjusting parameters, or 0 if `ask` is `FALSE` or the user chooses to use different parameters for each spectrum.
#' @noRd
get_adjno <- function(spectra, sfr, wshw, ask) {
    if (!ask || length(spectra) == 1) {
        return(0)
    }
    same_param <- get_yn_input("Use same parameters for all spectra?")
    if (!same_param) {
        return(0)
    }
    namestr <- paste(seq_along(spectra), names(spectra), sep = ": ", collapse = ", ")
    prompt <- sprintf("Number of spectrum for adjusting parameters? (%s)", namestr)
    get_num_input(prompt, min = 1, max = length(spectra), int = TRUE)
}

#' @noRd
#' @title Confirm and add SFR to each spectrum
#' @description Confirms the signal free region (SFR) of each spectrum, converts it to different units and adds it to the respective spectrum object in the list of spectra.
#' @param ss A list of spectra.
#' @param sfr A numeric vector of length 2, giving the left and right boundaries of the signal free region in ppm, or a list of such vectors (one for each spectrum).
#' @param ask A logical value indicating whether to ask the user to confirm the SFR(s).
#' @param adjno The spectrum number for adjusting the SFR. If 0, the user will be asked to confirm/adjust the SFR for each spectrum. If > 0, the SFR of the spectrum with this number will be used for all spectra.
#' @return A list of spectra with the SFR added to each spectrum.
get_sfrs <- function(ss, sfr, ask, adjno) {
    stopifnot(is.numeric(adjno) && length(adjno) == 1)
    n <- length(ss)
    # styler: off
    sfrs <- if (is_list_of_nums(sfr, n, 2)) sfr
        else if (is_num(sfr, 2)) rep(list(sfr), n)
        else stop(sprintf("Argument `sfr` should be either\n- a numeric vector of length 2, giving the left and right boundaries of the signal free region in ppm or\n- a list of length %d of such vectors (one for each spectrum)", n))
    if (ask && adjno == 0) { # Different SFR for each spectrum.
        sfrs <- lapply(seq_len(n), function(i) confirm_sfr(ss[[i]], sfrs[[i]]))
    }
    if (ask && adjno != 0) { # Same SFR for each spectrum.
        sfr_adjno <- confirm_sfr(ss[[adjno]], sfrs[[adjno]])
        sfrs <- rep(list(sfr_adjno), n)
    }
    names(sfrs) <- names(ss)
    sfrs
    # styler: on
}

#' @title Confirm and add WSR to each spectrum
#' @description Confirms the water signal region (WSR) of each spectrum, converts it to different units and adds it to the respective spectrum object in the list of spectra.
#' @param ss A list of spectra.
#' @param wshw A numeric value giving the half width of the water artefact in ppm, or a vector of such values (one for each spectrum).
#' @param ask A logical value indicating whether to ask the user to confirm the WSHW(s).
#' @param adjno The spectrum number for adjusting the WSHW. If 0, the user will be asked to confirm/adjust the WSHW for each spectrum. If > 0, the WSHW of the spectrum with this number will be used for all spectra.
#' @return A list of ss with the WSR added to each spectrum.
#' @noRd
get_wsrs <- function(ss, wshw, ask, adjno) {
    n <- length(ss)
    # Check initial, user provided WSHW(s). Convert to list if necessary.
    wshws <- if (is_num(wshw, 1)) {
        rep(list(wshw), n)
    } else if (is_num(wshw, n) || is_list_of_nums(wshw, n, 1)) {
        as.list(wshw)
    } else {
        stop(sprintf("Argument `wshw` should be either\n- a single value, giving the half width of the water artefact in ppm or\n- a vector of length %d of such values (one for each spectrum)", n))
    }
    # Ask user to confirm/adjust initial WSHW(s).
    if (ask) {
        wshws <- if (adjno == 0) {
            lapply(seq_len(n), function(i) confirm_wshw(ss[[i]], wshws[[i]]))
        } else if (adjno > 0) {
            wshw_adjno <- confirm_wshw(ss[[adjno]], wshws[[adjno]])
            rep(list(wshw_adjno), n)
        }
    }
    # Calculate SFR in different units for each spectrum
    wsrs <- if (adjno == 0) { # Different WSR for each spectrum.
        mapply(convert_wsr, ss, wshws, SIMPLIFY = FALSE)
    } else { # Same WSR for each spectrum.
        wsr_adjno <- convert_wsr(ss[[adjno]], wshws[[adjno]])
        lapply(ss, function(s) wsr_adjno)
    }
    wsrs
}

# Helpers for get_sfrs / get_wsrs #####

#' @title Confirm signal free region (SFR)
#' @description Repeatedly asks the user to refine and/or confirm the borders of the SFR.
#' @return A vector containing the final left and right borders of the SFR in ppm.
#' @noRd
confirm_sfr <- function(x, sfr = c(11.44494, -1.8828)) {
    plot_sfr(x, sfr[1], sfr[2])
    sfr_ok <- get_yn_input("Signal free region correctly selected?")
    while (!sfr_ok) {
        sfr[1] <- get_num_input("Choose another left border: [e.g. 12]", min = x$ppm_min, max = x$ppm_max)
        sfr[2] <- get_num_input("Choose another right border: [e.g. -2]", min = x$ppm_min, max = x$ppm_max)
        plot_sfr(spec, sfr[1], sfr[2])
        sfr_ok <- get_yn_input("Signal free region correctly selected?")
    }
    sfr
}

#' @title Confirm water signal half width (WSHW)
#' @description Repeatly asks the user to refine and/or confirm the half width of the water artefact.
#' @param spec A list representing the spectrum.
#' @param wshw The initial half width in ppm.
#' @return The confirmed half width of the water artefact in ppm units.
#' @noRd
confirm_wshw <- function(spec, wshw) {
    plot_ws(spec, wshw)
    ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    while (!ws_ok) {
        wshw <- get_num_input("Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154]")
        plot_ws(spec, wshw)
        ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    }
    wshw
}

#' @title Convert SFR from PPM to DP/SDP
#' @description Converts the sigal free region (SFR) from parts per million (PPM) to data points (DP) and scaled data points (SDP). Note that the conversion to DP/SDP is spectrum specific and that the conversion is slightly off (by 1-2 data points) to maintain backwards compatibility with the old MetaboDecon1D function. For details see `CHECK-2: signal free region calculation` in `TODOS.md`.
#' @param x A list representing the spectrum.
#' @param sfr A vector of length 2 containing the left and right borders of the SFR in ppm.
#' @return The spectrum with added SFR info.
#' @noRd
convert_sfr <- function(x, sfr) {
    within(list(), {
        left_ppm <- sfr[1]
        right_ppm <- sfr[2]
        left_dp <- (x$n + 1) - (x$ppm_max - left_ppm) / x$ppm_nstep
        left_sdp <- left_dp / x$sf[1] # nolint: object_usage_linter
        right_dp <- (x$n + 1) - (x$ppm_max - right_ppm) / x$ppm_nstep
        right_sdp <- right_dp / x$sf[2] # nolint: object_usage_linter
    })
}

#' @title Convert the WSR from PPM to DP
#' @description Converts the water signal half width (WSHW) from parts per million (PPM) to data points (DP) and scaled data points (SDP) and calculates the corresponding left and right boundaries from it. Note that the conversion to DP/SDP is spectrum specific and that the conversion is slightly off (by 1-2 data points) to maintain backwards compatibility with the old MetaboDecon1D function. For details see `CHECK-3: water signal calculation` in `TODOS.md`.
#' @param x A list representing the spectrum.
#' @param wshw The half width of the water artefact in ppm.
#' @return The spectrum with added WSR info.
#' @noRd
convert_wsr <- function(x, wshw, version = 1) {
    if (!is_spectrum(x) && !is_gspec(x)) {
        stop("Input must be of class `spectrum` or `gspec`, not ", class(x))
    } else if (version == 1) {
        g <- as_gspec(x)
        hwidth_ppm <- wshw
        hwidth_dp <- hwidth_ppm / g$ppm_nstep # half width in dp
        center_dp <- g$n / 2 # center line in dp
        right_dp <- center_dp + hwidth_dp # right border in dp
        left_dp <- center_dp - hwidth_dp # left border in dp
        if (left_dp <= 1 || right_dp >= g$n) stop("WSR is out of range")
        center_ppm <- g$ppm[center_dp] # center in ppm # nolint: object_usage_linter.
        right_ppm <- g$ppm[right_dp] # right border in ppm # nolint: object_usage_linter.
        left_ppm <- g$ppm[left_dp] # left border in ppm # nolint: object_usage_linter.
    } else if (version == 2) {
        s <- if (is_gspec(x)) as_spectrum(x)
        center_ppm <- x$cs[1] - x$cs[length(x$cs)]
        right_ppm <- center_ppm + wshw
        left_ppm <- center_ppm - wshw
        if (left_ppm < min(s$cs) || right_ppm > max(s$cs)) stop("WSR is out of range")
    } else {
        stop("Argument `x` should be a spectrum or a list of spectra.")
    }
    wsr <- c(left_ppm, right_ppm)
}
