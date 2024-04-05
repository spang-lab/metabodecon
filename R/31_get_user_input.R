# Private Helpers of [generate_lorentz_curves()] #####

#' @title Get the spectrum number for adjusting parameters
#' @description Asks the user which spectrum should be used for adjusting signal free region (SFR) and water signal half width (WSHW) and returns the corresponding spectrum number. If `ask` is `FALSE` or the user chooses to use different parameters for each spectrum, 0 is returned.
#' @param spectra A list of spectra.
#' @param sfr A numeric value or a list specifying the SFR for each spectrum.
#' @param wshw A numeric value or a list specifying the WSHW for each spectrum.
#' @param ask A logical value indicating whether to ask the user to confirm the spectrum number. Default is TRUE.
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

#' @title Confirm and add SFR to each spectrum
#' @description Confirms the signal free region (SFR) of each spectrum, converts it to different units and adds it to the respective spectrum object in the list of spectra.
#' @param spectra A list of spectra.
#' @param sfr A numeric vector of length 2, giving the left and right boundaries of the signal free region in ppm, or a list of such vectors (one for each spectrum).
#' @param ask A logical value indicating whether to ask the user to confirm the SFR(s).
#' @param adjno The spectrum number for adjusting the SFR. If 0, the user will be asked to confirm/adjust the SFR for each spectrum. If > 0, the SFR of the spectrum with this number will be used for all spectra.
#' @return A list of spectra with the SFR added to each spectrum.
#' @noRd
get_sfrs <- function(spectra, sfr, ask, adjno) {
    n <- length(spectra)
    stopifnot(is.numeric(adjno) && length(adjno) == 1)

    # Check initial, user provided SFR(s). Convert to list if necessary.
    sfrs <- if (is_num(sfr, 2)) {
        rep(list(sfr), n)
    } else if (is_list_of_nums(sfr, n, 2)) {
        sfr
    } else {
        stop(sprintf("Argument `sfr` should be either\n- a numeric vector of length 2, giving the left and right boundaries of the signal free region in ppm or\n- a list of length %d of such vectors (one for each spectrum)", n))
    }

    # Ask user to confirm/adjust initial SFR(s) if ask is TRUE..
    if (ask) {
        if (adjno == 0) { # Different SFR for each spectrum.
            sfrs <- lapply(seq_len(n), function(i) confirm_sfr(spectra[[i]], sfrs[[i]]))
        } else { # Same SFR for each spectrum.
            sfr_adjno <- confirm_sfr(spectra[[adjno]], sfrs[[adjno]])
            sfrs <- rep(list(sfr_adjno), n)
        }
    }

    # Add SFR in different units to each spectrum.
    if (adjno == 0) { # Different SFR for each spectrum.
        for (i in seq_len(n)) {
            spectra[[i]]$sfr <- convert_sfr(spectra[[i]], sfrs[[i]])
        }
    } else { # Same SFR for each spectrum.
        sfr_adjno <- convert_sfr(spectra[[adjno]], sfrs[[adjno]])
        for (i in seq_len(n)) spectra[[i]]$sfr <- sfr_adjno
    }
    spectra
}

#' @title Confirm and add WSR to each spectrum
#' @description Confirms the water signal region (WSR) of each spectrum, converts it to different units and adds it to the respective spectrum object in the list of spectra.
#' @param spectra A list of spectra.
#' @param wshw A numeric value giving the half width of the water artefact in ppm, or a vector of such values (one for each spectrum).
#' @param ask A logical value indicating whether to ask the user to confirm the WSHW(s).
#' @param adjno The spectrum number for adjusting the WSHW. If 0, the user will be asked to confirm/adjust the WSHW for each spectrum. If > 0, the WSHW of the spectrum with this number will be used for all spectra.
#' @return A list of spectra with the WSR added to each spectrum.
#' @noRd
get_wsrs <- function(spectra, wshw, ask, adjno) {
    n <- length(spectra)

    # Check initial, user provided WSHW(s). Convert to list if necessary.
    if (is_num(wshw, 1)) {
        wshws <- rep(list(wshw), n)
    } else if (is_num(wshw, n) || is_list_of_nums(wshw, n, 1)) {
        wshws <- as.list(wshw)
    } else {
        stop(sprintf("Argument `wshw` should be either\n- a single value, giving the half width of the water artefact in ppm or\n- a vector of length %d of such values (one for each spectrum)", n))
    }

    # Ask user to confirm/adjust initial WSHW(s).
    if (ask) {
        if (adjno == 0) {
            wshws <- lapply(seq_len(n), function(i) confirm_wshw(spectra[[i]], wshws[[i]]))
        } else if (adjno > 0) {
            wshw_adjno <- confirm_wshw(spectra[[adjno]], wshws[[adjno]])
            wshws <- rep(list(wshw_adjno), n)
        }
    }

    # Calculate WSR from WSHW for each spectrum.
    if (adjno == 0) { # Different WSR for each spectrum.
        for (i in seq_len(n)) {
            spectra[[i]]$wsr <- convert_wsr(spectra[[i]], wshws[[i]])
        }
    } else { # Same WSR for each spectrum.
        wsr_adjno <- convert_wsr(spectra[[adjno]], wshws[[adjno]])
        for (i in seq_len(n)) spectra[[i]]$wsr <- wsr_adjno
    }

    spectra
}

# Helpers for get_sfrs / get_wsrs #####

#' @title Confirm signal free region (SFR)
#' @description Repeatedly asks the user to refine and/or confirm the borders of the SFR.
#' @param spec A list containing the spectrum details including 'n', 'ppm_max', 'ppm_nstep', and 'sfx'.
#' @param sfr A vector of length 2 containing the initial left and right borders of the SFR in ppm.
#' @return A list containing the final left and right borders of the SFR in both ppm and dp units.
#' @noRd
confirm_sfr <- function(spec, sfr = c(11.44494, -1.8828)) {
    plot_sfr(spec, sfr[1], sfr[2])
    sfr_ok <- get_yn_input("Signal free region correctly selected?")
    while (!sfr_ok) {
        sfr[1] <- get_num_input("Choose another left border: [e.g. 12]", min = spec$ppm_min, max = spec$ppm_max)
        sfr[2] <- get_num_input("Choose another right border: [e.g. -2]", min = spec$ppm_min, max = spec$ppm_max)
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
#' @param spec A list representing the spectrum.
#' @param sfr A vector of length 2 containing the left and right borders of the SFR in ppm.
#' @return The spectrum with added SFR info.
#' @noRd
convert_sfr <- function(spec, sfr) {
    within(list(), {
        left_ppm <- sfr[1]
        right_ppm <- sfr[2]
        left_dp <- (spec$n + 1) - (spec$ppm_max - left_ppm) / spec$ppm_nstep
        left_sdp <- left_dp / spec$sfx # nolint: object_usage_linter
        right_dp <- (spec$n + 1) - (spec$ppm_max - right_ppm) / spec$ppm_nstep
        right_sdp <- right_dp / spec$sfx # nolint: object_usage_linter
    })
}

#' @title Convert the WSR from PPM to DP
#' @description Converts the water signal half width (WSHW) from parts per million (PPM) to data points (DP) and scaled data points (SDP) and calculates the corresponding left and right boundaries from it. Note that the conversion to DP/SDP is spectrum specific and that the conversion is slightly off (by 1-2 data points) to maintain backwards compatibility with the old MetaboDecon1D function. For details see `CHECK-3: water signal calculation` in `TODOS.md`.
#' @param spec A list representing the spectrum.
#' @param wshw The half width of the water artefact in ppm.
#' @return The spectrum with added WSR info.
#' @noRd
convert_wsr <- function(spec, wshw) {
    within(list(), {
        hwidth_ppm <- wshw
        hwidth_dp <- hwidth_ppm / spec$ppm_nstep # half width in dp
        center_dp <- spec$n / 2 # center line in dp
        right_dp <- center_dp + hwidth_dp # right border in dp
        left_dp <- center_dp - hwidth_dp # left border in dp
        center_ppm <- spec$ppm[center_dp] # center in ppm # nolint: object_usage_linter.
        right_ppm <- spec$ppm[right_dp] # right border in ppm # nolint: object_usage_linter.
        left_ppm <- spec$ppm[left_dp] # left border in ppm # nolint: object_usage_linter.
    })
}
