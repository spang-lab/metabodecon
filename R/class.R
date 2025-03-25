# S3 Public #####

#' @name print_methods
#' @rdname print_methods
#'
#' @title S3 Methods for Printing Metabodecon Objects
#'
#' @description
#' S3 Methods for printing metabodecon objects as described in the [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/).
#'
#' @param x
#' The object to print.
#'
#' @param name
#' Logical. If TRUE, the name of the object is printed before the object.
#'
#' @param ...
#' Not used. Only accepted to comply with generic [base::print()].
#'
#' @return
#' NULL, called for side effect of printing to the standard output device.
#'
#' @examples
#' print(sim[[1]])
#' print(sim[[1]], name = TRUE)
#' print(sim)
#' decon <- deconvolute(sim[[1]], sfr = c(3.55, 3.35))
#' print(decon)
NULL

#' @export
#' @rdname print_methods
print.spectrum <- function(x, name = FALSE, ...) {
    namestr <- if (name) paste0(x$meta$name %||% "NULL", ": ") else ""
    fmt <- "%sspectrum object (%d dp, %.1f to %.1f ppm)\n"
    catf(fmt, namestr, length(x$cs), max(x$cs), min(x$cs))
}

#' @export
#' @rdname print_methods
print.decon1 <- function(x, name = FALSE, ...) {
    ppm <- x$x_values_ppm
    n <- length(ppm)
    name <- if (name) paste0(x$filename %||% "NULL", ": ") else ""
    fmt <- "%sdecon1 object (%d dp, %.1f to %.1f ppm, %d peaks)\n"
    catf(fmt, name, n, max(ppm), min(ppm), length(x$A))
}

#' @export
#' @rdname print_methods
print.decon2 <- function(x, name = FALSE, ...) {
    name <- if (name) paste0(x$meta$name %||% "NULL", ": ") else ""
    fmt <- "%sdecon2 object (%d dp, %.1f to %.1f ppm, %d peaks)\n"
    catf(fmt, name, length(x$cs), max(x$cs), min(x$cs), length(x$lcpar$A))
}

#' @export
#' @rdname print_methods
print.align <- function(x, name = FALSE, ...) {
    name <- if (name) paste0(x$meta$name %||% "NULL", ": ") else ""
    fmt <- "%salign object (%d dp, %.1f to %.1f ppm, %d peaks)\n"
    catf(fmt, name, length(x$cs), max(x$cs), min(x$cs), length(x$lcpar$A))
}

#' @export
#' @rdname print_methods
print.spectra <- function(x, ...) {
    msg <- "spectra object consisting of %d spectrum objects:\n"
    catf(msg, length(x, ...))
    nams <- get_names(x, ...)
    msg <- "%s (%d datapoints from %.2f - %.2f ppm)\n"
    mapply(x, ..., nams, FUN = function(x, nam) {
        catf(msg, nam, length(x$si), min(x$cs), max(x$cs))
    })
    invisible(NULL)
}

#' @export
#' @rdname print_methods
print.decons1 <- function(x, ...) {
    catf("decons1 object with %s decon1 elements\n", length(x))
    invisible(sapply(x, print, name = TRUE))
}

#' @export
#' @rdname print_methods
print.decons2 <- function(x, ...) {
    catf("decons2 object with %s decon2 elements\n", length(x))
    invisible(sapply(x, print, name = TRUE))
}

#' @export
#' @rdname print_methods
print.aligns <- function(x, ...) {
    catf("aligns object with %s align elements\n", length(x))
    invisible(sapply(x, print, name = TRUE))
}

# S3 Private #####

#' @export
print.ispec <- function(x, name = FALSE, ...) {
    name <- {
        if (isTRUE(name)) paste0(get_name(x) %||% "NULL", ": ")
        else if (is.character(name)) paste0(name, ": ")
        else ""
    }
    fmt <- "%sispec object (%d dp, %.1f to %.1f ppm)\n"
    catf(fmt, name, length(x$ppm), max(x$ppm), min(x$ppm))
}

#' @export
print.idecon <- function(x, name = FALSE, ...) {
    name <- {
        if (isTRUE(name)) paste0(get_name(x) %||% "NULL", ": ")
        else if (is.character(name)) paste0(name, ": ")
        else ""
    }
    fmt <- "%sidecon object (%d dp, %.1f to %.1f ppm, %d peaks)\n"
    catf(fmt, name, length(x$ppm), max(x$ppm), min(x$ppm), length(x$lcr$A))
}

#' @export
print.rdecon <- function(x, name = FALSE, ...) {
    name <- {
        if (isTRUE(name)) paste0(get_name(x) %||% "NULL", ": ")
        else if (is.character(name)) paste0(name, ": ")
        else ""
    }
    catf("%srdecon object\n", name)
}

#' @export
print.ispecs <- function(x, ...) {
    catf("ispecs object with %s ispec elements\n", length(x))
    nams <- get_names(x)
    invisible(mapply(print, x, nams))
}

#' @export
#' @rdname print_methods
print.idecons <- function(x, ...) {
    catf("idecons object with %s idecon elements\n", length(x))
    nams <- get_names(x)
    invisible(mapply(print, x, nams))
}

#' @export
#' @rdname print_methods
print.rdecons <- function(x, ...) {
    catf("rdecons object with %s rdecon elements\n", length(x))
    nams <- get_names(x)
    invisible(mapply(print, x, nams))
}

`[.collection` <- function(x, i, ...) {
    result <- NextMethod("[")
    class(result) <- class(x)
    result
}

#' @export
`[.spectra` <- `[.collection`

#' @export
`[.decons0` <- `[.collection`

#' @export
`[.decons1` <- `[.collection`

#' @export
`[.decons2` <- `[.collection`

#' @export
`[.aligns` <- `[.collection`

#' @export
`[.ispecs` <- `[.collection`

#' @export
`[.idecons` <- `[.collection`

#' @export
`[.rdecons` <- `[.collection`

# Checks (Public) #####

#' @export
#'
#' @name is_metabodecon_class
#'
#' @title Is an Object from a Metabodecon Class?
#'
#' @description
#' Check if an object is an instance of a specific 'Metabodecon Class'. See
#' [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html) for a
#' list of classes.
#'
#' @param x
#' The object to check.
#'
#' @param check_class
#' Logical indicating whether to check the class of the object.
#'
#' @param check_contents
#' Logical indicating whether to check the contents of the object.
#'
#' @param check_child_classes
#' Logical indicating whether to check the class of each element of the object.
#'
#' @return
#' TRUE if the object is an instance of the specified class, otherwise FALSE.
#'
#' @examples
#' ss <- sim[1:2]
#' dd <- deconvolute(ss, sfr = c(3.55, 3.35))
#' aa <- align(dd)
#' s1 <- sim[[1]]
#' d1 <- dd[[1]]
#' a1 <- aa[[1]]
#'
#' is_spectrum(s1) # TRUE
#' is_spectrum(s1, check_contents = TRUE) # TRUE
#' is_decon0(d1) # FALSE
#' is_decon1(d1) # FALSE
#' is_decon2(d1) # TRUE
#' is_align(a1) # TRUE
#'
#' is_spectra(ss) # TRUE
#' is_decons0(dd) # FALSE
#' is_decons1(dd) # FALSE
#' is_decons2(dd) # TRUE
#' is_aligns(aa) # TRUE
#'
is_spectrum <- function(x,
                        check_class = TRUE,
                        check_contents = FALSE) {
    # styler: off
    if (check_class && !inherits(x, "spectrum")) return(FALSE)
    if (!check_contents) return(TRUE)
    if (!is.list(x)) return(FALSE)
    mandatory <- c("si", "cs")
    if (!all(mandatory %in% names(x))) return(FALSE)
    # styler: on
    return(TRUE)
}

#' @export
#' @rdname is_metabodecon_class
is_decon0 <- function(x) {
    is.list(x) && all(decon0_members_mandatory %in% names(x)) && !is_decon1(x)
}

#' @export
#' @rdname is_metabodecon_class
is_decon1 <- function(x) inherits(x, "decon1")

#' @export
#' @rdname is_metabodecon_class
is_decon2 <- function(x) inherits(x, "decon2")

#' @export
#' @rdname is_metabodecon_class
is_align <- function(x) inherits(x, "align")

#' @export
#' @rdname is_metabodecon_class
is_spectra <- function(x,
                       check_class = TRUE,
                       check_contents = FALSE,
                       check_child_classes = FALSE) {
    # styler: off
    if (check_class && !inherits(x, "spectra")) return(FALSE)
    if (check_child_classes && !all(sapply(x, is_spectrum))) return(FALSE)
    if (!check_contents) return(TRUE)
    if (!is.list(x)) return(FALSE)
    if (!all(sapply(x, is_spectrum, check_contents = TRUE))) return(FALSE)
    # styler: on
    return(TRUE)
}

#' @export
#' @rdname is_metabodecon_class
is_decons0 <- function(x) all(sapply(x, is_decon0))

#' @export
#' @rdname is_metabodecon_class
is_decons1 <- function(x) inherits(x, "decons1")

#' @export
#' @rdname is_metabodecon_class
is_decons2 <- function(x) inherits(x, "decons2")

#' @export
#' @rdname is_metabodecon_class
is_aligns <- function(x) inherits(x, "aligns")

# Checks (Private) #####

is_spectrum_or_spectra <- function(x) is_spectrum(x) || is_spectra(x)
is_ispec <- function(x) inherits(x, "ispec")
is_idecon <- function(x) inherits(x, "idecon")
is_rdecon <- function(x) inherits(x, "rdecon")
is_ispecs <- function(x) inherits(x, "ispecs")
is_idecons <- function(x) inherits(x, "idecons")
is_rdecons <- function(x) inherits(x, "rdecons")

# Convert (Public) #####

#' @export
#'
#' @name as_metabodecon_object
#' @rdname as_metabodecon_object
#'
#' @title Convert to a Metabodecon Object
#'
#' @description Convert a object to a Metabodecon object.
#'
#' @param x
#' The object to convert.
#'
#' @param sf
#' Scale factor used during Only required if `x` is a decon0 object.
#'
#' @param sfs
#' List of scale factors. Only required if `x` is a list of decon0 objects.
#'
#' @param spectrum,spectra
#' The `spectrum`/`spectra` object corresponding to `x` as returned by
#' [read_spectrum()] / [read_spectra]. Only required if `x` is a decon0 object.
#'
#' @param sfr,sfrs
#' `sfr` should be a vector specifying the borders of the signal free region.
#' `sfrs` should be a list of such vectors. Only required if `x` is a `decon0`
#' object where element `signal_free_region` is missing (or a `decons0` objected
#' containing such `decon0` objects).
#'
#' @param wshw,wshws
#' `wshw` should specify the half width of the water signal region. `wshws`
#' should be a list of such values. Only required if `x` is a `decon0` object
#' where element `range_water_signal_ppm` is missing (or a `decons0` objected
#' containing such `decon0` objects).
#'
#' @param bwc
#' Level of backwards compatibility. If `bwc == 0`, bug fixes introduced after
#' version 0.2.2 of Metabodecon are not used. If `bwc == 1`, new features
#' introduced after version 0.2.2 of Metabodecon (e.g. faster algorithms) are
#' not used. If `bwc == 2`, all bug fixes and features introduced after version
#' 0.2.2 are used. Support for `bwc == 0` will be removed in 'metabodecon v2.0'.
#'
#' @param optional
#' Logical. If `TRUE`, the two optional elements `signal_free_region` and
#' `range_water_signal_ppm` are included in the returned `decon0` object.
#'
#' @param nworkers
#' Number of workers for parallel processing.
#'
#' @return An object of the specified class.
#'
#' @examples
#' dirpath <- metabodecon_file("sim_subset")
#' spectra <- read_spectra(dirpath)
#' spectrum <- spectra[[1]]
#' decons1 <- generate_lorentz_curves_sim(spectra)
#' decon1 <- generate_lorentz_curves_sim(spectrum)
#' decon2 <- as_decon2(decon1)
as_spectrum <- function(x, sf = c(1e3, 1e6)) {
    if (is_spectrum(x)) {
        return(x)
    } else if (is_decon1(x)) {
        cs <- x$x_values_ppm
        si <- x$y_values_raw %||% (x$y_values * sf[2])
        name <- x$filename
        fq <- x$x_values_hz
        meta <- named(name, fq)
        obj <- named(cs, si, meta)
        return(structure(obj, class = "spectrum"))
    } else {
        msg <- "Converting %s to spectrum is not suppoorted"
        msg <- sprintf(msg, class(x)[1])
        stop(msg)
    }
}

#' @export
#' @rdname as_metabodecon_object
as_decon0 <- function(x,
                      sf = NULL,
                      spectrum = NULL,
                      optional = TRUE) {
    if (is_decon0(x)) return(x)
    y <- as_decon1(x)
    y <- unclass(y)
    y[if (optional) decon0_members else decon0_members_mandatory]
}

#' @export
#' @rdname as_metabodecon_object
as_decon1 <- function(x,
                      sf = c(1e3, 1e6),
                      spectrum = NULL,
                      sfr = NULL,
                      wshw = NULL,
                      bwc = 2) {
    if (is_decon0(x)) {
        if (is.null(sf)) stop("Please provide `sf`")
        if (is.null(spectrum)) stop("Please provide `spectrum`")
        # Define some shorthands
        fq <- spectrum$meta$fq
        si <- spectrum$si
        ssp <- as.numeric(x$spectrum_superposition)
        ppm <- x$x_values_ppm
        sdp <- x$x_values
        dp <- round(x$x_values * sf[1])
        ppm_nstep <- diff(range(ppm)) / (length(ppm))
        y <- x
        # Append optional elements if missing
        if (is.null(x[["signal_free_region"]])) {
            if (is.null(sfr)) stop("Please provide `sfr`")
            y[["signal_free_region"]] <- sfr_in_sdp_bwc(sfr, ppm, sf)
        }
        if (is.null(x[["range_water_signal_ppm"]])) {
            if (is.null(wshw)) stop("Please provide `wshw`")
            y[["range_water_signal_ppm"]] <- wshw
        }
        # Make sure elements are in correct order
        y <- y[decon0_members]
        # Calculate decon1 elements
        y$y_values_raw <- si
        y$x_values_hz <- fq
        y$mse_normed_raw <- mean(((si / sum(si)) - (ssp / sum(ssp)))^2)
        y$signal_free_region_ppm <- sfr %||% sfr_in_ppm_bwc(x[["signal_free_region"]], sdp, ppm)
        y$x_0_hz <- convert_pos(x$x_0, sdp, fq)
        y$x_0_dp <- convert_pos(x$x_0, sdp, dp)
        y$x_0_ppm <- convert_pos(x$x_0, sdp, ppm)
        y$A_hz <- convert_width(x$A, sdp, fq)
        y$A_dp <- convert_width(x$A, sdp, dp)
        y$A_ppm <- convert_width(x$A, sdp, ppm)
        y$lambda_hz <- convert_width(x$lambda, sdp, fq)
        y$lambda_dp <- convert_width(x$lambda, sdp, dp)
        y$lambda_ppm <- convert_width(x$lambda, sdp, ppm)
        class(y) <- "decon1"
    } else if (is_decon1(x)) {
        y <- x
    } else if (is_decon2(x)) {
        # Helper vars
        cs <- x$cs
        si <- x$si
        n <- length(si)
        dpn <- (n - 1):0
        sdp <- dpn / 1e3
        fq <- x$meta$fq
        cs_step <- width(cs) / (n - 1)
        dpn_step <- 1
        fq_step <- if (!is.null(fq)) width(fq) / (n - 1)
        sdp_step <- dpn_step / 1e3
        x0_ppm <- x$lcpar$x0
        A_raw_ppm <- x$lcpar$A
        lambda_ppm <- x$lcpar$lambda
        x0_dp <- convert_pos(x0_ppm, cs, dpn)
        x0_sdp <- convert_pos(x0_ppm, cs, sdp)
        x0_hz <- if (!is.null(fq)) convert_pos(x0_ppm, cs, fq)
        A_raw_dp <- A_raw_ppm * (dpn_step / cs_step)
        A_raw_sdp <- A_raw_ppm * (sdp_step / cs_step)
        A_raw_hz <- if (!is.null(fq)) A_raw_ppm * (fq_step / cs_step)
        A_sc_ppm <- A_raw_ppm / 1e6
        A_sc_dp <- A_raw_dp / 1e6
        A_sc_sdp <- A_raw_sdp / 1e6
        A_sc_hz <- if (!is.null(fq)) A_raw_hz / 1e6
        lambda_dp <- convert_width(lambda_ppm, cs, dpn)
        lambda_sdp <- convert_width(lambda_ppm, cs, sdp)
        lambda_hz <- if (!is.null(fq)) abs(convert_width(lambda_ppm, cs, fq))
        limits_sdp <- range(sdp)
        if (bwc < 2) limits_sdp[2] <- limits_sdp[2] + sdp_step
        integrals <- t(lorentz_int(x0_sdp, A_sc_sdp, lambda_sdp, limits = limits_sdp))
        # Outputs
        y <- structure(class = "decon1", .Data = list())
        y$number_of_files <- 1
        y$filename <- x$meta$name
        y$x_values <- seq.int(length(x$cs) - 1, 0, -1) / sf[1]
        y$x_values_ppm <- x$cs
        y$y_values <- x$sit$sm / 1e6
        y$spectrum_superposition <- t(x$sit$sup / 1e6)
        y$mse_normed <- x$mse$smnorm
        y$index_peak_triplets_middle <- x$peak$center
        y$index_peak_triplets_left <- x$peak$right # decon[01] has left and right inverted
        y$index_peak_triplets_right <- x$peak$left # decon[01] has left and right inverted
        y$peak_triplets_middle <- x$cs[x$peak$center]
        y$peak_triplets_left <- x$cs[x$peak$right] # decon[01] has left and right inverted
        y$peak_triplets_right <- x$cs[x$peak$left] # decon[01] has left and right inverted
        sdp <- ((length(x$cs) - 1):0) / sf[1]
        y$integrals <- integrals
        y$signal_free_region <- sfr_in_sdp_bwc(x$args$sfr, x$cs, sf)
        y$range_water_signal_ppm <- x$args$wshw
        y$A <- -A_sc_sdp
        y$lambda <- -lambda_sdp
        y$x_0 <- x0_sdp
        y$y_values_raw <- x$si
        y$x_values_hz <- if (!is.null(fq)) x$meta$fq
        y$mse_normed_raw <- x$mse$norm
        y$signal_free_region_ppm <- x$args$sfr
        y$x_0_hz <- if (!is.null(fq)) x0_hz
        y$x_0_dp <- x0_dp
        y$x_0_ppm <- x0_ppm
        y$A_hz <- if (!is.null(fq)) (A_sc_hz)
        y$A_dp <- -A_sc_dp
        y$A_ppm <- -A_sc_ppm
        y$lambda_hz <- if (!is.null(fq)) (lambda_hz)
        y$lambda_dp <- -lambda_dp
        y$lambda_ppm <- -lambda_ppm
    } else if (is_idecon(x)) {
        y <- structure(class = "decon1", .Data = list())
        y$number_of_files <- 1
        y$filename <- x$name
        y$x_values <- x$sdp
        y$x_values_ppm <- x$ppm
        y$y_values <- x$y_smooth
        y$spectrum_superposition <- t(s <- lorentz_sup(x = x$sdp, lcpar = x$lcr))
        y$mse_normed <- mean(((x$y_smooth / sum(x$y_smooth)) - (s / sum(s)))^2)
        y$index_peak_triplets_middle <- as.numeric(x$peak$center[x$peak$high])
        y$index_peak_triplets_left <- as.numeric(x$peak$right[x$peak$high]) # decon[01] has left and right inverted
        y$index_peak_triplets_right <- as.numeric(x$peak$left[x$peak$high]) # decon[01] has left and right inverted
        y$peak_triplets_middle <- x$ppm[x$peak$center[x$peak$high]]
        y$peak_triplets_left <- x$ppm[x$peak$right[x$peak$high]] # decon[01] has left and right inverted
        y$peak_triplets_right <- x$ppm[x$peak$left[x$peak$high]] # decon[01] has left and right inverted
        y$integrals <- t(x$lcr$integrals)
        y$signal_free_region <- sfr_in_sdp_bwc(x$args$sfr, x$ppm, sf)
        y$range_water_signal_ppm <- x$args$wshw
        y$A <- x$lcr$A
        y$lambda <- x$lcr$lambda
        y$x_0 <- x$lcr$w
        y$y_values_raw <- x$y_raw
        y$x_values_hz <- x$hz
        y$mse_normed_raw <- mean(((x$y_raw / sum(x$y_raw)) - (s / sum(s)))^2)
        y$signal_free_region_ppm <- x$args$sfr
        y$x_0_hz <- convert_pos(x$lcr$w, x$sdp, x$hz)
        y$x_0_dp <- convert_pos(x$lcr$w, x$sdp, x$dp)
        y$x_0_ppm <- convert_pos(x$lcr$w, x$sdp, x$ppm)
        y$A_hz <- convert_width(x$lcr$A, x$sdp, x$hz)
        y$A_dp <- convert_width(x$lcr$A, x$sdp, x$dp)
        y$A_ppm <- convert_width(x$lcr$A, x$sdp, x$ppm)
        y$lambda_hz <- convert_width(x$lcr$lambda, x$sdp, x$hz)
        y$lambda_dp <- convert_width(x$lcr$lambda, x$sdp, x$dp)
        y$lambda_ppm <- convert_width(x$lcr$lambda, x$sdp, x$ppm)
    } else {
        stopf("Converting %s to decon1 is not supported", class(x)[1])
    }
    y
}

#' @export
#' @rdname as_metabodecon_object
as_decon2 <- function(x, sf = c(1e3, 1e6), spectrum = NULL, sfr = NULL, wshw = NULL, bwc = 2) {
    if (is_decon0(x)) as_decon2.decon1(x = as_decon1(x, sf, spectrum, sfr, wshw, bwc))
    else if (is_decon1(x)) as_decon2.decon1(x)
    else if (is_decon2(x)) x
    else if (is_idecon(x)) as_decon2.idecon(x)
    else if (is_rdecon(x)) as_decon2.rdecon(x)
    else stop(sprintf("Converting %s to decon2 is not supported", class(x)[1]))
}

as_decon2.decon1 <- function(x, ...) {
    cs <- x$x_values_ppm
    si <- x$y_values_raw
    meta <- list(
        name = x$filename,
        fq = x$x_values_hz
    )
    args <- list(
        nfit = NA, smopts = NA, delta = NA,
        sfr = sfr_in_ppm_bwc(x$signal_free_region, x$x_values, x$x_values_ppm),
        wshw = x$range_water_signal_ppm,
        ask = NA, force = NA, verbose = NA, bwc = NA, nworkers = NA
    )
    sit <- data.frame(
        wsrm = NA, nvrm = NA,
        sm = x$y_values * 1e6,
        sup = x$spectrum_superposition[1, ] * 1e6
    )
    peak <- data.frame(
        left = x$index_peak_triplets_right, # decon[01] has left and right inverted
        center = x$index_peak_triplets_middle,
        right = x$index_peak_triplets_left
    )
    lcpar <- data.frame(
        x0 = x$x_0_ppm,
        A = -(x$A_ppm * 1e6),
        lambda = -(x$lambda_ppm)
    )
    mse <- list(
        raw = mse(si, sit$sup, normed = FALSE),
        norm = x$mse_normed_raw,
        sm = mse(sit$sm, sit$sup, normed = FALSE),
        smnorm = x$mse_normed
    )
    obj <- named(cs, si, meta, args, sit, peak, lcpar, mse)
    class(obj) <- "decon2"
    obj
}

as_decon2.idecon <- function(x, ...) {
    cs <- x$ppm
    si <- x$y_raw
    meta <- x$meta
    args <- x$args
    lcpar <- data.frame(
        x0 = convert_pos(x$lcr$w, x$sdp, x$ppm),
        A = -convert_width(x$lcr$A, x$sdp, x$ppm) * 1e6,
        lambda = -convert_width(x$lcr$lambda, x$sdp, x$ppm)
    )
    sit <- data.frame(
        wsrm = x$y_nows * 1e6,
        nvrm = x$y_pos * 1e6,
        sm = x$y_smooth * 1e6,
        sup = lorentz_sup(cs, lcpar = lcpar)
    )
    peak <- data.frame(
        left = x$peak$left[x$peak$high],
        center = x$peak$center[x$peak$high],
        right = x$peak$right[x$peak$high]
    )
    mse <- list(
        raw = mse(si, sit$sup, normed = FALSE),
        norm = mse(si, sit$sup, normed = TRUE),
        sm = mse(sit$sm, sit$sup, normed = FALSE),
        smnorm = mse(sit$sm, sit$sup, normed = TRUE)
    )
    obj <- named(cs, si, meta, args, sit, peak, lcpar, mse)
    class(obj) <- "decon2"
    obj
}

as_decon2.rdecon <- function(x, ...) {
    assert(is_rdecon(x))
    cs <- x$mdrb_spectrum$chemical_shifts()
    si <- x$mdrb_spectrum$intensities()
    meta <- x$spectrum$meta
    args <- x$args
    sup <- x$mdrb_decon$superposition_vec(cs)
    lcpar <- as.data.frame(x$mdrb_decon$lorentzians())
    mse_raw <- x$mdrb_decon$mse()
    mse_norm <- mse_raw / sum(sup)
    wsrm <- si # Rust backend doesn't remove the water signal
    nvrm <- si # Rust backend doesn't remove negative values
    spec <- list(y_pos = nvrm)
    reps <- args$smopts[[1]]
    size <- args$smopts[[2]]
    sm <- smooth_signals(spec, reps, size, verbose = FALSE)$y_smooth
    sup = sup
    sit <- named(wsrm, nvrm, sm, sup)
    mse <- list(raw = mse_raw, norm = mse_norm, sm = NA, smnorm = NA)
    peak <- get_peak(lcpar$x0, cs) # Should be provided directly by Rust backend in future versions
    obj <- named(cs, si, meta, args, sit, peak, lcpar, mse)
    class(obj) <- "decon2"
    obj
}

#' @export
#' @rdname as_metabodecon_object
#' @inheritParams read_spectra
as_spectra <- function(x,
                       file_format = "bruker",
                       expno = 10,
                       procno = 10,
                       raw = FALSE,
                       silent = TRUE,
                       force = FALSE) {
    if (is_spectrum(x)) {
        xx <- structure(list(x), class = "spectra")
        xx <- set_names(xx, get_names(xx))
    } else if (all(sapply(x, is_spectrum))) {
        xx <- structure(x, class = "spectra")
        xx <- set_names(xx, get_names(xx))
    } else if (is.character(x) && file.exists(x)) {
        xx <- read_spectra(x, file_format, expno, procno, raw, silent, force)
    } else {
        stop("Input must be a path, spectrum or list of spectrum objects")
    }
    xx
}

#' @export
#' @rdname as_metabodecon_object
as_decons0 <- function(x,
                       sfs = list(c(1e3, 1e6)),
                       spectra = list(NULL),
                       nworkers = 1) {
    if (is_decons0(x)) {
        return(x)
    } else if (is_decons1(x) || is_decons2(x) || is_idecons(x)) {
        decons0 <- mcmapply(as_decon0, x, sfs, spectra, nw = nworkers)
    } else if (is.list(x) && all(sapply(x, is_decon0))) {
        decons0 <- x
    } else {
        stop(paste(
            "Input must be a list of decon0 objects or a single object",
            "of type decons0, decons1, decons2 or idecons."
        ))
    }
    # Don't set names or class for decons0, as the original MetaboDecon1D
    # objects didn't have names or classes as well and we want to stay backwards
    # compatible. If someone wants to have names, they can use `decons1` or
    # `decons2` instead.
    n <- length(decons0)
    for (i in seq_len(n)) decons0[[i]]$number_of_files <- n
    decons0
}

#' @export
#' @rdname as_metabodecon_object
as_decons1 <- function(x,
                       sfs = list(c(1e3, 1e6)),
                       spectra = list(NULL),
                       sfrs = list(NULL),
                       wshws = list(NULL),
                       bwc = 2,
                       nworkers = 1) {
    if (is_decons1(x)) {
        return(x)
    } else if (is_decons0(x) || is_decons2(x) || is_idecons(x)) {
        decons1 <- mcmapply(as_decon1, x, sfs, spectra, sfrs, wshws, bwc, nw = nworkers)
    } else if (is.list(x) && all(sapply(x, is_decon1))) {
        decons1 <- x
    } else {
        stop(paste(
            "Input must be a list of decon1 objects or a single object",
            "of type decons0, decons1, decons2 or idecons."
        ))
    }
    names(decons1) <- get_names(x)
    class(decons1) <- "decons1"
    n <- length(decons1)
    for (i in seq_len(n)) decons1[[i]]$number_of_files <- n
    decons1
}

#' @export
#' @rdname as_metabodecon_object
as_decons2 <- function(x,
                       sfs = list(c(1e3, 1e6)),
                       spectra = list(NULL),
                       sfrs = list(NULL),
                       wshws = list(NULL),
                       bwc = 2,
                       nworkers = 1) {
    if (is_decons2(x)) {
        return(x)
    } else if (is_decons0(x) || is_decons1(x) || is_idecons(x)) {
        decons2 <- mcmapply(as_decon2, x, sfs, spectra, sfrs, wshws, bwc, nw = nworkers)
    } else if (is.list(x) && all(sapply(x, is_decon2))) {
        decons2 <- x
    } else {
        stop(paste(
            "Input must be a list of decon2 objects or a single object",
            "of type decons0, decons1, decons2 or idecons."
        ))
    }
    names(decons2) <- get_names(x)
    class(decons2) <- "decons2"
    decons2
}

# Convert (Private) #####

as_ispec <- function(x, sf = c(1e3, 1e6)) {
    if (is_ispec(x)) {
        return(x)
    }
    if (is_char(x)) x <- read_spectrum(x)
    s <- if (is_spectrum(x)) x else as_spectrum(x)
    y_raw <- s$si # Raw signal intensities
    y_scaled <- y_raw / sf[2] # Scaled signal intensities
    n <- length(y_raw) # Number of data points
    dp <- seq(n - 1, 0, -1) # Data point numbers
    sdp <- seq((n - 1) / sf[1], 0, -1 / sf[1]) # Scaled data point numbers [^1]
    ppm <- s$cs # Parts per million
    hz <- s$meta$fq # Frequency in Hz
    ppm_range <- diff(range(s$cs)) # Range of the chemical shifts in ppm.
    ppm_max <- max(s$cs) # Maximum chemical shift in ppm.
    ppm_min <- min(s$cs) # Minimum chemical shift in ppm.
    ppm_step <- ppm_range / (n - 1) # Step size calculated correctly.
    ppm_nstep <- ppm_range / n # Wrong, but backwards compatible [^2].
    name <- s$meta$name # Name of the spectrum
    meta <- s$meta # Other Metadata to the spectrum
    g <- locals()[ispec_members]
    structure(g, class = "ispec")
    # [^1]: Same as `dp / sf[1]`, but with slight numeric differences, so we
    #       stick with the old calculation method for backwards compatibility.
    # [^2]: Example: ppm = 11, 23, 35, 47 ==> ppm_step == 12, ppm_nstep ~= 10.6
    #       (not really useful, but we need it for backwards compatibility with
    #       MetaboDecon1D results)
}

as_idecon <- function(x) {
    if (is_idecon(x)) {
        x
    } else if (is_rdecon(x)) {
        y_raw <- x$spectrum$si
        ppm <- x$spectrum$cs
        hz <- x$spectrum$meta$fq
        name <- x$spectrum$meta$name
        meta <- x$spectrum$meta
        args <- x$args
        lcpar <- x$mdrb_decon$lorentzians()
        y_sup <- x$mdrb_decon$superposition_vec(ppm)
        n <- length(y_raw)
        dp <- seq(n - 1, 0, -1)
        sdp <- dp / 1e3
        sf <- c(1e3, 1e6)
        ppm_range <- width(ppm)
        ppm_max <- max(ppm)
        ppm_min <- min(ppm)
        ppm_step <- ppm_range / (n - 1)
        ppm_nstep <- ppm_range / n
        y_scaled <- y_raw / 1e6
        y_nows <- y_scaled
        y_pos <- abs(y_scaled)
        y_smooth <- smooth_signals(named(y_pos), args$smopts[[1]], args$smopts[[2]], verbose = FALSE)$y_smooth
        d <- (calc_second_derivative(y_smooth))
        A <- (-convert_width(lcpar$A, ppm, sdp) / 1e6)
        lambda <- (-convert_width(lcpar$lambda, ppm, sdp))
        w <- convert_pos(lcpar$x0, ppm, sdp)
        integrals <- (-A) * pi
        lcr <- named(A, lambda, w, integrals)
        center <- round(convert_pos(lcpar$x0, ppm, 1:n))
        npeaks <- length(center)
        left <- center - 1
        center <- center
        right <- center + 1
        score <- rep(999, npeaks)
        high <- rep(TRUE, npeaks)
        region <- rep("norm", npeaks)
        peak <- data.frame(left, center, right, score, high, region)
        Z <- lci <- lca <- NA
        structure(locals(without = "x")[idecon_members], class = "idecon")
    } else if (all(idecon_members %in% names(x))) {
        structure(x, class = "idecon")
    } else {
        stop("Input must have all elements listed in `idecon_members`")
    }
}

as_rdecon <- function(x) {
    assert( # (1)
        is_spectrum(x$spectrum),
        is.list(x$args),
        typeof(x$mdrb_spectrum) == "externalptr",
        typeof(x$mdrb_deconvr) == "externalptr",
        typeof(x$mdrb_decon) == "externalptr",
        class(x$mdrb_spectrum) == "Spectrum",
        class(x$mdrb_deconvr) == "Deconvoluter",
        class(x$mdrb_decon) == "Deconvolution"
    )
    stopifnot(length(x) == 5) # (1)
    structure(x, class = "rdecon")
    # (1) This function is private, so in theory it can never be called with
    # invalid arguments, as all public functions validate their inputs first.
    # Therefore, using assert for type checking is correct, as assert-checks are
    # deactivated when the package is loaded via library(), i.e. the
    # "production" code will run faster.
    #
    # However, in practice, it's very easy to run into nasty problems as soon as
    # assertions are disabled (e.g. when calling this function with invalid
    # arguments during unit testing). To prevent such scenarios, we include a
    # very tiny, super-fast sanity check here, using stopifnot. This check
    # will always run, even in production code, and might us save a lot of
    # headaches in the future.
}

#' @export
#' @rdname as_metabodecon_object
as_ispecs <- function(x, sf = c(1e3, 1e6)) {
    if (is_ispecs(x)) {
        return(x)
    }
    gg <- if (is_ispec(x)) {
        list(x)
    } else if (is_spectrum(x)) {
        list(as_ispec(x))
    } else if (is_spectra(x)) {
        lapply(x, as_ispec, sf = sf)
    } else {
        stop("Input must be ispec, ispecs or spectra, not ", class(x))
    }
    gg <- structure(gg, class = "ispecs")
    gg <- set_names(gg, get_names(gg))
    gg
}

#' @export
#' @rdname as_metabodecon_object
as_idecons <- function(x) {
    if (is_idecons(x)) return(x)
    stopifnot(is.list(x), all(sapply(x, is_idecon)))
    structure(x, class = "idecons")
}

as_rdecons <- function(x) {
    if (is_rdecons(x)) return(x)
    assert(is.list(x), all(sapply(x, is_rdecon)))
    structure(x, class = "rdecons")
}

# Convert a list of singlets to a collection
as_collection <- function(x, cls) {
    assert(
        is.list(x),
        is_char(cls, 1, "(decon[0-2]|idecon|rdecon)"),
        cls == "decon0" || all(sapply(x, class) == cls)
    )
    decons <- switch(cls,
        "decon0" = as_decons0(x),
        "decon1" = as_decons1(x),
        "decon2" = as_decons2(x),
        "idecon" = as_idecons(x),
        "rdecon" = as_rdecons(x),
        stop("Unsupported class: ", cls)
    )
}

# Convert a metabodecon singlet to a metabodecon "v1.2+" singlet
as_v12_singlet <- function(obj) {
    if (is_spectrum(obj)) obj
    else if (is_ispec(obj)) as_spectrum(obj)
    else if (is_idecon(obj)) as_decon2(obj)
    else if (is_rdecon(obj)) as_decon2(obj)
    else if (is_decon0(obj)) stop("decon0 objects are not supported. Convert with as_decon2.")
    else if (is_decon1(obj)) as_decon2(obj)
    else if (is_decon2(obj)) obj
    else if (is_align(obj)) obj
    else stop(sprintf("Objects of class %s are not supported.", class(obj)))
}

# Convert a metabodecon collection to a metabodecon "v1.2+" collection
as_v12_collection <- function(obj) {
    if (is_spectra(obj)) obj
    else if (is_ispecs(obj)) as_spectra(obj)
    else if (is_idecons(obj)) as_decons2(obj)
    else if (is_decons0(obj)) stop("decons0 objects are not supported. Convert with as_decons2.")
    else if (is_decons1(obj)) as_decons2(obj)
    else if (is_decons2(obj)) obj
    else if (is_aligns(obj)) obj
    else stop(sprintf("Objects of class %s are not supported.", class(obj)))
}

# Constructors (Private) #####

new_rdecon <- function(spectrum, args, mdrb_spectrum, mdrb_deconvr, mdrb_decon) {
    x <- named(spectrum, args, mdrb_spectrum, mdrb_deconvr, mdrb_decon)
    as_rdecon(x)
}

# Getters (Private) #####

#' @noRd
#' @title Returns the name of an iterable.
#' @param x An iterable object, e.g. a single metabodecon object.
#' @param default Default name if no name is found.
#' @return The name of the object as string or whatever is given as `default`.
#' @examples
#' s1 <- list()
#' s2 <- list(name = "foo")
#' s3 <- list(name = "foo", meta = list(name = "bar"))
#' get_name(s1) # ""
#' get_name(s2) # "foo"
#' get_name(s3) # "bar"
get_name <- function(x, default = "") {
    (if (is.list(x)) x$meta$name %||% x$name) %||% default
}

#' @noRd
#' @title Returns the names of a metabodecon collection object.
#' @param x A metabodecon collection object.
#' @param default Default names if no names are found. Passed on to `get_default_names`.
#' @return A character vector of names.
#' @examples
#' s1 <- list()
#' s2 <- list(name = "foo")
#' s3 <- list(name = "foo", meta = list(name = "bar"))
#'
#' get_names(list(s1, s1)) # c("spectrum_1", "spectrum_2")
#' get_names(list(s1, myspec = s1)) # c("spectrum_1", "myspec")
#' get_names(list(s1, myspec = s2)) # c("spectrum_1", "foo")
#' get_names(list(s1, myspec = s3)) # c("spectrum_1", "bar")
get_names <- function(x, default = "spectrum_%d") {
    obj_names <- sapply(x, get_name, "")
    obj_names_empty <- obj_names == ""
    if (any(obj_names_empty)) {
        list_names <- names(x) %||% rep("", length(x))
        list_names_empty <- list_names == ""
        if (any(list_names_empty)) {
            default_names <- get_default_names(x, default)
            list_names[list_names_empty] <- default_names[list_names_empty]
        }
        obj_names[obj_names_empty] <- list_names[obj_names_empty]
    }
    names(obj_names) <- NULL
    obj_names
}

get_default_names <- function(x, default) {
    if (length(default) == 1 && grepl("%d", default)) {
        return(sprintf(default, seq_along(x)))
    }
    if (length(unique(default)) == length(x)) {
        return(default)
    }
    stop(paste(
        "Default names must be a single string with a `%d` placeholder",
        "or a character vector of unique spectrum names."
    ))
}

# Temporary Helper Function until metabodecon-rust/mdrb offer a method to
# retrieve peak selection results from the mdrb_decon element.
get_peak <- function(x0, cs, target = "decon2") {
    assert(
        is_num(x0),
        is_num(cs),
        is_char(target, 1, "(idecon|decon2)")
    )
    center <- round(convert_pos(x0, cs, seq_along(cs)))
    npeaks <- length(center)
    left <- center - 1
    center <- center
    right <- center + 1
    if (target == "idecon") {
        score <- rep(999, npeaks)
        high <- rep(TRUE, npeaks)
        region <- rep("norm", npeaks)
        peak <- data.frame(left, center, right, score, high, region)
    } else {
        peak <- data.frame(left, center, right)
    }
    peak
}

# Setters #####

set_names <- function(x, nams) {
    assert(is.list(x))
    has_names <- all(sapply(x, function(e) "name" %in% names(e)))
    has_meta_names <- all(sapply(x, function(e) "name" %in% names(e$meta)))
    names(x) <- nams
    if (has_names) for (i in seq_along(x)) x[[i]]$name <- nams[[i]]
    if (has_meta_names) for (i in seq_along(x)) x[[i]]$meta$name <- nams[[i]]
    x
}

# Members (Private) #####

spectrum_members <- c("cs", "si", "meta")

ispec_members <- c(
    "y_raw",
    "y_scaled",
    "n",
    "dp",
    "sdp",
    "sf",
    "ppm",
    "hz",
    "ppm_range",
    "ppm_max",
    "ppm_min",
    "ppm_step",
    "ppm_nstep",
    "name",
    "meta"
)

idecon_members <- c(
    ispec_members,
    "args",
    "y_nows",
    "y_pos",
    "Z",
    "y_smooth",
    "d",
    "peak",
    "lci",
    "lca",
    "lcr"
)

rdecon_members <- c(
    "spectrum",
    "args",
    "mdrb_spectrum",
    "mdrb_deconvr",
    "mdrb_decon"
)

decon0_members <- c(
    "number_of_files",
    "filename",
    "x_values",
    "x_values_ppm",
    "y_values",
    "spectrum_superposition",
    "mse_normed",
    "index_peak_triplets_middle",
    "index_peak_triplets_left",
    "index_peak_triplets_right",
    "peak_triplets_middle",
    "peak_triplets_left",
    "peak_triplets_right",
    "integrals",
    "signal_free_region",
    "range_water_signal_ppm",
    "A",
    "lambda",
    "x_0"
)

decon0_members_optional <- c(
    "signal_free_region",
    "range_water_signal_ppm"
)

decon0_members_mandatory <- setdiff(
    decon0_members,
    decon0_members_optional
)

decon1_members <- c(
    decon0_members,
    "y_values_raw",
    "x_values_hz",
    "mse_normed_raw",
    "signal_free_region_ppm",
    "x_0_hz",
    "x_0_dp",
    "x_0_ppm",
    "A_hz",
    "A_dp",
    "A_ppm",
    "lambda_hz",
    "lambda_dp",
    "lambda_ppm"
)

decon2_members <- c(
    "cs",
    "si",
    "meta",
    "args",
    "sit",
    "peak",
    "lcpar",
    "mse"
)

align_members <- decon2_members
