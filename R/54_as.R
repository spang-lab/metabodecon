# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Convert to Anything #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @export
#' @name convert
#' @title Convert to a Metabodecon Class
#' @description Convert an object to a Metabodecon class.
#' @param x The object to convert.
#' @param to The class to convert to. For details see [metabodecon_classes].
#' @param ... Additional arguments passed to the conversion function.
#' @return An object of the specified class.
#' @examples
#' dirpath <- metabodecon_file("sim_subset")
#' spectra <- read_spectra(dirpath)
#' spectrum <- spectra[[1]]
#' decons1 <- generate_lorentz_curves_sim(spectra)
#' decon1 <- generate_lorentz_curves_sim(spectrum)
#' decon2 <- as_decon2(decon1)
convert <- function(x, to, ...) {
    as_class <- get(paste0("as_", to), envir = environment(convert))
    as_class(x, ...)
}

# Convert to Singleton Object #####

#' @export
#' @rdname convert
as_spectrum <- function(x, sf = c(1e3, 1e6)) {
    if (is_spectrum(x)) {
        return(x)
    } else if (is_decon1(x)) {
        cs = x$x_values_ppm
        si = x$y_values_raw %||% (x$y_values * sf[2])
        name = x$filename
        fq = x$x_values_hz
        meta <- named(name, fq)
        obj <- named(cs, si, meta)
        return(structure(obj, class = "spectrum"))
    } else {
        msg <- "Converting %s to spectrum is not suppoorted"
        stop(sprintf(msg, class(x)[1]))
    }
}

#' @export
#' @rdname convert
#' @inheritParams as_decon1
as_ispec <- function(x, sf = c(1e3, 1e6)) {
    if (is_ispec(x)) return(x)
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

#' @export
#' @rdname convert
as_idecon <- function(x) {
    if (is_idecon(x)) {
        x
    } else if (all(idecon_members %in% names(x))) {
        structure(x, class = "idecon")
    } else {
        stop("Input must have all elements listed in `idecon_members`")
    }
}

#' @export
#' @rdname convert
#' @param x Object to convert.
#' @param sf Vector with two entries consisting of the factor to scale the
#'           x-axis and the factor to scale the y-axis.
as_decon0 <- function(x, sf = NULL) {
    if (is_decon0(x)) return(x)
    x <- as_decon1(x)
    x[decon0_members]
}

#' @export
#' @rdname convert
#' @param sf Vector with two entries consisting of the factor to scale the
#'           x-axis and the factor to scale the y-axis.
#' @param hz Frequency of each data point in Hz. Only required when converting
#'           from `decon0` to `decon1`.
#' @param si Raw signal intensities.
as_decon1 <- function(x,
                      sf = c(1e3, 1e6),
                      spectrum  = NULL) {
    if (is_decon0(x)) {
        if (is.null(sf)) stop("Please provide `sf`")
        if (is.null(spectrum)) stop("Please provide `spectrum`")
        fq <- spectrum$meta$fq
        si <- spectrum$si
        ssp <- as.numeric(x$spectrum_superposition)
        ppm <- x$x_values_ppm
        sdp <- x$x_values
        dp <- round(x$x_values * sf[1])
        y <- x
        y$y_values_raw <- si
        y$x_values_hz <- fq
        y$mse_normed_raw <- mean(((si / sum(si)) - (ssp / sum(ssp)))^2)
        y$x_0_hz <- convert_pos(x$x_0, sdp, fq)
        y$x_0_dp <- convert_pos(x$x_0, sdp, dp)
        y$x_0_ppm <- convert_pos(x$x_0, sdp, ppm)
        y$A_hz <- convert_width(x$A, sdp, fq)
        y$A_dp <- convert_width(x$A, sdp, dp)
        y$A_ppm <- convert_width(x$A, sdp, ppm)
        y$lambda_hz <- convert_width(x$lambda, sdp, fq)
        y$lambda_dp <- convert_width(x$lambda, sdp, dp)
        y$lambda_ppm <- convert_width(x$lambda, sdp, ppm)
    } else if (is_decon1(x)) {
        y <- x
    } else if (is_idecon(x)) {
        y <- structure(class = "decon1", .Data = list())
        y$number_of_files <- 1
        y$filename <- x$name
        y$x_values <- x$sdp
        y$x_values_ppm <- x$ppm
        y$y_values <- x$y_smooth
        y$spectrum_superposition <- s <- lorentz_sup(x = x$sdp, lcpar = x$lcr)
        y$mse_normed <- mean(((x$y_smooth / sum(x$y_smooth)) - (s / sum(s)))^2)
        y$index_peak_triplets_middle <- x$peak$center[x$peak$high]
        y$index_peak_triplets_left <- x$peak$right[x$peak$high]
        y$index_peak_triplets_right <- x$peak$left[x$peak$high]
        y$peak_triplets_middle <- x$ppm[x$peak$center[x$peak$high]]
        y$peak_triplets_left <- x$ppm[x$peak$right[x$peak$high]]
        y$peak_triplets_right <- x$ppm[x$peak$left[x$peak$high]]
        y$integrals <- x$lcr$integrals
        y$signal_free_region <- c(x$sfr$left_sdp, x$sfr$right_sdp)
        y$range_water_signal_ppm <- x$wsr$hwidth_ppm
        y$A <- x$lcr$A
        y$lambda <- x$lcr$lambda
        y$x_0 <- x$lcr$w
        # Fields only available in `decon1` (but not in `decon0`):
        y$y_values_raw <- x$y_raw
        y$x_values_hz <- x$hz
        y$mse_normed_raw <- mean(((x$y_raw / sum(x$y_raw)) - (s / sum(s)))^2)
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
        stop("Converting %s to decon1 is not supported", class(x)[1])
    }
    y
}

#' @export
#' @rdname convert
as_decon2 <- function(x) {
    if (is_decon2(x)) {
        return(x)
    } else if (is_decon1(x)) {
        cs <- x$x_values_ppm
        si <- x$y_values_raw
        meta <- list(
            name = x$filename,
            path = NULL,
            type = NULL,
            fq = x$x_values_hz,
            mfs = NULL,
            simpar = NULL
        )
        args <- list(
            nfit = NA,
            smopts = NA,
            delta = NA,
            sfr = x$signal_free_region,
            wsr = x$range_water_signal_ppm
        )
        sit <- list(
            wsrm = NA,
            nvrm = NA,
            sm = x$y_values,
            sup = x$spectrum_superposition,
            al = NULL
        )
        peak <- list(
            center = x$index_peak_triplets_middle,
            left = x$index_peak_triplets_left,
            right = x$index_peak_triplets_right
        )
        lcpar <- list(
            A = x$A,
            lambda = x$lambda,
            x0 = x$x_0
        )
        mse <- list(
            raw = NA,
            norm = x$mse_normed_raw,
            sm = NA,
            smnorm = x$msw_normed
        )
    } else if (is_idecon(x)) {
        cs <- x$ppm
        si <- x$y_raw
        meta <- x$meta
        args <- x$args
        sit <- list(
            wsrm = x$y_nows,
            nvrm = x$y_pos,
            smooth = x$y_smooth,
            sup = NULL
        )
        peak <- x$peak
        lcpar <- x$lcr
        mse <- list(
            raw = NULL,
            normed = NULL,
            sm = NULL,
            smnorm = NULL
        )
    } else {
        stop(sprintf("Converting %s to decon2 is not supported", class(x)[1]))
    }
    obj <- named(cs, si, meta, args, sit, peak, lcpar, mse)
    class(obj) <- "decon2"
    obj
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Convert to Collection Object #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @export
#' @rdname convert
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
#' @rdname convert
as_ispecs <- function(x, sf = c(1e3, 1e6)) {
    if (is_ispecs(x)) return(x)
    gg <- if (is_ispec(x)) {
        list(x)
    } else if (is_spectra(x)) {
        lapply(x, as_ispec, sf = sf)
    } else {
        errmsg <- "Input must be ispec, ispecs or spectra, not "
        stop(errmsg, class(x))
    }
    gg <- structure(gg, class = "ispecs")
    gg <- set_names(gg, get_names(gg))
    gg
}

#' @export
#' @rdname convert
as_idecons <- function(x) {
    if (is_idecons(x)) return(x)
    stopifnot(is.list(x))
    stopifnot(all(sapply(x, is_idecon)))
    structure(x, class = "idecons")
}

#' @export
#' @rdname convert
as_decons1 <- function(x, sf = c(1e3, 1e6), spectra = NULL) {
    if (is_decons1(x)) return(x)
    stopifnot(is_decons0(x) || is_idecons(x))
    decons1 <- structure(vector("list", length(x)), class = "decons1")
    names(decons1) <- names(x)
    for (i in seq_along(x)) {
        decons1[[i]] <- as_decon1(x[[i]], sf = sf, spectrum = spectra[[i]])
    }
    for (i in seq_along(decons1)) decons1[[i]]$number_of_files <- length(decons1)
    decons1
}

#' @export
#' @rdname convert
as_decons2 <- function(x) {
    if (is_decons2(x)) return(x)
    decons2 <- lapply(x, as_decon2)
    names(decons2) <- names(x)
    class(decons2) <- "decons2"
    decons2
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Convert Unit #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @noRd
#' @description
#' Converts a vector of chemical shifts given in ppm to a vector of frequencies
#' in Hz.
#' @param cs Vector of chemical shifts in ppm.
#' @param fqref Frequency of the reference molecule in Hz.
as_frequency <- function(cs, fqref) {
    fqmax <- fqref - (min(cs) * 1e-6 * fqref) # Highest frequency in Hz
    fqmin <- fqref - (max(cs) * 1e-6 * fqref) # Lowest frequency in Hz
    fq <- seq(fqmin, fqmax, length.out = length(cs))
    fq <- fq[rev(order(cs))] # (1)
    # (1) Sort frequencies in descending order of chemical shifts, as the
    #     highest chemical shift corresponds to the lowest frequency.
    fq
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Helpers #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

set_names <- function(x, nams) {
    if (!is.list(x)) stop("Input must be a list.")
    has_names <- all(sapply(x, function(e) "name" %in% names(e)))
    has_meta_names <- all(sapply(x, function(e) "name" %in% names(e$meta)))
    names(x) <- nams
    if (has_names) for (i in seq_along(x)) x[[i]]$name <- nams[[i]]
    if (has_meta_names) for (i in seq_along(x)) x[[i]]$meta$name <- nams[[i]]
    x
}

get_names <- function(x, default = "spectrum_%d") {
    stopifnot(class(x)[1] %in% c("idecons", "ispecs", "spectra"))
    dn <- get_default_names(x, default)
    en <- names(x) # Element name
    sn <- sapply(x, function(s) s$meta$name %||% s$name) # Spectrum name
    sapply(seq_along(x), function(i) sn[i] %||% en[i] %||% dn[i])
}

get_default_names <- function(x, default) {
    if (length(default) == 1 && grepl("%d", default))
        return(sprintf(default, seq_along(x)))
    if (length(unique(default)) == length(x))
        return(default)
    stop("Default names must be a single string with a `%d` placeholder or a character vector of same length as the spectra object.")
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Class Members #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

spectrum_members <- c(
    "cs",
    "si",
    "meta"
)

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
    "sf",
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

decon1_members <- c(
    decon0_members,
    "y_values_raw",
    "x_values_hz",
    "mse_normed_raw",
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
