# Convert to anything #####

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

# Singular classes #####

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
as_gspec <- function(x, sf = c(1e3, 1e6)) {
    if (is_gspec(x)) return(x)
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
    g <- locals(without = "x")
    structure(g, class = "gspec")
    # [^1]: Same as `dp / sf[1]`, but with slight numeric differences, so we stick with the old calculation method for backwards compatibility.
    # [^2]: Example: ppm = 11, 23, 35, 47 ==> ppm_step == 12, ppm_nstep ~= 10.6 (not really useful, but we need it for backwards compatibility with MetaboDecon1D results)
}

#' @export
#' @rdname convert
as_gdecon <- function(x) {
    if (is_gdecon(x)) return(x)
    stopifnot(is_gspec(x))
    stopifnot(all(gdecon_members %in% names(x)))
    gdecon <- structure(x, class = "gdecon")
}

#' @export
#' @rdname convert
#' @param sf Vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis.
as_decon1 <- function(x, sf = NULL) {
    if (is_decon1(x)) return(x)
    x <- as_gdecon(x)
    structure(class = "decon1", .Data = list(
        number_of_files = 1,
        filename = x$name,
        x_values = x$sdp,
        x_values_ppm = x$ppm,
        y_values = x$y_smooth,
        spectrum_superposition = s <- lorentz_sup(
            x = x$sdp,
            x0 = x$lcr$w,
            A = x$lcr$A,
            lambda = x$lcr$lambda
        ),
        mse_normed = mean(((x$y_smooth / sum(x$y_smooth)) - (s / sum(s)))^2),
        index_peak_triplets_middle = x$peak$center[x$peak$high],
        index_peak_triplets_left = x$peak$right[x$peak$high],
        index_peak_triplets_right = x$peak$left[x$peak$high],
        peak_triplets_middle = x$ppm[x$peak$center[x$peak$high]],
        peak_triplets_left = x$ppm[x$peak$right[x$peak$high]],
        peak_triplets_right = x$ppm[x$peak$left[x$peak$high]],
        integrals = x$lcr$integrals,
        signal_free_region = c(x$sfr$left_sdp, x$sfr$right_sdp),
        range_water_signal_ppm = x$wsr$hwidth_ppm,
        A = x$lcr$A,
        lambda = x$lcr$lambda,
        x_0 = x$lcr$w,
        y_values_raw = x$y_raw, # Following fields only available in `decon1`, but not `decon0` (since v1.2.0)
        x_values_hz = x$hz,
        mse_normed_raw = mean(((x$y_raw / sum(x$y_raw)) - (s / sum(s)))^2),
        x_0_hz = convert_pos(x$lcr$w, x$sdp, x$hz),
        x_0_dp = convert_pos(x$lcr$w, x$sdp, x$dp),
        x_0_ppm = convert_pos(x$lcr$w, x$sdp, x$ppm),
        A_hz = convert_width(x$lcr$A, x$sdp, x$hz),
        A_dp = convert_width(x$lcr$A, x$sdp, x$dp),
        A_ppm = convert_width(x$lcr$A, x$sdp, x$ppm),
        lambda_hz = convert_width(x$lcr$lambda, x$sdp, x$hz),
        lambda_dp = convert_width(x$lcr$lambda, x$sdp, x$dp),
        lambda_ppm = convert_width(x$lcr$lambda, x$sdp, x$ppm)
    ))
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
            lcpt = NULL
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
    } else if (is_gdecon(x)) {
        x <- as_decon1(x)
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
            normed = NULL
        )
    }
    obj <- named(cs, si, meta, args, sit, peak, lcpar, mse)
    structure(obj, class = "decon2")
    # TODO: check whether all members of decon2 are really defined using
    # is_decon2(check_content = TRUE)
    obj
}

# Public Plural #####

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
as_gspecs <- function(x, sf = c(1e3, 1e6)) {
    if (is_gspecs(x)) return(x)
    gg <- if (is_gspec(x)) {
        list(x)
    } else if (is_spectra(x)) {
        lapply(x, as_gspec, sf = sf)
    } else {
        errmsg <- "Input must be gspec, gspecs or spectra, not "
        stop(errmsg, class(x))
    }
    gg <- structure(gg, class = "gspecs")
    gg <- set_names(gg, get_names(gg))
    gg
}

#' @export
#' @rdname convert
as_gdecons <- function(x) {
    if (is_gdecons(x)) return(x)
    stopifnot(is.list(x))
    stopifnot(all(sapply(x, is_gdecon)))
    structure(x, class = "gdecons")
}

#' @export
#' @rdname convert
as_decons1 <- function(x, sf = c(1e3, 1e6)) {
    if (is_decons1(x)) return(x)
    decons1 <- lapply(x, as_decon1)
    names(decons1) <- names(x)
    class(decons1) <- "decons1"
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

# Helpers #####

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
    stopifnot(class(x)[1] %in% c("gdecons", "gspecs", "spectra"))
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

# Class Members #####

spectrum_members <- c(
    "cs",
    "si",
    "meta"
)

gspec_members <- c(
    "y_raw",
    "y_scaled",
    "n",
    "dp",
    "sdp",
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

gdecon_members <- c(
    gspec_members,
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
