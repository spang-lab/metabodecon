# Classes #####

metabodecon_classes <- NULL

#' @name metabodecon_classes
#'
#' @title Metabodecon Classes
#'
#' @description
#' Metabodecon introduces a set of classes to highlight the presence of  certain
#' elements in corresponding objects.
#'
#' The order of elements may vary between  different  versions  of  Metabodecon,
#' thus elements should always be accessed by name, for example, using `x$si` or
#' `x[["cs"]]`. A short description of each class is given in the listing below.
#'
#' -  `spectrum`: One NMR spectrum
#' -  `decon0`: One deconvoluted NMR spectrum stored in [MetaboDecon1D()] format
#' -  `decon1`: One deconvoluted NMR spectrum stored in [generate_lorentz_curves()] format
#' -  `decon2`: One deconvoluted NMR spectrum stored in [deconvolute()] format
#' -  `align`: One aligned NMR spectrum
#'
#' The classes mentioned above represent individual objects, such  as  a  single
#' spectrum, deconvolution,  or  alignment.  However,  it  is  often  useful  to
#' describe collections  of  these  objects,  such  as  a  list  of  spectra  or
#' deconvolutions.  Therefore,  for  each  individual  class,  a   corresponding
#' "collection"  class  is  provided.  These  collection  classes   are   named:
#' `spectra`, `decons0`, `decons1`, `decons2`, and `aligns`.
#'
#' More details can be found in Metabodecon's online documentation at
#' [spang-lab.github.io/metabodecon/articles/Metabodecon-Classes](
#' https://spang-lab.github.io/metabodecon/articles/Classes.html).
NULL

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

# Methods #####

print_methods <- NULL

#' @name print_methods
#' @rdname print_methods
#'
#' @title S3 Methods for Printing Metabodecon Objects
#'
#' @description
#' S3 Methods for printing metabodecon objects as described in the  [Metabodecon
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
#' @examples
#' si <- c(1, 1, 3, 7, 8, 3, 8, 5, 2, 1)
#' cs_max <- 14.8
#' cs_width <- 20.0
#' fq_ref <- 600.25 * 1e6
#' fq_width <- 12005
#' spectrum <- read_spectrum()
#' print(spectrum)
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
print.ispec <- function(x, name = FALSE, ...) {
    # fmt <- "%sispec object (%d dp, %.1f to %.1f ppm)\n"
    # namestr <- if (name) paste0(x$name %||% "NULL", ": ") else ""
    # catf(fmt, namestr, length(x$ppm), max(x$ppm), min(x$ppm))
    str(x, 1)
}

#' @export
#' @rdname print_methods
print.idecon <- function(x, name = FALSE, ...) {
    # ppm <- x$ppm
    # n <- length(ppm)
    # name <- if (name) paste0(x$name %||% "NULL", ": ") else ""
    # fmt <- "%sidecon object (%d dp, %.1f to %.1f ppm, %d peaks)\n"
    # catf(fmt, name, n, max(ppm), min(ppm), length(x$A))
    str(x, 1)
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
print.ispecs <- function(x, ...) {
    # catf("ispecs object with %s ispec elements\n", length(x))
    # invisible(sapply(x, print, name = TRUE))
    str(x, 2, give.attr = FALSE)
}

#' @export
#' @rdname print_methods
print.idecons <- function(x, ...) {
    catf("idecons object with %s idecon elements\n", length(x))
    nams <- get_names(x)
    mapply(x, nams, FUN = function(xi, nam) {
        catf("%s: ", nam)
        print(xi, ...)
    })
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

#' @export
`[.spectra` <- function(x, i, ...) {
    result <- NextMethod("[")
    class(result) <- class(x)
    result
}

# Checks #####

is_metabodecon_class <- NULL

#' @export
#'
#' @name is_metabodecon_class
#'
#' @title Is an Object from a Metabodecon Class?
#'
#' @description
#' Check if an object is an instance of a specific [Metabodecon
#' Class](metabodecon_classes).
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
#' is_align(a1)  # TRUE
#'
#' is_spectra(ss) # TRUE
#' is_decons0(dd) # FALSE
#' is_decons1(dd) # FALSE
#' is_decons2(dd) # TRUE
#' is_aligns(aa)  # TRUE
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
is_ispec <- function(x) inherits(x, "ispec")

#' @export
#' @rdname is_metabodecon_class
is_idecon <- function(x) inherits(x, "idecon")

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
is_ispecs <- function(x) inherits(x, "ispecs")

#' @export
#' @rdname is_metabodecon_class
is_idecons <- function(x) inherits(x, "idecons")

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

# Convert #####

as_metabodecon_class <- NULL

#' @export
#'
#' @name as_metabodecon_class
#' @rdname as_metabodecon_class
#'
#' @title Convert to a Metabodecon Class
#'
#' @description Convert an object to a Metabodecon class.
#'
#' @param x
#' The object to convert.
#'
#' @param sf
#' Scale factor. Only required if `x` is a decon0 object.
#'
#' @param sfs List of scale factors. Only required if `x` is a list of decon0
#' objects.
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
#' @rdname as_metabodecon_class
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
#' @rdname as_metabodecon_class
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
#' @rdname as_metabodecon_class
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
#' @rdname as_metabodecon_class
as_decon1 <- function(x,
                      sf = c(1e3, 1e6),
                      spectrum  = NULL,
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
        dpn <- (n-1):0
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
        A_raw_dp <- A_raw_ppm * (dpn_step  / cs_step)
        A_raw_sdp <- A_raw_ppm * (sdp_step / cs_step)
        A_raw_hz  <- if (!is.null(fq)) A_raw_ppm * (fq_step  / cs_step)
        A_sc_ppm <- A_raw_ppm / 1e6
        A_sc_dp <- A_raw_dp  / 1e6
        A_sc_sdp <- A_raw_sdp / 1e6
        A_sc_hz <- if (!is.null(fq)) A_raw_hz  / 1e6
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
#' @rdname as_metabodecon_class
as_decon2 <- function(x,
                      sf = c(1e3, 1e6),
                      spectrum = NULL,
                      sfr = NULL,
                      wshw = NULL,
                      bwc = 2) {
    if (is_decon2(x)) {
        return(x)
    } else if (is_decon1(x)) {
        cs <- x$x_values_ppm
        si <- x$y_values_raw
        meta <- list(
            name = x$filename,
            fq = x$x_values_hz
        )
        args <- list(
            nfit = NA,
            smopts = NA,
            delta = NA,
            sfr = sfr_in_ppm_bwc(x$signal_free_region, x$x_values, x$x_values_ppm),
            wshw = x$range_water_signal_ppm,
            ask = NA,
            force = NA,
            verbose = NA,
            bwc = NA,
            nworkers = NA
        )
        sit <- data.frame(
            wsrm = NA,
            nvrm = NA,
            sm = x$y_values * 1e6,
            sup = x$spectrum_superposition[1, ] * 1e6
        )
        peak <- data.frame(
            left = x$index_peak_triplets_right, # decon[01] has left and right inverted
            center = x$index_peak_triplets_middle,
            right = x$index_peak_triplets_left  # decon[01] has left and right inverted
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
    } else if (is_decon0(x)) {
        x <- as_decon1(x, sf, spectrum, sfr, wshw, bwc)
        return(as_decon2(x))
    } else if (is_idecon(x)) {
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
    } else {
        stop(sprintf("Converting %s to decon2 is not supported", class(x)[1]))
    }
    obj <- named(cs, si, meta, args, sit, peak, lcpar, mse)
    class(obj) <- "decon2"
    obj
}

as_v2_obj <- function(obj) {
    if (is_spectrum(obj)) obj
    else if (is_ispec(obj)) as_spectrum(obj)
    else if (is_idecon(obj)) as_decon2(obj)
    else if (is_decon0(obj)) stop("decon0 objects are not supported. Convert with as_decon2.")
    else if (is_decon1(obj)) as_decon2(obj)
    else if (is_decon2(obj)) obj
    else if (is_align(obj)) obj
    else stopf("Objects of class %s are not supported.", class(obj))
}

as_v2_objs <- function(obj) {
    if (is_spectra(obj)) obj
    else if (is_ispecs(obj)) as_spectra(obj)
    else if (is_idecons(obj)) as_decons2(obj)
    else if (is_decons0(obj)) stop("decons0 objects are not supported. Convert with as_decons2.")
    else if (is_decons1(obj)) as_decons2(obj)
    else if (is_decons2(obj)) obj
    else if (is_aligns(obj)) obj
    else stopf("Objects of class %s are not supported.", class(obj))
}

#' @export
#' @rdname as_metabodecon_class
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
#' @rdname as_metabodecon_class
as_ispecs <- function(x, sf = c(1e3, 1e6)) {
    if (is_ispecs(x)) return(x)
    gg <- if (is_ispec(x)) list(x)
        else if (is_spectrum(x)) list(as_ispec(x))
        else if (is_spectra(x)) lapply(x, as_ispec, sf = sf)
        else stop("Input must be ispec, ispecs or spectra, not ", class(x))
    gg <- structure(gg, class = "ispecs")
    gg <- set_names(gg, get_names(gg))
    gg
}

#' @export
#' @rdname as_metabodecon_class
as_idecons <- function(x) {
    if (is_idecons(x)) return(x)
    stopifnot(is.list(x))
    stopifnot(all(sapply(x, is_idecon)))
    structure(x, class = "idecons")
}

#' @export
#' @rdname as_metabodecon_class
as_decons0 <- function(x,
                       sfs = list(c(1e3, 1e6)),
                       spectra = list(NULL),
                       nworkers = 1) {
    stopifnot(is_decons0(x) || is_decons1(x) || is_decons2(x) || is_idecons(x))
    if (is_decons0(x)) return(x)
    decons0 <- mcmapply(as_decon0, x, sfs, spectra, nw = nworkers)
    n <- length(decons0)
    for (i in seq_len(n)) decons0[[i]]$number_of_files <- n
    decons0
}

#' @export
#' @rdname as_metabodecon_class
as_decons1 <- function(x,
                       sfs = list(c(1e3, 1e6)),
                       spectra = list(NULL),
                       sfrs = list(NULL),
                       wshws = list(NULL),
                       bwc = 2,
                       nworkers = 1) {
    stopifnot(is_decons0(x) || is_decons1(x) || is_decons2(x) || is_idecons(x))
    if (is_decons1(x)) return(x)
    decons1 <- mcmapply(as_decon1, x, sfs, spectra, sfrs, wshws, bwc, nw = nworkers)
    class(decons1) <- "decons1"
    n <- length(decons1)
    for (i in seq_len(n)) decons1[[i]]$number_of_files <- n
    decons1
}

#' @export
#' @rdname as_metabodecon_class
as_decons2 <- function(x,
                       sfs = list(c(1e3, 1e6)),
                       spectra = list(NULL),
                       sfrs = list(NULL),
                       wshws = list(NULL),
                       bwc = 2,
                       nworkers = 1) {
    stopifnot(is_decons0(x) || is_decons1(x) || is_decons2(x) || is_idecons(x))
    if (is_decons2(x)) return(x)
    decons2 <- mcmapply(as_decon2, x, sfs, spectra, sfrs, wshws, bwc, nw = nworkers)
    names(decons2) <- get_names(x)
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
    stopifnot(is.list(x))
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
