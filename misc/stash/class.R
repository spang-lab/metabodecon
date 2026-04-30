# Stashed decon0 / decon1 / decons0 / decons1 class methods.
# Moved here when depr.R was removed; kept for reference.

# Print #####

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
print.decons1 <- function(x, ...) {
    catf("decons1 object with %s decon1 elements\n", length(x))
    invisible(sapply(x, print, name = TRUE))
}

# Subset #####

#' @export
`[.decons0` <- `[.collection`

#' @export
`[.decons1` <- `[.collection`

# Concat #####

#' @export
c.decon1 <- function(..., recursive = FALSE) {
    concat_collection_args(
        args = list(...),
        recursive = recursive,
        is_elem = is_decon1,
        is_coll = is_decons1,
        coll_class = "decons1",
        default_names = "decon1_%d",
        update_n_files = TRUE,
        err_msg = "All arguments to c.decon1 must be decon1 or decons1."
    )
}

#' @export
c.decons1 <- function(..., recursive = FALSE) {
    c.decon1(..., recursive = recursive)
}

# Format #####

#' @export
format.decon1 <- function(x, ...) {
    ppm <- x$x_values_ppm
    fmt <- "decon1 object (%d dp, %.1f to %.1f ppm, %d peaks)"
    sprintf(fmt, length(ppm), max(ppm), min(ppm), length(x$A))
}

#' @export
format.decons1 <- function(x, ...) {
    sprintf("decons1 object with %d decon1 elements", length(x))
}

# Summary #####

#' @export
summary.decon1 <- function(object, ...) {
    x <- object
    list(
        name = x$filename %||% NA_character_,
        n_dp = length(x$x_values_ppm),
        ppm_min = min(x$x_values_ppm),
        ppm_max = max(x$x_values_ppm),
        n_peaks = length(x$A),
        mse_normed = x$mse_normed
    )
}

#' @export
summary.decons1 <- function(object, ...) {
    summary_collection(object, summary.decon1)
}

# Checks (Public) #####

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
is_decons0 <- function(x) all(sapply(x, is_decon0))

#' @export
#' @rdname is_metabodecon_class
is_decons1 <- function(x) inherits(x, "decons1")

# Convert (Public) #####

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
                      spectrum = NULL,
                      sfr = NULL,
                      wshw = NULL,
                      bwc = 2) {
    if (is_decon0(x)) as_decon1.decon0(x, sf, spectrum, sfr, wshw, bwc)
    else if (is_decon1(x)) x
    else if (is_decon2(x)) as_decon1.decon2(x, sf, spectrum, sfr, wshw, bwc)
    else stop(sprintf("Converting %s to decon1 is not supported", class(x)[1]))
}

#' @export
#' @rdname as_metabodecon_class
as_decons0 <- function(x,
                       sfs = list(c(1e3, 1e6)),
                       spectra = list(NULL),
                       nworkers = 1) {
    if (is_decons0(x)) {
        return(x)
    } else if (is_decons1(x) || is_decons2(x)) {
        decons0 <- mcmapply(as_decon0, x, sfs, spectra, nw = nworkers)
    } else if (is.list(x) && all(sapply(x, is_decon0))) {
        decons0 <- x
    } else {
        stop(paste(
            "Input must be a list of decon0 objects or a single object",
            "of type decons0, decons1 or decons2."
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
#' @rdname as_metabodecon_class
as_decons1 <- function(x,
                       sfs = list(c(1e3, 1e6)),
                       spectra = list(NULL),
                       sfrs = list(NULL),
                       wshws = list(NULL),
                       bwc = 2,
                       nworkers = 1) {
    if (is_decons1(x)) {
        return(x)
    } else if (is_decons0(x) || is_decons2(x)) {
        decons1 <- mcmapply(as_decon1, x, sfs, spectra, sfrs, wshws, bwc, nw = nworkers)
    } else if (is.list(x) && all(sapply(x, is_decon1))) {
        decons1 <- x
    } else {
        stop(paste(
            "Input must be a list of decon1 objects or a single object",
            "of type decons0, decons1 or decons2."
        ))
    }
    names(decons1) <- get_names(x)
    class(decons1) <- "decons1"
    n <- length(decons1)
    for (i in seq_len(n)) decons1[[i]]$number_of_files <- n
    decons1
}

# Convert (Private) #####

as_decon1.decon0 <- function(x,
                            sf = c(1e3, 1e6),
                            spectrum = NULL,
                            sfr = NULL,
                            wshw = NULL,
                            bwc = 2) {
    if (is.null(sf)) stop("Please provide `sf`")
    if (is.null(spectrum)) stop("Please provide `spectrum`")
    # Define some shorthands
    fq <- spectrum$meta$fq
    si <- spectrum$si
    ssp <- as.numeric(x$spectrum_superposition)
    ppm <- x$x_values_ppm
    sdp <- x$x_values
    dp <- round(x$x_values * sf[1])
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
    y$mse_normed_raw <- mse(si, ssp, normed = TRUE)
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
    y
}

as_decon1.decon2 <- function(x, sf, spectrum, sfr, wshw, bwc) {
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
    limits_sdp <- NULL
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
    y
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
as_decon2.decon1 <- function(x, ...) {
    cs <- x$x_values_ppm
    si <- x$y_values_raw
    meta <- list(
        name = x$filename,
        fq = x$x_values_hz
    )
    args <- list(
        nfit = NA, smit = NA, smws = NA, delta = NA,
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
        raw  = mse(si, sit$sup, normed = FALSE),
        norm = x$mse_normed_raw,
        sm = mse(sit$sm, sit$sup, normed = FALSE),
        smnorm = x$mse_normed
    )
    obj <- named(cs, si, meta, args, sit, peak, lcpar, mse)
    class(obj) <- "decon2"
    obj
}

# Members (Private) #####

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

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_decon_obj <- function(x) {
    keys <- c(
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
    if (is.list(x) && all(keys %in% names(x))) TRUE else FALSE
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
is_decon_list <- function(x) {
    if (is.list(x) && all(sapply(x, is_decon_obj))) TRUE else FALSE
}
