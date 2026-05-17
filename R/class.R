
#' @name metabodecon-classes
#' @aliases spectrum spectra decon2 decons2 align aligns
#' @aliases is_spectrum is_spectra as_spectra as_decon2 as_decons2 get_names
#'
#' @title Metabodecon Classes and Helpers
#'
#' @description
#' Metabodecon represents NMR data using a small set of S3 classes connected
#' by **cumulative inheritance**. A raw spectrum has class `"spectrum"`. After
#' [metabodecon::deconvolute()] it gains class `"decon2"`, so its class vector becomes
#' `c("decon2", "spectrum")`. After [metabodecon::align()] it gains class `"align"`, with
#' class vector `c("align", "decon2", "spectrum")`. The corresponding
#' collection classes follow the same pattern.
#'
#' Every deconvoluted or aligned object is still a `spectrum` in the
#' [base::inherits()] sense, so S3 generic behavior for `spectrum` or
#' `spectra` also works at every stage. Element order may vary between
#' versions; always access fields by name, e.g. `x$si` or `x[["cs"]]`.
#' Elements marked optional may be absent or `NULL`.
#'
#' @usage
#' is_spectrum(x)
#' is_spectra(x)
#' as_spectra(x, ...)
#' as_decon2(x)
#' as_decons2(x)
#' get_names(x, default = "spectrum_\045d")
#'
#' @param x A metabodecon object, collection, list of objects, or path.
#' @param default Used by [metabodecon::get_names()] when no object names are present.
#' @param ... Parameters passed to [metabodecon::read_spectrum()] when `x` is a path.
#'
#' @return
#' `is_spectrum()` and `is_spectra()` return `TRUE` or `FALSE`.
#' The `as_*()` functions return an object of the requested class.
#' `get_names()` returns a character vector.
#'
#' @section Singlet classes:
#'
#' - `spectrum`: A single NMR spectrum. Class vector: `"spectrum"`.
#'   Constructed by [metabodecon::read_spectrum()], [metabodecon::make_spectrum()], or
#'   [metabodecon::simulate_spectrum()]. Carries the fields under *Always present
#'   (spectrum)* below.
#'
#' - `decon2`: A single deconvoluted NMR spectrum. Class vector:
#'   `c("decon2", "spectrum")`. Produced by [metabodecon::deconvolute()]. In addition to
#'   the `spectrum` fields, a `decon2` carries the *Added by deconvolute()*
#'   fields below.
#'
#' - `align`: A single deconvoluted NMR spectrum whose peak positions have
#'   been aligned across a collection. Class vector:
#'   `c("align", "decon2", "spectrum")`. Produced by [metabodecon::align()]. Carries
#'   everything a `decon2` does, plus the *Added by align()* fields below.
#'
#' @section Collection classes:
#'
#' For each singlet class there is a collection class that wraps a list of
#' those singlets:
#'
#' - `spectra`: List of `spectrum`. Class vector `"spectra"`.
#' - `decons2`: List of `decon2`. Class vector `c("decons2", "spectra")`.
#' - `aligns`: List of `align`. Class vector
#'   `c("aligns", "decons2", "spectra")`.
#'
#' Collections inherit from `"spectra"`, so generic methods written for
#' `spectra` also work on `decons2` and `aligns`. Constructed by
#' [metabodecon::read_spectra()] (returns `spectra`), [metabodecon::deconvolute()] when given a
#' `spectra` (returns `decons2`), and [metabodecon::align()] (returns `aligns`).
#' Concatenation follows the cumulative rule: the result class is the
#' most-general, least-specific class among the inputs. Mixing an `align` with
#' a plain `decon2` yields `decons2`; mixing any plain `spectrum` in yields
#' `spectra`.
#'
#' @section Always present (spectrum):
#'
#' 1. `cs`: Vector of chemical shifts in ppm. Same length as `si`.
#' 2. `si`: Vector of signal intensities (au). `si[i]` is the intensity at
#'    `cs[i]`.
#' 3. `meta`: Optional list of metadata, e.g. `name` (spectrum name), `path`
#'    (source path), `type` (experiment type), `fq` (signal frequencies in Hz),
#'    `mfs` (magnetic field strength), or `simpar` (true Lorentz-curve
#'    parameters for simulated spectra).
#'
#' @section Added by deconvolute():
#'
#' A `decon2` object additionally has:
#'
#' 1. `args`: List of deconvolution parameters used (`nfit`, `smit`, `smws`,
#'    `delta`, `sfr`, `igrs`, `npmax`, `use_rust`, `verbose`).
#' 2. `sit`: Data frame of signal intensities after transformations: `sm`
#'    (smoothed), `sup` (superposition of fitted Lorentz curves), and `supal`
#'    (superposition of aligned Lorentz curves, added by `align()`).
#' 3. `peak`: Data frame of peak triplets with columns `center`, `left`,
#'    `right`: integer indices into `cs`.
#' 4. `lcpar`: Data frame of Lorentz-curve parameters. Always carries
#'    `x0` (center in ppm), `A` (amplitude), `lambda` (half-width) and
#'    `pcide` (integer column index into `cssh` for `x0`). After
#'    [metabodecon::clupa()] also `x0al` / `pcial` (post-CluPA center
#'    and cssh index). After [metabodecon::snap_to_ref()] also `x0sn`
#'    / `pcisn` (post-RefPA center and cssh index, with `NA` for peaks
#'    snapped beyond `maxCombine`). `A` and `lambda` are preserved
#'    through every stage.
#'
#' @section Added by align():
#'
#' An `align` object has the same fields as `decon2`, but with the alignment
#' slots populated:
#'
#' - `lcpar$x0al`: Peak Centers after CluPA alignment in ppm
#' - `lcpar$pcial`: Peak Centers after CluPA alignment as `cssh` indices
#' - `lcpar$x0sn`: Peak Centers after RefPA snap in ppm (NA when snapped out)
#' - `lcpar$pcisn`: Peak Centers after RefPA snap as `cssh` indices (NA when snapped out)
#' - `sit$supal`: Signal Intensities of the superposition of aligned Lorentz curves
#'
#' @section Predicates:
#'
#' `is_spectrum()` and `is_spectra()` test inheritance from the base
#' metabodecon classes. Since `decon2` and `align` inherit from `spectrum`,
#' and `decons2` and `aligns` inherit from `spectra`, they satisfy these
#' checks. To test for a specific lifecycle stage, use [base::inherits()]
#' directly, e.g. `inherits(x, "decon2")` or `inherits(x, "aligns")`.
#'
#' @section Converters:
#' `as_spectra()` turns a path, `spectrum`, or list of `spectrum` objects into
#' a `spectra` collection. `as_decon2()` and `as_decons2()` are identity
#' converters that validate their input.
#'
#' @section Naming helpers:
#' `get_names()` returns collection names by checking each element's metadata,
#' each element's direct `name`, the list names, and finally generated default
#' names.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' s <- sim[[1]]
#' inherits(s, "spectrum")
#' is_spectrum(s)
#'
#' d <- deconvolute(s, sfr = c(3.55, 3.35))
#' class(d) # c("decon2", "spectrum")
#' inherits(d, "spectrum") # TRUE
#'
#' ds <- deconvolute(sim[1:3], sfr = c(3.55, 3.35))
#' class(ds)              # c("decons2", "spectra")
#' as_spectra(s)
#' get_names(list(s, myspec = s))
#'
NULL

# API #####

#' @export
is_spectrum <- function(x) inherits(x, "spectrum")

#' @export
is_spectra <- function(x) inherits(x, "spectra")

#' @export
as_spectra <- function(x, ...) {
    if (inherits(x, "spectra")) {
        x
    } else if (inherits(x, "spectrum")) {
        xx <- structure(list(x), class = "spectra")
        set_names(xx, get_names(xx))
    } else if (is.list(x) && all(sapply(x, inherits, "spectrum"))) {
        xx <- structure(x, class = "spectra")
        set_names(xx, get_names(xx))
    } else if (is.character(x) && file.exists(x)) {
        read_spectra(x, ...)
    } else {
        stop("Input must be a path, spectrum, or list of spectrum objects.")
    }
}

#' @export
as_decon2 <- function(x) {
    if (inherits(x, "decon2")) x
    else stop(sprintf("Cannot convert %s to decon2.", class(x)[1]))
}

#' @export
as_decons2 <- function(x) {
    if (inherits(x, "decons2")) return(x)
    if (is.list(x) && all(sapply(x, inherits, "decon2"))) {
        out <- structure(x, class = c("decons2", "spectra"))
        return(set_names(out, get_names(out)))
    }
    stop("Input must be a list of decon2 objects or a decons2 object.")
}

#' @export
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

#' @export
print.spectrum <- function(x, name = FALSE, ...) {
    cat(format(x, name = name), "\n", sep = "")
    invisible(NULL)
}

#' @export
print.spectra <- function(x, ...) {
    sg <- if (length(x)) class(x[[1]])[1] else "spectrum"
    catf("%s object with %d %s elements:\n", class(x)[1], length(x), sg)
    invisible(sapply(x, print, name = TRUE))
}

#' @export
format.spectrum <- function(x, name = FALSE, ...) {
    nam <- {
        if (isTRUE(name)) paste0(get_name(x, "NULL"), ": ")
        else if (is.character(name)) paste0(name, ": ")
        else ""
    }
    np <- if (inherits(x, "decon2")) sprintf(", %d peaks", length(x$lcpar$A)) else ""
    fmt <- "%s%s object (%d dp, %.1f to %.1f ppm%s)"
    sprintf(fmt, nam, class(x)[1], length(x$cs), max(x$cs), min(x$cs), np)
}

#' @export
format.spectra <- function(x, ...) {
    sg <- if (length(x)) class(x[[1]])[1] else "spectrum"
    sprintf("%s object with %d %s elements", class(x)[1], length(x), sg)
}

#' @export
summary.spectrum <- function(object, ...) {
    x <- object
    base <- list(
        name = get_name(x, NA_character_),
        n_dp = length(x$cs),
        ppm_min = min(x$cs),
        ppm_max = max(x$cs)
    )
    if (inherits(x, "decon2"))
        c(base, list(n_peaks = length(x$lcpar$A)))
    else
        c(base, list(si_min = min(x$si), si_max = max(x$si)))
}

#' @export
summary.spectra <- function(object, ...) {
    rows <- lapply(object, function(e) as.data.frame(summary(e)))
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

#' @export
`[.spectra` <- function(x, i, ...) {
    result <- NextMethod("[")
    class(result) <- class(x)
    result
}

#' @export
c.spectrum <- function(..., recursive = FALSE) {
    elems <- list()
    for (a in list(...)) {
        if (is.null(a)) next
        if (inherits(a, "spectra")) elems <- c(elems, unclass(a))
        else if (inherits(a, "spectrum")) elems <- c(elems, list(a))
        else if (is.list(a) && all(sapply(a, inherits, "spectrum")))
            elems <- c(elems, a)
        else stop("All arguments must be spectrum or spectra.", call. = FALSE)
    }
    if (all(sapply(elems, inherits, "align")))
        cls <- c("aligns", "decons2", "spectra")
    else if (all(sapply(elems, inherits, "decon2")))
        cls <- c("decons2", "spectra")
    else
        cls <- "spectra"
    sg <- switch(cls[1], spectra = "spectrum", decons2 = "decon2", aligns = "align")
    out <- structure(elems, class = cls)
    set_names(out, get_names(out, default = paste0(sg, "_%d")))
}

#' @export
c.spectra <- c.spectrum

# Private #####

get_name <- function(x, default = "") {
    (if (is.list(x)) x$meta$name %||% x$name) %||% default
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

get_peak <- function(x0, cs) {
    center <- round(convert_pos(x0, cs, seq_along(cs)))
    data.frame(left = center - 1, center = center, right = center + 1)
}

set_names <- function(x, nams) {
    has_names <- all(sapply(x, function(e) "name" %in% names(e)))
    has_meta_names <- all(sapply(x, function(e) "name" %in% names(e$meta)))
    names(x) <- nams
    if (has_names) for (i in seq_along(x)) x[[i]]$name <- nams[[i]]
    if (has_meta_names) for (i in seq_along(x)) x[[i]]$meta$name <- nams[[i]]
    x
}

decon2_members <- c("cs", "cssh", "si", "meta", "args", "sit", "peak", "lcpar")
