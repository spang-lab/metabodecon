# Class Documentation #####

#' @name metabodecon-classes
#' @aliases spectrum spectra decon2 decons2 align aligns
#'
#' @title Metabodecon Class Hierarchy
#'
#' @description
#' Metabodecon represents NMR data using a small set of S3 classes connected
#' by **cumulative inheritance**. A raw spectrum has class `"spectrum"`. After
#' [deconvolute()] it gains the class `"decon2"` (so its class vector becomes
#' `c("decon2", "spectrum")`). After [align()] it gains the class `"align"`
#' (class vector `c("align", "decon2", "spectrum")`). The corresponding
#' collection classes follow the same pattern.
#'
#' Because every deconvoluted/aligned object is still a `spectrum` (in the
#' `inherits()` sense), generic functions defined on `spectrum`/`spectra`
#' (such as `print`, `format`, `summary`, `c`, `plot`) keep working at every
#' stage. The label printed by `print()`/`format()` reflects the most-specific
#' class, e.g. `"align object (...)"`.
#'
#' Element order in an object may vary between versions; always access fields
#' by name (`x$si`, `x[["cs"]]`). Elements marked *optional* may be absent or
#' `NULL`.
#'
#' @section Singlet classes:
#'
#' \describe{
#' \item{`spectrum`}{A single NMR spectrum. Class vector: `"spectrum"`.
#'   Constructed by [read_spectrum()], [make_spectrum()], or
#'   [simulate_spectrum()]. Carries the fields under
#'   *Always present (spectrum)* below.}
#' \item{`decon2`}{A single deconvoluted NMR spectrum. Class vector:
#'   `c("decon2", "spectrum")`. Produced by [deconvolute()]. In addition to
#'   the `spectrum` fields, a `decon2` carries the *Added by deconvolute()*
#'   fields below.}
#' \item{`align`}{A single deconvoluted NMR spectrum whose peak positions
#'   have been aligned across a collection. Class vector:
#'   `c("align", "decon2", "spectrum")`. Produced by [align()]. Carries
#'   everything a `decon2` does, plus the *Added by align()* fields below.}
#' }
#'
#' @section Collection classes:
#'
#' For each singlet class there is a collection class that wraps a list of
#' those singlets:
#'
#' \describe{
#' \item{`spectra`}{List of `spectrum`. Class vector `"spectra"`.}
#' \item{`decons2`}{List of `decon2`. Class vector `c("decons2", "spectra")`.}
#' \item{`aligns`}{List of `align`. Class vector
#'   `c("aligns", "decons2", "spectra")`.}
#' }
#'
#' Collections inherit from `"spectra"`, so generic methods written for
#' `spectra` also work on `decons2` and `aligns`. Constructed by
#' [read_spectra()] (returns `spectra`), [deconvolute()] when given a
#' `spectra` (returns `decons2`), and [align()] (returns `aligns`).
#' Concatenation with `c()` follows the cumulative rule: the result class is
#' the most-general (least-specific) class present among the inputs. Mixing
#' an `align` with a plain `decon2` yields a `decons2`; mixing any plain
#' `spectrum` in yields a `spectra`.
#'
#' @section Always present (spectrum):
#'
#' \enumerate{
#' \item `cs`: Vector of chemical shifts in ppm. Same length as `si`.
#' \item `si`: Vector of signal intensities (au). `si[i]` is the intensity
#'   at `cs[i]`.
#' \item `meta`: Optional list of metadata, e.g.:
#'   \itemize{
#'     \item `name`: Name of the spectrum, e.g. `"Blood 1"`.
#'     \item `path`: Path to the source file/folder.
#'     \item `type`: Experiment type, e.g. `"H1 CPMG"` or `"H1 NOESY"`.
#'     \item `fq`: Signal frequencies in Hz (same length as `si`/`cs`).
#'     \item `mfs`: Magnetic field strength in Tesla.
#'     \item `simpar`: True Lorentz-curve parameters (simulated spectra only).
#'   }
#' }
#'
#' @section Added by deconvolute():
#'
#' A `decon2` object additionally has:
#'
#' \enumerate{
#'   \setcounter{enumi}{3}
#' \item `args`: List of deconvolution parameters used (`nfit`, `smit`,
#'   `smws`, `delta`, `sfr`, `igrs`, `npmax`, `use_rust`, `verbose`).
#' \item `sit`: Data frame of signal intensities after transformations:
#'   `sm` (smoothed), `sup` (superposition of fitted Lorentz curves), and
#'   `supal` (superposition of *aligned* Lorentz curves, added by `align()`).
#' \item `peak`: Data frame of peak triplets with columns `center`, `left`,
#'   `right`: integer indices into `cs`.
#' \item `lcpar`: Data frame of Lorentz-curve parameters with columns `A`
#'   (amplitude), `lambda` (half-width), `x0` (center, in `cs` units), and
#'   `x0al`/`pcial` (aligned center and integer index into `cs`, added by
#'   `align()`).
#' \item `mse`: List of mean-squared errors: `raw` (between `si` and
#'   `sit$sup`), `norm` (`raw` divided by `sum(sit$sup)`), `sm` (between
#'   `sit$sm` and `sit$sup`), `smnorm` (`sm` divided by `sum(sit$sup)`).
#' }
#'
#' @section Added by align():
#'
#' An `align` object has the same fields as `decon2`, but with the
#' alignment slots populated: `lcpar$x0al`, `lcpar$pcial`, `sit$supal`.
#'
#' @section Methods, predicates, and converters:
#'
#' Methods defined for `spectrum`/`spectra` (and inherited by all
#' subclasses): [print()][print.spectrum], `format()`, `summary()`,
#' `c()`, `[`.
#'
#' Predicates: [is_spectrum()], [is_spectra()]. For lifecycle-specific
#' checks use `inherits(x, "decon2")`, `inherits(x, "aligns")`, etc.
#'
#' Converters: [as_spectra()] turns a path or list of `spectrum` into a
#' `spectra`. [as_decon2()] / [as_decons2()] are identity converters that
#' validate their input.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' s <- sim[[1]]
#' inherits(s, "spectrum")
#'
#' d <- deconvolute(s, sfr = c(3.55, 3.35))
#' class(d)               # c("decon2", "spectrum")
#' inherits(d, "spectrum") # TRUE
#'
#' ds <- deconvolute(sim[1:3], sfr = c(3.55, 3.35))
#' class(ds)              # c("decons2", "spectra")
NULL


# Print #####

#' @export
#'
#' @name print.spectrum
#' @rdname print.spectrum
#'
#' @title Print Method for spectrum and spectra Objects
#'
#' @description
#' S3 print methods for the base metabodecon classes. Subclasses (`decon2`,
#' `align`, `decons2`, `aligns`) inherit these methods; the printed label
#' reflects the most-specific class. See [metabodecon-classes].
#'
#' @param x The object to print.
#' @param name Logical or string. If `TRUE`, prepend the object's name. If a
#' string, prepend that string.
#' @param ... Unused. Accepted to comply with [base::print()].
#'
#' @return `NULL`, invisibly. Called for the side effect of printing.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' print(sim[[1]])
#' print(sim[[1]], name = TRUE)
#' print(sim)
#' print(deconvolute(sim[[1]], sfr = c(3.55, 3.35)))
print.spectrum <- function(x, name = FALSE, ...) {
    cat(format(x, name = name), "\n", sep = "")
    invisible(NULL)
}

#' @export
#' @rdname print.spectrum
print.spectra <- function(x, ...) {
    sg <- if (length(x)) class(x[[1]])[1] else "spectrum"
    catf("%s object with %d %s elements:\n", class(x)[1], length(x), sg)
    invisible(sapply(x, print, name = TRUE))
}

# Format #####

#' @export
#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' @noRd
format.spectra <- function(x, ...) {
    sg <- if (length(x)) class(x[[1]])[1] else "spectrum"
    sprintf("%s object with %d %s elements", class(x)[1], length(x), sg)
}

# Summary #####

#' @export
#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
summary.spectrum <- function(object, ...) {
    x <- object
    base <- list(
        name = get_name(x, NA_character_),
        n_dp = length(x$cs),
        ppm_min = min(x$cs),
        ppm_max = max(x$cs)
    )
    if (inherits(x, "decon2"))
        c(base, list(n_peaks = length(x$lcpar$A), mse_norm = x$mse$norm))
    else
        c(base, list(si_min = min(x$si), si_max = max(x$si)))
}

#' @export
#' @noRd
summary.spectra <- function(object, ...) {
    rows <- lapply(object, function(e) as.data.frame(summary(e)))
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

# Subset #####

#' @export
#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
`[.spectra` <- function(x, i, ...) {
    result <- NextMethod("[")
    class(result) <- class(x)
    result
}

# Concat #####

#' @export
#' @noRd
#' @title Concatenate spectrum/spectra Objects
#' @description
#' Combines any mix of `spectrum`-family singlets, `spectra`-family
#' collections, and lists of singlets into a single collection. The output
#' class chain is the most-general (least-specific) common class among the
#' inputs: a mix of `align` and plain `decon2` yields `c("decons2",
#' "spectra")`; any plain `spectrum` in the mix yields just `"spectra"`.
#' @author 2024-2025 Tobias Schmidt: initial version.
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
#' @noRd
c.spectra <- c.spectrum

# Predicates #####

#' @export
#'
#' @name is_spectrum
#' @rdname is_spectrum
#'
#' @title Is an Object a spectrum or spectra?
#'
#' @description
#' Check if an object inherits from one of the base metabodecon classes
#' (`spectrum` or `spectra`). Since deconvoluted (`decon2`) and aligned
#' (`align`) objects inherit from `spectrum` (and `decons2`/`aligns` from
#' `spectra`), they also satisfy these checks. To test for a specific
#' lifecycle stage, use [base::inherits()] directly, e.g.
#' `inherits(x, "decon2")` or `inherits(x, "aligns")`. See
#' [metabodecon-classes].
#'
#' @param x The object to check.
#'
#' @return `TRUE` if the object inherits from the named class, else `FALSE`.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' is_spectrum(sim[[1]])  # TRUE
#' is_spectra(sim[1:2])   # TRUE
#'
#' d <- deconvolute(sim[[1]], sfr = c(3.55, 3.35))
#' is_spectrum(d)              # TRUE (decon2 inherits from spectrum)
#' inherits(d, "decon2")       # TRUE
is_spectrum <- function(x) inherits(x, "spectrum")

#' @export
#' @rdname is_spectrum
is_spectra <- function(x) inherits(x, "spectra")

# Converters #####

#' @export
#'
#' @name as_spectra
#' @rdname as_spectra
#'
#' @title Convert to a Metabodecon Object
#'
#' @description
#' Identity-or-validate converters between metabodecon classes. See
#' [metabodecon-classes] for the class hierarchy.
#'
#' @param x
#' The object to convert. For [as_spectra()], either a `spectrum`,
#' a list of `spectrum`, or a path passed to [read_spectra()].
#'
#' @param file_format,expno,procno,raw,silent,force
#' Passed to [read_spectra()] when `x` is a path.
#'
#' @return An object of the requested class.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' as_spectra(sim[[1]])
#' as_decon2(deconvolute(sim[[1]], sfr = c(3.55, 3.35)))
as_spectra <- function(x,
                       file_format = "bruker",
                       expno = 10,
                       procno = 10,
                       raw = FALSE,
                       silent = TRUE,
                       force = FALSE) {
    if (inherits(x, "spectra")) {
        x
    } else if (inherits(x, "spectrum")) {
        xx <- structure(list(x), class = "spectra")
        set_names(xx, get_names(xx))
    } else if (is.list(x) && all(sapply(x, inherits, "spectrum"))) {
        xx <- structure(x, class = "spectra")
        set_names(xx, get_names(xx))
    } else if (is.character(x) && file.exists(x)) {
        read_spectra(x, file_format, expno, procno, raw, silent, force)
    } else {
        stop("Input must be a path, spectrum, or list of spectrum objects.")
    }
}

#' @export
#' @rdname as_spectra
as_decon2 <- function(x) {
    if (inherits(x, "decon2")) x
    else stop(sprintf("Cannot convert %s to decon2.", class(x)[1]))
}

#' @export
#' @rdname as_spectra
as_decons2 <- function(x) {
    if (inherits(x, "decons2")) return(x)
    if (is.list(x) && all(sapply(x, inherits, "decon2"))) {
        out <- structure(x, class = c("decons2", "spectra"))
        return(set_names(out, get_names(out)))
    }
    stop("Input must be a list of decon2 objects or a decons2 object.")
}

# Getters (Private) #####

#' @noRd
#' @title Returns the name of an iterable.
#' @param x An iterable object, e.g. a single metabodecon object.
#' @param default Default name if no name is found.
#' @return The name of the object as string or whatever is given as `default`.
#' @author 2024-2025 Tobias Schmidt: initial version.
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

#' @export
#' @title Returns the names of a metabodecon collection object.
#' @param x A metabodecon collection object.
#' @param default Default names if no names are found. Passed on to `get_default_names`.
#' @return A character vector of names.
#' @author 2024-2025 Tobias Schmidt: initial version.
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

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
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

#' @noRd
#' @description
#' Generates dummy peak selection results based on peak centers and chemical
#' shifts, in the format expected by `decon2` objects. Necessary because the
#' Rust backend does not yet expose the true peak selection results.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
get_peak <- function(x0, cs) {
    center <- round(convert_pos(x0, cs, seq_along(cs)))
    data.frame(left = center - 1, center = center, right = center + 1)
}

# Setters #####

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
set_names <- function(x, nams) {
    has_names <- all(sapply(x, function(e) "name" %in% names(e)))
    has_meta_names <- all(sapply(x, function(e) "name" %in% names(e$meta)))
    names(x) <- nams
    if (has_names) for (i in seq_along(x)) x[[i]]$name <- nams[[i]]
    if (has_meta_names) for (i in seq_along(x)) x[[i]]$meta$name <- nams[[i]]
    x
}

# Members (Private) #####

decon2_members <- c("cs", "si", "meta", "args", "sit", "peak", "lcpar", "mse")
