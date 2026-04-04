# Exported Main #####

#' @export
#'
#' @title Align Spectra
#'
#' @description
#' Align signals across  a  list  of  deconvoluted  spectra  using  the  'CluPA'
#' algorithm from the 'speaq' package, described  in  Beirnaert  et  al.  (2018)
#' <doi:10.1371/journal.pcbi.1006018>     and     Vu     et      al.      (2011)
#' <doi:10.1186/1471-2105-12-405> plus the additional peak combination described
#' in [metabodecon::combine_peaks()].
#'
#' @param x
#' An object of  type  `decons1`  or  `decons2`  as  described  in  [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).   To
#' align `decons0` objects (as returned by the now deprecated
#' [metabodecon::MetaboDecon1D]), you can use [metabodecon::as_decons2()] to
#' convert it to a `decons2` object first.
#'
#' @param maxShift
#' Maximum number of points along the "ppm-axis" a value can  be  moved  by  the
#' 'speaq' package. 50 is a suitable starting value for plasma  spectra  with  a
#' digital resolution of 128K. Note that this parameter has to  be  individually
#' optimized  depending  on  the  type  of  analyzed  spectra  and  the  digital
#' resolution. For urine which is more prone to chemical shift  variations  this
#' value most probably has to be increased. Passed  as  argument  `maxShift`  to
#' [metabodecon::speaq_align()].
#'
#' @param maxCombine
#' Amount of adjacent columns which may be combined for improving
#' the alignment during the CluPA step. We recommend setting this
#' to 0 and instead relying on the peak snapping implemented  in
#' [metabodecon::get_si_mat()]. Since version 1.7.0  the  default
#' is 0. Non-zero values are deprecated and support  will  be
#' removed in a future version.
#'
#' @param verbose
#' Whether to print additional information during the alignment process.
#'
#' @param nworkers
#' Number of parallel workers for the alignment. Default is 1 (no parallelism).
#'
#' @param ref
#' Optional reference spectrum of type `align` or `decon2`. When supplied,
#' all spectra in `x` are aligned towards this reference. The reference is
#' prepended to `x` internally and removed from the result. If `NULL`
#' (default), the reference is chosen automatically via
#' `speaq::findRef()`.
#'
#' @param use_speaq
#' If `FALSE` (default), alignment uses metabodecon's  own  built-in
#' implementation.   If  `TRUE`,  the  external  'speaq'  package  is
#' used  instead.   The  'speaq'  backend  is  no  longer  recommended
#' since  version  1.7.0  and  will  be  removed  in  a  future  version.
#' When `use_speaq = TRUE`, the packages 'speaq', 'MassSpecWavelet'
#' and 'impute' must be installed (see `install_deps`).
#'
#' @param install_deps
#' Only used when `use_speaq = TRUE`.  'speaq' relies on the
#' 'MassSpecWavelet' and 'impute' packages. Both, 'MassSpecWavelet' and 'impute'
#' are   not   available    on    CRAN,    but    can    be    installed    from
#' [Bioconductor](https://www.bioconductor.org/)                              or
#' [R-Universe](https://r-universe.dev/). If `install_deps=TRUE`, these packages
#' will  be  automatically  installed  from  R-Universe   without   asking   for
#' confirmation. If `install_deps=NULL` (default), the user will  be  asked  for
#' confirmation  before  installing  missing   dependencies.   If   asking   for
#' confirmation is not possible or `install_deps=FALSE`, the function will raise
#' an error if the packages are not installed.
#'
#' @return
#' An object of type `align` as described in [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' if (interactive()) {
#'     # Example requires an interactive R session, because in case of missing
#'     # dependencies the user will be asked for confirmation to install them.
#'     decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
#'     aligned <- align(decons)
#' }
align <- function(x,
                  maxShift = 50,
                  maxCombine = 0,
                  verbose = TRUE,
                  use_speaq = FALSE,
                  install_deps = NULL,
                  nworkers = 1,
                  ref = NULL) {

    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())

    # If an external reference is supplied, prepend it to the
    # input so speaq aligns everything towards it.
    has_ext_ref <- !is.null(ref)
    if (has_ext_ref) {
        ref_d2 <- ref
        if (is_align(ref_d2)) class(ref_d2) <- "decon2"
        x <- c(ref_d2, x)
    }

    # Check and convert inputs
    xx <- as_decons2(x)
    stopifnot(length(xx) > 1)
    ndps <- sapply(xx, function(xi) length(xi$cs))
    if (length(unique(ndps)) != 1) {
        fmt <- paste(
            "All spectra must have the same number of data points",
            "but spectrum 1 has %d",
            "and spectrum %d has %d."
        )
        msg <- sprintf(fmt, ndps[1], i, ndps[i])
        stop(msg, call. = FALSE)
    }

    # Check for required packages (only needed for speaq backend)
    if (use_speaq) {
        pkgvec <- c("MassSpecWavelet", "impute")
        if (isTRUE(install_deps)) bioc_install(pkgvec, ask = FALSE, verbose = FALSE)
        if (is.null(install_deps)) bioc_install(pkgvec, ask = TRUE, verbose = FALSE)
        is_installed <- sapply(pkgvec, requireNamespace, quietly = TRUE)
        if (any(!is_installed)) {
            pkgvec_missing <- pkgvec[!is_installed]
            pkgstr_missing <- paste(pkgvec_missing, collapse = ", ")
            msg <- paste(
                "The following required packages are missing: %s.",
                "Please install them manually from Bioconductor",
                "or R-Universe or try again with",
                "`install_deps = TRUE`"
            )
            stop(sprintf(msg, pkgstr_missing))
        }
    }

    # Start alignment
    nworkers_str <- if (nworkers == 1) "1 worker" else sprintf("%d workers", nworkers)
    logf("Starting alignment of %d deconvoluted spectra with %s", length(xx), nworkers_str)
    starttime <- Sys.time()

    # Do initial alignment. The result
    # object contains element `new_peakList` which contains the "peak center
    # indices after alignment" (pciaa). The indices are given as continuous
    # numbers. E.g. a value of 1044.28 means that the aligned peak center is
    # between the datapoint 1044 and 1045.
    backend <- if (use_speaq) "speaq" else "built-in"
    logf("Performing %s alignment with maxShift = %d", backend, maxShift)
    X <- get_sup_mat(xx)
    peakList <- lapply(xx, get_peak_indices)
    find_ref_fn <- if (use_speaq) speaq::findRef else find_ref
    refInd <- if (has_ext_ref) 1L else find_ref_fn(peakList)$refInd
    if (nworkers == 1) {
        obj <- dohCluster(
            X = X, peakList = peakList, refInd = refInd,
            maxShift = maxShift, verbose = verbose,
            use_speaq = use_speaq
        )
    } else {
        # Split seq_len(now(X)) `nworkers` submatrices, where refInd is always the first row.
        # Example: if nworkers == 3, nrow(X) == 10 and refInd == 4, the submatrices would be:
        # X[c(4,1,2,3), ], X[c(4,5,6,7), ], X[c(4,8,9,10), ]
        # Then we can call dohCluster() in parallel on each submatrix with
        # refInd == 1 and then combine the results. For parallel execution we use mclapply,
        # as is done in `deconvolute_spectra()`.
        idx <- setdiff(seq_len(nrow(X)), refInd)
        k <- min(nworkers, length(idx))
        grp <- cut(seq_along(idx), breaks = k, labels = FALSE)
        chunks <- lapply(split(idx, grp), function(ids) c(refInd, ids))
        XX <- lapply(chunks, function(rows) X[rows, , drop = FALSE])
        peakLists <- lapply(chunks, function(rows) peakList[rows])
        nw_apply <- min(nworkers, length(chunks))
        objs <- mcmapply(
            nw_apply, dohCluster, XX, peakLists,
            refInd = 1, maxShift = maxShift,
            verbose = verbose, use_speaq = use_speaq
        )
        Y <- X
        new_peakList <- peakList
        for (j in seq_along(chunks)) {
            rows <- chunks[[j]]
            rows_no_ref <- rows[rows != refInd]
            if (length(rows_no_ref) == 0) next
            Y[rows_no_ref, ] <- objs[[j]]$Y[-1, , drop = FALSE]
            new_peakList[rows_no_ref] <- objs[[j]]$new_peakList[-1]
        }
        obj <- list(Y = Y, new_peakList = new_peakList)
    }
    pciaa <- obj$new_peakList

    # Discretize the continous peak center indices from above by rounding. Then
    # prepare the "signal-matrix-after-alignment" (smat). A non-zero value at
    # position i,j indicates that spectrum i has a peak centered at datapoint j.
    # The value of the signal equals the area parameter of the peak. A value of
    # zero at position i,j indicates, that spectrum i has no peak centered at
    # datapoint j. Example:
    #
    # pciaa[[1]] == c(143.2, 120.0, 548.9); xx[[1]]$lcpar$A == c(10, 20, 30)
    # pciaa[[2]] == c(143.4, 122.4, 548.7); xx[[2]]$lcpar$A == c(40, 50, 60)
    # pciaa[[3]] == c(141.9, 123.1, 548.6); xx[[3]]$lcpar$A == c(70, 80, 90)
    #
    # becomes
    #
    # i |... | 142 | 143 | ... | 120 | 121 | 122 | 123 | ... | 548 |
    # --|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
    # A |... |   0 |  10 | ... |  20 |   0 |   0 |     | ... |  30 |
    # B |... |   0 |  40 | ... |   0 |   0 |  40 |     | ... |  60 |
    # C |... |  70 |   0 | ... |   0 |   0 |  80 |     | ... |  90 |
    #
    logf("Discretizing peak center indices")
    smat <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    for (i in seq_len(nrow(X))) {
        d <- round(pciaa[[i]])
        A <- xx[[i]]$lcpar$A
        if (length(d) == 0) stop(sprintf("No peaks found in spectrum %d", i))
        dups <- which(duplicated(d))
        while (length(dups) > 0) {
            fmt <- paste(
                "%d peak center indicies led to duplicate positions.",
                "Original positions: %s"
            )
            origstr <- paste(pciaa[[i]][dups], collapse = ", ")
            logf(fmt, length(dups), origstr)
            d[dups] <- d[dups] + 1
            dups <- which(duplicated(d))
        }
        smat[i, d] <- A
    }

    # Combine signals ACROSS spectra. Example: the matrix from above would
    # become:
    #
    # i |... | 142 | 143 | ... | 120 | 121 | 122 | 123 | ... | 548 |
    # --|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
    # A |... |   0 |  10 | ... |   0 |   0 |  20 |     | ... |  30 |
    # B |... |   0 |  40 | ... |   0 |   0 |  40 |     | ... |  60 |
    # C |... |   0 |  70 | ... |   0 |   0 |  80 |     | ... |  90 |
    #
    logf("Combining peaks across spectra with maxCombine = %d", maxCombine)
    cmat <- if (maxCombine >= 1) combine_peaks(smat, maxCombine)$long else smat

    # Create `align` objects from the `decon2` objects:
    # 1. Store the new peak centers in their respective slot `$lcpar$x0_al`
    # 2. Calculate new signal intensities as superposition of lorentz curves,
    #    (using the updated peak centers) and store them in `$sit$supal`.
    # 3. Store the new signal intensities as integrals in `$sit$al`.
    logf("Creating aligns object")
    cs <- xx[[1]]$cs
    n <- length(cs)
    for (i in seq_along(xx)) {
        al <- rep(0, n)
        pciac <- which(cmat[i, ] != 0)  # Peak center indices after combine_peaks
        lcpar <- xx[[i]]$lcpar          # Lorentzian curve parameters
        if (length(pciac) != length(lcpar$A)) {
            msg <- paste(
                "Number of peak centers after alignment (%d)",
                "does not match number of peaks",
                "in original spectrum (%d)."
            )
            stop(sprintf(msg, length(pciac), length(lcpar$A)))
        }
        x0_al <- cs[pciac]              # Chemical shifts of peak centers
        al[pciac] <- lcpar$A * pi       # SIs as integrals of aligned lorentzians
        xx[[i]]$lcpar$x0_al <- x0_al
        xx[[i]]$sit$al <- al
        xx[[i]]$sit$supal <- lorentz_sup(cs, x0_al, lcpar$A, lcpar$lambda)
        class(xx[[i]]) <- "align"
    }
    aligns <- structure(xx, class = "aligns")

    # If an external reference was prepended, strip it from
    # the result so the caller only gets the aligned input.
    if (has_ext_ref) {
        aligns <- aligns[-1]
        class(aligns) <- "aligns"
    }

    duration <- format(round(Sys.time() - starttime, 3))
    logf("Finished alignment in %s", duration)


    aligns
}

# Exported Helpers #####

#' @export
#' @title Extract Matrix of aligned Signal Intensities
#'
#' @description
#' Extracts a peak-area matrix from aligned spectra. Rows are
#' chemical-shift positions, columns are spectra. With
#' `maxSnap = FALSE` the raw aligned integral vectors are returned.
#' With `maxSnap = TRUE` peaks are snapped to the peaks of a
#' reference spectrum, reducing the row count to the number
#' of reference peaks.
#'
#' @param x
#' An object of type `aligns`.
#'
#' @param maxSnap
#' Controls peak snapping. `FALSE` or `0` (default): off.
#' `TRUE` or `1`: snap within one half-width. Any positive
#' number scales the radius, e.g. `maxSnap = 2` allows two
#' half-widths. See 'Details'.
#'
#' @param ref
#' A single `align` or `decon2` object whose peaks define
#' the rows of the output matrix. If `NULL` (default), the
#' reference is auto-detected from `x` via
#' `speaq::findRef()`. Ignored when `maxSnap = 0`.
#'
#' @param drop_zero
#' If `TRUE`, rows where all values are zero are removed.
#'
#' @details
#' When `maxSnap > 0`, each peak in every spectrum is mapped to
#' the nearest reference peak. A peak is kept only if the
#' distance is at most `maxSnap * lambda_ref / dp` data points,
#' where `lambda_ref` is the half-width of the corresponding
#' reference peak and `dp` is the chemical-shift step size.
#' Areas of peaks that map to the same reference peak are
#' summed.
#'
#' @return
#' A numeric matrix with chemical shifts as rownames and
#' spectrum names as colnames.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' if (interactive()) {
#'     decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
#'     aligns <- align(decons, maxCombine = 0)
#'     si_mat_0 <- get_si_mat(aligns)                  # raw
#'     si_mat_1 <- get_si_mat(aligns, maxSnap = 1)     # 1x hw
#'     si_mat_2 <- get_si_mat(aligns, maxSnap = 2)     # 2x hw
#' }
get_si_mat <- function(x, maxSnap = 0, ref = NULL,
                       drop_zero = FALSE) {
    stopifnot(is_aligns(x))
    cs <- x[[1]]$cs
    if (maxSnap == 0) {
        mat <- sapply(x, function(xi) xi$sit$al)
        rownames(mat) <- cs
    } else {
        if (is.null(ref)) {
            pl <- lapply(x, get_peak_indices)
            idx <- find_ref(pl)$refInd
            ref <- x[[idx]]
        }
        ref_x0 <- if (is_align(ref)) ref$lcpar$x0_al else ref$lcpar$x0
        ref_cs <- ref_x0
        ref_idx <- match(ref_cs, cs)
        dp <- abs(cs[2] - cs[1])
        thresh <- maxSnap * ref$lcpar$lambda / dp
        nr <- length(ref_idx)
        ns <- length(x)
        mat <- matrix(0, nrow = nr, ncol = ns)
        for (s in seq_len(ns)) {
            peak_idx <- match(x[[s]]$lcpar$x0_al, cs)
            A <- x[[s]]$lcpar$A
            for (p in seq_along(peak_idx)) {
                dists <- abs(peak_idx[p] - ref_idx)
                best <- which.min(dists)
                if (dists[best] <= thresh[best]) {
                    mat[best, s] <- mat[best, s] + A[p] * pi
                }
            }
        }
        rownames(mat) <- ref_cs
    }
    colnames(mat) <- get_names(x)
    if (drop_zero) {
        mat <- mat[rowSums(mat != 0) > 0, , drop = FALSE]
    }
    mat
}

# Exported Helpers (Deprecated) #####

#' @export
#'
#' @title Get PPM Range covered by Spectra
#'
#' @description
#' Helper function of `align()`. Should not be called directly by the user.
#'
#' Returns the ppm range across all peaks of the provided deconvoluted spectra.
#'
#' Direct usage of this function has been deprecated with metabodecon version
#' 1.4.3 and will be removed with metabodecon version 2.0.0.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param spectrum_data
#' A list of deconvoluted spectra as returned by [metabodecon::generate_lorentz_curves()].
#'
#' @param full_range
#' If TRUE, the full range of the spectra is returned. If FALSE, only the range
#' from the lowest to the highest peak center is returned.
#'
#' @return
#' A vector containing the lowest and highest ppm value over all peaks of the
#' provided deconvoluted spectra.
#'
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
#'
#' @examples
#' spectrum_data <- generate_lorentz_curves(
#'     data_path = sim[1:2],
#'     nfit = 3,
#'     sfr = c(3.55, 3.35),
#'     wshw = 0,
#'     ask = FALSE,
#'     verbose = FALSE
#' )
#' ppm_rng <- get_ppm_range(spectrum_data)
#' print(ppm_rng)
get_ppm_range <- function(spectrum_data, full_range = FALSE) {
    if (is_decons0(spectrum_data) || is_decons1(spectrum_data)) {
        if (full_range) xx <- lapply(spectrum_data, function(s) s$x_values_ppm)
        else            xx <- lapply(spectrum_data, function(s) s$peak_triplets_middle)
    } else if (is_decons2(spectrum_data) || is_aligns(spectrum_data)) {
        if (full_range) xx <- lapply(spectrum_data, function(s) s$cs)
        else            xx <- lapply(spectrum_data, function(s) s$lcpar$x0)
    } else {
        stop(paste0(
            "spectrum_data must be a decons0, decons1, decons2 or aligns object.\n",
            "See help('metabodecon_classes') for details."
        ))
    }
    x_min <- min(sapply(xx, min))
    x_max <- max(sapply(xx, max))
    c(x_min, x_max)
}

#' @export
#'
#' @title Generate Feature Matrix.
#'
#' @description
#' Helper function of `align()`. Should not be called directly by the user.
#'
#' Generates a list of elements required by [metabodecon::speaq_align()].
#' See 'Value' for a detailed description of the list elements.
#'
#' Direct usage of this function has been deprecated with metabodecon version
#' 1.4.3 and will be removed with metabodecon version 2.0.0.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param data_path A list of deconvoluted spectra as returned by
#' `generate_lorentz_curves()`. In older versions, this could also be the path
#' passed to `generate_lorentz_curves()`, but this is deprecated and will
#' trigger a warning. See 'Details' for more information.
#'
#' @param ppm_range The ppm range over which your signals are distributed.
#'
#' @param si_size_real_spectrum Number of data points in your spectra.
#'
#' @param scale_factor_x The x scale factor used during the deconvolution.
#'
#' @param warn Whether to print a warning in case a file path is passed to
#' `data_path` instead of a list of deconvoluted spectra.
#'
#' @details Before version 1.2 of metabodecon, the deconvolution functions
#' `generate_lorentz_curves` and `MetaboDecon1D` wrote their output partially as
#' txt files to their input folder. Back then, `gen_feat_mat()` used those txt
#' files as input to generate the feature matrix. Since version 1.2 these txt
#' files are no longer created by default, to prevent accidental modifications
#' of the input folders. Therefore, the recommended way to pass the required
#' information to `gen_feat_mat()` is to directly pass the output of
#' `generate_lorentz_curves()` to `gen_feat_mat()`. However, to stay backwards
#' compatible, the name of parameter `data_path` was not changed and passing an
#' actual path to `data_path` is still possible, but will result in a warning
#' (unless `warn` is set to `FALSE`).
#'
#' @return A list with the following elements:
#'
#' `data_matrix`: A data.frame where each row corresponds to one spectrum and
#'                each column to one data point, i.e. for 10 input spectra
#'                with 131072 data points each `data_matrix` would have
#'                dimensions 10 x 131072.
#'
#' `peakList`: A list of vectors, where each vector contains the indices of
#'             the peaks in the corresponding spectrum. The indices increase
#'             from left to right, i.e. the smallest index corresponds to the
#'             highest ppm value, as the ppm values decrease from left to
#'             right.
#'
#' `w`: A list of vectors where each vector contains the "position parameter"
#'      of the peaks in the corresponding spectrum.
#'
#' `A`: A list of vectors where each vector contains the "area parameter" of
#'      the peaks in the corresponding spectrum.
#'
#' `lambda`: A list of vectors where each vector contains the "width
#'           parameter" of the peaks in the corresponding spectrum.
#'
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
#'
#' @examples
#' sim_subset <- metabodecon_file("sim_subset")
#' decons <- generate_lorentz_curves_sim(sim_subset)
#' obj <- gen_feat_mat(decons)
#' str(obj, 2, give.attr = FALSE)
gen_feat_mat <- function(data_path,
                         ppm_range = get_ppm_range(data_path),
                         si_size_real_spectrum = length(data_path$y_values),
                         scale_factor_x = 1000,
                         warn = TRUE) {

    D <- get_decon_params(data_path, warn, check = TRUE)
    X <- do.call(rbind, D$spectrum_superposition)
    P <- lapply(seq_along(D$A), function(i) ncol(X) - D$w[[i]] * scale_factor_x) # (1)
    names(P) <- names(D$A)
    # (1) We get the indices in SDPs, which decrease from left to right, but we
    # need them as normal indices for speaq (i.e. increasing from left to
    # right). For example if we have n = 8 datapoints, and the third and the
    # seventh point from left are peaks, we want the numbers 3 and 7, but we get
    # 0.005 and 0.001. To fix this, we first revert the scaling and then
    # subtract from n = 8, i.e.
    #
    # 3 = 8 - (0.005 * 1000)
    # 7 = 8 - (0.001 * 1000)
    #
    # _____  _____  Peak_  _____  _____  _____  Peak_  _____
    # 0.007  0.006  0.005  0.004  0.003  0.002  0.001  0.000
    # 1      2      3      4      5      6      7      8
    #
    # In the code line above, ncol(X) corresponds to 8, scale_factor_x to 1000
    # and D$w[[i]] to c(0.005, 0.001).
    list(data_matrix = X, peakList = P, w = D$w, A = D$A, lambda = D$lambda)
}

#' @export
#'
#' @title Align Signals using 'speaq'
#'
#' @description
#' Helper function of `align()`. Should not be called directly by the user.
#'
#' Performs signal alignment across the individual spectra using the 'speaq'
#' package (Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S, Pieters L, et al.
#' (2018) speaq 2.0: A complete workflow for high-throughput 1D NMRspectra
#' processing and quantification. PLoS Comput Biol 14(3): e1006018.
#' https://www.doi.org/10.1371/journal.pcbi.1006018). The spectra deconvolution
#' process yields the signals of all spectra. Due to slight changes in
#' measurement conditions, e.g. pH variations, signal positions may vary
#' slightly across spectra. As a consequence, prior to further analysis signals
#' belonging to the same compound have to be aligned across spectra. This is the
#' purpose of the 'speaq' package.
#'
#' Direct usage of this function has been deprecated with metabodecon version
#' 1.4.3 and will be removed with metabodecon version 2.0.0.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param feat Output of `gen_feat_mat()`.
#'
#' @param maxShift Maximum number of points along the "ppm-axis" which a value
#' can be moved by speaq package e.g. 50. 50 is a suitable starting value for
#' plasma spectra with a digital resolution of 128K. Note that this parameter
#' has to be individually optimized depending on the type of analyzed spectra
#' and the digital resolution. For urine which is more prone to chemical shift
#' variations this value most probably has to be increased.
#'
#' @param spectrum_data Output of `generate_lorentz_curves()`.
#'
#' @param si_size_real_spectrum Number of real data points in your original
#' spectra.
#'
#' @param verbose Whether to print additional information during the alignment
#' process.
#'
#' @param show Whether to plot the original and aligned spectra.
#'
#' @param mfrow Layout to use for the plot. Passed on to `par()`. Use `mfrow =
#' NULL` if the plot layout should not be changed.
#'
#' @return
#' A matrix containing the integral values of the spectra after alignment.
#'
#' There is one row per spectrum and one column per ppm value. The entry at
#' position `i, j` holds the integral value of the signal from spectrum `i` that
#' has its center at position `j` after alignment by speaq. If there is no
#' signal with center `j` in spectrum `i`, entry `i, j` is set to NA. The column
#' names of the matrix are the ppm values of the original spectra.
#'
#' Example return matrix:
#'
#' ```
#'     ...  3.59  3.55  3.57  3.56  3.55  3.54  3.53
#'   .----------------------------------------------> PPM
#' 1 | NA   NA    0.20  NA    NA    NA    0.25  NA
#' 2 | NA   NA    0.15  NA    NA    NA    0.13  NA
#' 3 | NA   NA    NA    0.2   NA    NA    0.18  NA
#' SpNr
#' ```
#'
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
#'
#' @examples
#' deps <- c("MassSpecWavelet", "impute")
#' deps_installed <- sapply(deps, requireNamespace, quietly = TRUE)
#' if (all(deps_installed)) {
#'     # 'speaq' requires 'MassSpecWavelet' and 'impute' to be installed
#'     sim_subset <- metabodecon_file("bruker/sim_subset")
#'     spectrum_data <- generate_lorentz_curves_sim(sim_subset)
#'     feat <- gen_feat_mat(spectrum_data)
#'     maxShift <- 200
#'     M <- speaq_align(feat, maxShift, spectrum_data, show = TRUE)
#'     str(M)
#' }
speaq_align <- function(feat = gen_feat_mat(spectrum_data),
                        maxShift = 50,
                        spectrum_data,
                        si_size_real_spectrum = length(spectrum_data[[1]]$y_values),
                        verbose = TRUE,
                        show = FALSE,
                        mfrow = c(2, 1)) {
    Y <- feat$data_matrix
    nsp <- nrow(Y) # Number of spectra
    ndp <- ncol(Y) # Number of data points
    upci_list <- feat$peakList # Unaligned peak center indices
    idx_ref <- find_ref(upci_list)$refInd # Reference spectrum
    clust_obj <- dohCluster(Y, upci_list, idx_ref, maxShift, verbose)
    apci_list <- clust_obj$new_peakList # Aligned peak center indices
    ppm <- spectrum_data[[1]]$x_values_ppm
    M <- matrix(nrow = nsp, ncol = ndp, dimnames = list(NULL, ppm))
    for (i in seq_len(nsp)) {
        M[i, round(apci_list[[i]])] <- feat$A[[i]] * (-pi)
    }
    if (show) {
        opar <- par(mfrow = mfrow, mar = c(5.1, 4.1, 2.1, 0.1))
        on.exit(par(opar), add = TRUE)
        plot_si_mat(feat$data_matrix, main = "Original Spectra")
        plot_si_mat(clust_obj$Y, main = "Aligned Spectra")
    }
    return(M)
}

#' @export
#'
#' @title Combine Peaks
#'
#' @description
#'
#' Helper function of `align()`. Should not be called directly by the user.
#'
#' Even after calling [metabodecon::speaq_align()], the alignment of individual signals is
#' not always perfect, as 'speaq' performs a segment-wise alignment i.e. groups
#' of signals are aligned. For further improvements, partly filled neighboring
#' columns are merged. See 'Details' for an illustrative example.
#'
#' Direct usage of this function has been deprecated with metabodecon version
#' 1.4.3 and will be removed with metabodecon version 2.0.0.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param shifted_mat The matrix returned by `speaq_align()`.
#'
#' @param range Amount of adjacent columns which are permitted to be used for
#' improving the alignment.
#'
#' @param lower_bound Minimum amount of non-zero elements per column to trigger
#' the alignment improvement.
#'
#' @param spectrum_data The list of deconvoluted spectra as returned by
#' `generate_lorentz_curves()` that was used to generate `shifted_mat`. No
#' longer required since version 1.2 of Metabodecon.
#'
#' @param data_path If not NULL, the returned dataframes `long` and `short` are
#' written to `data_path` as "aligned_res_long.csv" and "aligned_res_short.csv".
#'
#' @return A list containing two data frames `long` and `short`. The first data
#' frame contains one column for each data point in the original spectrum. The
#' second data frame contains only columns where at least one entry is non-zero.
#'
#' @details
#'
#' Example of what the function does:
#'
#' ```txt
#' |            | 1    | 2    | 3    | 4    | 5    |
#' |----------- |------|------|------|------|------|
#' | Spectrum 1 | 0.13 | 0    | 0    | 0.11 | 0    |
#' | Spectrum 2 | 0    | 0.88 | 0    | 0.12 | 0    |
#' | Spectrum 3 | 0.07 | 0.56 | 0.30 | 0    | 0    |
#' | Spectrum 4 | 0.08 | 0    | 0.07 | 0    | 0.07 |
#' | Spectrum 5 | 0.04 | 0    | 0    | 0.04 | 0    |
#' ```
#'
#' becomes
#'
#' ```txt
#' |            | 1    | 2    | 3    | 4    | 5    |
#' |----------- |------|------|------|------|------|
#' | Spectrum 1 | 0.13 | 0    | 0    | 0.11 | 0    |
#' | Spectrum 2 | 0    | 0.88 | 0    | 0.12 | 0    |
#' | Spectrum 3 | 0.07 | 0.56 | 0    | 0.30 | 0    |
#' | Spectrum 4 | 0.08 | 0    | 0    | 0.07 | 0.07 |
#' | Spectrum 5 | 0.04 | 0    | 0    | 0.04 | 0    |
#' ```
#'
#' I.e.
#'
#' 1. Column 1 and 2 get NOT merged, because they have a common non-zero entry.
#' 2. Column 3 and 4 get merged, because they are in `range` of each other and
#'    have no common non-zero entries.
#' 3. Column 4 and 5 get NOT merged, because it is more beneficial to merge
#'    column 3 and 4, as they have more mergeable entries and after merging
#'    column 3 and 4, column 4 and 5 have a common non-zero entry.
#'
#' @author
#'
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
#'
#' @examples
#' deps <- c("MassSpecWavelet", "impute")
#' deps_installed <- sapply(deps, requireNamespace, quietly = TRUE)
#' if (all(deps_installed)) {
#'     # 'speaq' requires 'MassSpecWavelet' and 'impute' to be installed
#'     sim_subset <- metabodecon_file("bruker/sim_subset")
#'     spectrum_data <- generate_lorentz_curves_sim(sim_subset)
#'     shifted_mat <- speaq_align(spectrum_data = spectrum_data, verbose = FALSE)
#'     range <- 5
#'     lower_bound <- 1
#'     obj <- combine_peaks(shifted_mat, range, lower_bound)
#'     str(obj)
#' }
combine_peaks <- function(shifted_mat,
                          range = 5,
                          lower_bound = 1,
                          spectrum_data = NULL,
                          data_path = NULL) {
    M <- replace(shifted_mat, is.na(shifted_mat), 0)
    # U[i,j] is TRUE if M[i,j] is nonzero, else FALSE.
    U <- M != 0
    # uu[j] gives the number of nonzero elements in M[,j].
    uu <- colSums(U)
    nc <- ncol(M)
    for (i in (nrow(M) - 1):lower_bound) {
        for (j in which(uu == i)) {
            if (uu[j] == 0) next
            from <- max(1, j - range)
            to <- min(nc, j + range)
            nn <- seq(from, to)
            nn <- nn[nn != j]
            if (length(nn) == 0) next
            mj <- M[, j]
            uj <- U[, j]
            repeat {
                nn <- nn[uu[nn] > 0]
                if (length(nn) == 0) break
                cc <- combine_scores(U, uu, j, nn, uj = uj)
                if (max(cc) == 0) break
                n <- nn[which.max(cc)]
                mj <- mj + M[, n]
                uj <- uj | U[, n]
                uu[j] <- uu[j] + uu[n]
                M[, n] <- 0
                U[, n] <- FALSE
                uu[n] <- 0
                nn <- nn[nn != n]
                if (length(nn) == 0) break
            }
            M[, j] <- mj
            U[, j] <- uj
        }
    }
    if (!is.null(spectrum_data)) colnames(M) <- spectrum_data[[1]]$x_values_ppm
    S <- M[, uu > 0]
    if (!is.null(data_path)) {
        long.csv <- file.path(data_path, "aligned_res_long.csv")
        short.csv <- file.path(data_path, "aligned_res_short.csv")
        utils::write.csv2(M, file = long.csv)
        utils::write.csv2(S, file = short.csv)
    }
    list("short" = S, "long" = M)
}

#' @export
#'
#' @title Cluster Based Peak Alignment
#'
#' @description
#'
#' Helper function of `align()`. Should not be called directly by the user.
#'
#' Rewrite of `speaq::dohCluster()`, compatible with the data format returned by
#' 'generate_lorentz_curves()' and 'gen_feat_mat()'. The function name
#' "dohCluster" comes from "Do Hierarchical Clustering" which is part of the
#' Alignment algorithm proposed by Vu et al. (2011) in
#' <doi:10.1186/1471-2105-12-405>.
#'
#' Direct usage of this function has been deprecated with metabodecon version
#' 1.4.3 and will be removed with metabodecon version 2.0.0.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param X Dataframe of signal intensities from all spectra as returned by
#' [metabodecon::gen_feat_mat()].
#'
#' @param peakList List of peak indices as returned [metabodecon::gen_feat_mat()].
#'
#' @param refInd Number of the reference spectrum i.e. the spectrum to which all
#' signals will be aligned to.
#'
#' @param maxShift Maximum number of points a value can be moved.
#'
#' @param verbose Whether to print additional information during the alignment
#' process.
#'
#' @return
#' A list containing two data frames `Y` and `new_peakList`. The first one
#' contains the aligned spectra, the second one contains the aligned signals of
#' each spectrum.
#'
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
#'
#' @examples
#' deps <- c("MassSpecWavelet", "impute")
#' deps_installed <- sapply(deps, requireNamespace, quietly = TRUE)
#' if (all(deps_installed)) {
#'     # 'speaq' requires 'MassSpecWavelet' and 'impute' to be installed
#'     sim_subset <- metabodecon_file("bruker/sim_subset")
#'     decons <- generate_lorentz_curves_sim(sim_subset)
#'     feat <- gen_feat_mat(decons)
#'     refObj <- speaq::findRef(feat$peakList)
#'     hclObj <- dohCluster(
#'         X = feat$data_matrix,
#'         peakList = feat$peakList,
#'         refInd = refObj$refInd,
#'         maxShift = 100,
#'         verbose = TRUE
#'     )
#'     str(hclObj, 1)
#' }
dohCluster <- function(X,
                       peakList,
                       refInd = 0,
                       maxShift = 100,
                       verbose = TRUE,
                       use_speaq = FALSE) {
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())
    Y <- X
    peakListNew <- peakList
    startTime <- proc.time()
    nspec <- nrow(X)
    logf("Running dohCluster with maxShift = %d on %d spectra", maxShift, nspec)
    refSpec <- Y[refInd, ]
    for (tarInd in seq_len(nspec)) {
        if (tarInd != refInd) {
            logf("Aligning spectrum %d/%d", tarInd, nspec)
            targetSpec <- Y[tarInd, ]
            n_ref <- length(peakList[[refInd]])
            n_tar <- length(peakList[[tarInd]])
            myPeakList <- c(peakList[[refInd]], peakList[[tarInd]])
            myPeakLabel <- c(rep(1, n_ref), rep(0, n_tar))
            startP <- 1
            endP <- length(targetSpec)
            res <- if (use_speaq) {
                speaq::hClustAlign(
                    refSpec, targetSpec,
                    myPeakList, myPeakLabel,
                    startP, endP,
                    maxShift = maxShift,
                    acceptLostPeak = FALSE
                )
            } else {
                hclust_align(
                    refSpec, targetSpec,
                    myPeakList, myPeakLabel,
                    startP, endP,
                    maxShift = maxShift
                )
            }
            Y[tarInd, ] <- res$tarSpec
            n_res <- length(res$peakList)
            n_end <- min(n_ref + n_tar, n_res)
            peakListNew[[tarInd]] <- res$peakList[
                (n_ref + 1):n_end
            ]
        }
    }
    elapsed <- (proc.time() - startTime)[3]
    logf("Finished dohCluster in %.1f s", elapsed)
    list("Y" = Y, "new_peakList" = peakListNew)
}

# Private Helpers #####

#' @noRd
#'
#' @description
#' Calculates a "combine score" for a set of columns `nn` of a
#' Matrix `M`. The score describes how "beneficial" it is to merge column `n`
#' (from `nn`) into column `j`. A column `n` is considered "combinable" if there
#' is no row where columns `n` and `j` both have a non-zero element. If a column
#' is not combinable, its combine score is 0. If a column is combinable, its
#' combine score is the amount of non-zero elements in the column. This function
#' calculates the combine score for all columns `n` in `nn` and then returns the
#' calculated combine scores as vector.
#'
#' @param U Matrix describing nonzero entries of `M.` `U[i,j]` should be `TRUE`
#' if `M[i,j]` is nonzero, else `FALSE.` I.e. for any matrix `M` you can
#' generate `U` as `U <- M != 0`.
#'
#' @param uu Vector describing the amount of nonzero elements in each column of
#' `M`. I.e. for any matrix `M` you can generate `uu` as `uu <- colSums(M !=
#' 0)`. Called `uu` because it denotes the amount of elements that are *unequal
#' zero*.
#'
#' @param j Index of the column where the neighbor columns `nn` should be
#' combined into.
#'
#' @param nn Indices of the columns to calculate the combine score for. Called
#' `nn` because this should be a set of *neighbor columns*.
#'
#' @details
#' Since we only need to know whether an element is nonzero in order to
#' calculate the combine scores, the function takes the matrix `U` and the
#' vector `uu` as input instead of their common ancestor `M`. This is much more
#' efficient if the function is called multiple times, as the conversion must
#' not be done multiple times.
#'
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
#'
#' @examples
#' M <- rbind(
#'     c(2, 0, 2, 2, 0),
#'     c(2, 1, 0, 2, 0),
#'     c(0, 1, 0, 2, 0),
#'     c(0, 0, 3, 0, 1)
#' )
#' U <- M != 0
#' uu <- colSums(U) # 2 2 2 3 1
#'
#' cc <- combine_scores(U, uu, j = 2, nn = c(1, 3))
#' cc[1] == 0  # M[,1] and M[, 2] are not combinable
#' cc[2] == 2  # M[,3] and M[, 2] are combinable and M[, 3] has two nonzero elements
#'
#' cc <- combine_scores(U, uu, j = 1, nn = 2:5)
#' cc[1] == 0  # M[, 2] and M[, 1] are not combinable
#' cc[2] == 0  # M[, 3] and M[, 1] are not combinable
#' cc[3] == 0  # M[, 4] and M[, 1] are not combinable
#' cc[4] == 1  # M[, 5] and M[, 1] are combinable and M[, 5] has one nonzero element
combine_scores <- function(U, uu, j, nn, uj = NULL) {
    nn <- nn[nn >= 1 & nn <= ncol(U)]
    if (length(nn) == 0) return(numeric(0))
    if (is.null(uj)) uj <- U[, j]

    # A neighbor is combinable if it has no shared nonzero row with column j.
    overlaps <- .colSums(U[, nn, drop = FALSE] & uj, nrow(U), length(nn))
    cc <- uu[nn]
    cc[overlaps > 0] <- 0
    unname(cc)
}



#' @noRd
#'
#' @description
#' Helper function of [metabodecon::gen_feat_mat()] to extract the deconvolution parameters
#' from `data_path`, where `data_path` can be a `decon1` or `decons1` object or
#' a folder containing a `"parameters.txt"` and `"approximated_spectrum.txt"`
#' file, as created when calling `MetaboDecon1D()` before version 1.2.
#'
#' @param data_path
#' A list of deconvoluted spectra as returned by [metabodecon::generate_lorentz_curves()] or
#' a path to a folder containing ".* parameters.txt" and ".*
#' approximated_spectrum.txt" files.
#'
#' @param warn
#' (logical) Whether to print warning in case a file path is provided instead of
#' a list of deconvoluted spectra.
#'
#' @param check
#' (logical) Whether to sanity check the deconvolution parameters before
#' returning them.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @return
#' A list containing the deconvolution parameters, i.e. `w`, `lambda`, `A`, and
#' `spectrum_superposition`.
get_decon_params <- function(data_path, warn = TRUE, check = TRUE) {
    dd <- data_path
    if (is.list(dd)) {
        dd <- if (is_decon_list(dd)) {
            dd
        } else if (is_decon_obj(dd)) {
            list(dd)
        } else {
            stop("data_path must be a list of deconvoluted spectra.")
        }
        w <- lapply(dd, function(d) d$x_0)
        lambda <- lapply(dd, function(d) d$lambda)
        A <- lapply(dd, function(d) d$A)
        spectrum_superposition <- lapply(dd, function(d) d$spectrum_superposition)
        params <- named(w, lambda, A, spectrum_superposition)
    } else {
        if (!file.exists(dd)) {
            stop(dd, " does not exist.")
        } else if (warn) {
            warning(
                "You have provided a path to `gen_feat_mat()`. Since",
                "metabodecon v1.2 it is recommended to provide the output of",
                "`generate_lorentz_curves()` directly to speed up",
                "computations. For details see section 'Details' after",
                "calling `help('gen_feat_mat')`."
            )
        }
        params <- read_decon_params(dd)
    }
    if (check) check_decon_params(params)
    params
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
read_decon_params <- function(data_path) {
    par_txt <- dir(data_path, "(.*) parameters.txt", full.names = TRUE)
    spc_txt <- dir(data_path, "(.*) approximated_spectrum.txt", full.names = TRUE)
    par_nam <- sub(" parameters.txt", "", basename(par_txt))
    spc_nam <- sub(" approximated_spectrum.txt", "", basename(spc_txt))
    if (length(par_txt) != length(spc_txt)) stop("Number of parameter files and spectrum files differs.")
    if (length(par_txt) == 0) stop("No parameter files found in the given directory.")
    mapply(par_nam, spc_nam, FUN = function(p, s) {
        if (p != s) warning(sprintf("Mismatched file names: %s and %s\n", p, s))
    })
    par_lst <- lapply(par_txt, function(file) {
        data <- as.matrix(data.table::fread(file, header = FALSE))
        rownames(data) <- data[, 1]
        data[, -1]
    })
    names(par_lst) <- par_nam
    w <- sapply(par_lst, function(obj) as.numeric(obj["w_new", ]), simplify = FALSE)
    lambda <- sapply(par_lst, function(obj) as.numeric(obj["lambda_new", ]), simplify = FALSE)
    A <- sapply(par_lst, function(obj) as.numeric(obj["A_new", ]), simplify = FALSE)
    spectrum_superposition <- lapply(spc_txt, function(file) {
        as.vector(unlist(data.table::fread(file, header = FALSE)))[-1]
    })
    names(spectrum_superposition) <- spc_nam
    named(w, lambda, A, spectrum_superposition)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
check_decon_params <- function(params) {
    nulls <- unlist(sapply(params, function(pp) sapply(pp, is.null)))
    if (any(nulls)) warning("Detected missing params: ", names(nulls)[nulls])
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
rm_zero_width_peaks <- function(params) {
    for (i in seq_len(params$A)) {
        not_zero <- params$w[[i]] != 0
        params$lambda[[i]] <- params$lambda[[i]][not_zero]
        params$w[[i]] <- params$w[[i]][not_zero]
        params$A[[i]] <- params$A[[i]][not_zero]
    }
    params
}

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

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
get_peak_indices <- function(decon2) {
    x0 <- decon2$lcpar$x0
    cs <- decon2$cs
    convert_pos(x0, cs, seq_along(cs))
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
get_sup_mat <- function(decons2) {
    do.call(rbind, lapply(decons2, function(d) d$sit$sup))
}

#' @noRd
#' @title Install Bioconductor packages from R-Universe
#'
#' @description
#' Installs `pkgs` from [R-Universe](https://r-universe.dev/). (By installing
#' from R Universe instead of Bioconductor, we can skip the additional
#' dependency of BiocManager.)
#'
#' @param pkgs A character vector of package names to install.
#' @param ask Whether to ask the user for confirmation before installing the
#' packages.
#'
#' @details
#' The function works as follows:
#'
#' | req | ask  | ia   | action                                     |
#' | --- | ---- | ---- | ------------------------------------------ |
#' | F   | TF   | TF   | Return                                     |
#' | T   | T    | T    | Ask user and install or return accordingly |
#' | T   | T    | F    | Stop with error                            |
#' | T   | F    | TF   | Install                                    |
#'
#' - req = Is installation required?
#' - ask = Is `ask` argument TRUE?
#' - ia = Interactive session?
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' pkgs <- c("MassSpecWavelet", "impute")
#' evalwith(answers = "n", bioc_install(pkgs))
bioc_install <- function(pkgs, ask = TRUE, verbose = TRUE) {
    is_installed <- sapply(pkgs, requireNamespace, quietly = TRUE)
    pkgs_missing <- pkgs[!is_installed]
    if (length(pkgs_missing) == 0) {
        if (verbose) logf("All requested packages are already installed.\n")
        return(invisible())
    }
    if (ask) {
        if (!interactive()) stop("Cannot ask for confirmation in non-interactive mode.")
        pkgs_missing_str <- paste(pkgs_missing, collapse = " and ")
        pkgs_word <- if (length(pkgs_missing) == 1) "package" else "packages"
        msg <- "Proceeding will install the following %s: %s. Continue?"
        msg <- sprintf(msg, pkgs_word, pkgs_missing_str)
        if (isFALSE(get_yn_input(msg))) return()
    }
    repos <- c("https://bioc.r-universe.dev", "https://cloud.r-project.org")
    utils::install.packages(pkgs_missing, repos = repos, keep_outputs = tmpdir("bioc_install"))
}
