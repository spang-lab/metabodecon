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
#' in [combine_peaks()].
#'
#' @param x
#' An object of  type  `decons1`  or  `decons2`  as  described  in  [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).   To
#' align `decons0` objects (as returned by the now deprecated  [MetaboDecon1D]),
#' you can use [as_decons2()] to convert it to a `decons2` object first.
#'
#' @param maxShift
#' Maximum number of points along the "ppm-axis" a value can  be  moved  by  the
#' 'speaq' package. 50 is a suitable starting value for plasma  spectra  with  a
#' digital resolution of 128K. Note that this parameter has to  be  individually
#' optimized  depending  on  the  type  of  analyzed  spectra  and  the  digital
#' resolution. For urine which is more prone to chemical shift  variations  this
#' value most probably has to be increased. Passed  as  argument  `maxShift`  to
#' [speaq_align()].
#'
#' @param maxCombine
#' Amount of adjacent columns which may be combined for improving the alignment.
#' Passed as argument `range` to [combine_peaks()].
#'
#' @param verbose
#' Whether to print additional information during the alignment process.
#'
#' @param install_deps
#' Alignment  relies  on  the  'speaq'  package,  which  itself  relies  on  the
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
align <- function(x, maxShift = 50, maxCombine = 5, verbose = TRUE, install_deps = NULL) {

    # Check for required packages
    pkgvec <- c("MassSpecWavelet", "impute")
    if (isTRUE(install_deps)) bioc_install(pkgvec, ask = FALSE, verbose = verbose)
    if (is.null(install_deps)) bioc_install(pkgvec, ask = TRUE, verbose = verbose)
    is_installed <- sapply(pkgvec, requireNamespace, quietly = TRUE)
    if (any(!is_installed)) {
        pkgvec_missing <- pkgvec[!is_installed]
        pkgstr_missing <- paste(pkgvec_missing, collapse = ", ")
        msg <- paste(
            "The following required packages are missing: %s.",
            "Please install them manually from Bioconductor or R-Universe",
            "or try again with `install_deps = TRUE`"
        )
        stop(sprintf(msg, pkgstr_missing))
    }

    # Check and convert inputs
    xx <- as_decons2(x)
    stopifnot(length(xx) > 1)
    for (i in 2:length(xx)) {
        if (!is_equal(xx[[i-1]]$cs, xx[[i]]$cs)) {
            stop("Chemical shifts must be equal across all spectra.")
        }
    }

    # Do initial alignment using speaq. Parameter `acceptLostPeak` must be FALSE
    # so we can backtrack which peak has shifted to which position. The result
    # object contains element `new_peakList` which contains the "peak center
    # indices after alignment" (pciaa). The indices are given as continuous
    # numbers. E.g. a value of 1044.28 means that the aligned peak center is
    # between the datapoint 1044 and 1045.
    obj <- dohCluster(
        X <- get_sup_mat(xx),
        peakList <- lapply(xx, get_peak_indices),
        refInd = speaq::findRef(peakList)$refInd,
        maxShift = maxShift,
        acceptLostPeak = FALSE,
        verbose = verbose
    )
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
    smat <- matrix(nrow = nrow(X), ncol = ncol(X))
    for (i in seq_len(nrow(X))) smat[i, round(pciaa[[i]])] <- xx[[i]]$lcpar$A

    # Combine signals ACROSS spectra. Example: the matrix from above would
    # become:
    #
    # i |... | 142 | 143 | ... | 120 | 121 | 122 | 123 | ... | 548 |
    # --|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
    # A |... |   0 |  10 | ... |   0 |   0 |  20 |     | ... |  30 |
    # B |... |   0 |  40 | ... |   0 |   0 |  40 |     | ... |  60 |
    # C |... |   0 |  70 | ... |   0 |   0 |  80 |     | ... |  90 |
    #
    cmat <- combine_peaks(smat, maxCombine)$long

    # Create `align` objects from the `decon2` objects:
    # 1. Store the new peak centers in their respective slot `$lcpar$x0_al`
    # 2. Calculate new signal intensities as superposition of lorentz curves,
    #    (using the updated peak centers) and store them in `$sit$supal`.
    # 3. Store the new signal intensities as integrals in `$sit$al`.
    cs <- xx[[1]]$cs
    n <- length(cs)
    for (i in seq_along(xx)) {
        al <- rep(0, n)
        pciac <- which(cmat[i, ] != 0)  # Peak center indices after combine_peaks
        lcpar <- xx[[i]]$lcpar          # Lorentzian curve parameters
        x0_al <- cs[pciac]              # Chemical shifts of peak centers
        al[pciac] <- lcpar$A * pi       # SIs as integrals of aligned lorentzians
        xx[[i]]$lcpar$x0_al <- x0_al
        xx[[i]]$sit$al <- al
        xx[[i]]$sit$supal <- lorentz_sup(cs, x0_al, lcpar$A, lcpar$lambda)
        class(xx[[i]]) <- "align"
    }
    aligns <- structure(xx, class = "aligns")
    aligns
}

# Exported Helpers #####

#' @export
#' @title Extract Matrix of aligned Signal Intensities
#'
#' @description
#' Takes an object of type `aligns`, i.e., a list of deconvoluted and aligned
#' spectra, extracts the vector of aligned signal integrals for each spectrum
#' and returns them as a matrix with datapoints in rows and spectra in columns.
#'
#' @param x
#' An object of type `aligns`.
#'
#' @return
#' A matrix of aligned signal intensities.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
#' aligns <- align(decons)
#' si_mat <- get_si_mat(aligns) # 2048 x 2 matrix (2048 datapoints, 2 spectra)
get_si_mat <- function(x) {
    stopifnot(is_aligns(x))
    sapply(x, function(xi) xi$sit$al)
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
#' A list of deconvoluted spectra as returned by [generate_lorentz_curves()].
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
#' Generates a list of elements required by [speaq_align()].
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
#' sim_subset <- metabodecon_file("bruker/sim_subset")
#' spectrum_data <- generate_lorentz_curves_sim(sim_subset)
#' feat <- gen_feat_mat(spectrum_data)
#' maxShift <- 200
#' M <- speaq_align(feat, maxShift, spectrum_data, show = TRUE)
#' str(M)
speaq_align <- function(feat = gen_feat_mat(spectrum_data),
                        maxShift = 50,
                        spectrum_data,
                        si_size_real_spectrum = length(spectrum_data[[1]]$y_values),
                        verbose = TRUE,
                        show = FALSE,
                        mfrow = c(2, 1)) {
    acceptLostPeak <- FALSE # Must be FALSE so we can assign A and lambda to the respective peaks
    Y <- feat$data_matrix
    nsp <- nrow(Y) # Number of spectra
    ndp <- ncol(Y) # Number of data points
    upci_list <- feat$peakList # Unaligned peak center indices
    idx_ref <- speaq::findRef(upci_list)$refInd # Reference spectrum
    clust_obj <- dohCluster(Y, upci_list, idx_ref, maxShift, acceptLostPeak, verbose) # list(Y, new_peakList)
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
#' Even after calling [speaq_align()], the alignment of individual signals is
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
#' 2021-2024 Wolfram Gronwald: initial version.\cr 2024-2025 Tobias Schmidt:
#' refactored initial version.
#'
#' @examples
#'
#' sim_subset <- metabodecon_file("bruker/sim_subset")
#' spectrum_data <- generate_lorentz_curves_sim(sim_subset)
#' shifted_mat <- speaq_align(spectrum_data = spectrum_data, verbose = FALSE)
#' range <- 5
#' lower_bound <- 1
#' obj <- combine_peaks(shifted_mat, range, lower_bound)
#' str(obj)
combine_peaks <- function(shifted_mat,
                          range = 5,
                          lower_bound = 1,
                          spectrum_data = NULL,
                          data_path = NULL) {
    M <- replace(shifted_mat, is.na(shifted_mat), 0)
    U <- M != 0 # Unequal zero matrix. U[i,j] is TRUE if M[i,j] is nonzero, else FALSE.
    uu <- colSums(U) # Unequal zero vector. u[j] gives the amount of nonzero elements in M[,j].
    for (i in (nrow(M) - 1):lower_bound) {
        for (j in which(uu == i)) {
            if (uu[j] == 0) next
            nn <- c((j - range):(j - 1), (j + 1):(j + range)) # Neighbors of j.
            cc <- combine_scores(U, uu, j, nn)
            while (!all(cc == 0)) {
                n <- nn[which.max(cc)]
                M[, j] <- M[, j] + M[, n]
                U[, j] <- U[, j] | U[, n]
                uu[j] <- uu[j] + uu[n]
                M[, n] <- 0
                U[, n] <- FALSE
                uu[n] <- 0
                cc <- combine_scores(U, uu, j, nn)
            }
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
#' [gen_feat_mat()].
#'
#' @param peakList List of peak indices as returned [gen_feat_mat()].
#'
#' @param refInd Number of the reference spectrum i.e. the spectrum to which all
#' signals will be aligned to.
#'
#' @param maxShift Maximum number of points a value can be moved.
#'
#' @param acceptLostPeak Whether to allow the the alignment algorithm to ignore
#' peaks that cannot easily be aligned with the reference spectrum.
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
#' sim_subset <- metabodecon_file("bruker/sim_subset")
#' decons <- generate_lorentz_curves_sim(sim_subset)
#' feat <- gen_feat_mat(decons)
#' refObj <- speaq::findRef(feat$peakList)
#' hclObj <- dohCluster(
#'      X = feat$data_matrix,
#'      peakList = feat$peakList,
#'      refInd = refObj$refInd,
#'      maxShift = 100,
#'      acceptLostPeak = TRUE,
#'      verbose = TRUE
#' )
#' str(hclObj, 1)
dohCluster <- function(X,
                       peakList,
                       refInd = 0,
                       maxShift = 100,
                       acceptLostPeak = TRUE,
                       verbose = TRUE) {
    res <- if (is.null(maxShift)) {
        if (verbose) {
            cat("\n --------------------------------")
            cat("\n maxShift=NULL, thus dohCluster will automatically detect the optimal value of maxShift.")
            cat("\n --------------------------------\n")
        }
        dohCluster_withoutMaxShift(X, peakList, refInd, acceptLostPeak, verbose)
    } else {
        if (verbose) {
            cat("\n --------------------------------")
            cat("\n dohCluster will run with maxShift=", maxShift)
            cat("\n If you want dohCluster to detect the optimal maxShift automatically,")
            cat("\n use dohCluster(..., maxShift = NULL, ...)")
            cat("\n --------------------------------\n")
        }
        dohCluster_withMaxShift(X, peakList, refInd, maxShift, acceptLostPeak, verbose)
    }
    return(res)
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
combine_scores <- function(U, uu, j, nn) {
    sapply(seq_along(nn), function(k) {
        n <- nn[k] # Index of neighboring column
        combinable <- sum(U[, n] & U[, j]) == 0
        if (combinable) uu[n] else 0
    })
}

#' @noRd
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
dohCluster_withoutMaxShift <- function(X,
                                       peakList,
                                       refInd = 0,
                                       acceptLostPeak = TRUE,
                                       verbose = TRUE) {
    if (verbose) {
        startTime <- proc.time()
    }
    maxShift_ladder <- 2^(c(1:trunc(log2(ncol(X) / 2))))
    bestCor <- -1
    corVec <- NULL
    bestY <- NULL
    bestMaxShift <- 0
    for (maxShift_val in maxShift_ladder) {
        if (verbose) {
            cat("\n maxShift=", maxShift_val)
        }
        Y <- X
        peakListNew <- peakList
        refSpec <- Y[refInd, ]
        for (tarInd in seq_len(nrow(X))) {
            if (tarInd != refInd) {
                targetSpec <- Y[tarInd, ]
                myPeakList <- c(peakList[[refInd]], peakList[[tarInd]])
                myPeakLabel <- c(
                    rep(1, length(peakList[[refInd]])),
                    rep(0, length(peakList[[tarInd]]))
                )
                startP <- 1
                endP <- length(targetSpec)
                res <- speaq::hClustAlign(
                    refSpec,
                    targetSpec,
                    myPeakList,
                    myPeakLabel,
                    startP,
                    endP,
                    maxShift = maxShift_val,
                    acceptLostPeak = acceptLostPeak
                )
                Y[tarInd, ] <- res$tarSpec
                if (length(myPeakList) > length(res$peakList)) {
                    peakListNew[[tarInd]] <- res$peakList[(length(peakList[[refInd]]) +
                        1):length(res$peakList)]
                } else {
                    peakListNew[[tarInd]] <- res$peakList[(length(peakList[[refInd]]) +
                        1):length(myPeakList)]
                }
            }
        }
        Z <- stats::cor(t(Y))
        newCor <- stats::median(Z[lower.tri(Z)])
        corVec <- c(corVec, newCor)
        if (verbose) {
            cat(
                "\n Median Pearson correlation coefficent:",
                newCor, ", the best result:", bestCor
            )
        }
        if (newCor > bestCor) {
            bestCor <- newCor
            bestY <- Y
            bestMaxShift <- maxShift_val
        }
    }
    if (verbose) {
        cat(
            "\nOptimal maxShift=", bestMaxShift, "with median Pearson correlation of aligned spectra=",
            bestCor
        )
        plot(log2(maxShift_ladder), corVec,
            type = "b",
            xlab = "log2(maxShift)", ylab = "Median Pearson correlation coefficent",
            main = paste("Optimal maxShift=", bestMaxShift,
                " (red star) \n with median Pearson correlation coefficent of ",
                round(bestCor, 6),
                sep = ""
            )
        )
        graphics::points(log2(bestMaxShift), bestCor, col = "red", pch = 8, cex = 2)
    }
    if (verbose) {
        endTime <- proc.time()
        cat("\n Alignment time: ", (endTime[3] - startTime[3]) / 60, " minutes")
    }
    # Added by me at 31.08.21
    return_list <- list("Y" = bestY, "new_peakList" = peakListNew)
    return(return_list)
}

#' @noRd
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
dohCluster_withMaxShift <- function(X,
                                    peakList,
                                    refInd = 0,
                                    maxShift = 100,
                                    acceptLostPeak = TRUE,
                                    verbose = TRUE) {
    Y <- X
    peakListNew <- peakList
    if (verbose) {
        startTime <- proc.time()
    }
    refSpec <- Y[refInd, ]
    for (tarInd in seq_len(nrow(X))) {
        if (tarInd != refInd) {
            if (verbose) {
                cat("\n aligning spectrum ", tarInd)
            }
            targetSpec <- Y[tarInd, ]
            myPeakList <- c(peakList[[refInd]], peakList[[tarInd]])
            myPeakLabel <- double(length(myPeakList))
            for (i in seq_along(peakList[[refInd]])) myPeakLabel[i] <- 1
            startP <- 1
            endP <- length(targetSpec)
            res <- speaq::hClustAlign(
                refSpec,
                targetSpec,
                myPeakList,
                myPeakLabel,
                startP,
                endP,
                maxShift = maxShift,
                acceptLostPeak = acceptLostPeak
            )
            Y[tarInd, ] <- res$tarSpec
            if (length(myPeakList) > length(res$peakList)) {
                peakListNew[[tarInd]] <- res$peakList[(length(peakList[[refInd]]) +
                    1):length(res$peakList)]
            } else {
                peakListNew[[tarInd]] <- res$peakList[(length(peakList[[refInd]]) +
                    1):length(myPeakList)]
            }
        }
    }
    if (verbose) {
        Z <- stats::cor(t(Y))
        newCor <- stats::median(Z[lower.tri(Z)])
        cat(
            "\n Median pearson correlation of aligned spectra:",
            newCor
        )
        endTime <- proc.time()
        cat(
            "\n Alignment time: ", (endTime[3] - startTime[3]) / 60,
            " minutes"
        )
    }

    # Added modifications:
    return_list <- list("Y" = Y, "new_peakList" = peakListNew)
    return(return_list)
}

#' @noRd
#'
#' @description
#' Helper function of [gen_feat_mat()] to extract the deconvolution parameters
#' from `data_path`, where `data_path` can be a `decon1` or `decons1` object or
#' a folder containing a `"parameters.txt"` and `"approximated_spectrum.txt"`
#' file, as created when calling `MetaboDecon1D()` before version 1.2.
#'
#' @param data_path
#' A list of deconvoluted spectra as returned by [generate_lorentz_curves()] or
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
        dd <- if (is_decon_list(dd)) dd else if (is_decon_obj(dd)) list(dd) else stop("data_path must be a list of deconvoluted spectra.")
        w <- lapply(dd, function(d) d$x_0)
        lambda <- lapply(dd, function(d) d$lambda)
        A <- lapply(dd, function(d) d$A)
        spectrum_superposition <- lapply(dd, function(d) d$spectrum_superposition)
        params <- named(w, lambda, A, spectrum_superposition)
    } else {
        if (!file.exists(dd)) stop(dd, " does not exist.") else if (warn) warning("You have provided a path to `gen_feat_mat()`. Since metabodecon v1.2 it is recommended to provide the output of `generate_lorentz_curves()` directly to speed up computations. For details see section 'Details' after calling `help('gen_feat_mat')`.")
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
    mapply(par_nam, spc_nam, FUN = function(p, s) if (p != s) warning(sprintf("Mismatched file names: %s and %s\n", p, s)) )
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
#' | F   | [TF] | [TF] | Return                                     |
#' | T   | T    | T    | Ask user and install or return accordingly |
#' | T   | T    | F    | Stop with error                            |
#' | T   | F    | [TF] | Install                                    |
#'
#' req = Installation required (because requested packages are not yet installed
#' on the system)
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
