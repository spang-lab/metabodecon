# Exported Main #####

#' @export
#'
#' @title Align Spectra
#'
#' @description
#' Align signals across a list of deconvoluted spectra using the 'CluPA'
#' algorithm from the 'speaq' package, described in Beirnaert et al. (2018)
#' <doi:10.1371/journal.pcbi.1006018> and Vu et al. (2011)
#' <doi:10.1186/1471-2105-12-405> plus the additional peak combination described
#' in [combine_peaks()].
#'
#' @param decons An object of type `MD1Decons` or `GLCDecons` as returned by
#' [MetaboDecon1D()] or [generate_lorentz_curves()].
#'
#' @param maxShift Maximum number of points along the "ppm-axis" a value can be
#' moved by the 'speaq' package. 50 is a suitable starting value for plasma
#' spectra with a digital resolution of 128K. Note that this parameter has to be
#' individually optimized depending on the type of analyzed spectra and the
#' digital resolution. For urine which is more prone to chemical shift
#' variations this value most probably has to be increased. Passed as argument
#' `maxShift` to [speaq_align()].
#'
#' @param maxCombine Amount of adjacent columns which may be combined for
#' improving the alignment. Passed as argument `range` to [combine_peaks()].
#'
#' @return
#' An object of type `AlignedSpectra` as described in [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
#' @examples
#' sim_dir <- metabodecon_file("bruker/sim")
#' spectra <- read_spectra(sim_dir)
#' decons <- generate_lorentz_curves_sim(spectra)
#' aligned <- align_spectra(decons)
#' str(aligned)
align_spectra <- function(decons, maxShift = 50, maxCombine = 5) {
    smat <- speaq_align(spectrum_data = decons, maxShift = 50)
    tmp <- combine_peaks(
        shifted_mat = smat,
        range = maxCombine,
        spectrum_data = decons
    )
}

# Exported Helpers #####

#' @export
#'
#' @title Get PPM Range covered by Spectra
#'
#' @description
#' Returns the ppm range across all peaks of the provided deconvoluted spectra.
#'
#' @param spectrum_data A list of deconvoluted spectra as returned by
#' [generate_lorentz_curves()].
#'
#' @param show Whether to plot the ppm range on the spectrum plot.
#'
#' @param mar The margins of the plot. Passed on to `par()`.
#'
#' @return
#' A vector containing the lowest and highest ppm value over all peaks of the
#' provided deconvoluted spectra.
#'
#' @author
#' Initial version from Wolfram Gronwald.
#' Refactored by Tobias Schmidt in 2024.
#'
#' @examples
#' sim_subset <- metabodecon_file("sim_subset")
#' spectrum_data = generate_lorentz_curves_sim(sim_subset)
#' ppm_rng <- get_ppm_range(spectrum_data, show = TRUE)
#' print(ppm_rng)
get_ppm_range <- function(spectrum_data = generate_lorentz_curves_sim(),
                          show = FALSE,
                          mar = c(4.1, 4.1, 1.1, 0.1)) {
    msg <- "spectrum_data must be a list of deconvoluted spectra."
    ss <- spectrum_data
    ss <- if (is_decon_list(ss)) ss else if (is_decon_obj(ss)) list(ss) else stop(msg)
    ppm_min <- min(sapply(ss, function(s) min(s$peak_triplets_middle)))
    ppm_max <- max(sapply(ss, function(s) max(s$peak_triplets_middle)))
    ppm_rng <- c(ppm_min, ppm_max)
    if (show) plot_spectra(ss, mar = mar, peak_rng = ppm_rng)
    ppm_rng
}

#' @export
#'
#' @title Generate Feature Matrix.
#'
#' @description Generate a feature matrix.
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
#' @details Before version 1.2 of 'metabodecon', the deconvolution functions
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
#' `w`: A list of vectors where each vector contains the "position paramater"
#'      of the peaks in the corresponding spectrum.
#'
#' `A`: A list of vectors where each vector contains the "area parameter" of
#'      the peaks in the corresponding spectrum.
#'
#' `lambda`: A list of vectors where each vector contains the "width
#'           parameter" of the peaks in the corresponding spectrum.
#'
#' @author
#' Initial version from Wolfram Gronwald.
#' Refactored by Tobias Schmidt in 2024.
#'
#' @examples
#' sim_subset <- metabodecon_file("sim_subset")
#' decons <- generate_lorentz_curves_sim(sim_subset)
#' obj <- gen_feat_mat(decons)
#' str(obj, 2, give.attr = FALSE)
gen_feat_mat <- function(data_path = generate_lorentz_curves_sim(),
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
#' postion `i, j` holds the integral value of the signal from spectrum `i` that
#' has its center at position `j` after alignment by speaq. If there is no
#' signal with center `j` in spectrum `i`, entry `i, j` is set to NA. The column
#' names of the matrix are the ppm values of the original spectra.
#'
#' Example return matrix:
#'
#' ```
#'     ...  3.59  3.58  3.57  3.56  3.55  3.54  3.53
#'   .----------------------------------------------> PPM
#' 1 | NA   NA    0.20  NA    NA    NA    0.25  NA
#' 2 | NA   NA    0.15  NA    NA    NA    0.13  NA
#' 3 | NA   NA    NA    0.2   NA    NA    0.18  NA
#' SpNr
#' ```
#'
#' @author
#' Initial version from Wolfram Gronwald.
#' Refactored by Tobias Schmidt in 2024.
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
                        spectrum_data = generate_lorentz_curves_sim(),
                        si_size_real_spectrum = length(spectrum_data[[1]]$y_values),
                        verbose = TRUE,
                        show = FALSE,
                        mfrow = c(2, 1)) {
    acceptLostPeak <- FALSE # Must be FALSE so we can assign A and lambda to the respective peaks
    Y <- feat$data_matrix
    s <- nrow(Y) # Number of spectra
    n <- ncol(Y) # Number of data points
    U <- feat$peakList # Unaligned peak centers indices (PCIs) in increasing order (i.e. the leftmost value index 1).
    r <- speaq::findRef(U)$refInd # Reference spectrum
    C <- dohCluster(Y, U, r, maxShift, acceptLostPeak, verbose) # list(Y, new_peakList)
    Y2 <- C$Y
    Q <- C$new_peakList # Aligned PCIs in increasing order.
    P <- lapply(Q, function(p) n - p) # Aligned PCIs in decreasing order.
    integrals <- mapply(SIMPLIFY = FALSE, feat$A, feat$lambda, P, FUN = function(A, l, p) {
        integral <- A * (atan((n - p) / l) - atan((0 - p) / l))
        integral[p == 0] <- NA
        integral
    })
    ppm <- spectrum_data[[1]]$x_values_ppm
    M <- matrix(nrow = s, ncol = n, dimnames = list(NULL, ppm))
    for (i in seq_len(s)) M[i, round(Q[[i]]) - 1] <- integrals[[i]]
    if (show) {
        opar <- par(mfrow = mfrow, mar = c(5.1, 4.1, 2.1, 0.1))
        on.exit(par(opar), add = TRUE)
        plot_si_mat(Y = feat$data_matrix, main = "Original Spectra")
        plot_si_mat(Y = C$Y, main = "Aligned Spectra")
    }
    return(M)
}

#' @export
#'
#' @title Combine Peaks
#'
#' @description Even after calling `speaq_align()`, the alignment of individual
#'
#' signals is not always perfect, as 'speaq' performs a segment-wise alignment
#' i.e. groups of signals are aligned. For further improvements, partly filled
#' neighboring columns are merged.
#'
#' @param shifted_mat The matrix returned by `speaq_align()`.
#'
#' @param range Amount of adjacent columns which are permitted to be used for
#'
#' improving the alignment.
#'
#' @param lower_bound Minimum amount of non-zero elements per column to trigger
#'
#' the alignment improvement.
#'
#' @param spectrum_data The list of deconvoluted spectra as returned by
#'
#' `generate_lorentz_curves()` that was used to generate `shifted_mat`. No
#' longer required since version 1.2 of Metabodecon.
#'
#' @param data_path If not NULL, the returned dataframes `long` and `short` are
#'
#' written to `data_path` as "aligned_res_long.csv" and "aligned_res_short.csv".
#'
#' @return A list containing two data frames `long` and `short`. The first data
#'
#' frame contains one one column for each data point in the original spectrum.
#' The second data frame contains only columns where at least one entry is
#' non-zero.
#'
#' @details
#'
#' Example of what the function does:
#'
#' ```txt
#' |            | 3.56 | 3.54 | 3.51 | 3.51 | 3.50 |
#' |----------- |------|------|------|------|------|
#' | Spectrum 1 | 0.13 | 0    | 0.11 | 0    | 0    |
#' | Spectrum 2 | 0.13 | 0    | 0.12 | 0    | 0    |
#' | Spectrum 3 | 0.07 | 0    | 0    | 0    | 0    |
#' | Spectrum 4 | 0.08 | 0    | 0    | 0.07 | 0    |
#' | Spectrum 5 | 0.04 | 0    | 0.04 | 0    | 0    |
#' ```
#'
#' becomes
#'
#' ```txt
#' |            | 3.56 | 3.54 | 3.51 | 3.50 |
#' |----------- |------|------|------|------|
#' | Spectrum 1 | 0.13 | 0    | 0.11 | 0    |
#' | Spectrum 2 | 0.13 | 0    | 0.12 | 0    |
#' | Spectrum 3 | 0.07 | 0    | 0    | 0    |
#' | Spectrum 4 | 0.08 | 0    | 0.07 | 0    |
#' | Spectrum 5 | 0.04 | 0    | 0.04 | 0    |
#' ```
#'
#' I.e. column 3 and 4 get merged, because they are in `range` of each other
#' and have no common non-zero entries.
#'
#' @author Initial version from Wolfram Gronwald. Refactored by Tobias Schmidt
#'
#' in 2024.
#'
#' @examples
#'
#' spectrum_data <- generate_lorentz_curves_sim("bruker/sim")
#' shifted_mat <- speaq_align(spectrum_data = spectrum_data, verbose = FALSE)
#' range <- 5
#' lower_bound <- 1
#' obj <- combine_peaks(shifted_mat, range, lower_bound)
#' str(obj)
combine_peaks <- function(shifted_mat = speaq_align(),
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


# Private Helpers #####

#' @export
#'
#' @title Cluster Based Peak Alignment
#'
#' @description
#' Rewrite of `speaq::dohCluster()`, compatible with the data format returned by
#' 'generate_lorentz_curves()' and 'gen_feat_mat()'. The function name
#' "dohCluster" comes from "Do Hierarchical Clustering" which is part of the
#' Alignment algorithm proposed by Vu et al. (2011) in
#' <doi:10.1186/1471-2105-12-405>.
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
#' Initial version from Wolfram Gronwald.
#' Refactored by Tobias Schmidt in 2024.
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

#' @author Wolfram Gronwald.
#'
#' @noRd
#'
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
        cat("\n Alignment time: ", (endTime[3] - startTime[3]) / 60," minutes")
    }
    # Added by me at 31.08.21
    return_list <- list("Y" = bestY, "new_peakList" = peakListNew)
    return(return_list)
}

#' @author Wolfram Gronwald.
#'
#' @noRd
#'
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
#' @description Helper function of [gen_feat_mat()] to extract the deconvolution
#' parameters from `data_path`, where `data_path` can be:
#'
#' 1. A single devonvolved spectrum as obtained by calling
#'    `generate_lorentz_curves(spectrum)`
#' 2. A list of deconvoluted spectra as obtained by calling
#'    `generate_lorentz_curves(spectra)`
#' 3. A folder containing ".* parameters.txt" and ".* approximated_spectrum.txt"
#'    files, as created by when calling `MetaboDecon1D()` before version 1.2.
#'
#' @param data_path A list of deconvoluted spectra as returned by
#' [generate_lorentz_curves()] or a path to a folder containing ".*
#' parameters.txt" and ".* approximated_spectrum.txt" files.
#'
#' @param warn (logical) Whether to print warning in case a file path is
#' provided instead of a list of deconvoluted spectra.
#'
#' @param check (logical) Whether to sanity check the deconvolution parameters
#' before returning them.
#'
#' @author Tobias Schmidt
#'
#' @return A list containing the deconvolution parameters, i.e. `w`, `lambda`,
#' `A`, and `spectrum_superposition`.
#'
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
        if (!file.exists(dd)) stop(dd, " does not exist.") else if (warn) warning("You have provided a path to `gen_feat_mat()`. Since metabodecon v1.2 it is recommended to provide the output of `generate_lorentz_curves()` directly instead to speed up the computations. For details see section 'Details' after calling `help('gen_feat_mat')`.")
        params <- read_decon_params(dd)
    }
    if (check) check_decon_params(params)
    params
}

#' @author Tobias Schmidt
#'
#' @noRd
#'
read_decon_params <- function(data_path) {
    par_txt <- dir(data_path, "(.*) parameters.txt", full.names = TRUE)
    spc_txt <- dir(data_path, "(.*) approximated_spectrum.txt", full.names = TRUE)
    par_nam <- sub(" parameters.txt", "", basename(par_txt))
    spc_nam <- sub(" approximated_spectrum.txt", "", basename(spc_txt))
    if (length(par_txt) != length(spc_txt)) stop("Number of parameter files and spectrum files differs.")
    if (length(par_txt) == 0) stop("No parameter files found in the given directory.")
    tmp <- mapply(par_nam, spc_nam, FUN = function(p, s) if (p != s) warning(sprintf("Mismatched file names: %s and %s\n", p, s)) )
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

#' @author Tobias Schmidt
#'
#' @noRd
#'
check_decon_params <- function(params) {
    nulls <- unlist(sapply(params, function(pp) sapply(pp, is.null)))
    if (any(nulls)) warning("Detected missing params: ", names(nulls)[nulls])
}

#' @author Tobias Schmidt
#'
#' @noRd
#'
rm_zero_width_peaks <- function(params) {
    for (i in seq_len(params$A)) {
        not_zero <- params$w[[i]] != 0
        params$lambda[[i]] <- params$lambda[[i]][not_zero]
        params$w[[i]] <- params$w[[i]][not_zero]
        params$A[[i]] <- params$A[[i]][not_zero]
    }
    params
}

#' @author Tobias Schmidt
#'
#' @noRd
#'
is_decon_obj <- function(x) {
    keys <- c(
        "number_of_files", "filename", "x_values", "x_values_ppm",
        "y_values", "spectrum_superposition", "mse_normed", "index_peak_triplets_middle",
        "index_peak_triplets_left", "index_peak_triplets_right", "peak_triplets_middle",
        "peak_triplets_left", "peak_triplets_right", "integrals", "signal_free_region",
        "range_water_signal_ppm", "A", "lambda", "x_0"
    )
    if (is.list(x) && all(keys %in% names(x))) TRUE else FALSE
}

#' @author Tobias Schmidt
#'
#' @noRd
#'
is_decon_list <- function(x) {
    if (is.list(x) && all(sapply(x, is_decon_obj))) TRUE else FALSE
}

# Deprecated #####

#' @author Initial version written by from Wolfram Gronwald as part of
#' `gen_feat_mat`. Moved into seperate function in 2024 by Tobias Schmidt.
#'
#' @noRd
#'
read_decon_params_original <- function(data_path) {
    files <- list.files(data_path, ".txt", full.names = TRUE)
    num_spectra <- length(files) / 2
    spectrum_superposition <- vector(mode = "list", length = num_spectra)
    data_info <- vector(mode = "list", length = num_spectra)
    numerator_info <- 1
    numerator_spectrum <- 1
    for (i in 1:length(files)) {
        if (grepl(" parameters", files[i])) {
            data_info[[numerator_info]] <- as.matrix(data.table::fread(files[i], header = FALSE))
            numerator_info <- numerator_info + 1
        }
        if (grepl(" approximated_spectrum", files[i])) {
            spectrum_superposition[[numerator_spectrum]] <- as.vector(unlist(data.table::fread(files[i], header = FALSE)))[-1]
            numerator_spectrum <- numerator_spectrum + 1
        }
    }
    w <- vector(mode = "list", length = num_spectra)
    lambda <- vector(mode = "list", length = num_spectra)
    A <- vector(mode = "list", length = num_spectra)
    for (i in 1:num_spectra) {
        w[[i]] <- as.numeric(data_info[[i]][1, ][-1])
        lambda[[i]] <- as.numeric(data_info[[i]][2, ][-1])
        A[[i]] <- as.numeric(data_info[[i]][3, ][-1])
   }
    return(named(w, lambda, A, spectrum_superposition))
}

#' @noRd
#'
#' @title Align Signals using 'speaq'
#'
#' @description Performs signal alignment across the individual spectra using
#' the 'speaq' package (Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S,
#' Pieters L, et al. (2018) speaq 2.0: A complete workflow for high-throughput
#' 1D NMRspectra processing and quantification. PLoS Comput Biol 14(3):
#' e1006018. https://www.doi.org/10.1371/journal.pcbi.1006018). The spectra
#' deconvolution process yields the signals of all spectra. Due to slight
#' changes in measurement conditions, e.g. pH variations, signal positions may
#' vary slightly across spectra. As a consequence, prior to further analysis
#' signals belonging to the same compound have to be aligned across spectra.
#' This is the purpose of the 'speaq' package.
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
#' @return A matrix containing the aligned integral values of all spectra. Each
#' row contains the data of each spectrum and each column corresponds to one
#' data point. Each entry corresponds to the integral of a deconvoluted signal
#' with the signal center at this specific position after alignment by speaq.
#'
#' @author Wolfram Gronwald
#'
#' @examples
#'
#' spectrum_data <- generate_lorentz_curves_sim("bruker/sim")
#' feat <- gen_feat_mat(spectrum_data)
#' maxShift <- 100
#' si_size_real_spectrum <- length(spectrum_data[[1]]$y_values)
#' after_speaq_mat <- speaq_align(feat, maxShift, spectrum_data, si_size_real_spectrum)
speaq_align_original <- function(feat = gen_feat_mat(spectrum_data),
                           maxShift = 50,
                           spectrum_data = generate_lorentz_curves_sim(),
                           si_size_real_spectrum = length(spectrum_data[[1]]$y_values)) {

    # Find reference spectrum, i.e the spectrum to which all others will be aligned to
    resFindRef <- speaq::findRef(feat$peakList)
    refInd <- resFindRef$refInd

    # Use overwritten CluPA function for multiple spectra
    data_result_aligned <- dohCluster(
        X = feat$data_matrix,
        peakList = feat$peakList,
        refInd = refInd,
        maxShift = maxShift,
        acceptLostPeak = FALSE,
        verbose = TRUE
    )
    # > str(data_result_aligned, 2)
    # List of 2
    # $ Y: num [1:16, 1:1309] 0.013 0.024 ... (matrix of aligned spectra)
    # $ new_peakList: List of 16              (aligned peaks of each spectrum)
    # ..$ sim_01: num [1:13] 350 416 457 ...
    # ..$ sim_02: num [1:11] 416 458 542 ...
    # ...

    # Aligned spectra values
    data_matrix_aligned <- data_result_aligned$Y
    # new peak list values
    new_peakList <- data_result_aligned$new_peakList
    # Recalculate peakList for the original format
    num_spectra <- dim(data_matrix_aligned)[1]
    peakList_result <- vector(mode = "list", length = num_spectra)
    for (i in 1:length(new_peakList)) {
        peakList_result[[i]] <- (dim(feat$data_matrix)[2]) - new_peakList[[i]]
    }

    # Generate feature matrix and integral matrix by using distance matrix
    # message("Warning: Generation of feature matrix and integral matrix have a high time duration (about 4h!)")
    # Calculate integral results for each positions
    integral_results <- vector(mode = "list", length = num_spectra)
    # Calculate integral_results
    for (spec in 1:length(peakList_result)) {
        for (ent in 1:length(peakList_result[[spec]])) {
            # Calculate integral value for each entry which is unequal 0
            if (peakList_result[[spec]][ent] != 0) {
                integral_results[[spec]][ent] <- feat$A[[spec]][ent] * (atan((-peakList_result[[spec]][ent] + si_size_real_spectrum) / feat$lambda[[spec]][ent]) - atan((-peakList_result[[spec]][ent]) / feat$lambda[[spec]][ent]))
            }
        }
    }

    # Generate a matrix containing the integrals of all spectra at their positions after alignment by speaq
    after_speaq_mat <- matrix(nrow = num_spectra, ncol = si_size_real_spectrum)
    # set ppm values for data points of matrix based on original data
    # Note signals will be shifted across data points, data points itself will keep their positions.
    colnames(after_speaq_mat) <- spectrum_data[[1]]$x_values_ppm
    mid_spec <- round(si_size_real_spectrum / 2)
    for (i in 1:num_spectra)
    {
        for (j in 1:length(peakList_result[[i]]))
        {
            k <- round(peakList_result[[i]][j])
            k1 <- mid_spec - (k - mid_spec) # mirror imaging of data at mid data point
            after_speaq_mat[i, k1] <- integral_results[[i]][j]
        }
    }
    return(after_speaq_mat)
}
