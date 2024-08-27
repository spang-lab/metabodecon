# Exported Main #####

#' @export
#' @title Get PPM Range covered by Spectra
#' @description Returns the ppm range across all peaks of the provided deconvoluted spectra.
#' @param spectrum_data A list of deconvoluted spectra as returned by [generate_lorentz_curves()].
#' @param show Whether to plot the ppm range on the spectrum plot.
#' @param mar The margins of the plot. Passed on to `par()`.
#' @return A vector containing the lowest and highest ppm value over all peaks of the provided deconvoluted spectra.
#' @author Initial version from Wolfram Gronwald, 2021. Refactored by Tobias Schmidt in 2024.
#' @examples
#' spectrum_data = glc_sim()
#' ppm_rng <- get_ppm_range(spectrum_data, show = TRUE)
#' print(ppm_rng)
get_ppm_range <- function(spectrum_data = glc_sim(),
                          show = FALSE,
                          mar = c(4.1, 4.1, 1.1, 0.1)) {
    msg <- "spectrum_data must be a list of deconvoluted spectra."
    ss <- spectrum_data
    ss <- if (is_decon_list(ss)) ss else if (is_decon_obj(ss)) list(ss) else stop(msg)
    ppm_min <- min(sapply(ss, function(s) min(s$peak_triplets_middle)))
    ppm_max <- max(sapply(ss, function(s) max(s$peak_triplets_middle)))
    ppm_rng <- c(ppm_min, ppm_max)
    if (show) plot_spectra(ss, mar = mar, ppm_rng = ppm_rng)
    ppm_rng
}

#' @export
#' @title Generate Feature Matrix.
#' @description Generate a feature matrix.
#' @param data_path A list of deconvoluted spectra as returned by `generate_lorentz_curves()`. In older versions, this could also be the path passed to `generate_lorentz_curves()`, but this is deprecated and will trigger a warning. See 'Details' for more information.
#' @param ppm_range The ppm range over which your signals are distributed.
#' @param si_size_real_spectrum Number of data points in your spectra.
#' @param scale_factor_x The x scale factor used during the deconvolution.
#' @param warn Whether to print a warning in case a file path is passed to `data_path` instead of a list of deconvoluted spectra.
#' @details Before version 1.2 of 'metabodecon', the deconvolution functions `generate_lorentz_curves` and `MetaboDecon1D` wrote their output partially as txt files to their input folder. Back then, `gen_feat_mat()` used those txt files as input to generate the feature matrix. Since version 1.2 these txt files are no longer created by default, to prevent accidental modifications of the input folders. Therefore, the recommended way to pass the required information to `gen_feat_mat()` is to directly pass the output of `generate_lorentz_curves()` to `gen_feat_mat()`. However, to stay backwards compatible, the name of parameter `data_path` was not changed and passing an actual path to `data_path` is still possible, but will result in a warning (unless `warn` is set to `FALSE`).
#' @return A list with the following elements:
#' - `data_matrix`: A data.frame where each row corresponds to one spectrum and each column to one data point, i.e. for 10 input spectra with 131072 data points each `data_matrix` would have dimensions 10 x 131072.
#' - `peakList`: A list of vectors, where each vector contains the indices of the peaks in the corresponding spectrum. The indices increase from left to right, i.e. the smallest index corresponds to the highest ppm value, as the ppm values decrease from left to right.
#' - `w`: A list of vectors where each vector contains the "position paramater" of the peaks in the corresponding spectrum.
#' - `A`: A list of vectors where each vector contains the "area parameter" of the peaks in the corresponding spectrum.
#' - `lambda`: A list of vectors where each vector contains the "width parameter" of the peaks in the corresponding spectrum.
#' @examples
#' sim <- metabodecon_file("bruker/sim_subset")
#' decons <- generate_lorentz_curves(
#'      sim,
#'      sfr = c(3.58, 3.42), wshw = 0, delta = 0.1,
#'      ask = FALSE, verbose = FALSE
#' )
#' obj <- gen_feat_mat(decons)
#' str(obj, 2, give.attr = FALSE)
gen_feat_mat <- function(data_path = generate_lorentz_curves(),
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
    # _____  _____  PPPPP  _____  _____  _____  PPPPP  _____
    # 0.007  0.006  0.005  0.004  0.003  0.002  0.001  0.000
    # 1      2      3      4      5      6      7      8
    #
    # In the code line above, ncol(X) corresponds to 8, scale_factor_x to 1000
    # and D$w[[i]] to c(0.005, 0.001).
    list(data_matrix = X, peakList = P, w = D$w, A = D$A, lambda = D$lambda)
}

#' @export
#' @title Align Signals using 'speaq'
#' @description Performs signal alignment across the individual spectra using the 'speaq' package (Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S, Pieters L, et al. (2018) speaq 2.0: A complete workflow for high-throughput 1D NMRspectra processing and quantification. PLoS Comput Biol 14(3): e1006018. https://www.doi.org/10.1371/journal.pcbi.1006018). The spectra deconvolution process yields the signals of all spectra. Due to slight changes in measurement conditions, e.g. pH variations, signal positions may vary slightly across spectra. As a consequence, prior to further analysis signals belonging to the same compound have to be aligned across spectra. This is the purpose of the 'speaq' package.
#' @param feat Output of `gen_feat_mat()`.
#' @param maxShift Maximum number of points along the "ppm-axis" which a value can be moved by speaq package e.g. 50. 50 is a suitable starting value for plasma spectra with a digital resolution of 128K. Note that this parameter has to be individually optimized depending on the type of analyzed spectra and the digital resolution. For urine which is more prone to chemical shift variations this value most probably has to be increased.
#' @param spectrum_data Output of `generate_lorentz_curves()`.
#' @param si_size_real_spectrum Number of real data points in your original spectra.
#' @param verbose Whether to print additional information during the alignment process.
#' @param show Whether to plot the original and aligned spectra.
#' @param mfcol Layout to use for the plot. Passed on to `par()`. Use `mfcol = NULL` if the plot layout should not be changed.
#' @return A matrix containing the integral values of the spectra after alignment.
#' There is one row per spectrum and one column per ppm value.
#' The entry at postion `i, j` holds the integral value of the signal from spectrum `i` that has its center at position `j` after alignment by speaq.
#' If there is no signal with center `j` in spectrum `i`, entry `i, j` is set to NA.
#' The column names of the matrix are the ppm values of the original spectra.
#'
#' Example return matrix:
#'
#' ```
#'      ...  3.59  3.58  3.57  3.56  3.55  3.54  3.53
#'    .----------------------------------------------> PPM
#'  1 | NA   NA    0.20  NA    NA    NA    0.25  NA
#'  2 | NA   NA    0.15  NA    NA    NA    0.13  NA
#'  3 | NA   NA    NA    0.2   NA    NA    0.18  NA
#'  SpNr
#' ```
#'
#' Spectrum 3 has a signal with an integral value of 0.2 at
#' position 3.57 after alignment. Spectrum 1 and 2 have no
#' signal at this position after alignment.
#' ```
#'
#' @examples
#' spectrum_data <- glc_sim("bruker/sim")
#' feat <- gen_feat_mat(spectrum_data)
#' maxShift <- 200
#' M <- speaq_align(feat, maxShift, spectrum_data, show = TRUE)
#' str(M)
speaq_align <- function(feat = gen_feat_mat(spectrum_data),
                        maxShift = 50,
                        spectrum_data = glc_sim(),
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

#' @noRd
#' @title Combine Peaks
#' @description Even after calling `speaq_align()`, the alignment of individual signals is not always perfect, as 'speaq' performs a segment-wise alignment i.e. groups of signals are aligned. For further improvements, partly filled neighboring columns are merged.
#' @param shifted_mat The matrix returned by `speaq_align()`.
#' @param range Amount of adjacent columns which are permitted to be used for improving the alignment.
#' @param lower_bound Amount of columns that need to be skipped (e.g. because they contain rownames instead of values). Only modify in case of errors.
#' @param spectrum_data The list of deconvoluted spectra as returned by `generate_lorentz_curves()`.
#' @param data_path If not NULL, the returned dataframes `long` and `short` are written to `data_path` as "aligned_res_long.csv" and "aligned_res_short.csv".
#' @return A list containing two data frames `long` and `short`. The first data frame contains one one column for each data point in the original spectrum. The second data frame contains only columns where at least one entry is non-zero.
#' @details
#' Example:
#'
#' ```txt
#'               3.56   3.54  3.51  3.51   3.50
#' Spectrum 1  | 0.13   0     0.11  0      0
#' Spectrum 2  | 0.13   0     0.12  0      0
#' Spectrum 3  | 0.07   0     0     0      0
#' Spectrum 4  | 0.08   0     0     0.07   0
#' Spectrum 5  | 0.04   0     0.04  0      0
#'
#' becomes
#'
#'               3.54   3.54  3.51  3.50
#' Spectrum 1  | 0.13   0     0.11  0
#' Spectrum 2  | 0.13   0     0.12  0
#' Spectrum 3  | 0.07   0     0     0
#' Spectrum 4  | 0.08   0     0.07  0
#' Spectrum 5  | 0.04   0     0.04  0
#'
#' I.e. column 3 and 4 get merged, because they are adjacent
#' and have no common non-zero entries.
#' ```
#'
#' @examples
#' shifted_mat <- speaq_align(spectrum_data = glc_sim("bruker/sim"))
#' M2 <- combine_peaks(shifted_mat, spectrum_data = glc_sim("bruker/sim"))
combine_peaks <- function(shifted_mat = speaq_align(spectrum_data = spectrum_data, show = TRUE),
                          range = 5,
                          lower_bound = 1,
                          spectrum_data = glc_sim("bruker/sim"),
                          data_path = NULL) {
    S <- replace(shifted_mat, is.na(shifted_mat), 0)
    nz <- numeric(ncol(S))
    for (i in 1:ncol(S)) nz[i] <- length(which(S[, i] != 0))
    bla <- S != 0
    for (i in (nrow(S) - 1):lower_bound) {
        for (j in which(nz == i)) {
            if (nz[j] != 0) {
                while (TRUE) {
                    bli <- rep(0, 2 * range)
                    l <- 0
                    search_range <- c((j - range):(j - 1), (j + 1):(j + range))
                    for (k in search_range) {
                        l <- l + 1
                        if (all(which(bla[, k]) %in% which(!bla[, j]))) {
                            bli[l] <- nz[k]
                        }
                    }
                    if (all(bli == 0)) {
                        break
                    } else {
                        index_max <- search_range[which.max(bli)]

                        S[, j] <- S[, j] + S[, index_max]
                        bla[, j] <- bla[, j] | bla[, index_max]
                        nz[j] <- nz[j] + nz[index_max]
                        S[, index_max] <- 0
                        bla[, index_max] <- F
                        nz[index_max] <- 0
                    }
                }
            }
        }
    }
    colnames(S) <- as.character(round(spectrum_data[[1]]$x_values_ppm, digits = 4))
    # removal of columns containing only zeros
    kick.out <- apply(S, 2, function(z) {
        all(z == 0)
    })
    shifted_mat_no_na_s <- S[, !kick.out]
    # utils::write.csv2(S, file = file.path(data_path, "aligned_res_long.csv"))
    # utils::write.csv2(shifted_mat_no_na_s, file = file.path(data_path, "aligned_res_short.csv"))
    return_list <- list("short" = shifted_mat_no_na_s, "long" = S)
    return(return_list)
}

# Private Helpers #####

#' @export
#' @title Cluster Based Peak Alignment
#' @description Rewrite of `speaq::dohCluster()`, compatible with the data format returned by 'generate_lorentz_curves()' and 'gen_feat_mat()'. The function name "dohCluster" comes from "Do Hierarchical Clustering" which is part of the Alignment algorithm proposed by Vu et al. (2011) in <doi:10.1186/1471-2105-12-405>.
#' @param X Dataframe of signal intensities from all spectra as returned by [gen_feat_mat()].
#' @param peakList List of peak indices as returned [gen_feat_mat()].
#' @param refInd Number of the reference spectrum i.e. the spectrum to which all signals will be aligned to.
#' @param maxShift Maximum number of points a value can be moved.
#' @param acceptLostPeak Whether to allow the the alignment algorithm to ignore peaks that cannot easily be aligned with the reference spectrum.
#' @param verbose Whether to print additional information during the alignment process.
#' @return A list containing two data frames `Y` and `new_peakList`. The first one contains the aligned spectra, the second one contains the aligned signals of each spectrum.
#' @examples
#' decons <- glc_sim()
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

#' @noRd
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

#' @noRd
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
#' @description Helper function of [gen_feat_mat()] to extract the deconvolution parameters from `data_path`, where `data_path` can be:
#' 1. A single devonvolved spectrum as obtained by calling `generate_lorentz_curves(spectrum)`
#' 2. A list of deconvoluted spectra as obtained by calling `generate_lorentz_curves(spectra)`
#' 3. A folder containing ".* parameters.txt" and ".* approximated_spectrum.txt" files, as created by when calling `MetaboDecon1D(dd)`
#' @param data_path A list of deconvoluted spectra as returned by [generate_lorentz_curves()] or a path to a folder containing ".* parameters.txt" and ".* approximated_spectrum.txt" files.
#' @param warn (logical) Whether to print warning in case a file path is provided instead of a list of deconvoluted spectra.
#' @param check (logical) Whether to sanity check the deconvolution parameters before returning them.
#' @return A list containing the deconvolution parameters, i.e. `w`, `lambda`, `A`, and `spectrum_superposition`.
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
        params <- read_decon_params_v2(dd)
    }
    if (check) check_decon_params(params)
    params
}

#' @noRd
read_decon_params_v2 <- function(data_path) {
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

#' @noRd
check_decon_params <- function(params) {
    nulls <- unlist(sapply(params, function(pp) sapply(pp, is.null)))
    if (any(nulls)) warning("Detected missing params: ", names(nulls)[nulls])
}

#' @noRd
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

#' @noRd
is_decon_list <- function(x) {
    if (is.list(x) && all(sapply(x, is_decon_obj))) TRUE else FALSE
}

plot_spectra <- function(ss = glc_sim(),
                         mar = c(4.1, 4.1, 1.1, 0.1),
                         peak_rng = get_ppm_range(ss, show = FALSE)) {
    if (is_decon_obj(ss)) ss <- list(ss)
    if (is_decon_list(ss)) {
        xrng <- range(c(sapply(ss, function(s) s$x_values_ppm)))
        ymax <- max(sapply(ss, function(s) max(s$y_values)))
    } else if (is.data.frame(ss)) {
        stop("ss must be a list of deconvoluted spectra.")
    }
    a <- peak_rng[1]
    b <- peak_rng[2]
    w <- (b - a) / 4
    y8 <- ymax * 0.8
    cols <- rainbow(length(ss))
    ltxt <- paste("Spectrum", 1:length(ss))
    opar <- par(mar = mar)
    on.exit(par(opar))
    plot(
        x = NA,
        type = "n",
        xlim = xrng[2:1],
        ylim = c(0, ymax),
        xlab = "Chemical Shift [ppm]",
        ylab = "Signal Intensity [au]"
    )
    abline(v = peak_rng, lty = 2)
    for (i in 1:length(ss)) {
        lines(
            x = ss[[i]]$x_values_ppm,
            y = ss[[i]]$y_values,
            col = cols[i]
        )
    }
    arrows(
        x0 = c(a + w, b - w),
        x1 = c(a, b),
        y0 = y8,
        y1 = y8,
        length = 0.1,
        lty = 2,
        col = "black"
    )
    text(
        x = mean(c(a, b)),
        y = y8,
        labels = "ppm range"
    )
    mtext(
        text = round(c(a, b), 4),
        side = 3,
        line = 0,
        at = c(a, b)
    )
    legend(
        x = "topright",
        legend = ltxt,
        col = cols,
        lty = 1
    )
}

plot_si_mat <- function(Y = glc_sim("bruker/sim"),
                        lgdcex = "auto",
                        main = NULL,
                        mar = par("mar")) {
    n <- nrow(Y) # Number of Spectra
    p <- ncol(Y) # Number of Datapoints
    cols <- rainbow(n) # Colors
    dpis <- seq_len(p) # Datapoint Indices
    spis <- seq_len(n) # Spectrum Indices
    ltxt <- paste("Spectrum", spis)
    xlab <- "Datapoint Number"
    ylab <- "Signal Intensity"
    xlim <- c(1, p)
    ylim <- c(0, max(Y))
    args <- named(x = NA, type = "n", xlim, ylim, xlab, ylab)
    opar <- par(mar = mar)
    on.exit(par(opar))
    do.call(plot, args)
    for (i in spis) lines(x = dpis, y = Y[i, ], col = cols[i])
    if (lgdcex == "auto") lgdcex <- 1 / max(1, log(n, base = 8))
    legend(x = "topright", legend = ltxt, col = cols, lty = 1, cex = lgdcex)
    if (!is.null(main)) title(main = main)
}

# Deprecated #####

#' @noRd
read_decon_params_v1 <- function(data_path) {
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
#' @title Align Signals using 'speaq'
#' @description Performs signal alignment across the individual spectra using the 'speaq' package (Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S, Pieters L, et al. (2018) speaq 2.0: A complete workflow for high-throughput 1D NMRspectra processing and quantification. PLoS Comput Biol 14(3): e1006018. https://www.doi.org/10.1371/journal.pcbi.1006018). The spectra deconvolution process yields the signals of all spectra. Due to slight changes in measurement conditions, e.g. pH variations, signal positions may vary slightly across spectra. As a consequence, prior to further analysis signals belonging to the same compound have to be aligned across spectra. This is the purpose of the 'speaq' package.
#' @param feat Output of `gen_feat_mat()`.
#' @param maxShift Maximum number of points along the "ppm-axis" which a value can be moved by speaq package e.g. 50. 50 is a suitable starting value for plasma spectra with a digital resolution of 128K. Note that this parameter has to be individually optimized depending on the type of analyzed spectra and the digital resolution. For urine which is more prone to chemical shift variations this value most probably has to be increased.
#' @param spectrum_data Output of `generate_lorentz_curves()`.
#' @param si_size_real_spectrum Number of real data points in your original spectra.
#' @return A matrix containing the aligned integral values of all spectra. Each row contains the data of each spectrum and each column corresponds to one data point. Each entry corresponds to the integral of a deconvoluted signal with the signal center at this specific position after alignment by speaq.
#' @examples
#' spectrum_data <- glc_sim("bruker/sim")
#' feat <- gen_feat_mat(spectrum_data)
#' maxShift <- 100
#' si_size_real_spectrum <- length(spectrum_data[[1]]$y_values)
#' after_speaq_mat <- speaq_align(feat, maxShift, spectrum_data, si_size_real_spectrum)
speaq_align_v1 <- function(feat = gen_feat_mat(spectrum_data),
                           maxShift = 50,
                           spectrum_data = glc_sim(),
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

# Work in Progress #####

plot_alignment <- function(Y, U, Q, C, mfrow) {
    if (!is.null(mfrow)) {
        opar <- par(mfrow = mfrow, mar = c(0, 2, 0, 0), oma = c(4.1, 2.1, 0, 0))
        on.exit(par(opar), add = TRUE)
    }
    for (i in seq_len(s)) {
        plot(x = seq_len(ncol(Y)), y = Y[i, ],
            type = "l", lty = 1, col = "darkgrey",
            xlim = c(1, ncol(Y)), ylim = c(0, max(Y)),
            xaxt = if (i == s) "s" else "n",
            xlab = "Datapoint Number", ylab = "Signal Intensity"
        )
        abline(v = U[[i]], col = transp("darkgrey", 0.5), lty = 2)
        abline(v = Q[[i]], col = transp("blue", 0.5), lty = 2)
        lines(x = seq_len(ncol(Y)), y = C$Y[i, ], col = "blue", lty = 1)
    }
    mtext("Datapoint Number", side = 1, outer = TRUE, line = 3)
    mtext("Signal Intensity", side = 2, outer = TRUE, line = 1)
}

#' @export
#' @title Combine Peaks
#' @description Even after calling `speaq_align()`, the alignment of individual signals is not always perfect, as 'speaq' performs a segment-wise alignment i.e. groups of signals are aligned. For further improvements, partly filled neighboring columns are merged.
#' @param shifted_mat The matrix returned by `speaq_align()`.
#' @param range Amount of adjacent columns which are permitted to be used for improving the alignment.
#' @param lower_bound Minimum amount of non-zero elements per column to trigger the alignment improvement.
#' @param spectrum_data The list of deconvoluted spectra as returned by `generate_lorentz_curves()`.
#' @param data_path If not NULL, the returned dataframes `long` and `short` are written to `data_path` as "aligned_res_long.csv" and "aligned_res_short.csv".
#' @return A list containing two data frames `long` and `short`. The first one contains one for one column for each data point in the original spectrum. The second one contains only columns where at least one entry is non-zero. The returned data.frames are additionally written to disk as .csv files `aligned_res_long.csv` and `aligned_res_short.csv`.
#' @examples
#' shifted_mat <- speaq_align(spectrum_data = glc_sim("bruker/sim"))
#' spectrum_data <- glc_sim("bruker/sim")
#' system.time(M1 <- combine_peaks_v1(shifted_mat, spectrum_data = spectrum_data))
#' system.time(M2 <- combine_peaks(shifted_mat, spectrum_data = spectrum_data))
#' all.equal(M1, M2)
combine_peaks_v2 <- function(shifted_mat = speaq_align(spectrum_data = spectrum_data, show = TRUE),
                             lower_bound = 1,
                             range = 5,
                             spectrum_data = glc_sim("bruker/sim"),
                             data_path = NULL) {
    M <- replace(shifted_mat, is.na(shifted_mat), 0) # Shifted matrix with NAs replaced by 0
    s <- nrow(M) # Number of spectra
    nz <- apply(M, 2, function(x) sum(x != 0)) # Number of Non-zero Entries (NZs) per column
    nzIDs <- which(nz >= lower_bound) # Columns IDs (CIDs) with at least `lower_bound` NZs.
    nzIDs <- nzIDs[order(nz[nzIDs], decreasing = TRUE)] # Sort column IDs with most NZs first
    for (j in nzIDs) {
        if (nz[j] == 0) next # Can happen if column j has already been absorbed by another column
        while (TRUE) {
            nbIDs <- c((j - range):(j - 1), (j + 1):(j + range)) # Indices of Neighbour columns
            nbIDs <- nbIDs[nbIDs > 0 & nbIDs <= ncol(M)] # Remove out-of-bounds indices
            cvs <- nz[nbIDs] # Number of Combinable Values (CVs) per column
            if (all(cvs %in% c(0, s))) break # Stop early if neighbour cols are empty or completely filled already
            for (i in seq_along(cvs)[cvs > 0 & cvs < s]) {
                n <- nbIDs[i] # Index of neighbour column
                if (any(M[, n] != 0 & M[, j] != 0)) cvs[n] <- 0
                # Example:
                # c(0, 0, 3, 4, 0) -> M[, n]
                # c(2, 3, 0, 2, 1) -> M[, j]
                # c(F, F, F, T, F) -> M[, n] != 0 & M[, j] != 0
                # Not combinable because there is a common non-zero element at index 4
            }
            if (all(cvs %in% c(0, nrow(M)))) break
            x <- nbIDs[which.max(cvs)] # Index of column to combine with current column
            if (any(M[, x] == 0.03336328)) browser()
            M[, j] <- M[, j] + M[, x]
            M[, x] <- 0
            nz[x] <- 0
        }
    }
    ppm_rounded <- round(spectrum_data[[1]]$x_values_ppm, 4)
    colnames(M) <- as.character(ppm_rounded)
    nz <- apply(M, 2, function(x) sum(x != 0))
    nzIDs <- which(nz >= lower_bound)
    S <- M[, nzIDs] # Shifted matrix without zero-only-columns
    # utils::write.csv2(M2, file = file.path(data_path, "aligned_res_long.csv"))
    # utils::write.csv2(M, file = file.path(data_path, "aligned_res_short.csv"))
    list("short" = S, "long" = M)
}
