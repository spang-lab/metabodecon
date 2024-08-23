# Exported Main #####

#' @export
#' @title Get ppm range
#' @author Wolfram Gronwald, Tobias Schmidt
#' @description Returns the ppm range across all peaks of the provided deconvoluted spectra.
#' @param spectrum_data A list of deconvoluted spectra as returned by [generate_lorentz_curves()].
#' @param show Whether to plot the ppm range on the spectrum plot.
#' @param mar The margins of the plot. Passed on to `par()`.
#' @return A vector containing the lowest and highest ppm value over all peaks of the provided deconvoluted spectra.
#' @examples
#' sim <- metabodecon_file("bruker/sim_subset")
#' decons <- generate_lorentz_curves(
#'     sim,
#'     sfr = c(3.58, 3.42), wshw = 0,
#'     ask = FALSE, verbose = FALSE
#' )
#' ppm_rng <- get_ppm_range(decons, show = TRUE)
#' print(ppm_rng)
get_ppm_range <- function(spectrum_data,
                          show = FALSE,
                          mar = c(4.1, 4.1, 1.1, 0.1)) {
    msg <- "spectrum_data must be a list of deconvoluted spectra."
    ss <- spectrum_data
    ss <- if (is_decon_list(ss)) ss else if (is_decon_obj(ss)) list(ss) else stop(msg)
    ppm_min <- min(sapply(ss, function(s) min(s$peak_triplets_middle)))
    ppm_max <- max(sapply(ss, function(s) max(s$peak_triplets_middle)))
    ppm_rng <- c(ppm_min, ppm_max)
    if (show) {
        ymax <- max(sapply(ss, function(s) max(s$y_values)))
        xrng <- range(c(sapply(ss, function(s) s$x_values_ppm)))
        y0.8 <- ymax * 0.8
        a <- ppm_rng[1]
        b <- ppm_rng[2]
        cols <- rainbow(length(ss))
        alen <- (b - a) / 4
        lgdtxt <- paste("Spectrum", 1:length(ss))
        op <- par(mar = mar)
        on.exit(par(op))
        plot(NA, type = "n", xlim = xrng[2:1], ylim = c(0, ymax), ylab = "Signal Intensity [au]", xlab = "Chemical Shift [ppm]")
        abline(v = ppm_rng, lty = 2)
        for (i in 1:length(ss)) lines(ss[[i]]$x_values_ppm, ss[[i]]$y_values, col = cols[i])
        arrows(
            x0 = c(a + alen, b - alen),
            x1 = c(a, b),
            y0 = y0.8,
            y1 = y0.8,
            length = 0.1,
            lty = 2,
            col = "black"
        )
        text(mean(c(a, b)), y0.8, "ppm range")
        mtext(round(c(a, b), 4), side = 3, line = 0, at = c(a, b))
        legend("topright", legend = lgdtxt, col = cols, lty = 1)
    }
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
    P <- lapply(seq_along(D$A), function(i) ncol(X) - D$w[[i]] * scale_factor_x)
    names(P) <- names(D$A)
    # We get the indices in SDP, which decreases from right to left, but we need
    # them as normal indices (increasing from right to left). For example if we
    # have n = 8 datapoints, and the third and the seventh point from left are
    # peaks, we want the numbers 3 and 7, but we get 0.005 and 0.001. To fix
    # this, we first revert the scaling and then subtract from n = 8, i.e.
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
#' @return A matrix containing the aligned integral values of all spectra. Each row contains the data of each spectrum and each column corresponds to one data point. Each entry corresponds to the integral of a deconvoluted signal with the signal center at this specific position after alignment by speaq.
#' @examples
#' \dontrun{
#' after_speaq_mat <- speaq_align(feat, maxShift)
#' }
speaq_align <- function(feat = gen_feat_mat(spectrum_data),
                        maxShift = 50,
                        spectrum_data,
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
    # $ Y: num [1:16, 1:1309] 0.01356 0.02418 ... (matrix of aligned spectra)
    # $ new_peakList: List of 16                  (aligned peaks of each spectrum)
    # ..$ sim_01: num [1:13] 350 416 457 583 683 ...
    # ..$ sim_02: num [1:11] 416 458 542 583 683 ...
    # ...

    # CONTINUE HERE
    # TODO: plot `feat$data_matrix` and `data_result_aligned$Y, i.e. spectra before and after alignment

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
                integral_results[[spec]][ent] <- feat$A[[spec]][ent] * (atan((-peakList_result[[spec]][ent] + 131072) / feat$lambda[[spec]][ent]) - atan((-peakList_result[[spec]][ent]) / feat$lambda[[spec]][ent]))
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

#' @export
#' @title Combine Peaks
#' @description Even after alignment of spectra, alignment of individual signals is not always perfect, as a segment-wise alignment is performed i.e. groups of signals are aligned. For further improvements, partly filled neighboring columns are merged.
#' @param shifted_mat (matrix) the matrix obtained after alignment by speaq.
#' @param range (positive integer) amount of adjacent columns which are permitted to be used for improving the alignment e.g. 5.
#' @param lower_bound (positive integer) amount of columns that need to be skipped (f.e. because they contain rownames instead of values, only modify in case of errors) default=1.
#' @param spectrum_data (data frame) The output generated by the function generate_lorentz_curves
#' @param data_path (string) Path to the parent folder where the original spectra are stored. After deconvolution this folder also contains for each spectrum two .txt files which contain for each spectrum the spectrum approximated from all deconvoluted signals and a parameter file that contains all numerical values of the deconvolution
#' @return A list containing two data frames `long` and `short`. The first one contains one for one column for each data point in the original spectrum. The second one contains only columns where at least one entry is non-zero. The returned data.frames are additionally written to disk as .csv files `aligned_res_long.csv` and `aligned_res_short.csv`.
#' @examples
#' \dontrun{
#' aligned_res <- combine_peaks(after_speaq_mat, range, lower_bound, spectrum_data, data_path)
#' }
combine_peaks <- function(shifted_mat,
                          range,
                          lower_bound,
                          spectrum_data,
                          data_path) {
    shifted_mat_no_na <- replace(is.na(shifted_mat), 0)
    col_entries <- numeric(ncol(shifted_mat_no_na))
    for (i in 1:ncol(shifted_mat_no_na)) {
        col_entries[i] <- length(which(shifted_mat_no_na[, i] != 0))
    }
    bla <- shifted_mat_no_na != 0
    shifted_mat_no_na <- shifted_mat_no_na
    for (i in (nrow(shifted_mat_no_na) - 1):lower_bound) {
        print(i)
        indices <- which(col_entries == i)
        for (j in indices) {
            if (col_entries[j] != 0) {
                while (T) {
                    bli <- rep(0, 2 * range)
                    l <- 0
                    search_range <- c((j - range):(j - 1), (j + 1):(j + range))
                    for (k in search_range) {
                        l <- l + 1
                        if (all(which(bla[, k]) %in% which(!bla[, j]))) {
                            bli[l] <- col_entries[k]
                        }
                    }
                    if (all(bli == 0)) {
                        break
                    } else {
                        index_max <- search_range[which.max(bli)]

                        shifted_mat_no_na[, j] <- shifted_mat_no_na[, j] + shifted_mat_no_na[, index_max]
                        bla[, j] <- bla[, j] | bla[, index_max]
                        col_entries[j] <- col_entries[j] + col_entries[index_max]
                        shifted_mat_no_na[, index_max] <- 0
                        bla[, index_max] <- F
                        col_entries[index_max] <- 0
                    }
                }
            }
        }
    }
    colnames(shifted_mat_no_na) <- as.character(round(spectrum_data[[1]]$x_values_ppm, digits = 4))
    # removal of columns containing only zeros
    kick.out <- apply(shifted_mat_no_na, 2, function(z) {
        all(z == 0)
    })
    shifted_mat_no_na_s <- shifted_mat_no_na[, !kick.out]
    utils::write.csv2(shifted_mat_no_na, file = file.path(data_path, "aligned_res_long.csv"))
    utils::write.csv2(shifted_mat_no_na_s, file = file.path(data_path, "aligned_res_short.csv"))
    return_list <- list("short" = shifted_mat_no_na_s, "long" = shifted_mat_no_na)
    return(return_list)
}

#' @export
#' @title dohCluster
#' @description Helper of [speaq_align()].
#' @param X `data_matrix` as returned by [gen_feat_mat()]
#' @param peakList `peakList` as returned by function [gen_feat_mat()]
#' @param refInd Number of the reference spectrum i.e. the spectrum to which all signals will be aligned to.
#' @param maxShift Maximum number of points along the "ppm-axis" which a value can be moved (e.g. 50).
#' @param acceptLostPeak (logic) default is TRUE
#' @param verbose (logic) default is TRUE
#' @return A list containing two data frames `Y` and `new_peakList`. The first one contains the aligned spectra, the second one contains the aligned signals of each spectrum.
#' @details The function dohCluster of the 'speaq' package has been rewritten to be compatible with the data generated by 'metabodecon' and the function 'gen_feat_mat' and to return a new peakList of aligned spectra.
#' @examples
#' \dontrun{
#' res <- dohCluster(X, peakList, refInd = 0, maxShift = 100, acceptLostPeak = TRUE, verbose = TRUE)
#' }
dohCluster <- function(X, peakList, refInd = 0, maxShift = 100, acceptLostPeak = TRUE, verbose = TRUE) {
    if (is.null(maxShift)) {
        if (verbose) {
            cat("\n --------------------------------")
            cat("\n maxShift=NULL, thus dohCluster will automatically detect the optimal value of maxShift.")
            cat("\n --------------------------------\n")
        }
        res <- dohCluster_withoutMaxShift(X = X, peakList = peakList, refInd = refInd, acceptLostPeak = acceptLostPeak, verbose = verbose)
    } else {
        if (verbose) {
            cat("\n --------------------------------")
            cat("\n dohCluster will run with maxShift=", maxShift)
            cat("\n If you want dohCluster to detect the optimal maxShift automatically,")
            cat("\n use dohCluster(..., maxShift = NULL, ...)")
            cat("\n --------------------------------\n")
        }
        res <- dohCluster_withMaxShift(X = X, peakList = peakList, refInd = refInd, maxShift = maxShift, acceptLostPeak = acceptLostPeak, verbose = verbose)
    }
    return(res)
}

# Private Helpers #####

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
    tmp <- mapply(par_nam, spc_nam, FUN = function(p, s) if (p != s) warning(sprintf("Mismatched file names: %s and %s", p, s)) )
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
