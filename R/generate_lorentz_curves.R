# Develop API Functions #####

#' @title Generate Lorentz Curves from NMR Spectra
#' @description Deconvolutes NMR spetra and generates a Lorentz curve for each detected signal within a spectra.
#' @param data_path (string) Path to the folder where the original spectra are stored. After deconvolution this folder contains two additional .txt files for each spectrum which contain the spectrum approximated from all deconvoluted signals and a parameter file that contains all numerical values of the deconvolution.
#' @param file_format (string) Format of the spectra files.
#' @param make_rds (bool) Store results as a rds file on disk? Should be set to TRUE if many spectra are evaluated to decrease computation time.
#' @param number_iterations (int) Number of iterations for the approximation of the parameters for the Lorentz curves.
#' @param range_water_signal_ppm (float) Half width of the water artefact in ppm.
#' @param signal_free_region (float) Row vector with two entries consisting of the ppm positions for the left and right border of the signal free region of the spectrum.
#' @param smoothing_param (int) Row vector with two entries consisting of the number of smoothing repeats for the whole spectrum and the number of data points (uneven) for the mean calculation.
#' @param delta (float) Threshold value to distinguish between signal and noise.
#' @param scale_factor (int) Row vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis.
#' @param ask (bool) Whether to ask for user input during the deconvolution process. If set to FALSE, the provided default values will be used.
#' @details First, an automated curvature based signal selection is performed. Each signal is represented by 3 data points to allow the determination of initial Lorentz curves. These Lorentz curves are then iteratively adjusted to optimally approximate the measured spectrum. For each spectrum two text files will be created in the parent folder i.e. the folder given in data path. The spectrum approximated from all deconvoluted signals and a parameter file that contains all numerical values of the deconvolution. Furthermore, the numerical values of the deconvolution are also stored in a data_frame.
#' @details Shall replace `generate_lorentz_curves` as soon as implementation is finished.
#' @noRd
generate_lorentz_curves_v1 <- function(data_path,
                                       file_format = c("bruker", "jcampdx"),
                                       make_rds = FALSE,
                                       number_iterations = 10,
                                       range_water_signal_ppm = 0.1527692,
                                       signal_free_region = c(11.44494, -1.8828),
                                       smoothing_param = c(2, 5),
                                       delta = 6.4,
                                       scale_factor = c(1000, 1000000),
                                       ask = TRUE) {
    # Check arguments
    file_format <- match.arg(file_format)

    # Switch to data directory
    data_path <- normalizePath(data_path)
    owd <- getwd()
    setwd(data_path)
    on.exit(setwd(owd))

    # Get input files
    if (file_format == "jcampdx") {
        files <- dir(data_path, pattern = "\\.dx$") # `.dx` files inside `data_path`
        spectroscopy_value <- NULL
        processing_value <- NULL
    } else if (file_format == "bruker") {
        files <- list.dirs(data_path, recursive = FALSE, full.names = FALSE) # folders inside `data_path`
        spectroscopy_value <- readline(prompt = "What is the name of the subfolder of your filepath: (e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10) ")
        processing_value <- readline(prompt = "What is the name of the subsubsubfolder of your filepath: (e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10/pdata/10) ")
    }

    # Reorder files in case user wants to use one spectrum to determine parameters for all spectra
    same_parameter <- get_yn_input("Do you want to use the same parameters (signal_free_region, range_water_signal_ppm) for all spectra?")
    if (same_parameter) {
        print(files)
        number <- get_num_input("Choose number of file which is used to adjust all parameters: [e.g. 1] ", min = 1, max = length(files), int = TRUE)
        message(paste("The selected file to adjust all parameters for all spectra is: ", files[number]))
        files <- c(files[number], files[-number])
    }

    # Do actual deconvolution
    spectrum_data <- list()
    for (i in seq_along(files)) {
        name <- files[i] # bruker: `urine_2`, jcampdx: `urine_2.dx`
        filepath <- switch(file_format, # see [FAQ](../vignettes/FAQ.Rmd#file-structure) for example file structures
            "bruker" = paste(data_path, name, spectroscopy_value, sep = "/"),
            "jcampdx" = data_path,
            stop("Invalid file format")
        )
        x <- deconvolute_spectrum(filepath, name, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, current_filenumber = i, number_of_files = length(files))
        spectrum_data[[name]] <- x
        # Save `range_water_signal` and `signal_free_region` for next loop passage as those might have been adjusted interactively by the user
        range_water_signal_ppm <- x$range_water_signal_ppm
        signal_free_region <- x$signal_free_region
    }

    # Save results
    if (make_rds) {
        saveRDS(object = spectrum_data, file = file.path(data_path, "spectrum_data.rds"))
    }
    return(spectrum_data)
}

# Exported API Functions #####

#' @export
#' @title generate_lorentz_curves
#' @description This function starts the actual deconvolution of your spectra
#' and generates for each detected signal a Lorentz curve
#' @param data_path (string) path to the parent folder where the original
#' spectra are stored. After deconvolution this folder also contains for each
#' spectrum two .txt files which contain for each spectrum the spectrum
#' approximated from all deconvoluted signals and a parameter file that
#' contains all numerical values of the deconvolution
#' @param file_format (string) default is bruker, the other supported format
#' is jcamp-dx
#' @param make_rds (logic) if you would like to store your results as a rds
#' file, default is set to FALSE, should be set to true to save your results if
#' many spectra are evaluated and therefore computing time increases
#' @examples \dontrun{
#' data_path <- c("C:/example_data")
#' spectrum_data <- generate_lorentz_curves(
#'     data_path = data_path,
#'     file_format = "bruker",
#'     make_rds = FALSE
#' )
#' }
#' @details First, an automated curvature based signal selection is performed.
#' Each signal is represented by 3 data points to allow the determination of
#' initial Lorentz curves. These Lorentz curves are then iteratively adjusted
#' to optimally approximate the measured spectrum. For each spectrum two text
#' files will be created in the parent folder i.e. the folder given in data
#' path. The spectrum approximated from all deconvoluted signals and a
#' parameter file that contains all numerical values of the deconvolution.
#' Furthermore, the numerical values of the deconvolution are also stored in a
#' data_frame.
generate_lorentz_curves <- function(data_path,
                                    file_format = "bruker",
                                    make_rds = FALSE) {
    spectrum_data <- MetaboDecon1D(data_path, file_format = file_format)
    if (make_rds) {
        saveRDS(object = spectrum_data, file = file.path(data_path, "spectrum_data.rds"))
    }
    return(spectrum_data)
}

# Exported Helpers #####

#' @export
#' @title get_ppm_range
#' @description function to get ppm range of your spectra, this is required for
#' spectra alignment following deconvolution
#' @param spectrum_data Data frame generated from the function
#' generate_lorentz_curves, contains results from the deconvolution of your
#' spectra
#' @examples \dontrun{
#' ppm_range <- get_ppm_range(spectrum_data = spectrum_data)
#' }
#' @details determines after deconvolution the signal with the highest and the
#' signal with the lowest ppm value.
get_ppm_range <- function(spectrum_data) {
    ppm_scale_min <- min(spectrum_data[[1]]$peak_triplets_middle)
    ppm_scale_max <- max(spectrum_data[[1]]$peak_triplets_middle)

    for (entry in spectrum_data) {
        if (ppm_scale_min > min(entry$peak_triplets_middle)) {
            ppm_scale_min <- min(entry$peak_triplets_middle)
        }
        if (ppm_scale_max < max(entry$peak_triplets_middle)) {
            ppm_scale_max <- max(entry$peak_triplets_middle)
        }
    }
    return(c(ppm_scale_min, ppm_scale_max))
}

#' @export
#' @title gen_feat_mat
#' @description function to generate feature matrix
#' @param data_path (string): path to the parent folder of where the original
#' spectra are stored. After deconvolution this folder also contains for each
#' spectrum two .txt files which contain for each spectrum the spectrum
#' approximated from all deconvoluted signals and a parameter file that
#' contains all numerical values of the deconvolution
#' @param ppm_range (numeric) this is the result from the function
#' get_ppm_range i.e the ppm range over which your signals are distributed
#' @param si_size_real_spectrum (positive integer) number of real data points
#' in your original spectra, e.g. 128k = 131072 data points
#' @param scale_factor_x  (positive integer): A factor which is used to avoid
#' rounding errors e.g. 1000
#' @examples \dontrun{
#' feat <- gen_feat_mat(
#'     data_path = data_path,
#'     ppm_range = ppm_range,
#'     si_size_real_spectrum = si_size_real_spectrum,
#'     scale_factor_x = scale_factor
#' )
#' }
#' @details The output of this function is a data frame containing a matrix of
#' all integral values found in your spectra (data_matrix). Here, each row
#' corresponds to one spectrum and each column to one data point of the
#' spectra, for example 128k data points in each spectrum correspond to 128K
#' columns in the data matrix. Furthermore, a list of all signals(peakList),
#' and individual parameters of the signals are stored in the data structure.
gen_feat_mat <- function(data_path,
                         ppm_range,
                         si_size_real_spectrum,
                         scale_factor_x) {
    # Set highest/lowest ppm value and range
    message("Loading ppm range")
    ppm_highest_value <- max(ppm_range)
    ppm_lowest_value <- min(ppm_range)
    ppm_range <- ppm_highest_value - ppm_lowest_value

    # Import data sets
    message("Importing data")
    your_path_spectra <- data_path
    # Be careful: `data_path` should only contain txt files of MetaboDecon1D
    files <- list.files(your_path_spectra, ".txt")
    num_spectra <- length(files) / 2

    # Import spectra datasets
    message("Importing spectra")
    data_spectrum <- vector(mode = "list", length = num_spectra)
    data_info <- vector(mode = "list", length = num_spectra)
    numerator_info <- 1
    numerator_spectrum <- 1
    for (i in 1:length(files)) {
        if (grepl(" parameters", files[i])) {
            data_info[[numerator_info]] <- as.matrix(data.table::fread(files[i], header = FALSE))
            numerator_info <- numerator_info + 1
        }

        if (grepl(" approximated_spectrum", files[i])) {
            data_spectrum[[numerator_spectrum]] <- as.vector(unlist(data.table::fread(files[i], header = FALSE)))[-1]
            numerator_spectrum <- numerator_spectrum + 1
        }
    }
    message("Converting parameters to list")

    # Save parameters in list
    w <- vector(mode = "list", length = num_spectra)
    lambda <- vector(mode = "list", length = num_spectra)
    A <- vector(mode = "list", length = num_spectra)
    noise_threshold <- vector(mode = "list", length = num_spectra)
    spectrum_approx <- vector(mode = "list", length = num_spectra)
    vec_num_Lorentz_curves <- numeric(num_spectra)

    for (i in 1:num_spectra) {
        w[[i]] <- as.numeric(data_info[[i]][1, ][-1])
        lambda[[i]] <- as.numeric(data_info[[i]][2, ][-1])
        A[[i]] <- as.numeric(data_info[[i]][3, ][-1])
        noise_threshold[[i]] <- as.numeric(data_info[[i]][4, ][2])
        spectrum_approx[[i]] <- data_spectrum[[i]]
        vec_num_Lorentz_curves[i] <- length(w[[i]])
    }

    message("Checking for irregular values")
    # Remove zeros
    for (i in 1:num_spectra) {
        w_values <- w[[i]]
        A_values <- A[[i]]
        lambda_values <- lambda[[i]]

        # Find position of 0's
        zero_positions <- which(w[[i]] == 0)

        # If zero_positions is filled
        if (length(zero_positions) > 0) {
            w[[i]] <- w_values[-zero_positions]
            A[[i]] <- A_values[-zero_positions]
            lambda[[i]] <- lambda_values[-zero_positions]
        } else {
            w[[i]] <- w_values
            A[[i]] <- A_values
            lambda[[i]] <- lambda_values
        }
    }

    # Generate x values
    spectrum_x <- seq(length(spectrum_approx[[1]]), 1, -1)
    # Generate x ppm values
    spectrum_x_ppm <- seq(ppm_highest_value, ppm_lowest_value, -ppm_range / (si_size_real_spectrum - 1))

    # Generate data matrix
    data_matrix <- matrix(nrow = 0, ncol = length(data_spectrum[[1]]))
    for (i in 1:num_spectra) {
        data_matrix <- rbind(data_matrix, data_spectrum[[i]])
    }

    # Generate peakList
    # Here the peaks are substracted of real size of spectra to compensate
    # problems of shifting of speaq package. This is changed later again
    peakList <- vector(mode = "list", length = num_spectra)
    for (i in 1:num_spectra) {
        peakList[[i]] <- (dim(data_matrix)[2]) - (w[[i]] * scale_factor_x)
    }
    return_list <- list("data_matrix" = data_matrix, "peakList" = peakList, "w" = w, "A" = A, "lambda" = lambda)
    return(return_list)
}

#' @export
#' @title dohCluster
#' @description function for the speaq package. The function dohCluster will be
#' automatically called by the function speaq_align. The speaq package is used
#' to perform the signal alignment across the individual spectra. For speaq
#' please cite: (Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S, Pieters L,
#' et al. (2018) speaq 2.0: A complete workflow for high-throughput 1D
#' NMRspectra processing and quantification. PLoS Comput Biol 14(3): e1006018.
#' https://www.doi.org/10.1371/journal.pcbi.1006018. The spectra deconvolution
#' process yields the signals of all spectra. Due to slight changes in
#' measurement conditions, e.g. pH variations, signal positions may vary
#' slightly across spectra. As a consequence, prior to further analysis signals
#' belonging to the same compound have to be aligned across spectra. This is
#' the purpose of the speaq package.
#' @param X `data_matrix` as returned by [gen_feat_mat()]
#' @param peakList `peakList` as returned by function [gen_feat_mat()]
#' @param refInd (positive integer) the number of the reference spectrum i.e.
#' the spectrum to which all signals will be aligned to. This number will be
#' automatically determined by the function speaq_align which then calls
#' dohCluster
#' @param maxShift (positive integer) maximum number of points along the
#' "ppm-axis" which a value can be moved by speaq package e.g. 50
#' @param acceptLostPeak (logic) default is TRUE
#' @param verbose (logic) default is TRUE
#' @details The function dohCluster of the speaq package has been rewritten to
#' be compatible with the data generated by MetaboDecon1D and the function
#' gen_feat_mat and to return a new peakList of aligned spectra.
#' Overwrite original dohCluster function of speaq package. Function is able to
#' return the new peakList of aligned spectra.
dohCluster <- function(X, peakList, refInd = 0, maxShift = 100, acceptLostPeak = TRUE, verbose = TRUE) {
    .withMaxShift <- function(X, peakList, refInd = 0, maxShift = 100,
                              acceptLostPeak = TRUE, verbose = TRUE) {
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
    .withoutMaxShift <- function(X, peakList, refInd = 0, acceptLostPeak = TRUE,
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
            cat(
                "\n Alignment time: ", (endTime[3] - startTime[3]) / 60,
                " minutes"
            )
        }
        # Added by me at 31.08.21
        return_list <- list("Y" = bestY, "new_peakList" = peakListNew)
        return(return_list)
    }
    if (is.null(maxShift)) {
        if (verbose) {
            cat("\n --------------------------------")
            cat("\n maxShift=NULL, thus CluPA will automatically detect the optimal value of maxShift.")
            cat("\n --------------------------------\n")
        }
        res <- .withoutMaxShift(
            X = X, peakList = peakList, refInd = refInd,
            acceptLostPeak = acceptLostPeak, verbose = verbose
        )
    } else {
        if (verbose) {
            cat("\n --------------------------------")
            cat("\n CluPA will run with maxShift=", maxShift)
            cat("\n If you want CluPA automatically detect the optimal maxShift,")
            cat("\n let try with dohCluster(...,maxShift=NULL,..)")
            cat("\n --------------------------------\n")
        }
        res <- .withMaxShift(
            X = X, peakList = peakList, refInd = refInd,
            maxShift = maxShift, acceptLostPeak = acceptLostPeak,
            verbose = verbose
        )
    }
    return(res)
}

#' @export
#' @title speaq_align
#' @description Alignment by speaq. The speaq package is used to perform the
#' signal alignment across the individual spectra. For speaq please cite:
#' (Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S, Pieters L, et al. (2018)
#' speaq 2.0: A complete workflow for high-throughput 1D NMRspectra processing
#' and quantification. PLoS Comput Biol 14(3): e1006018.
#' https://www.doi.org/10.1371/journal.pcbi.1006018. The spectra deconvolution
#' process yields the signals of all spectra. Due to slight changes in
#' measurement conditions, e.g. pH variations, signal positions may vary
#' slightly across spectra. As a consequence, prior to further analysis signals
#' belonging to the same compound have to be aligned across spectra. This is
#' the purpose of the speaq package.
#' @param feat (data frame) Contains data of the deconvoluted signals of all
#' analyzed spectra before alignment.
#' @param maxShift (positive integer) maximum number of points along the
#' "ppm-axis" which a value can be moved by speaq package e.g. 50. 50 is a
#' suitable starting value for plasma spectra with a digital resolution of
#' 128K. Note that this parameter has to be individually optimized depending on
#' the type of analyzed spectra and the digital resolution. For urine which is
#' more prone to chemical shift variations this value most probably has to be
#' increased.
#' @param spectrum_data (data frame) The output generated by the function
#' generate_lorentz_curves.
#' @param si_size_real_spectrum (positive integer) number of real data points
#' in your original spectra, e.g. 128k = 131072 data points
#' @examples \dontrun{
#' after_speaq_mat <- speaq_align(feat, maxShift)
#' }
#' @details The output of speaq_align is a data matrix of aligned integral
#' values. Each row contains the data of each spectrum and each column
#' corresponds to one data point. Each entry corresponds to the integral of a
#' deconvoluted signal with the signal center at this specific position after
#' alignment by speaq.
speaq_align <- function(feat, maxShift, spectrum_data, si_size_real_spectrum) {
    # First identify reference spectrum
    resFindRef <- speaq::findRef(feat$peakList)
    # Save index of reference spectrum i.e the spectrum which all others will be aligned to
    refInd <- resFindRef$refInd
    # Use overwritten CluPA function for multiple spectra
    data_result_aligned <- dohCluster(feat$data_matrix, peakList = feat$peakList, refInd = refInd, maxShift = maxShift, acceptLostPeak = FALSE, verbose = TRUE)
    # data_result_aligned$Y contains matrix of aligned spectra each row corresponds to one spectrum
    # data_result_aligned$new_peakList contains the aligned signals of each spectrum

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
#' @title combine_peaks
#' @description Even after alignment of spectra, alignment of individual
#' signals is not always perfect, as a segment-wise alignment is performed i.e.
#' groups of signals are aligned. For further improvements, partly filled
#' neighboring columns are merged.
#' @param shifted_mat (matrix) the matrix obtained after alignment by speaq.
#' @param range (positive integer) amount of adjacent columns which are
#' permitted to be used for improving the alignment e.g. 5
#' @param lower_bound (positive integer) amount of columns that need to be
#' skipped (f.e. because they contain rownames instead of values, only modify
#' in case of errors) default=1.
#' @param spectrum_data (data frame) The output generated by the function
#' generate_lorentz_curves
#' @param data_path (string) Path to the parent folder where the original
#' spectra are stored. After deconvolution this folder also contains for each
#' spectrum two .txt files which contain for each spectrum the spectrum
#' approximated from all deconvoluted signals and a parameter file that
#' contains all numerical values of the deconvolution
#' @examples \dontrun{
#' aligned_res <- combine_peaks(
#'     after_speaq_mat,
#'     range,
#'     lower_bound,
#'     spectrum_data,
#'     data_path
#' )
#' }
#' @details As result two .csv files will be generated. One that contain all
#' columns for each data point in the original spectrum one column
#' (aligned_res_long.csv) and one where all columns that contain only zeros
#' (aligned_res_short.csv) are omitted. Both files will be stored in the
#' directory specified by data_path.
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

# Private Helpers: deconvolute_spectrum #####

#' @title Deconvolute one single spectrum
#' @description Deconvolute one single spectrum
#' @param path Path to file or folder containing the spectra files.
#' @param type Format of the spectra files. Either `"bruker"` or `"jcampdx"`.
#' @param procno Processing value for the file. E.g. `"10"`. Called `procno` in the Bruker TopSpin Manual.
#' @param expno Spectroscopy value for the file. E.g. `"10"`. Called `expno` in the Bruker TopSpin Manual.
#' @param nfit Number of iterations for the approximation of the parameters for the Lorentz curves.
#' @param wshw Half width of the water artefact in ppm.
#' @param sfr Row vector with two entries consisting of the ppm positions for the left and right border of the signal free region of the spectrum.
#' @param smopts Row vector with two entries consisting of the number of smoothing repeats for the whole spectrum and the number of data points (uneven) for the mean calculation.
#' @param delta Threshold value to distinguish between signal and noise.
#' @param sf Row vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis.
#' @param ask Whether the function should ask the user to confirm the signal free region and the water signal. Must be TRUE if either sfr or wshw are not given.
#' @param filno Current file number. Only used for progress prints.
#' @param nfils Total number of files. Only used for progress prints.
#' @param bwc Use the old, slightly incorrect method for calculating the signal free region and water signal to maintain backwards compatibility with MetaboDecon1D results? For details see `Check: ...` issues in `TODOS.md`.
#' @return A list containing the deconvoluted spectrum data.
#' @examples \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1")
#' type <- "bruker"
#' deconvolute_spectrum_v2(path, type)
#' }
#' @noRd
deconvolute_spectrum_v2 <- function(path = file.path(download_example_datasets(), "bruker/urine/urine_1"),
                                    type = "bruker",
                                    expno = 10,
                                    procno = 10,
                                    nfit = 10,
                                    wshw = 0.1527692,
                                    sfr = c(11.44494, -1.8828),
                                    smopts = c(2, 5),
                                    delta = 6.4,
                                    sf = c(1e3, 1e6),
                                    filno = 1,
                                    nfils = 1,
                                    ask = TRUE,
                                    bwc = TRUE) {
    type <- match.arg(type, c("bruker", "jcampdx"))
    spec <- read_spectrum(path, type, sf, expno, procno)
    # vcomp(n1 <- debugenv$spectrum_length, n <- n2 <- spec$n)
    # vcomp(v1 <- debugenv$spectrum_x, v2 <- spec$sdp)
    # vcomp(v1 <- debugenv$spectrum_x_ppm, v2 <- spec$ppm)
    # vcomp(v1 <- debugenv$spectrum_y_raw, v2 <- spec$Y$raw)
    # vcomp(v1 <- debugenv$spectrum_y_scaled, v2 <- spec$Y$scaled)

    spec <- determine_signal_free_region(spec, sfr, ask)
    spec <- determine_water_signal(spec, hwidth_ppm = wshw, bwc, ask)
    # vcomp(v1 <- debugenv$signal_free_region_left, v2 <- spec$sfr$left_sdp)
    # vcomp(v1 <- debugenv$signal_free_region_right, v2 <- spec$sfr$right_sdp)
    # vcomp(v1 <- debugenv$water_signal_left, v2 <- spec$ws$left_dp)
    # vcomp(v1 <- debugenv$water_signal_right, v2 <- spec$ws$right_dp)

    spec <- remove_water_signal(spec, bwc)
    spec <- remove_negative_signals(spec)
    # vcomp(v1 <- debugenv$spectrum_y_no_ws, v2 <- spec$Y$nows)
    # vcomp(v1 <- debugenv$spectrum_y_no_ws_no_neg, v2 <- spec$Y$pos)

    spec <- smooth_signals(spec, reps = smopts[1], k = smopts[2], bwc)
    # vcomp(v1 <- debugenv$spectrum_y_no_ws_no_neg_smoothed, v2 <- spec$Y$smooth)

    spec <- select_peaks_v2(spec, bwc)
    vcomp(n1 <- debugenv$spectrum_length, n <- n2 <- spec$n) # nolint
    vcomp(v1 <- debugenv$second_derivative[1, ], v2 <- spec$sdp[2:(n - 1)]) # nolint
    vcomp(v1 <- debugenv$second_derivative[2, ], v2 <- spec$d[2:(n - 1)]) # nolint
    vcomp(v1 <- debugenv$peaks_index, v2 <- as.integer(spec$peaks - 1)) # nolint
    vcomp(v1 <- debugenv$peaks_x, v2 <- spec$sdp[spec$peaks]) # nolint

    spec <- find_peak_borders(spec)
    vcomp(v1 <- debugenv$left_position[1, ], v2 <- spec$left - 1)
    vcomp(v1 <- debugenv$right_position[1, ], v2 <- spec$right - 1)

    spec <- find_left_positions_v0(spec)
    spec <- find_right_positions_v0(spec)
    vcomp(v1 <- debugenv$left_position[1, ], v2 <- spec$left - 1)
    vcomp(v1 <- debugenv$right_position[1, ], v2 <- spec$right - 1)

    spec <- find_right_positions_v1(spec)
    vcomp(v1 <- debugenv$right_position[1, ], v2 <- spec$sdp[spec$peaks + 1])

    plot_peaks(spec)

    spec <- find_left_positions_v1(spec)
    spec <- find_right_positions_v1(spec)


    spec <- get_peak_triplets(spec)
    spec <- get_preak_triplet_scores(spec)
    spec <- filter_peak_triplets(spec)
    # TODO

    spec <- calculate_initial_lorentz_curves(spec)
    spec <- refine_lorentz_curves(spec, nfit)
    spec <- calculate_lorentz_curve_integrals(spec)

    return_list <- list(
        filename = name,
        spectrum_x = spec$sdp,
        spectrum_x_ppm = spec$x_ppm,
        spectrum_y = spec$y,
        lorentz_curves = lorentz_curves_initial,
        mse_normed = mse_normed,
        spectrum_approx = spectrum_approx,
        index_peak_triplets_middle = index_peak_triplets_middle,
        index_peak_triplets_left = index_peak_triplets_left,
        index_peak_triplets_right = index_peak_triplets_right,
        peak_triplets_middle = peak_triplets_middle,
        peak_triplets_left = peak_triplets_left,
        peak_triplets_right = peak_triplets_right,
        integrals = integrals,
        sfr = c(sfrl_sdp, sfrr_sdp),
        wshwidth_ppm = ws$hwidth_ppm,
        A = A_new,
        lambda = lambda_new,
        w = w_new
    )
    return(return_list)
}

deconvolute_spectrum <- function(filepath,
                                 name,
                                 file_format,
                                 same_parameter,
                                 processing_value,
                                 number_iterations,
                                 range_water_signal_ppm,
                                 signal_free_region,
                                 smoothing_param,
                                 delta,
                                 scale_factor,
                                 current_filenumber,
                                 number_of_files) {
    msgf("Start deconvolution of %s:", name)

    x <- deconvolution(
        filepath,
        name,
        file_format,
        same_parameter,
        processing_value,
        number_iterations,
        range_water_signal_ppm,
        signal_free_region,
        smoothing_param,
        delta,
        scale_factor,
        current_filenumber
    )
    y <- list(
        "number_of_files" = number_of_files, # [1] add entry
        "filename" = x$filename,
        "x_values" = x$spectrum_x, # [3] rename
        "x_values_ppm" = x$spectrum_x_ppm, # [4] rename
        "y_values" = x$spectrum_y, # [5] rename
        "spectrum_superposition" = x$spectrum_approx, # [6] rename
        "mse_normed" = x$mse_normed,
        "index_peak_triplets_middle" = x$index_peak_triplets_middle,
        "index_peak_triplets_left" = x$index_peak_triplets_left,
        "index_peak_triplets_right" = x$index_peak_triplets_right,
        "peak_triplets_middle" = x$peak_triplets_middle,
        "peak_triplets_left" = x$peak_triplets_left,
        "peak_triplets_right" = x$peak_triplets_right,
        "integrals" = x$integrals,
        "signal_free_region" = x$signal_free_region %||% signal_free_region,
        "range_water_signal_ppm" = x$range_water_signal_ppm %||% range_water_signal_ppm,
        "A" = x$A,
        "lambda" = x$lambda,
        "x_0" = x$w # [19] rename
    )
    return(y)
}

# Private Helpers: read_spectrum #####

#' @title Load Spectrum
#' @description Loads a single spectrum file and returns the spectrum data in ppm.
#' @param path The path of the file containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"`.
#' @param type The type of the spectrum file. E.g. `"bruker"` or `"jcampdx"`.
#' @param sf A vector of two elements to scale the x and y axis values.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @return A list containing the spectrum data.
#' @details For details about `procno` and `expno` see section [File Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure) in the metabodecon FAQ.
#' @examples
#' \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1.dx")
#' spectrum_data <- load_spectrum(path, type = "jcampdx")
#' str(spectrum_data, 1)
#' }
#' @noRd
read_spectrum <- function(path, type = "bruker", sf = c(1e3, 1e6), expno = 10, procno = 10) {
    sfx <- sf[1]
    sfy <- sf[2]
    switch(type,
        "bruker" = load_bruker_spectrum_v1(path, sfx, sfy, expno, procno),
        "jcampdx" = load_jcampdx_spectrum_v1(path, sfx, sfy)
    )
}

#' @title Load single JCAMPDX Spectrum
#' @description Loads a single JCAMPDX spectrum file and returns the spectrum data in ppm.
#' @param path The path of the file containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"`.
#' @param scale_factor A vector of two elements to scale the x and y axis values. Default is c(1, 1).
#' @return A list containing the spectrum data.
#' \dontrun{
#' # Example Usage (took 30s on the development machine)
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1.dx")
#' spectrum_data <- load_jcampdx_spectrum(path)
#' str(spectrum_data, 1)
#' }
#' @noRd
load_jcampdx_spectrum_v1 <- function(path, sfx = 1e3, sfy = 1e6) {
    # List of 5
    # $ dataGuide   :'data.frame':   3 obs. of  3 variables:
    # ..$ Format   : chr [1:3] "metadata" "XRR" "XII"
    # ..$ FirstLine: num [1:3] 1 2181 10054
    # ..$ LastLine : num [1:3] 2180 10052 18572
    # $ metadata    : chr [1:2180] "##TITLE=GCKD_U_FA_1152_06_03_2020" ...
    # $ commentLines: int 10053
    # $ real        :'data.frame':   131072 obs. of  2 variables:
    # ..$ x: num [1:131072] 12019 12019 12019 12019 12019 ...
    # ..$ y: num [1:131072] 1265 1003 105 -937 -1062 ...
    # $ imaginary   :'data.frame':   131072 obs. of  2 variables:
    # ..$ x: num [1:131072] 12019 12019 12019 12019 12019 ...
    # ..$ y: num [1:131072] -35831 -36561 -36864 -36185 -34777 ...
    data <- readJDX::readJDX(file = path, SOFC = TRUE, debug = 0) # reading urine_1.dx (~1MB) takes ~30s on machine r31
    real <- data$real
    meta <- data$metadata
    n <- length(real$x) # number of data points (previously called `spectrum_length`)
    ppm_range <- as.numeric(sub("\\D+", "", meta[startsWith(meta, "##$SW=")]))
    ppm_max <- as.numeric(sub("\\D+", "", meta[startsWith(meta, "##$OFFSET=")]))
    ppm_min <- ppm_max - ppm_range
    ppm_step <- ppm_range / (n - 1) # Example: data points in ppm = 1.1, 2.3, 3.5, 4.7 --> ppm_step == 1.2
    ppm_nstep <- ppm_range / n # Don't really know what this is, but it's used in later calculations
    ppm <- seq(ppm_max, ppm_min, -ppm_step) # ppm (previously called `x_ppm`)
    dp <- seq(n - 1, 0, -1) # data points
    sdp <- seq((n - 1) / sfx, 0, -1 / sfx) # scaled data points (previously called `x`). Same as `dp / sfx`, but with slight numeric differences, so we stick with the old calculation method for backwards compatibility.
    return(list(
        Y = list(raw = real$y, scaled = real$y / sfy),
        n = n, sfx = sfx, sfy = sfy, # misc
        dp = dp, sdp = sdp, ppm = ppm, # x-axis
        ppm_min = ppm_min, ppm_max = ppm_max, ppm_range = ppm_range, ppm_step = ppm_step, ppm_nstep = ppm_nstep # additional ppm info
        # , length = n, x_ppm = ppm, x = sdp, ppm_highest_value = ppm_max, ppm_lowest_value = ppm_min # backwards compatible names
    ))

    # TODO: return as few values as possible, e.g.
    # >>> return(list(ppm = ppm, si = ss))
    # This should be sufficient to calculate everything else with super simple functions like
    # >>> max(ppm)
    # >>> seq(length(ppm) - 1, 0, -1)
    # etc.
}

#' @title Load single JCAMPDX or Bruker Spectrum
#' @description Loads a single JCAMPDX spectrum file and returns the spectrum data in ppm.
#' @param path The path of the file containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"`.
#' @param scale_factor A vector of two elements to scale the x and y axis values. Default is c(1, 1).
#' @return A list containing the spectrum data.
#' \dontrun{
#' # Example Usage (took 30s on the development machine)
#' xds_path <- download_example_datasets()
#' urine_1_dx <- file.path(xds_path, "jcampdx/urine/urine_1.dx")
#' system.time(spectrum_data <- load_jcampdx_spectrum(urine_1_dx, c(2, 2)))
#' str(spectrum_data, 1)
#' }
#' @details Never published and superseded by `load_jcampdx_spectrum_v1()`.
#' @noRd
load_jcampdx_spectrum <- function(path, scale_factor = c(1e3, 1e6)) {
    # Import data using readJDX package
    data <- readJDX::readJDX(file = path, SOFC = TRUE, debug = 0) # reading urine_1.dx (~1MB) takes ~30s on machine r31

    # Scale factors for x and y axis
    factor_x <- scale_factor[1]
    factor_y <- scale_factor[2]

    # Get the length of the spectrum
    spectrum_length <- length(data[[4]]$x)
    spectrum_length_m1 <- spectrum_length - 1

    # Scale the x and y axis values
    spectrum_x <- seq((spectrum_length_m1 / factor_x), 0, -1 / factor_x)
    spectrum_y <- (data[[4]]$y) / factor_y

    # Extract spectral width in ppm
    ppm_range_index <- which(startsWith(data[[2]], "##$SW="))
    ppm_range <- as.numeric(sub("\\D+", "", data[[2]][ppm_range_index]))

    # Extract highest and lowest ppm values
    ppm_highest_index <- which(startsWith(data[[2]], "##$OFFSET="))
    ppm_highest_value <- as.numeric(sub("\\D+", "", data[[2]][ppm_highest_index]))
    ppm_lowest_value <- ppm_highest_value - ppm_range

    # Generate ppm sequence for x axis
    spectrum_x_ppm <- seq(ppm_highest_value, ppm_lowest_value, (-ppm_range / spectrum_length_m1))

    # Return the spectrum data
    return(list(
        x = spectrum_x,
        y = spectrum_y,
        x_ppm = spectrum_x_ppm,
        length = spectrum_length,
        ppm_range = ppm_range,
        ppm_highest_value = ppm_highest_value,
        ppm_lowest_value = ppm_lowest_value
    ))
}

#' @title Load single Bruker Spectrum
#' @description Loads a single Bruker spectrum and returns the spectrum data in ppm.
#' @param path The path of the directory holding the spectrum data. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param sfx Scaling factor for x axis in datapoints. E.g. `"1000"`. Only relevant for plotting.
#' @param sfy Scaling factor for y axis (signal strength). E.g. `"100000"`. Only relevant for plotting.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @return A list containing the spectrum data.
#' @details For details about `procno` and `expno` see section [File Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure) in the metabodecon FAQ.
#' @examples
#' \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "bruker/urine/urine_1")
#' spectrum_data <- load_bruker_spectrum_v1(path) # takes 3s on machine r31
#' str(spectrum_data, 1)
#' }
#' @noRd
load_bruker_spectrum_v1 <- function(path = file.path(download_example_datasets(), "bruker/urine/urine_1"),
                                    sfx = 1e3, sfy = 1e6,
                                    expno = 10, procno = 10) {
    acqus <- readLines(file.path(path, expno, "acqus"))
    procs <- readLines(file.path(path, expno, "pdata", procno, "procs"))
    ppm_range <- as.numeric(sub("\\D+", "", acqus[startsWith(acqus, "##$SW=")]))
    y <- read_1r_file(path = path, expno = expno, procno = procno, procs = procs) # signal strength (could also be called signal intensity)
    n <- length(y)
    ppm_max <- as.numeric(sub("\\D+", "", procs[startsWith(procs, "##$OFFSET=")]))
    ppm_min <- ppm_max - ppm_range
    ppm_step <- ppm_range / (n - 1) # Example: data points in ppm = 1.1, 2.3, 3.5, 4.7 --> ppm_step == 1.2
    ppm_nstep <- ppm_range / n # Not really useful, but we need it for backwards compatibility with MetaboDecon1D results
    ppm <- seq(ppm_max, ppm_max - ppm_range, by = -ppm_step) # parts per million
    dp <- seq(n - 1, 0, -1) # data points
    sdp <- seq((n - 1) / sfx, 0, -1 / sfx) # scaled data points (previously called `x`). Same as `dp / sfx`, but with slight numeric differences, so we stick with the old calculation method for backwards compatibility.
    return(list(
        Y = list(raw = y, scaled = y / sfy), # y-axis
        n = n, sfx = sfx, sfy = sfy, # misc
        dp = dp, sdp = sdp, ppm = ppm, # x-axis
        ppm_min = ppm_min, ppm_max = ppm_max, ppm_range = ppm_range, ppm_step = ppm_step, ppm_nstep = ppm_nstep # additional ppm info
        # , length = n, x_ppm = ppm, x = sdp, ppm_highest_value = ppm_max, ppm_lowest_value = ppm_min # backwards compatible names
    ))
}

#' @title Load single Bruker Spectrum
#' @description Loads a single Bruker spectrum and returns the spectrum data in ppm.
#' @param path The path of the directory holding the spectrum data. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param sf A vector of two elements to scale the x and y axis values. Default is c(1, 1).
#' @param proc The processing value for the file. E.g. `"10"`.
#' @param spec The spectroscopy value for the file. E.g. `"10"`.
#' @return A list containing the spectrum data.
#' @details For details about `proc` and `spec` see section [File Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure) in the metabodecon FAQ.
#' @examples
#' \dontrun{
#' xds_path <- download_example_datasets()
#' urine_1 <- file.path(xds_path, "bruker/urine/urine_1")
#' system.time(spectrum_data <- load_bruker_spectrum(urine_1)) # takes 3s on machine r31
#' str(spectrum_data, 1)
#' }
#' @details Never published and superseded by `load_bruker_spectrum_v1()`.
#' @noRd
load_bruker_spectrum <- function(path, sf = c(1e3, 1e6), spec = 10, proc = 10) {
    # Get file paths of all necessary documents
    acqus_file <- file.path(path, spec, "acqus")
    spec_file <- file.path(path, spec, paste("pdata", proc, "1r", sep = "/"))
    procs_file <- file.path(path, spec, paste("pdata", proc, "procs", sep = "/"))

    # Read content of files
    acqs_content <- readLines(acqus_file)
    procs_content <- readLines(procs_file)

    # Load necessary meta data of acqa content
    index_spectral_width <- which(startsWith(acqs_content, "##$SW="))
    ppm_range <- as.numeric(sub("\\D+", "", acqs_content[index_spectral_width]))

    # Load necessary meta data of procs content
    index_integer_type <- which(startsWith(procs_content, "##$DTYPP="))
    index_highest_ppm <- which(startsWith(procs_content, "##$OFFSET="))
    index_data_points <- which(startsWith(procs_content, "##$SI="))

    # Save digit value
    ppm_highest_value <- as.numeric(sub("\\D+", "", procs_content[index_highest_ppm]))
    data_points <- as.numeric(sub("\\D+", "", procs_content[index_data_points]))
    integer_type <- as.numeric(sub("\\D+", "", procs_content[index_integer_type]))

    # Establish connection , rb = read binary
    to_read <- file(spec_file, "rb")

    # If integer_type = 0, then size of each int is 4, else it is 8
    size_int <- ifelse(integer_type == 0, 4, 8)

    # Calculate lowest ppm value
    ppm_lowest_value <- ppm_highest_value - ppm_range

    # Read binary spectrum
    spectrum_y <- readBin(to_read, what = "int", size = size_int, n = data_points, signed = TRUE, endian = "little")
    close(to_read)

    # Choose factor to scale axis values
    factor_x <- sf[1]
    factor_y <- sf[2]

    # Calculate length of spectrum
    spectrum_length <- length(spectrum_y)

    # Calculate x axis
    spectrum_x_ppm <- seq(ppm_highest_value, ppm_highest_value - ppm_range, by = -(ppm_range / (spectrum_length - 1)))

    # Calculate spectrum_x and recalculate spectum_y
    spectrum_x <- seq((spectrum_length - 1) / factor_x, 0, -0.001)
    spectrum_y <- spectrum_y / factor_y

    # Return the spectrum data
    return(list(
        x = spectrum_x,
        y = spectrum_y,
        x_ppm = spectrum_x_ppm,
        length = spectrum_length,
        ppm_range = ppm_range,
        ppm_highest_value = ppm_highest_value,
        ppm_lowest_value = ppm_lowest_value
    ))
}

#' @title Read Bruker TopSpin spectrum from 1r file
#' @param path The path of the directory holding the spectrum data. E.g. `"example_datasets/bruker/urine/urine_1/"`.
#' @param procno The processing number for the file. E.g. `"10"`.
#' @param expno The experiment number for the file. E.g. `"10"`.
#' @param procs The content of the `procs` file. E.g. `readLines(file.path(path, expno, "pdata", procno, "procs"))`.
#' @examples \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "bruker/urine/urine_1")
#' procs <- readLines(file.path(path, 10, "pdata", 10, "procs"))
#' y <- read_1r_file(path, 10, 10, procs)
#' }
#' @return The signals read from the file as numeric vector
read_1r_file <- function(path, expno, procno, procs) {
    # TODO: see issue `Check: check use of DTYPP in load_spectrum` in TODOS.R
    n <- as.numeric(sub("\\D+", "", procs[startsWith(procs, "##$SI=")]))
    int_type <- as.numeric(sub("\\D+", "", procs[startsWith(procs, "##$DTYPP=")]))
    int_size <- if (int_type == 0) 4 else 8
    path_1r <- file.path(path, expno, "pdata", procno, "1r")
    spec_stream <- file(path_1r, "rb")
    on.exit(close(spec_stream), add = TRUE)
    y <- readBin(spec_stream, what = "int", size = int_size, n = n, signed = TRUE, endian = "little")
    return(y)
}

# Private Helpers: signal_free_region #####

#' @title Determine Signal Free Region
#' @description This function determines the signal free region (SFR) of a given spectrum. It asks the user to confirm the left and right borders of the SFR, and allows them to adjust these borders if necessary. The function returns a list containing the left and right borders in both ppm and data points (dp), as well as the scaled data points (sdp).
#' @param spec A list representing the spectrum, which should include the minimum and maximum ppm (`$ppm_min` and `$ppm_max`), and the scaling factor (`$sfx`).
#' @param sfr Initial values for the left and right borders of the SFR in ppm. If not provided, the function will ask the user to select the borders.
#' @param ask Logical. If TRUE, the function will ask the user to confirm or adjust the borders of the SFR. Default is TRUE.
#' @param bwc Use the old, slightly incorrect method for conversion from ppm to data points to maintain backwards compatibility with MetaboDecon1D results? For details see issue `Check: ppm to dp conversion` in TODOS.md
#' @return A list containing the left and right borders of the SFR in ppm (`$left_ppm` and `$right_ppm`), data points (`$left_dp` and `$right_dp`), and scaled data points (`$left_sdp` and `$right_sdp`).
#' @noRd
determine_signal_free_region <- function(spec, sfr = NULL, ask = TRUE, bwc = TRUE) {
    left_ppm <- sfr[1]
    right_ppm <- sfr[2]
    if (is.null(sfr) && isFALSE(ask)) {
        stop("No signal free region (SFR) provided and `ask` is FALSE. Please provide the SFR or set `ask` to TRUE.")
    }
    if (ask) {
        plot_sfr(spec, left_ppm, right_ppm)
        sfr_ok <- get_yn_input("Signal free region borders correct selected? (Area left and right of the green lines)")
        while (!sfr_ok) {
            left_ppm <- get_num_input("Choose another left border: [e.g. 12]", min = spec$ppm_min, max = spec$ppm_max)
            right_ppm <- get_num_input("Choose another right border: [e.g. -2]", min = spec$ppm_min, max = spec$ppm_max)
            plot_sfr(spec, left_ppm, right_ppm)
            sfr_ok <- get_yn_input("Signal free region borders correct selected? (Area left and right of the green lines)")
        }
    }
    left_dp <- ppm_to_dp(left_ppm, spec, bwc)
    left_sdp <- left_dp / spec$sfx
    right_dp <- ppm_to_dp(right_ppm, spec, bwc)
    right_sdp <- right_dp / spec$sfx
    spec$sfr <- list(left_ppm = left_ppm, right_ppm = right_ppm, left_sdp = left_sdp, right_sdp = right_sdp, left_dp = left_dp, right_dp = right_dp)
    spec
}

#' @title Plot Signal Free Region
#' @description Draws the signal free region as green vertical lines into the given spectrum.
#' @param spec A list representing the spectrum as returned by [load_jcampdx_spectrum()] or [load_bruker_spectrum()].
#' @param left_ppm The left border of the signal free region in ppm.
#' @param right_ppm The right border of the signal free region in ppm.
#' @return NULL. Called for side effect of plotting the signal free region.
plot_sfr <- function(spec, left_ppm, right_ppm) {
    plot(
        x = spec$ppm,
        y = spec$Y$scaled,
        type = "l",
        xlab = "[ppm]",
        ylab = "Intensity [a.u.]",
        xlim = c(spec$ppm_max, spec$ppm_min)
    )
    graphics::abline(v = c(left_ppm, right_ppm), col = "green")
}

# Private Helpers: water_signal #####

#' @title Calculate water signal parameters
#' @description Calculates water signal parameters for a given spectrum.
#' @param spec A list representing the spectrum.
#' @param hwidth_ppm The half width in ppm. Default is wshw.
#' @param bwc Use the old, slightly incorrect methods for calculating water signal values to maintain backwards compatibility with MetaboDecon1D results? For details see issue `Check: water signal calculation` in `TODOS.md`.
#' @return List of parameters including half width in dp and ppm, center line in dp and ppm and right and left borders in dp and ppm.
determine_water_signal <- function(spec, hwidth_ppm, bwc = TRUE, ask = TRUE) {
    if (ask) {
        plot_ws(spec, hwidth_ppm)
        ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
        while (!ws_ok) {
            hwidth_ppm <- get_num_input("Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154]")
            plot_ws(spec, hwidth_ppm)
            ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
        }
    }
    hwidth_dp <- if (bwc) hwidth_ppm / spec$ppm_nstep else hwidth_ppm / spec$ppm_step # half width in dp
    center_dp <- if (bwc) spec$n / 2 else (spec$n - 1) / 2 # center line in dp
    right_dp <- center_dp + hwidth_dp # right border in dp
    left_dp <- center_dp - hwidth_dp # left border in dp
    center_ppm <- if (bwc) spec$ppm[center_dp] else dp_to_ppm(center_dp, spec) # center in ppm
    right_ppm <- if (bwc) spec$ppm[right_dp] else dp_to_ppm(right_dp, spec) # right border in ppm
    left_ppm <- if (bwc) spec$ppm[left_dp] else dp_to_ppm(left_dp, spec) # left border in ppm
    spec$ws <- list(center_dp = center_dp, hwidth_dp = hwidth_dp, left_dp = left_dp, right_dp = right_dp, center_ppm = center_ppm, hwidth_ppm = hwidth_ppm, left_ppm = left_ppm, right_ppm = right_ppm)
    spec
}

remove_water_signal <- function(spec, bwc = TRUE) {
    if (is.null(spec$ws)) {
        stop("No water signal parameters found. Please call `determine_water_signal()`.")
    }
    y <- spec$Y$scaled
    if (bwc) {
        left <- spec$ws$left_dp
        right <- spec$ws$right_dp
        y[right:left] <- 0.01 / spec$sfy # (1)
        # Order, i.e. `right:left` instead of `left:right`, is important here, because `right` and `left` are floats in the backwards compatible case. Example: `right <- 3.3; left <- 1.4` ==> `right:left == c(3.3, 2.3)` and `left:right == c(1.4, 2.4)`.
    } else {
        ppm <- spec$ppm
        ws_min <- spec$ws$center_ppm - spec$ws$hwidth_ppm
        ws_max <- spec$ws$center_ppm + spec$ws$hwidth_ppm
        ws_idx <- which(ppm >= ws_min & ppm <= ws_max)
        y[wsidx] <- 0
    }
    spec$Y$nows <- y
    spec
}

#' @title Plot Water Signal
#' @description Draws the water signal as red vertical lines into the given spectrum.
#' @param spec A list representing the spec as returned by [load_jcampdx_spectrum()] or [load_bruker_spectrum()].
#' @param hwidth_ppm The half width of the water signal in ppm.
#' @return NULL. Called for side effect of plotting the water signal.
plot_ws <- function(spec, hwidth_ppm) {
    center_ppm <- (spec$ppm_max + spec$ppm_min) / 2
    plot(
        spec$ppm,
        spec$y,
        type = "l",
        xlab = "[ppm]",
        ylab = "Intensity [a.u.]",
        xlim = c(center_ppm + 2 * hwidth_ppm, center_ppm - 2 * hwidth_ppm)
    )
    graphics::abline(v = center_ppm + hwidth_ppm, col = "red")
    graphics::abline(v = center_ppm - hwidth_ppm, col = "red")
}

# Private Helpers: remove_negative_signals #####

remove_negative_signals <- function(spec) {
    if (is.null(spec$Y$nows)) stop("Water signal not removed yet. Please call `remove_water_signal()` first.")
    spec$Y$pos <- abs(spec$Y$nows)
    spec
}

# Private Helpers: smooth #####

#' @inherit smooth_signals_v1
#' @param bwc Maintain backwards compatibility with MetaboDecon1D results by using the old and slow method for smoothing ([smooth_signals_v1()])? If FALSE, the new and fast method ([smooth_signals_v2()]) is used instead.
#' @noRd
smooth_signals <- function(spec, reps = 2, k = 5, bwc = TRUE) {
    if (bwc) {
        smooth_signals_v1(spec, reps, k)
    } else {
        smooth_signals_v2(spec, reps, k)
    }
}

#' @inherit smooth_signals_v1
#' @details New and fast version for smoothing of signals. Implements the same algorithm as `smooth_signal_v1` using different R functions (e.g. [stats::filter()]), causing a massive speedup but also numeric differences compared to the old version.
#' @noRd
smooth_signals_v2 <- function(spec, reps = 2, k = 5) {
    if (k %% 2 == 0) stop("k must be odd")

    Z <- vector("list", length = reps)
    y <- spec$Y$pos
    n <- length(y)

    for (i in 1:reps) {
        filter <- rep(1 / k, k)
        z <- stats::filter(y, filter, sides = 2) # (1)
        q <- (k - 1) / 2 # (2)
        for (j in seq_len(q)) {
            z[j] <- mean(y[1:(q + j)]) # (3)
            z[n - j + 1] <- mean(y[(n - q - j + 1):n]) # (4)
        }
        y <- Z[[i]] <- as.numeric(z)
        # Calling (1) gives NAs at both sides of vector, as there are not enough values for the moving average. The number of NAs at each side is given by (2). Example: if n==100 and k==5, then q==2, so z[1]==NA, z[2]==NA, z[99]==NA and z[100]==NA. To stay backwards compatible, these values must be filled with the mean of the values that are available. To do so, we iterate from 1:q, i.e. j==1 and j==2 and set
        # z[1]   <- mean(y[1:3])    # 3 == 2+1 == q+j            # (3)
        # z[2]   <- mean(y[1:4])    # 4 == 2+2 == q+j            # (3)
        # z[100] <- mean(y[98:100]) # 98 == 100-2-1+1 == n-q-j+1 # (4)
        # z[99]  <- mean(y[97:100]) # 97 == 100-2-2+1 == n-q-j+1 # (4)
        # Note: we could also think of leaving the NAs as they are, which would be more correct I think and even faster, but would break compatibility with the old version completely. So not even `all.equal(v1, v2)` would be TRUE anymore.
    }

    spec$Z <- Z
    spec$Y$smooth <- Z[[reps]]
    spec
}

#' @title Smooth signal intensities using a moving average
#' @description This function smooths signal intensities by applying a [moving average](https://en.wikipedia.org/wiki/Moving_average) filter with a window size of k.
#' @param spec A list representing the spectrum, which should include the scaled signal intensities, after removal of the water artefact and negative values (`spec$Y$pos`).
#' @param reps The number of times to apply the moving average.
#' @param k The number of points within the moving average window. Must be odd, so the smoothed point is in the middle of the window.
#' @return A numeric vector of the smoothed values.
#' @details Old and slow version producing the same results as the implementation within `deconvolution` from `MetaboDecon1D_deconvolution.R`.
#' @noRd
smooth_signals_v1 <- function(spec, reps = 2, k = 5) {
    if (k %% 2 == 0) stop("k must be odd")
    Z <- vector("list", length = reps)
    y <- spec$Y$pos
    n <- length(y)
    for (i in 1:reps) {
        z <- y
        for (j in 1:(n)) {
            left_border <- j - floor(k / 2)
            right_border <- j + floor(k / 2)
            if (left_border <= 0) {
                left_border <- 1
                z[j] <- (1 / right_border) * sum(y[left_border:right_border])
            } else if (right_border >= n) {
                right_border <- n
                z[j] <- (1 / (right_border - left_border + 1)) * sum(y[left_border:right_border])
            } else {
                z[j] <- (1 / k) * sum(y[left_border:right_border])
            }
        }
        y <- Z[[i]] <- as.numeric(z)
    }
    spec$Z <- Z
    spec$Y$smooth <- Z[[reps]]
    spec
}

# Private Helpers: peak_selection #####


#' Plot peaks of a spectrum
#'
#' This function plots the peaks of a spectrum, including the smoothed and scaled signal intensity and the second derivative.
#' It also allows for the specification of peak positions and the option to draw vertical lines at these positions.
#'
#' @param spec A data frame containing the spectrum data. It should have columns 'ppm', 'Y', 'ip', 'ip_left', 'ip_right', and 'd'.
#' @param ppm A vector of length 2 specifying the range of ppm values to consider for the plot. Default is c(3.402, 3.437).
#' @param dp A vector specifying the positions of the peaks. If NULL (default), the function will determine the peak positions based on the 'ppm' range.
#' @param vlines A logical value indicating whether to draw vertical lines at the peak positions. Default is FALSE.
#' @return A data frame with columns 'x' (ppm values), 'y' (smoothed and scaled signal intensity), 'd' (second derivative), 'is_ip' (whether the position is a peak), and 'is_ip_left' (whether the position is to the left of a peak).
#' @examples \dontrun{
#' plot_peaks(spec, ppm = c(3.402, 3.437), dp = NULL, vlines = FALSE) # region from 3.402 to 3.437 ppm
#' plot_peaks(spec, dp = 1:200, vlines = FALSE) # first 200 data points
#' }
#' @noRd
plot_peaks <- function(spec, ppm = c(3.402, 3.437), dp = NULL, vlines = FALSE) {
    if (is.null(dp)) dp <- which(spec$ppm > min(ppm) & spec$ppm < max(ppm))
    x <- spec$ppm[dp]
    y <- spec$Y$smooth[dp]
    ip <- which(dp %in% spec$peaks)
    ip_left <- which(dp %in% spec$ip_left) %||% numeric()
    ip_right <- which(dp %in% spec$ip_right) %||% numeric()
    d <- spec$d[dp]
    withr::with_par(list(mfrow = c(2, 1), mar = c(0, 6, 4, 2), las = 1), {
        plot(x, y, type = "l", xlab = "ppm", ylab = "", xaxt = "n", xlim = c(max(x), min(x)))
        mtext("smoothed and scaled signal intensity", side = 2, line = 5, las = 0)
        points(x, y, type = "p")
        points(x[ip], y[ip], col = "red", pch = 19)
        points(x[ip_left], y[ip_left], col = "blue", pch = 19)
        points(x[ip_right], y[ip_right], col = "cyan4", pch = 19)
        if (vlines) {
            abline(v = x[ip], col = "red", lty = 2, lwd = 0.5)
            abline(v = x[ip_left], col = "blue", lty = 2, lwd = 0.5)
            abline(v = x[ip_right], col = "cyan4", lty = 2, lwd = 0.5)
        }
        axis(3, at = x, labels = dp)
        legend("topright", legend = c("peaks", "left", "right"), col = c("red", "blue", "cyan4"), pch = 19)
        withr::with_par(list(mar = c(5, 6, 0, 2)), {
            plot(x, d, type = "l", xlab = "ppm", ylab = "", xlim = c(max(x), min(x)))
            mtext("second derivative", side = 2, line = 5, las = 0)
            points(x, d, type = "p")
            points(x[ip], d[ip], col = "red", pch = 19)
            points(x[ip_left], d[ip_left], col = "blue", pch = 19)
            points(x[ip_right], d[ip_right], col = "cyan4", pch = 19)
            if (vlines) {
                abline(v = x[ip], col = "red", lty = 2, lwd = 0.5)
                abline(v = x[ip_left], col = "blue", lty = 2, lwd = 0.5)
                abline(v = x[ip_right], col = "cyan4", lty = 2, lwd = 0.5)
            }
            abline(h = 0, col = "black", lty = 2, lwd = 0.5)
        })
    })
    df <- data.frame(x = x, y = y, d = d, is_ip = dp %in% spec$peaks, is_ip_left = dp %in% spec$ip_left)
    return(df)
}

select_peaks_v2 <- function(spec, bwc = TRUE) {
    d <- calculate_second_derivative(y = spec$Y$smooth, bwc = bwc)
    m <- length(d)
    d0 <- d[1:(m - 2)]
    d1 <- d[2:(m - 1)]
    d2 <- d[3:(m - 0)]
    spec$peaks <- which(d1 < 0 & d1 <= d0 & d1 < d2) + 1 # TODO: why <= instad of < once?
    spec$d <- d
    return(spec)
}


#' @inherit select_peaks_v0
#' @param details This is the same as `select_peaks_v0`, but using different variable names within the code to make it more readable.
select_peaks_v1 <- function(spec) {
    y <- spec$ssss
    x <- spec$sdp
    n <- spec$n
    # Calculate second derivative of smoothed and scaled data points
    d <- matrix(nrow = 2, ncol = n - 2)
    for (i in 2:n - 1) { # Hack: this should be 2:(n-1), because we want to iterate from 2 to (n-1), but the actual code iterates from 1 to (n-1), so we have `d[1, 0] <- x[1]` in the first iteration. It's not too bad, because R silently ignores this error, so the code works. But again, it's very hard to reason about the code.
        d[1, i - 1] <- x[i]
        d[2, i - 1] <- y[i - 1] + y[i + 1] - 2 * y[i]
    }
    # Find all peaks, i.e. local minima of second derivative (in fact a local minima of the second derivative is not necessarily a peak, but we will call it a peak here for simplicity)
    ipx <- c()
    ipi <- c()
    for (i in 2:(ncol(d) - 1)) {
        if (d[2, i] < 0) {
            if ((d[2, i] <= d[2, i - 1]) & (d[2, i] < d[2, i + 1])) {
                # Add local minima to peak list
                ipx <- c(ipx, d[1, i])
                ipi <- c(ipi, i)
            }
        }
    }
    return(list(d2 = d[2, ], sdp = ipx, idx = ipi))
}

#' @title Select peaks
#' @description This function calculates the second derivative of a spec and finds all local minima of the derivate, i.e. peaks. It returns a list containing the second derivative, the x-values of the local minima, and their indices.
#' @param spec A list with elements `sdp` and `ssss` representing the x-values in scaled data points and the y values in "smoothed and scaled signal strength" respectively.
#' @param details this is the code copy pasted from the original `deconvolution` function of `MetaboDecon1D` with only the inputs renamed to `spec$sdp` and `spec$ssss`.
select_peaks_v0 <- function(spec) {
    # Calculate second derivative of spec
    second_derivative <- matrix(nrow = 2, ncol = length(spec$sdp) - 2)
    for (i in 2:length(spec$sdp) - 1) {
        second_derivative[1, i - 1] <- spec$sdp[i]
        second_derivative[2, i - 1] <- spec$ssss[i - 1] + spec$ssss[i + 1] - 2 * spec$ssss[i]
    }

    # Find all local minima of second derivative
    peaks_x <- c()
    peaks_index <- c()
    second_derivative_border <- ncol(second_derivative) - 1
    for (i in 2:second_derivative_border) {
        if (second_derivative[2, i] < 0) {
            if ((second_derivative[2, i] <= second_derivative[2, i - 1]) & (second_derivative[2, i] < second_derivative[2, i + 1])) {
                # Add local minima to peak list
                peaks_x <- c(peaks_x, second_derivative[1, i])
                peaks_index <- c(peaks_index, i)
            }
        }
    }
    return(list(
        d2 = second_derivative[2, ],
        sdp = peaks_x,
        idx = peaks_index
    ))
}


# Combination of `find_left_positions_v0` and `find_right_positions_v0`
find_peak_borders <- function(spec) {
    d <- spec$d
    peaks <- spec$peaks
    spec$left <- spec$right <- rep(NA, length(peaks))
    for (i in seq_along(peaks)) {
        j <- peaks[i] + 1
        while (j < length(d) && !is_left_border(j, d)) j <- j + 1
        spec$left[i] <- j
        j <- peaks[i] - 1
        while (j > 1 && !is_right_border(j, d)) j <- j - 1
        spec$right[i] <- j
    }
    spec
}

is_left_border <- function(j, d) {
    c1 <- d[j] > d[j - 1]
    c2 <- d[j] >= d[j + 1]
    c3 <- d[j] < 0
    c4 <- d[j + 1] >= 0
    return((c1 && c2) || (c1 && c3 && c4))
}

is_right_border <- function(j, d) {
    c1 <- d[j] > d[j + 1]
    c2 <- d[j] >= d[j - 1]
    c3 <- d[j] < 0
    c4 <- d[j - 1] >= 0
    return((c1 && c2) || (c1 && c3 && c4))
}

find_left_positions_v0 <- function(spec) {
    # ToSc: in fact this function searches points to the right of the peak, not to the left. So I think the comments from the original `deconvolution` function from `MetaboDecon1D.R` are wrong. However, for reference we will keep the wrong naming here and combine `find_left_positions_v0` and `find_right_positions_v0` into one function called `find_peak_borders`. This makes the code more readable and as a side effect solves the naming problem.
    d <- spec$d
    ip <- spec$peaks
    spec$left <- rep(NA, length(ip))
    for (i in seq_along(ip)) {
        j <- ip[i] + 1
        while (j < length(d)) {
            c1 <- d[j] > d[j - 1]
            c2 <- d[j] >= d[j + 1]
            c3 <- d[j] < 0
            c4 <- d[j + 1] >= 0
            if ((c1 && c2) || (c1 && c3 && c4)) {
                spec$left[i] <- j
                break
            }
            j <- j + 1
        }
    }
    spec
}

find_right_positions_v0 <- function(spec) {
    # ToSc: in fact this function searches points to the left of the peak, not to the right. So I think the comments from the original `deconvolution` function from `MetaboDecon1D.R` are wrong. However, for reference we will keep the wrong naming here and combine `find_right_positions_v0` and `find_left_positions_v0` into one function called `find_peak_borders`. This makes the code more readable and as a side effect solves the naming problem.
    d <- spec$d
    ip <- spec$peaks
    spec$right <- rep(NA, length(ip))
    for (i in seq_along(ip)) {
        j <- ip[i] - 1
        while (j >= 2) {
            c1 <- d[j] > d[j + 1]
            c2 <- d[j] >= d[j - 1]
            c3 <- d[j] < 0
            c4 <- d[j - 1] >= 0
            if (((c1) && (c2)) || (c1 && c3 && c4)) {
                spec$right[i] <- j
                break
            }
            j <- j - 1
        }
    }
    spec
}


calculate_second_derivative <- function(y, bwc) {
    n <- length(y)
    if (bwc) {
        x <- c(NA, y[-n]) # x[i] == y[i-1]
        z <- c(y[-1], NA) # z[i] == y[i+1]
        d <- x + z - 2 * y
    } else {
        # Using diff is almost equivalent to the above, but due to numeric instabilities, it sometimes gives slightly different results (e.g. -5.51600000000001e-05 instead of -5.51599999999998e-05). Since we do `<=` and `<` comparisons further below this is not backwards compatible. Nonetheless, we keep it here, because it would make the code slightly more readable (runtime is almost the same).
        d <- c(NA, diff(y, differences = 2), NA)
    }
}

get_peak_triplets <- function(spec) {
    # Check borders of peak triplets
    # If NA values are available, remove corresponding peak triplet
    for (i in length(left_position):1) {
        if (is.na(left_position[i]) | (is.na(right_position[i]))) {
            ip$sdp <- ip$sdp[-i]
            spec$peaks <- spec$peaks[-i]
            left_position <- left_position[-i]
            right_position <- right_position[-i]
        }
    }
}

get_preak_triplet_scores <- function(spec) {
    # Calculate peak triplet score to distinguish between signal and noise
    scores <- matrix(nrow = 1, ncol = length(ip))
    scores_left <- matrix(nrow = 1, ncol = length(ip))
    scores_right <- matrix(nrow = 1, ncol = length(ip))
    for (i in 1:length(ip)) {
        # Calculate left score
        left_score <- 0
        for (j in spec$peaks[i]:left_position[i]) {
            left_score <- sum(left_score, abs(d2[2, j]))
        }
        scores_left[i] <- left_score
        # Calculate right score
        right_score <- 0
        for (k in right_position[i]:spec$peaks[i]) {
            right_score <- sum(right_score, abs(d2[2, k]))
        }
        scores_right[i] <- right_score
        # Save minimum score
        scores[i] <- min(left_score, right_score)
    }
}

filter_peak_triplets <- function(spec) {
    index_left <- which(spec$sdp[spec$peaks + 1] >= sfrl_sdp)
    index_right <- which(spec$sdp[spec$peaks + 1] <= sfrr_sdp)
    mean_score <- mean(c(scores[index_left], scores[index_right]))
    sd_score <- stats::sd(c(scores[index_left], scores[index_right]))
    filtered_peaks <- c()
    filtered_left_position <- c()
    filtered_right_position <- c()
    save_scores <- c()
    for (i in 1:length(ip)) {
        if (scores[i] >= mean_score + delta * sd_score) {
            # Save peak position
            filtered_peaks <- c(filtered_peaks, spec$peaks[i])
            # Save left position
            filtered_left_position <- c(filtered_left_position, left_position[i])
            # Save right position
            filtered_right_position <- c(filtered_right_position, right_position[i])
            # Save value of scores of filtered peaks
            save_scores <- c(save_scores, scores[i])
        }
    }
}

calculate_initial_lorentz_curves <- function(spec) {
    # Calculate parameters w, lambda and A for the initial lorentz curves
    for (i in 1:length(filtered_peaks)) {
        # Calculate position of peak triplets
        w_1 <- c(w_1, spec$sdp[filtered_left_position[i] + 1])
        w_2 <- c(w_2, spec$sdp[filtered_peaks[i] + 1])
        w_3 <- c(w_3, spec$sdp[filtered_right_position[i] + 1])

        # Calculate intensity of peak triplets
        y_1 <- c(y_1, spec$y[filtered_left_position[i] + 1])
        y_2 <- c(y_2, spec$y[filtered_peaks[i] + 1])
        y_3 <- c(y_3, spec$y[filtered_right_position[i] + 1])

        # Calculate mirrored points if necesccary
        # For ascending shoulders
        if ((y_1[i] < y_2[i]) & (y_2[i] < y_3[i])) {
            w_3[i] <- 2 * w_2[i] - w_1[i]
            y_3[i] <- y_1[i]
        }
        # For descending shoulders
        if ((y_1[i] > y_2[i]) & (y_2[i] > y_3[i])) {
            w_1[i] <- 2 * w_2[i] - w_3[i]
            y_1[i] <- y_3[i]
        }

        # Move triplet to zero position
        w_delta[i] <- w_1[i]
        w_1[i] <- w_1[i] - w_delta[i]
        w_2[i] <- w_2[i] - w_delta[i]
        w_3[i] <- w_3[i] - w_delta[i]

        # Calculate difference of position of peak triplets
        w_1_2 <- c(w_1_2, w_1[i] - w_2[i])
        w_1_3 <- c(w_1_3, w_1[i] - w_3[i])
        w_2_3 <- c(w_2_3, w_2[i] - w_3[i])

        # Calculate difference of intensity values of peak triplets
        y_1_2 <- c(y_1_2, y_1[i] - y_2[i])
        y_1_3 <- c(y_1_3, y_1[i] - y_3[i])
        y_2_3 <- c(y_2_3, y_2[i] - y_3[i])

        # Calculate w for each peak triplet
        w_result <- (w_1[i]^2 * y_1[i] * y_2_3[i] +
            w_3[i]^2 * y_3[i] * y_1_2[i] +
            w_2[i]^2 * y_2[i] * (-y_1_3[i])
        ) / (
            2 * w_1_2[i] * y_1[i] * y_2[i] -
                2 * (w_1_3[i] * y_1[i] + (-w_2_3[i]) * y_2[i]) * y_3[i]
        )
        w_result <- w_result + w_delta[i]
        w <- c(w, w_result)
        # Wenn y Werte nach der Hoehenanpassung 0 werden, so ist w_new[i] NaN
        if (is.nan(w[i])) {
            w[i] <- 0
        }

        # Calculate lambda for each peak triplet
        lambda_result <- -((sqrt(abs((-w_2[i]^4 * y_2[i]^2 * y_1_3[i]^2 - w_1[i]^4 * y_1[i]^2 * y_2_3[i]^2 - w_3[i]^4 * y_1_2[i]^2 * y_3[i]^2 + 4 * w_2[i] * w_3[i]^3 * y_2[i] * ((-y_1[i]) + y_2[i]) * y_3[i]^2 + 4 * w_2[i]^3 * w_3[i] * y_2[i]^2 * y_3[i] * ((-y_1[i]) + y_3[i]) + 4 * w_1[i]^3 * y_1[i]^2 * y_2_3[i] * (w_2[i] * y_2[i] - w_3[i] * y_3[i]) + 4 * w_1[i] * y_1[i] * (w_2[i]^3 * y_2[i]^2 * y_1_3[i] - w_2[i] * w_3[i]^2 * y_2[i] * (y_1[i] + y_2[i] - 2 * y_3[i]) * y_3[i] + w_3[i]^3 * y_1_2[i] * y_3[i]^2 - w_2[i]^2 * w_3[i] * y_2[i] * y_3[i] * (y_1[i] - 2 * y_2[i] + y_3[i])) + 2 * w_2[i]^2 * w_3[i]^2 * y_2[i] * y_3[i] * (y_1[i]^2 - 3 * y_2[i] * y_3[i] + y_1[i] * (y_2[i] + y_3[i])) + 2 * w_1[i]^2 * y_1[i] * (-2 * w_2[i] * w_3[i] * y_2[i] * y_3[i] * (-2 * y_1[i] + y_2[i] + y_3[i]) + w_3[i]^2 * y_3[i] * (y_1[i] * (y_2[i] - 3 * y_3[i]) + y_2[i] * (y_2[i] + y_3[i])) + w_2[i]^2 * y_2[i] * (y_1[i] * (-3 * y_2[i] + y_3[i]) + y_3[i] * (y_2[i] + y_3[i])))))))) / (2 * sqrt((w_1[i] * y_1[i] * y_2_3[i] + w_3[i] * y_1_2[i] * y_3[i] + w_2[i] * y_2[i] * ((-y_1[i]) + y_3[i]))^2))
        # If y and w are 0, then 0/0=NaN
        if (is.nan(lambda_result)) {
            lambda_result <- 0
        }
        lambda <- c(lambda, lambda_result)

        # Calculate scaling factor A for each peak triplet
        A_result <- (-4 * w_1_2[i] * w_1_3[i] * w_2_3[i] * y_1[i] * y_2[i] * y_3[i] * (w_1[i] * y_1[i] * y_2_3[i] + w_3[i] * y_3[i] * y_1_2[i] + w_2[i] * y_2[i] * (-y_1_3[i])) * lambda[i]) / (w_1_2[i]^4 * y_1[i]^2 * y_2[i]^2 - 2 * w_1_2[i]^2 * y_1[i] * y_2[i] * (w_1_3[i]^2 * y_1[i] + w_2_3[i]^2 * y_2[i]) * y_3[i] + (w_1_3[i]^2 * y_1[i] - w_2_3[i]^2 * y_2[i])^2 * y_3[i]^2)
        # If y and w are 0, then 0/0=NaN
        if (is.nan(A_result)) {
            A_result <- 0
        }
        A <- c(A, A_result)
    }

    # Calculate all initial lorentz curves
    lorentz_curves_initial <- matrix(nrow = length(filtered_peaks), ncol = length(spec$sdp))
    for (i in 1:length(filtered_peaks)) {
        # If A = 0, then the lorentz curve is a zero line
        if (A[i] == 0) {
            lorentz_curves_initial[i, ] <- 0
        } else {
            lorentz_curves_initial[i, ] <- abs(A[i] * (lambda[i] / (lambda[i]^2 + (spec$sdp - w[i])^2)))
        }
    }
}

refine_lorentz_curvs <- function(spec) {
    # Parameter approximation method
    w_1 <- c()
    w_2 <- c()
    w_3 <- c()
    y_1 <- c()
    y_2 <- c()
    y_3 <- c()
    w_1_2 <- c()
    w_1_3 <- c()
    w_2_3 <- c()
    y_1_2 <- c()
    y_1_3 <- c()
    y_2_3 <- c()
    w_delta <- c()
    w <- c()
    lambda <- c()
    A <- c()

    # Approximation of lorentz curves
    for (b in 1:nfit) {
        # Calculate new heights of peak triplets
        w_1_new <- c()
        w_2_new <- c()
        w_3_new <- c()
        y_1_new <- c()
        y_2_new <- c()
        y_3_new <- c()
        w_1_2_new <- c()
        w_1_3_new <- c()
        w_2_3_new <- c()
        y_1_2_new <- c()
        y_1_3_new <- c()
        y_2_3_new <- c()
        w_delta_new <- c()
        w_new <- c()
        lambda_new <- c()
        A_new <- c()
        sum_left <- c()
        sum_peaks <- c()
        sum_right <- c()
        proportion_left <- c()
        proportion_peaks <- c()
        proportion_right <- c()

        for (i in 1:length(filtered_peaks)) {
            # Calculate the position of the peak triplets
            w_1_new <- c(w_1_new, spec$sdp[filtered_left_position[i] + 1])
            w_2_new <- c(w_2_new, spec$sdp[filtered_peaks[i] + 1])
            w_3_new <- c(w_3_new, spec$sdp[filtered_right_position[i] + 1])

            # Calculate the sum of all lorentz curves for each data point
            sum_left[i] <- sum(lorentz_curves_initial[1:length(filtered_left_position), filtered_left_position[i] + 1])
            sum_peaks[i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), filtered_peaks[i] + 1])
            sum_right[i] <- sum(lorentz_curves_initial[1:length(filtered_right_position), filtered_right_position[i] + 1])

            # Calculate the proportion between original spectrum an the sum of the lorentz curves for each peak triplets position
            proportion_left[i] <- spec$y[filtered_left_position[i] + 1] / sum_left[i]
            proportion_peaks[i] <- spec$y[filtered_peaks[i] + 1] / sum_peaks[i]
            proportion_right[i] <- spec$y[filtered_right_position[i] + 1] / sum_right[i]

            # Calculate the new heights of the peak triplets
            y_1_new[i] <- lorentz_curves_initial[i, filtered_left_position[i] + 1] * proportion_left[i]
            y_2_new[i] <- lorentz_curves_initial[i, filtered_peaks[i] + 1] * proportion_peaks[i]
            y_3_new[i] <- lorentz_curves_initial[i, filtered_right_position[i] + 1] * proportion_right[i]

            # Calculate mirrored points if necesccary
            # For ascending shoulders
            if ((y_1_new[i] < y_2_new[i]) & (y_2_new[i] < y_3_new[i])) {
                w_3_new[i] <- 2 * w_2_new[i] - w_1_new[i]
                y_3_new[i] <- y_1_new[i]
            }
            # For descending shoulders
            if ((y_1_new[i] > y_2_new[i]) & (y_2_new[i] > y_3_new[i])) {
                w_1_new[i] <- 2 * w_2_new[i] - w_3_new[i]
                y_1_new[i] <- y_3_new[i]
            }

            # Move triplet to zero position
            w_delta_new[i] <- w_1_new[i]
            w_1_new[i] <- w_1_new[i] - w_delta_new[i]
            w_2_new[i] <- w_2_new[i] - w_delta_new[i]
            w_3_new[i] <- w_3_new[i] - w_delta_new[i]

            # Calculate difference of peak triplet positions
            w_1_2_new <- c(w_1_2_new, w_1_new[i] - w_2_new[i])
            w_1_3_new <- c(w_1_3_new, w_1_new[i] - w_3_new[i])
            w_2_3_new <- c(w_2_3_new, w_2_new[i] - w_3_new[i])

            # Calculate difference of new intensity values of peak triplets
            y_1_2_new <- c(y_1_2_new, y_1_new[i] - y_2_new[i])
            y_1_3_new <- c(y_1_3_new, y_1_new[i] - y_3_new[i])
            y_2_3_new <- c(y_2_3_new, y_2_new[i] - y_3_new[i])

            # Calculate w for each peak triplet
            w_result <- (w_1_new[i]^2 * y_1_new[i] * y_2_3_new[i] + w_3_new[i]^2 * y_3_new[i] * y_1_2_new[i] + w_2_new[i]^2 * y_2_new[i] * (-y_1_3_new[i])) / (2 * w_1_2_new[i] * y_1_new[i] * y_2_new[i] - 2 * (w_1_3_new[i] * y_1_new[i] + (-w_2_3_new[i]) * y_2_new[i]) * y_3_new[i])
            w_result <- w_result + w_delta_new[i]
            w_new <- c(w_new, w_result)

            # If y values are getting 0 after height adjustment, then w_new[i]=NaN
            if (is.nan(w_new[i])) {
                w_new[i] <- 0
            }

            # Calculate lambda for each peak triplet
            lambda_result <- -((sqrt(abs(((-w_2_new[i]^4 * y_2_new[i]^2 * y_1_3_new[i]^2 - w_1_new[i]^4 * y_1_new[i]^2 * y_2_3_new[i]^2 - w_3_new[i]^4 * y_1_2_new[i]^2 * y_3_new[i]^2 + 4 * w_2_new[i] * w_3_new[i]^3 * y_2_new[i] * ((-y_1_new[i]) + y_2_new[i]) * y_3_new[i]^2 + 4 * w_2_new[i]^3 * w_3_new[i] * y_2_new[i]^2 * y_3_new[i] * ((-y_1_new[i]) + y_3_new[i]) + 4 * w_1_new[i]^3 * y_1_new[i]^2 * y_2_3_new[i] * (w_2_new[i] * y_2_new[i] - w_3_new[i] * y_3_new[i]) + 4 * w_1_new[i] * y_1_new[i] * (w_2_new[i]^3 * y_2_new[i]^2 * y_1_3_new[i] - w_2_new[i] * w_3_new[i]^2 * y_2_new[i] * (y_1_new[i] + y_2_new[i] - 2 * y_3_new[i]) * y_3_new[i] + w_3_new[i]^3 * y_1_2_new[i] * y_3_new[i]^2 - w_2_new[i]^2 * w_3_new[i] * y_2_new[i] * y_3_new[i] * (y_1_new[i] - 2 * y_2_new[i] + y_3_new[i])) + 2 * w_2_new[i]^2 * w_3_new[i]^2 * y_2_new[i] * y_3_new[i] * (y_1_new[i]^2 - 3 * y_2_new[i] * y_3_new[i] + y_1_new[i] * (y_2_new[i] + y_3_new[i])) + 2 * w_1_new[i]^2 * y_1_new[i] * (-2 * w_2_new[i] * w_3_new[i] * y_2_new[i] * y_3_new[i] * (-2 * y_1_new[i] + y_2_new[i] + y_3_new[i]) + w_3_new[i]^2 * y_3_new[i] * (y_1_new[i] * (y_2_new[i] - 3 * y_3_new[i]) + y_2_new[i] * (y_2_new[i] + y_3_new[i])) + w_2_new[i]^2 * y_2_new[i] * (y_1_new[i] * (-3 * y_2_new[i] + y_3_new[i]) + y_3_new[i] * (y_2_new[i] + y_3_new[i]))))))))) / (2 * sqrt((w_1_new[i] * y_1_new[i] * y_2_3_new[i] + w_3_new[i] * y_1_2_new[i] * y_3_new[i] + w_2_new[i] * y_2_new[i] * ((-y_1_new[i]) + y_3_new[i]))^2))

            # If y and w are 0, then 0/0=NaN
            if (is.nan(lambda_result)) {
                lambda_result <- 0
            }
            lambda_new <- c(lambda_new, lambda_result)

            # Calculate scaling factor A for each peak triplet
            A_result <- (-4 * w_1_2_new[i] * w_1_3_new[i] * w_2_3_new[i] * y_1_new[i] * y_2_new[i] * y_3_new[i] * (w_1_new[i] * y_1_new[i] * y_2_3_new[i] + w_3_new[i] * y_3_new[i] * y_1_2_new[i] + w_2_new[i] * y_2_new[i] * (-y_1_3_new[i])) * lambda_new[i]) / (w_1_2_new[i]^4 * y_1_new[i]^2 * y_2_new[i]^2 - 2 * w_1_2_new[i]^2 * y_1_new[i] * y_2_new[i] * (w_1_3_new[i]^2 * y_1_new[i] + w_2_3_new[i]^2 * y_2_new[i]) * y_3_new[i] + (w_1_3_new[i]^2 * y_1_new[i] - w_2_3_new[i]^2 * y_2_new[i])^2 * y_3_new[i]^2)

            # If y and w are 0, then 0/0=NaN
            if (is.nan(A_result)) {
                A_result <- 0
            }
            A_new <- c(A_new, A_result)

            # Calculate new lorentz curves
            # If y values are zero, then lorentz curves should also be zero
            if ((w_new[i] == 0) | (lambda_new[i] == 0) | (A_new[i] == 0)) {
                lorentz_curves_initial[i, ] <- 0
            } else {
                lorentz_curves_initial[i, ] <- abs(A_new[i] * (lambda_new[i] / (lambda_new[i]^2 + (spec$sdp - w_new[i])^2)))
            }
        }

        # Calculate sum of lorentz curves
        spectrum_approx <- matrix(nrow = 1, ncol = length(spec$sdp))
        for (i in 1:length(spec$sdp)) {
            spectrum_approx[1, i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), i])
        }
        # ToSc: use vectorized functions, e.g.
        # spectrum_approx <- colSums(lorentz_curves_initial[1:length(filtered_peaks), ])

        # Standardize the spectra so that total area equals 1
        spectrum_y_normed <- spec$y / sum(spec$y)
        spectrum_approx_normed <- spectrum_approx / sum(spectrum_approx)

        # Calculate the difference between normed original spectrum and normed approximated spectrum
        difference_normed <- c()
        for (i in 1:length(spec$sdp)) {
            difference_normed[i] <- (spectrum_y_normed[i] - spectrum_approx_normed[i])^2
        }
        mse_normed <- (1 / length(difference_normed)) * sum(difference_normed)
        message(paste("\nNormed MSE value of iteration", b, "is: "))
        print(mse_normed)
    }
}

calculate_lorentz_curve_integrals <- function() {
    # Calculate the integrals for each lorentz curve
    integrals <- matrix(nrow = 1, ncol = length(lambda_new))
    for (i in 1:length(lambda_new)) {
        integrals[1, i] <- A_new[i] * (atan((-w_new[i] + (spec$length / sfx)) / lambda_new[i]) - atan((-w_new[i]) / lambda_new[i]))
    }

    # Save index of peak triplets
    index_peak_triplets_middle <- c()
    index_peak_triplets_left <- c()
    index_peak_triplets_right <- c()
    for (i in 1:length(filtered_peaks)) {
        index_peak_triplets_middle[i] <- filtered_peaks[i] + 1
        index_peak_triplets_left[i] <- filtered_left_position[i] + 1
        index_peak_triplets_right[i] <- filtered_right_position[i] + 1
    }

    # Save ppm x position of peak triplets
    peak_triplets_middle <- c()
    peak_triplets_left <- c()
    peak_triplets_right <- c()
    for (i in 1:length(filtered_peaks)) {
        peak_triplets_middle[i] <- spec$x_ppm[index_peak_triplets_middle[i]]
        peak_triplets_left[i] <- spec$x_ppm[index_peak_triplets_left[i]]
        peak_triplets_right[i] <- spec$x_ppm[index_peak_triplets_right[i]]
    }
}

init_debug_env <- function() {
    with(
        testdir = "deconvolute_spectrum/1",
        inputs = c(urine = "bruker/urine/urine_1"),
        output = "captured", message = "captured", plots = "plots.pdf",
        answers = c("y", "y"),
        expr = {
            set.seed(1234)
            deconvolute_spectrum(
                filepath = "urine_1/10",
                name = NULL,
                file_format = "bruker",
                same_parameter = FALSE,
                processing_value = 10,
                number_iterations = 1,
                range_water_signal_ppm = 0.1527692,
                signal_free_region = c(11.44494, -1.8828),
                smoothing_param = c(2, 5),
                delta = 6.4,
                scale_factor = c(1000, 1000000),
                current_filenumber = 1,
                number_of_files = 2
            )
        }
    )
}
