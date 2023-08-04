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
#' data_path=c("C:/example_data")
#' spectrum_data <- generate_lorentz_curves(
#'   data_path = data_path,
#'   file_format = "bruker",
#'   make_rds = FALSE)
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
  spectrum_data <- MetaboDecon1D(data_path, file_format = "bruker")
  if (make_rds) {
    saveRDS(object = spectrum_data, file = file.path(data_path, "spectrum_data.rds"))
  }
  return(spectrum_data)
}

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
#'   data_path = data_path,
#'   ppm_range = ppm_range,
#'   si_size_real_spectrum=si_size_real_spectrum,
#'   scale_factor_x=scale_factor
#' )
#' }
#' @details The output of this function is a data frame containig a matrix of
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
#' to perform the signal alinment across the individual spectra. For speaq
#' please cite: (Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S, Pieters L,
#' et al. (2018) speaq 2.0: A complete workflow for high-throughput 1D
#' NMRspectra processing and quantification. PLoS Comput Biol 14(3): e1006018.
#' https://doi.org/10.1371/journal.pcbi.1006018. The spectra deconvolution
#' process yields the signals of all spectra. Due to slight changes in
#' measurement conditions, e.g. pH variations, signal positions may vary
#' slightly across spectra. As a consequence, prior to further analysis signals
#' belonging to the same compound have to be aligned across spectra. This is
#' the purpose of the speaq package.
#' @param X From the data frame generated by the function gen_feat_mat we
#' require here the data_matix
#' @param peakList From the data frame generated by the function gen_feat_mat
#' we require her the peakList
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
#' signal alinment across the individual spectra. For speaq please cite:
#' (Beirnaert C, Meysman P, Vu TN, Hermans N, Apers S, Pieters L, et al. (2018)
#' speaq 2.0: A complete workflow for high-throughput 1D NMRspectra processing
#' and quantification. PLoS Comput Biol 14(3): e1006018.
#' https://doi.org/10.1371/journal.pcbi.1006018. The spectra deconvolution
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
#' after_speaq_mat<-speaq_align(feat, maxShift)
#' }
#' @details The output of speaq_align is a data matrix of aligned integral
#' values. Each row contains the data of each spectrum and each column
#' corresponds to one data point. Each entry corresponds to the integral of a
#' deconvoluted signal with the signal center at this specific positon after
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
#' @param range (positive integer) amount of adjecant columns which are
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
#'   after_speaq_mat,
#'   range,
#'   lower_bound,
#'   spectrum_data,
#'   data_path
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
  shifted_mat_no_na <-  replace(is.na(shifted_mat), 0)
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
