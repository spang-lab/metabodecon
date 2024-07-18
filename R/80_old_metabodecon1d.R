# Deconvolution #####

#' @export
#' @title Deconvolute 1D NMR spectrum
#' @description Automatic deconvolution of a 1D NMR spectrum into several Lorentz curves and the integration of them. The NMR file needs to be in Bruker format or jcamp-dx format.
#'
#' Deprecated since metabodecon v1.2.0. Please use [generate_lorentz_curves()] instead. See examples below for usage.
#'
#' `r lifecycle::badge("deprecated")`
#' @author Martina Haeckl, Tobias Schmidt
#' @param filepath Complete path of the file folder (Notice for Bruker format: filepath need to be the spectrum folder containing one or more different spectra (e.g."C:/Users/Username/Desktop/spectra_from_bruker"))
#' @param filename Name of the NMR file. (Notice for Bruker format: filename need to be the name of your spectrum which is also the name of the folder) (Default: filename = NA to analyze more spectra at once)
#' @param file_format Format (bruker or jcampdx) of the NMR file. (Default: file_format = "bruker")
#' @param number_iterations Number of iterations for the approximation of the parameters for the Lorentz curves (Default: number_iterations=10)
#' @param range_water_signal_ppm Half width of the water artefact in ppm (Default: range_water_signal=0.1527692 (e.g. for urine NMR spectra))
#' @param signal_free_region Row vector with two entries consisting of the ppm positions for the left and right border of the signal free region of the spectrum. (Default: signal_free_region=c(11.44494, -1.8828))
#' @param smoothing_param Row vector with two entries consisting of the number of smoothing repeats for the whole spectrum and the number of data points (uneven) for the mean calculation (Default: smoothing_param=c(2,5))
#' @param delta Defines the threshold value to distinguish between signal and noise (Default: delta=6.4)
#' @param scale_factor Row vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis (Default: scale_factor=c(1000,1000000))
#' @param debug Logical value to activate the debug mode (Default: debug=FALSE)
#' @return \loadmathjax
#' List containing
#' * `filename`: Name of the analyzed spectrum.
#' * `x_values`: Scaled datapoint numbers (SDP). Datapoints are numbered in descending order going from N to 0, where N equals the . Scaled data point numbers are obtained by dividing these numbers by the scale factor of the x-axis. I.e., for a spectrum with 131072 datapoints and a scale factor of 1000, the first scale datapoint has value 131.071 and the last one has value 0.
#' * `x_values_ppm`: The chemical shift of each datapoint in ppm (parts per million).
#' * `y_values`: The scaled signal intensity (SSI) of each datapoint. Obtained by reading the raw intensity values from the provided `data_path` as integers and dividing them scale factor of the y-axis.
#' * `spectrum_superposition`: Scaled signal intensity obtained by calculating the sum of all estimated Lorentz curves for each data point.
#' * `mse_normed`: Normalized mean squared error. Calculated as \mjeqn{\frac{1}{n} \sum_{i=1}^{n} (z_i - \hat{z}_i)^2}{1/n * sum((z_i - zhat_i)^2)} where \mjeqn{z_i}{z_i} is the normalized, smoothed signal intensity of data point i and \mjeqn{\hat{z}_i}{zhat_i} is the normalized superposition of Lorentz curves at data point i. Normalized in this context means that the vectors were scaled so the sum over all data points equals 1.
#' * `peak_triplets_middle`: Chemical shift of peak centers in ppm.
#' * `peak_triplets_left`: Chemical shift of left peak borders in ppm.
#' * `peak_triplets_right`: Chemical shift of right peak borders in ppm.
#' * `index_peak_triplets_middle`: Datapoint numbers of peak centers.
#' * `index_peak_triplets_left`: Datapoint numbers of left peak borders.
#' * `index_peak_triplets_right`: Datapoint numbers of right peak borders.
#' * `integrals`: Integrals of the Lorentz curves.
#' * `signal_free_region`: Borders of the signal free region of the spectrum in scaled datapoint numbers. Left of the first element and right of the second element no signals are expected.
#' * `range_water_signal_ppm`: Half width of the water signal in ppm. Potential signals in this region are ignored.
#' * `A`: Amplitude parameter of the Lorentz curves. Provided as negative number to maintain backwards compatibility with MetaboDecon1D. The area under the Lorentz curve is calculated as \mjeqn{A \cdot \pi}{A * pi}.
#' * `lambda`: Half width of the Lorentz curves in scaled data points. Provided as negative value to maintain backwards compatibility with MetaboDecon1D. Example: a value of -0.00525 corresponds to 5.25 data points. With a spectral width of 12019 Hz and 131072 data points this corresponds to a halfwidth at half height of approx. 0.48 Hz. The corresponding calculation is: (12019 Hz / 131071 dp) * 5.25 dp.
#' * `x_0`: Center of the Lorentz curves in scaled data points.
#'
#' \bold{Notice:} The parameters A, lambda and x_0 to calculate the Lorentz curves are saved in parameters.txt and the approximated spectrum is saved in approximated_spectrum.txt under the file path.
#' @details A list with parameters `A`, `lambda` and `x_0` required to calculate the Lorentz curves. The Lorentz curves could be calculated by using the function [calculate_lorentz_curves()]. This returns a matrix containing the generated and approximated Lorentz curves for each real peak of the spectrum. Each row of the matrix depicts one Lorentz curve. The Lorentz curves could be visualised and saved by using the function [plot_lorentz_curves_save_as_png()]. The superposition of all Lorentz curves, which reconstructs the original spectrum, could be visualised and saved with function [plot_spectrum_superposition_save_as_png()]. For the analytical calculation of the Lorentz curves peak triplets for each peak are used. To visualise these peak triplets and to illustrate the impact of the threshold `delta` the function `plot_triplets()` is available. The integral values for each generated Lorentz curves are saved in a vector.
#'
#' \bold{Notice}: It is feasible to load all spectra of a folder at once. Here the filename needs to be "NA" which is the default value. One selected spectrum could then be used to adjust the parameters (`signal_free_region` and `range_water_signal_ppm`) for the analysis of all spectra. Furthermore it is possible to adjust these parameters for each spectrum separate.
#' @import readJDX
#' @seealso [calculate_lorentz_curves()], [plot_triplets()], [plot_lorentz_curves_save_as_png()], [plot_spectrum_superposition_save_as_png()]
#' @references
#' Haeckl, M.; Tauber, P.; Schweda, F.; Zacharias, H.U.; Altenbuchinger, M.; Oefner, P.J.; Gronwald, W. An R-Package for the Deconvolution and Integration of 1D NMR Data: MetaboDecon1D. Metabolites 2021, 11, 452. https://doi.org/10.3390/metabo11070452
#' @examples
#' \donttest{
#' xds_path <- download_example_datasets()
#' urine <- file.path(xds_path, "bruker/urine")
#' urine_1 <- file.path(urine, "urine_1")
#'
#' \dontrun{
#' # Deprecated since metabodecon v1.2.0. Please use generate_lorentz_curves()
#' # instead. Shown below.
#' urine_1_deconv <- MetaboDecon1D(urine, "urine_1")
#' urine_all_deconv <- MetaboDecon1D(urine)
#' }
#'
#' urine_1_deconv <- generate_lorentz_curves(urine_1, ask = FALSE)
#' urine_all_deconv <- generate_lorentz_curves(urine, ask = FALSE)
#' }
MetaboDecon1D <- function(filepath,
                          filename = NA,
                          file_format = "bruker",
                          number_iterations = 10,
                          range_water_signal_ppm = 0.1527692,
                          signal_free_region = c(11.44494, -1.8828),
                          smoothing_param = c(2, 5),
                          delta = 6.4,
                          scale_factor = c(1000, 1000000),
                          debug = FALSE) {

  if (debug) logf("Starting MetaboDecon1D")
  store_results <- NULL

  example <- FALSE
  if (filepath == "load_example_path") {
    filepath <- system.file("extdata", package = "MetaboDecon1D", mustWork = TRUE)
    example <- TRUE
    owd <- setwd(filepath)
    on.exit(setwd(owd), add = TRUE)
  }

  # Check if filepath is a global file path (e.g. C:/) or local
  if (grepl("[ABCDEFGHIJKLMNOPQRSTUVWXYZ]+:+/", filepath)) {
    # filepath is a global file
    owd <- setwd(filepath)
    on.exit(setwd(owd), add = TRUE)
  } else {
    # Get current working directory and concat with current filepath and save as new working directory
    owd <- setwd(file.path(getwd(), filepath))
    on.exit(setwd(owd), add = TRUE)
    # Set afterwards filepath to global path
    filepath <- getwd()
  }

  # If filename is, as the default value, NA, then the all spectra of the filepath of the folder are analyzed
  if (is.na(filename)) {
    # Check which file_format is present
    if (file_format == "jcampdx") {
      # List files of filepath
      files_list <- list.files(filepath)

      # Save only file with .dx format
      files <- c()
      for (i in 1:length(files_list)) {
        if (endsWith(files_list[i], ".dx")) {
          files <- c(files, files_list[i])
        }
      }

      # Get number of files
      number_of_files <- length(files)

      # Ask User if he want to use same parameters for all spectra of the folder
      parameter_request <- readline(prompt = "Do you want to use the same parameters (signal_free_region, range_water_signal_ppm) for all spectra? (y/n) ")

      # Set parameter to TRUE or FALSE
      if (parameter_request == "y" | parameter_request == "n") {
        correct_input <- TRUE
      } else {
        correct_input <- FALSE
      }

      # Check if User input is correct or not
      while (correct_input == FALSE) {
        # Ask User if he want to use same parameters for all spectra of the folder
        message("Error. Please type only y or n.")
        parameter_request <- readline(prompt = "Do you want to use the same parameters (signal_free_region, range_water_signal_ppm) for all spectra? (y/n) ")

        if (parameter_request == "y" | parameter_request == "n") {
          correct_input <- TRUE
        } else {
          correct_input <- FALSE
        }
      }


      if (parameter_request == "y") {
        # Show User all files
        print(files)

        # Set variable to true
        same_parameter <- TRUE

        # Ask User which of the files should be used to adjust the parameters
        file_number <- readline(prompt = "Choose number of file which is used to adjust all parameters: [e.g. 1] ")

        number_of_files <- length(files)
        pattern <- paste("^[1-", number_of_files, "]", sep = "")

        # Check if input is a digit and smaller than number of files
        digit_true <- grepl(pattern, file_number)

        while (digit_true != TRUE) {
          # Ask User which of the files should be used to adjust the parameters
          message("Error. Please type only a digit which is smaller than number of available files.")
          file_number <- readline(prompt = "Choose number of file which is used to adjust all parameters: [e.g. 1] ")
          # Check if input is a digit
          digit_true <- grepl(pattern, file_number)
        }

        # Save as numeric
        file_number <- as.numeric(file_number)

        # Print which file the user selects
        message(paste("The selected file to adjust all parameters for all spectra is: ", files[file_number]))

        # Change order of files to analyze
        files_rearranged <- c()
        files_rearranged[1] <- files[file_number]
        files <- files[-file_number]
        files_rearranged <- c(files_rearranged, files)

        # Start deconvolution for each file
        return_list <- list()
        for (i in 1:length(files_rearranged)) {
          name <- files_rearranged[i]
          current_filenumber <- i

          print_text_1 <- "Start deconvolution of "
          print_text_2 <- ":"
          message(paste(print_text_1, files_rearranged[i], print_text_2, sep = ""))

          deconv_result <- deconvolution(filepath, name, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, current_filenumber, debug = debug, store_results = store_results)
          store_results <- deconv_result$store_results

          # Save return values in a list and return the list
          list_file <- list(
            "number_of_files" = number_of_files,
            "filename" = deconv_result$filename,
            "x_values" = deconv_result$spectrum_x,
            "x_values_ppm" = deconv_result$spectrum_x_ppm,
            "y_values" = deconv_result$spectrum_y,
            "spectrum_superposition" = deconv_result$spectrum_approx,
            "mse_normed" = deconv_result$mse_normed,
            "index_peak_triplets_middle" = deconv_result$index_peak_triplets_middle,
            "index_peak_triplets_left" = deconv_result$index_peak_triplets_left,
            "index_peak_triplets_right" = deconv_result$index_peak_triplets_right,
            "peak_triplets_middle" = deconv_result$peak_triplets_middle,
            "peak_triplets_left" = deconv_result$peak_triplets_left,
            "peak_triplets_right" = deconv_result$peak_triplets_right,
            "integrals" = deconv_result$integrals,
            "signal_free_region" = deconv_result$signal_free_region,
            "range_water_signal_ppm" = deconv_result$range_water_signal_ppm,
            "A" = deconv_result$A,
            "lambda" = deconv_result$lambda,
            "x_0" = deconv_result$w
          )
          # "lorentz_curves"=deconv_result$lorentz_curves,

          if (debug) {
            list_file$debuglist <- deconv_result$debuglist
          }

          # Save result list of current spectrum in a return list
          return_list[[paste0(files_rearranged[i])]] <- list_file

          # Save range_Water_signal and signal_free_region for next loop passage
          range_water_signal_ppm <- list_file$range_water_signal_ppm
          signal_free_region <- list_file$signal_free_region
        }
        return(return_list)
      }
      if (parameter_request == "n") {
        # User want to adjust parameters for each spectrum separately

        # Set variable to false
        same_parameter <- FALSE

        # Start deconvolution for each file
        return_list <- list()
        for (i in 1:length(files)) {
          name <- files[i]

          print_text_1 <- "Start deconvolution of "
          print_text_2 <- ":"
          message(paste(print_text_1, files[i], print_text_2, sep = ""))

          deconv_result <- deconvolution(filepath, name, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, debug = debug, store_results = store_results)
          store_results <- deconv_result$store_results

          # Save return values in a list and return the list
          list_file <- list("number_of_files" = number_of_files, "filename" = deconv_result$filename, "x_values" = deconv_result$spectrum_x, "x_values_ppm" = deconv_result$spectrum_x_ppm, "y_values" = deconv_result$spectrum_y, "spectrum_superposition" = deconv_result$spectrum_approx, "mse_normed" = deconv_result$mse_normed, "index_peak_triplets_middle" = deconv_result$index_peak_triplets_middle, "index_peak_triplets_left" = deconv_result$index_peak_triplets_left, "index_peak_triplets_right" = deconv_result$index_peak_triplets_right, "peak_triplets_middle" = deconv_result$peak_triplets_middle, "peak_triplets_left" = deconv_result$peak_triplets_left, "peak_triplets_right" = deconv_result$peak_triplets_right, "integrals" = deconv_result$integrals, "A" = deconv_result$A, "lambda" = deconv_result$lambda, "x_0" = deconv_result$w)
          # "lorentz_curves"=deconv_result$lorentz_curves,

          if (debug) {
            list_file$debuglist <- deconv_result$debuglist
          }

          # Save result list of current spectrum in a return list
          # return_list[[paste0("element", i)]] <- list_file
          return_list[[paste0(files[i])]] <- list_file
        }
        return(return_list)
      }
    }

    # Check which file_format is present
    if (file_format == "bruker") {
      # List files of filepath
      files_list <- list.files(filepath)

      # Check if files are only folders
      check_files <- dir.exists(files_list)

      # Delete file if it is not a folder
      files <- c()
      for (i in 1:length(check_files)) {
        if (check_files[i] == TRUE) {
          files <- c(files, files_list[i])
        }
      }

      # Get number of files
      number_of_files <- length(files)

      # If example is loaded use predefined values, else get some values from the user
      if (example == TRUE) {
        spectroscopy_value <- 10
        processing_value <- 10
      } else {
        spectroscopy_value <- readline(prompt = "What is the name of the subfolder of your filepath: \n[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10] ")
        processing_value <- readline(prompt = "What is the name of the subsubsubfolder of your filepath: \n[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10/pdata/10] ")
      }

      # Ask User if he want to use same parameters for all spectra of the folder
      parameter_request <- readline(prompt = "Do you want to use the same parameters (signal_free_region, range_water_signal_ppm) for all spectra? (y/n) ")

      # Set parameter to TRUE or FALSE
      if (parameter_request == "y" | parameter_request == "n") {
        correct_input <- TRUE
      } else {
        correct_input <- FALSE
      }

      # Check if User input is correct or not
      while (correct_input == FALSE) {
        # Ask User if he want to use same parameters for all spectra of the folder
        message("Error. Please type only y or n.")
        parameter_request <- readline(prompt = "Do you want to use the same parameters (signal_free_region, range_water_signal_ppm) for all spectra? (y/n) ")

        if (parameter_request == "y" | parameter_request == "n") {
          correct_input <- TRUE
        } else {
          correct_input <- FALSE
        }
      }



      if (parameter_request == "y") {
        # Show User all files
        print(files)

        # Set variable to true
        same_parameter <- TRUE

        # Ask User which of the files should be used to adjust the parameters
        file_number <- readline(prompt = "Choose number of file which is used to adjust all parameters: [e.g. 1] ")

        number_of_files <- length(files)
        pattern <- paste("^[1-", number_of_files, "]", sep = "")

        # Check if input is a digit and smaller than number of files
        digit_true <- grepl(pattern, file_number)

        while (digit_true != TRUE) {
          # Ask User which of the files should be used to adjust the parameters
          message("Error. Please type only a digit which is smaller than number of available files.")
          file_number <- readline(prompt = "Choose number of file which is used to adjust all parameters: [e.g. 1] ")
          # Check if input is a digit
          digit_true <- grepl(pattern, file_number)
        }

        # Save as numeric
        file_number <- as.numeric(file_number)


        # Print which file the user selects
        message(paste("The selected file to adjust all parameters for all spectra is: ", files[file_number]))

        # Change order of files to analyze
        files_rearranged <- c()
        files_rearranged[1] <- files[file_number]
        files <- files[-file_number]
        files_rearranged <- c(files_rearranged, files)


        # Start deconvolution for each file
        return_list <- list()
        for (i in 1:length(files_rearranged)) {
          name <- files_rearranged[i]
          current_filenumber <- i
          # Generate whole filepath of current folder
          filepath_completed <- paste(filepath, files_rearranged[i], spectroscopy_value, sep = "/")

          print_text_1 <- "Start deconvolution of "
          print_text_2 <- ":"
          message(paste(print_text_1, files_rearranged[i], print_text_2, sep = ""))

          deconv_result <- deconvolution(filepath_completed, name, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, current_filenumber, debug = debug, store_results = store_results)
          store_results <- deconv_result$store_results


          # Save return values in a list and return the list
          list_file <- list("number_of_files" = number_of_files, "filename" = deconv_result$filename, "x_values" = deconv_result$spectrum_x, "x_values_ppm" = deconv_result$spectrum_x_ppm, "y_values" = deconv_result$spectrum_y, "spectrum_superposition" = deconv_result$spectrum_approx, "mse_normed" = deconv_result$mse_normed, "index_peak_triplets_middle" = deconv_result$index_peak_triplets_middle, "index_peak_triplets_left" = deconv_result$index_peak_triplets_left, "index_peak_triplets_right" = deconv_result$index_peak_triplets_right, "peak_triplets_middle" = deconv_result$peak_triplets_middle, "peak_triplets_left" = deconv_result$peak_triplets_left, "peak_triplets_right" = deconv_result$peak_triplets_right, "integrals" = deconv_result$integrals, "signal_free_region" = deconv_result$signal_free_region, "range_water_signal_ppm" = deconv_result$range_water_signal_ppm, "A" = deconv_result$A, "lambda" = deconv_result$lambda, "x_0" = deconv_result$w)
          # "lorentz_curves"=deconv_result$lorentz_curves,

          if (debug) {
            list_file$debuglist <- deconv_result$debuglist
          }

          # Save result list of current spectrum in a return list
          return_list[[paste0(files_rearranged[i])]] <- list_file

          # Save range_Water_signal and signal_free_region for next loop passage
          range_water_signal_ppm <- list_file$range_water_signal_ppm
          signal_free_region <- list_file$signal_free_region
        }

        if (debug) {
          logf("Finished MetaboDecon1D")
        }

        return(return_list)
      }
      if (parameter_request == "n") {
        # User want to adjust parameters for each spectrum separately

        # Set variable to false
        same_parameter <- FALSE

        # Start deconvolution for each file
        return_list <- list()
        for (i in 1:length(files)) {
          name <- files[i]
          # Generate whole filepath of current folder
          filepath_completed <- paste(filepath, files[i], spectroscopy_value, sep = "/")

          print_text_1 <- "Start deconvolution of "
          print_text_2 <- ":"
          message(paste(print_text_1, files[i], print_text_2, sep = ""))

          deconv_result <- deconvolution(filepath_completed, name, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, debug = debug, store_results = store_results)
          store_results <- deconv_result$store_results


          # Save return values in a list and return the list
          list_file <- list("number_of_files" = number_of_files, "filename" = deconv_result$filename, "x_values" = deconv_result$spectrum_x, "x_values_ppm" = deconv_result$spectrum_x_ppm, "y_values" = deconv_result$spectrum_y, "spectrum_superposition" = deconv_result$spectrum_approx, "mse_normed" = deconv_result$mse_normed, "index_peak_triplets_middle" = deconv_result$index_peak_triplets_middle, "index_peak_triplets_left" = deconv_result$index_peak_triplets_left, "index_peak_triplets_right" = deconv_result$index_peak_triplets_right, "peak_triplets_middle" = deconv_result$peak_triplets_middle, "peak_triplets_left" = deconv_result$peak_triplets_left, "peak_triplets_right" = deconv_result$peak_triplets_right, "integrals" = deconv_result$integrals, "A" = deconv_result$A, "lambda" = deconv_result$lambda, "x_0" = deconv_result$w)
          # "lorentz_curves"=deconv_result$lorentz_curves,

          if (debug) {
            list_file$debuglist <- deconv_result$debuglist
          }

          # Save result list of current spectrum in a return list
          # return_list[[paste0("element", i)]] <- list_file
          return_list[[paste0(files[i])]] <- list_file
        }

        if (debug) {
          logf("Finished MetaboDecon1D")
        }

        return(return_list)
      }
    }



    # If a filename is given, thus unequal NA, only this file is analyzed
  } else {
    # Call deconvolution function

    # Set variable to false
    same_parameter <- FALSE

    # Get number of files
    number_of_files <- 1

    print_text_1 <- "Start deconvolution of "
    print_text_2 <- ":"
    message(paste(print_text_1, filename, print_text_2, sep = ""))

    # Check which file format is loaded
    if (file_format == "jcampdx") {
      deconv_result <- deconvolution(filepath, filename, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, debug = debug, store_results = store_results)
      store_results <- deconv_result$store_results


      # Save return values in a list and return the list
      return_list <- list("number_of_files" = number_of_files, "filename" = deconv_result$filename, "x_values" = deconv_result$spectrum_x, "x_values_ppm" = deconv_result$spectrum_x_ppm, "y_values" = deconv_result$spectrum_y, "spectrum_superposition" = deconv_result$spectrum_approx, "mse_normed" = deconv_result$mse_normed, "index_peak_triplets_middle" = deconv_result$index_peak_triplets_middle, "index_peak_triplets_left" = deconv_result$index_peak_triplets_left, "index_peak_triplets_right" = deconv_result$index_peak_triplets_right, "peak_triplets_middle" = deconv_result$peak_triplets_middle, "peak_triplets_left" = deconv_result$peak_triplets_left, "peak_triplets_right" = deconv_result$peak_triplets_right, "integrals" = deconv_result$integrals, "A" = deconv_result$A, "lambda" = deconv_result$lambda, "x_0" = deconv_result$w)
      # "lorentz_curves"=deconv_result$lorentz_curves,

      if (debug) {
        return_list$debuglist <- deconv_result$debuglist
        logf("Finished MetaboDecon1D")
      }

      return(return_list)
    }

    # Check which file format is loaded
    if (file_format == "bruker") {
      # If example is loaded use predefined values, else get some values from the user
      if (example == TRUE) {
        spectroscopy_value <- 10
        processing_value <- 10
      } else {
        spectroscopy_value <- readline(prompt = "What is the name of the subfolder of your filepath: \n[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10] ")
        processing_value <- readline(prompt = "What is the name of the subsubsubfolder of your filepath: \n[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10/pdata/10] ")
      }

      # Set variable to false
      same_parameter <- FALSE

      # Generate whole filepath of current folder
      filepath_completed <- paste(filepath, filename, spectroscopy_value, sep = "/")
      deconv_result <- deconvolution(filepath_completed, filename, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, debug = debug, store_results = store_results)
      store_results <- deconv_result$store_results

      # Save return values in a list and return the list
      return_list <- list("number_of_files" = number_of_files, "filename" = deconv_result$filename, "x_values" = deconv_result$spectrum_x, "x_values_ppm" = deconv_result$spectrum_x_ppm, "y_values" = deconv_result$spectrum_y, "spectrum_superposition" = deconv_result$spectrum_approx, "mse_normed" = deconv_result$mse_normed, "index_peak_triplets_middle" = deconv_result$index_peak_triplets_middle, "index_peak_triplets_left" = deconv_result$index_peak_triplets_left, "index_peak_triplets_right" = deconv_result$index_peak_triplets_right, "peak_triplets_middle" = deconv_result$peak_triplets_middle, "peak_triplets_left" = deconv_result$peak_triplets_left, "peak_triplets_right" = deconv_result$peak_triplets_right, "integrals" = deconv_result$integrals, "A" = deconv_result$A, "lambda" = deconv_result$lambda, "x_0" = deconv_result$w)

      if (debug) {
        return_list$debuglist <- deconv_result$debuglist
        logf("Finished MetaboDecon1D")
      }

      return(return_list)
    }
  }
}

#' @noRd
#' @details For examples on how to call see `tests/testthat/test-deconvolution_v11-123.R`
#' @author Martina Haeckl
deconvolution <- function(filepath,
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
                          debug = FALSE,
                          store_results = NULL) {
  if (debug) {
    logf("Starting deconvolution")
    args <- list(
      filepath = filepath,
      name = name,
      file_format = file_format,
      same_parameter = same_parameter,
      processing_value = tryCatch(processing_value, error = function(e) "missing"),
      number_iterations = number_iterations,
      range_water_signal_ppm = range_water_signal_ppm,
      signal_free_region = signal_free_region,
      smoothing_param = smoothing_param,
      delta = delta,
      scale_factor = scale_factor,
      current_filenumber = tryCatch(current_filenumber, error = function(e) "missing")
    )
    debuglist <- list(args = args)
  }

  # If loaded name is in JCAMP-DX format
  if (file_format == "jcampdx") {

    # Import data and install necessary package readJDX automatically
    data <- readJDX::readJDX(file = name, SOFC = TRUE, debug = 0)

    # Choose factor to scale axis values
    factor_x <- scale_factor[1]
    factor_y <- scale_factor[2]

    # Save data
    spectrum_length <- length(data[[4]]$x) - 1
    spectrum_x <- seq((spectrum_length / factor_x), 0, -1 / factor_x)
    spectrum_y <- (data[[4]]$y) / factor_y
    if (debug) {
      logf("Read raw data from %s", name)
      debuglist$data_read <- list(
        spectrum_y_raw = (data[[4]]$y),
        spectrum_y = spectrum_y
      )
    }

    # Calculate ppm x-axis
    # Get values of ppm range from data
    # Get spectral width
    for (i in 1:length(data[[2]])) {
      if (startsWith(data[[2]][i], "##$SW=")) {
        ppm_range_index <- i
      }
    }
    ppm_range <- as.numeric(sub("\\D+", "", data[[2]][ppm_range_index]))
    for (j in 1:length(data[[2]])) {
      if (startsWith(data[[2]][j], "##$OFFSET=")) {
        ppm_highest_index <- j
      }
    }
    ppm_highest_value <- as.numeric(sub("\\D+", "", data[[2]][ppm_highest_index]))
    ppm_lowest_value <- ppm_highest_value - ppm_range
    spectrum_x_ppm <- seq(ppm_highest_value, ppm_lowest_value, (-ppm_range / spectrum_length))
  }

  # If loaded name is in Bruker format
  if (file_format == "bruker") {
    # Get file paths of all necessary documents
    acqus_file <- file.path(filepath, "acqus")
    file_folders_spec <- paste("pdata", processing_value, "1r", sep = "/")
    spec_file <- file.path(filepath, file_folders_spec)
    file_folders_procs <- paste("pdata", processing_value, "procs", sep = "/")
    procs_file <- file.path(filepath, file_folders_procs)

    # Read content of files
    acqs_content <- readLines(acqus_file)
    procs_content <- readLines(procs_file)

    # Load necessary meta data of acqa content
    for (i in 1:length(acqs_content)) {
      if (startsWith(acqs_content[i], "##$SW=")) {
        index_spectral_width <- i
      }
    }

    # Save digit value
    ppm_range <- as.numeric(sub("\\D+", "", acqs_content[index_spectral_width]))

    # Load necessary meta data of procs content
    for (i in 1:length(procs_content)) {
      if (startsWith(procs_content[i], "##$DTYPP=")) {
        index_integer_type <- i
      }
      if (startsWith(procs_content[i], "##$OFFSET=")) {
        index_highest_ppm <- i
      }
      if (startsWith(procs_content[i], "##$SI=")) {
        index_data_points <- i
      }
    }
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
    factor_x <- scale_factor[1]
    factor_y <- scale_factor[2]

    # Calculate length of spectrum
    spectrum_length <- length(spectrum_y)

    # Calculate x axis
    spectrum_x_ppm <- seq(ppm_highest_value, ppm_highest_value - ppm_range, by = -(ppm_range / (spectrum_length - 1)))

    # Calculate spectrum_x and recalculate spectum_y
    spectrum_x <- seq((spectrum_length - 1) / factor_x, 0, -0.001)
    if (debug) {
      logf("Read raw data from %s", spec_file)
      debuglist$data_read <- list(
        spectrum_y_raw = spectrum_y,
        spectrum_y = spectrum_y / factor_y
      )
    }
    spectrum_y <- spectrum_y / factor_y
  }

  # Check if parameters are the same for all analyzed spectra
  if (same_parameter == FALSE) {
    # Calculate signal free region
    signal_free_region_left <- (spectrum_length + 1) - ((ppm_highest_value - signal_free_region[1]) / (ppm_range / spectrum_length))
    signal_free_region_right <- (spectrum_length + 1) - ((ppm_highest_value - signal_free_region[2]) / (ppm_range / spectrum_length))

    signal_free_region_left <- signal_free_region_left / factor_x
    signal_free_region_right <- signal_free_region_right / factor_x

    plot(spectrum_x_ppm, spectrum_y, type = "l", xlab = "[ppm]", ylab = "Intensity [a.u.]", xlim = rev(range(spectrum_x_ppm)))
    graphics::abline(v = signal_free_region[1], col = "green")
    graphics::abline(v = signal_free_region[2], col = "green")


    # Check for correct range of signal free region
    check_range_signal_free_region <- readline(prompt = "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): ")

    # Set parameter to TRUE or FALSE
    if (check_range_signal_free_region == "y" | check_range_signal_free_region == "n") {
      correct_input <- TRUE
    } else {
      correct_input <- FALSE
    }

    # Check if User input is correct or not
    while (correct_input == FALSE) {
      # Ask User if he want to use same parameters for all spectra of the folder
      message("Error. Please type only y or n.")
      check_range_signal_free_region <- readline(prompt = "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): ")

      if (check_range_signal_free_region == "y" | check_range_signal_free_region == "n") {
        correct_input <- TRUE
      } else {
        correct_input <- FALSE
      }
    }


    while (check_range_signal_free_region == "n") {
      signal_free_region_left_ppm <- readline(prompt = "Choose another left border: [e.g. 12] ")

      # Check if input is a digit
      digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", signal_free_region_left_ppm)

      while (digit_true != TRUE) {
        # Ask User which of the files should be used to adjust the parameters
        message("Error. Please only type a digit.")
        signal_free_region_left_ppm <- readline(prompt = "Choose another left border: [e.g. 12] ")
        # Check if input is a digit
        digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", signal_free_region_left_ppm)
      }
      # Save as numeric
      signal_free_region_left_ppm <- as.numeric(signal_free_region_left_ppm)
      signal_free_region_left <- ((spectrum_length + 1) - ((ppm_highest_value - signal_free_region_left_ppm) / (ppm_range / spectrum_length))) / factor_x

      signal_free_region_right_ppm <- readline(prompt = "Choose another right border: [e.g. -2] ")

      # Check if input is a digit
      digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", signal_free_region_right_ppm)

      while (digit_true != TRUE) {
        # Ask User which of the files should be used to adjust the parameters
        message("Error. Please only type a digit.")
        signal_free_region_right_ppm <- readline(prompt = "Choose another right border: [e.g. -2] ")
        # Check if input is a digit
        digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", signal_free_region_right_ppm)
      }
      # Save as numeric
      signal_free_region_right_ppm <- as.numeric(signal_free_region_right_ppm)
      signal_free_region_right <- ((spectrum_length + 1) - ((ppm_highest_value - signal_free_region_right_ppm) / (ppm_range / spectrum_length))) / factor_x

      plot(spectrum_x_ppm, spectrum_y, type = "l", xlab = "[ppm]", ylab = "Intensity [a.u.]", xlim = rev(range(spectrum_x_ppm)))
      graphics::abline(v = signal_free_region_left_ppm, col = "green")
      graphics::abline(v = signal_free_region_right_ppm, col = "green")

      check_range_signal_free_region <- readline(prompt = "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): ")

      # Set parameter to TRUE or FALSE
      if (check_range_signal_free_region == "y" | check_range_signal_free_region == "n") {
        correct_input <- TRUE
      } else {
        correct_input <- FALSE
      }

      # Check if User input is correct or not
      while (correct_input == FALSE) {
        # Ask User if he want to use same parameters for all spectra of the folder
        message("Error. Please type only y or n.")
        check_range_signal_free_region <- readline(prompt = "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): ")

        if (check_range_signal_free_region == "y" | check_range_signal_free_region == "n") {
          correct_input <- TRUE
        } else {
          correct_input <- FALSE
        }
      }
    }

    # Remove water signal
    water_signal_position <- length(spectrum_x) / 2
    water_signal_position_ppm <- spectrum_x_ppm[length(spectrum_x_ppm) / 2]

    # Recalculate ppm into data points
    range_water_signal <- range_water_signal_ppm / (ppm_range / spectrum_length)
    water_signal_left <- water_signal_position - range_water_signal
    water_signal_right <- water_signal_position + range_water_signal

    plot(spectrum_x_ppm, spectrum_y, type = "l", xlab = "[ppm]", ylab = "Intensity [a.u.]", xlim = rev(range((water_signal_position_ppm - 2 * range_water_signal_ppm), (water_signal_position_ppm + 2 * range_water_signal_ppm))))
    graphics::abline(v = spectrum_x_ppm[water_signal_left], col = "red")
    graphics::abline(v = spectrum_x_ppm[water_signal_right], col = "red")

    # Check for correct range of water artefact
    check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

    # Set parameter to TRUE or FALSE
    if (check_range_water_signal == "y" | check_range_water_signal == "n") {
      correct_input <- TRUE
    } else {
      correct_input <- FALSE
    }

    # Check if User input is correct or not
    while (correct_input == FALSE) {
      # Ask User if he want to use same parameters for all spectra of the folder
      message("Error. Please type only y or n.")
      check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

      if (check_range_water_signal == "y" | check_range_water_signal == "n") {
        correct_input <- TRUE
      } else {
        correct_input <- FALSE
      }
    }


    while (check_range_water_signal == "n") {
      range_water_signal_ppm <- readline(prompt = "Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154] ")

      # Check if input is a digit
      digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", range_water_signal_ppm)

      while (digit_true != TRUE) {
        # Ask User which of the files should be used to adjust the parameters
        message("Error. Please only type a digit.")
        range_water_signal_ppm <- readline(prompt = "Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154] ")
        # Check if input is a digit
        digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", range_water_signal_ppm)
      }
      # Save as numeric
      range_water_signal_ppm <- as.numeric(range_water_signal_ppm)



      # Remove water signal
      water_signal_position <- length(spectrum_x) / 2
      water_signal_position_ppm <- spectrum_x_ppm[length(spectrum_x_ppm) / 2]
      # Recalculate ppm into data points
      range_water_signal <- range_water_signal_ppm / (ppm_range / spectrum_length)
      water_signal_left <- water_signal_position - range_water_signal
      water_signal_right <- water_signal_position + range_water_signal

      plot(spectrum_x_ppm, spectrum_y, type = "l", xlab = "[ppm]", ylab = "Intensity [a.u.]", xlim = rev(range((water_signal_position_ppm - 2 * range_water_signal_ppm), (water_signal_position_ppm + 2 * range_water_signal_ppm))))
      graphics::abline(v = spectrum_x_ppm[water_signal_left], col = "red")
      graphics::abline(v = spectrum_x_ppm[water_signal_right], col = "red")

      check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

      # Set parameter to TRUE or FALSE
      if (check_range_water_signal == "y" | check_range_water_signal == "n") {
        correct_input <- TRUE
      } else {
        correct_input <- FALSE
      }

      # Check if User input is correct or not
      while (correct_input == FALSE) {
        # Ask User if he want to use same parameters for all spectra of the folder
        message("Error. Please type only y or n.")
        check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

        if (check_range_water_signal == "y" | check_range_water_signal == "n") {
          correct_input <- TRUE
        } else {
          correct_input <- FALSE
        }
      }
    }
  }

  # Check if parameters are the same for all analyzed spectra
  if (same_parameter == TRUE) {
    # Check if current file is the first file
    # If yes, this file is used to adjust the parameters for all spectra
    if (current_filenumber == 1) {
      # Calculate signal free region
      signal_free_region_left <- (spectrum_length + 1) - ((ppm_highest_value - signal_free_region[1]) / (ppm_range / spectrum_length))
      signal_free_region_right <- (spectrum_length + 1) - ((ppm_highest_value - signal_free_region[2]) / (ppm_range / spectrum_length))

      signal_free_region_left <- signal_free_region_left / factor_x
      signal_free_region_right <- signal_free_region_right / factor_x

      plot(spectrum_x_ppm, spectrum_y, type = "l", xlab = "[ppm]", ylab = "Intensity [a.u.]", xlim = rev(range(spectrum_x_ppm)))
      graphics::abline(v = signal_free_region[1], col = "green")
      graphics::abline(v = signal_free_region[2], col = "green")

      # Check for correct range of signal free region
      check_range_signal_free_region <- readline(prompt = "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): ")

      # Set parameter to TRUE or FALSE
      if (check_range_signal_free_region == "y" | check_range_signal_free_region == "n") {
        correct_input <- TRUE
      } else {
        correct_input <- FALSE
      }

      # Check if User input is correct or not
      while (correct_input == FALSE) {
        # Ask User if he want to use same parameters for all spectra of the folder
        message("Error. Please type only y or n.")
        check_range_signal_free_region <- readline(prompt = "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): ")

        if (check_range_signal_free_region == "y" | check_range_signal_free_region == "n") {
          correct_input <- TRUE
        } else {
          correct_input <- FALSE
        }
      }

      while (check_range_signal_free_region == "n") {
        signal_free_region_left_ppm <- readline(prompt = "Choose another left border: [e.g. 12] ")
        # Check if input is a digit
        digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", signal_free_region_left_ppm)

        while (digit_true != TRUE) {
          # Ask User which of the files should be used to adjust the parameters
          message("Error. Please only type a digit.")
          signal_free_region_left_ppm <- readline(prompt = "Choose another left border: [e.g. 12] ")
          # Check if input is a digit
          digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", signal_free_region_left_ppm)
        }
        # Save as numeric
        signal_free_region_left_ppm <- as.numeric(signal_free_region_left_ppm)
        signal_free_region_left <- ((spectrum_length + 1) - ((ppm_highest_value - signal_free_region_left_ppm) / (ppm_range / spectrum_length))) / factor_x

        signal_free_region_right_ppm <- readline(prompt = "Choose another right border: [e.g. -2] ")
        # Check if input is a digit
        digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", signal_free_region_right_ppm)

        while (digit_true != TRUE) {
          # Ask User which of the files should be used to adjust the parameters
          message("Error. Please only type a digit.")
          signal_free_region_right_ppm <- readline(prompt = "Choose another right border: [e.g. -2] ")
          # Check if input is a digit
          digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", signal_free_region_right_ppm)
        }
        # Save as numeric
        signal_free_region_right_ppm <- as.numeric(signal_free_region_right_ppm)
        signal_free_region_right <- ((spectrum_length + 1) - ((ppm_highest_value - signal_free_region_right_ppm) / (ppm_range / spectrum_length))) / factor_x

        plot(spectrum_x_ppm, spectrum_y, type = "l", xlab = "[ppm]", ylab = "Intensity [a.u.]", xlim = rev(range(spectrum_x_ppm)))
        graphics::abline(v = signal_free_region_left_ppm, col = "green")
        graphics::abline(v = signal_free_region_right_ppm, col = "green")

        check_range_signal_free_region <- readline(prompt = "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): ")
        # Set parameter to TRUE or FALSE
        if (check_range_signal_free_region == "y" | check_range_signal_free_region == "n") {
          correct_input <- TRUE
        } else {
          correct_input <- FALSE
        }

        # Check if User input is correct or not
        while (correct_input == FALSE) {
          # Ask User if he want to use same parameters for all spectra of the folder
          message("Error. Please type only y or n.")
          check_range_signal_free_region <- readline(prompt = "Signal free region borders correct selected? (Area left and right of the green lines) (y/n): ")

          if (check_range_signal_free_region == "y" | check_range_signal_free_region == "n") {
            correct_input <- TRUE
          } else {
            correct_input <- FALSE
          }
        }
      }

      # Save adjusted signal_free_region
      signal_free_region <- c(signal_free_region_left, signal_free_region_right)


      # Remove water signal
      water_signal_position <- length(spectrum_x) / 2
      water_signal_position_ppm <- spectrum_x_ppm[length(spectrum_x_ppm) / 2]
      # Recalculate ppm into data points
      range_water_signal <- range_water_signal_ppm / (ppm_range / spectrum_length)
      water_signal_left <- water_signal_position - range_water_signal
      water_signal_right <- water_signal_position + range_water_signal

      plot(spectrum_x_ppm, spectrum_y, type = "l", xlab = "[ppm]", ylab = "Intensity [a.u.]", xlim = rev(range((water_signal_position_ppm - 2 * range_water_signal_ppm), (water_signal_position_ppm + 2 * range_water_signal_ppm))))
      graphics::abline(v = spectrum_x_ppm[water_signal_left], col = "red")
      graphics::abline(v = spectrum_x_ppm[water_signal_right], col = "red")


      # Check for correct range of water artefact
      check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

      # Set parameter to TRUE or FALSE
      if (check_range_water_signal == "y" | check_range_water_signal == "n") {
        correct_input <- TRUE
      } else {
        correct_input <- FALSE
      }

      # Check if User input is correct or not
      while (correct_input == FALSE) {
        # Ask User if he want to use same parameters for all spectra of the folder
        message("Error. Please type only y or n.")
        check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

        if (check_range_water_signal == "y" | check_range_water_signal == "n") {
          correct_input <- TRUE
        } else {
          correct_input <- FALSE
        }
      }

      while (check_range_water_signal == "n") {
        range_water_signal_ppm <- readline(prompt = "Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154] ")

        # Check if input is a digit
        digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", range_water_signal_ppm)

        while (digit_true != TRUE) {
          # Ask User which of the files should be used to adjust the parameters
          message("Error. Please only type a digit.")
          range_water_signal_ppm <- readline(prompt = "Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154] ")
          # Check if input is a digit
          digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", range_water_signal_ppm)
        }
        # Save as numeric
        range_water_signal_ppm <- as.numeric(range_water_signal_ppm)


        # Remove water signal
        water_signal_position <- length(spectrum_x) / 2
        water_signal_position_ppm <- spectrum_x_ppm[length(spectrum_x_ppm) / 2]
        # Recalculate ppm into data points
        range_water_signal <- range_water_signal_ppm / (ppm_range / spectrum_length)
        water_signal_left <- water_signal_position - range_water_signal
        water_signal_right <- water_signal_position + range_water_signal

        plot(spectrum_x_ppm, spectrum_y, type = "l", xlab = "[ppm]", ylab = "Intensity [a.u.]", xlim = rev(range((water_signal_position_ppm - 2 * range_water_signal_ppm), (water_signal_position_ppm + 2 * range_water_signal_ppm))))
        graphics::abline(v = spectrum_x_ppm[water_signal_left], col = "red")
        graphics::abline(v = spectrum_x_ppm[water_signal_right], col = "red")

        check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")
        # Set parameter to TRUE or FALSE
        if (check_range_water_signal == "y" | check_range_water_signal == "n") {
          correct_input <- TRUE
        } else {
          correct_input <- FALSE
        }

        # Check if User input is correct or not
        while (correct_input == FALSE) {
          # Ask User if he want to use same parameters for all spectra of the folder
          message("Error. Please type only y or n.")
          check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

          if (check_range_water_signal == "y" | check_range_water_signal == "n") {
            correct_input <- TRUE
          } else {
            correct_input <- FALSE
          }
        }
      }

      # Save adjusted range_water_signal
      range_water_signal_ppm <- range_water_signal_ppm
    } else {
      # If current file is not the first file, parameters are already adjusted and only needs to be loaded
      signal_free_region_left <- signal_free_region[1]
      signal_free_region_right <- signal_free_region[2]

      # Recalculate ppm into data points
      water_signal_position <- length(spectrum_x) / 2
      water_signal_position_ppm <- spectrum_x_ppm[length(spectrum_x_ppm) / 2]
      range_water_signal <- range_water_signal_ppm / (ppm_range / spectrum_length)
      water_signal_left <- water_signal_position - range_water_signal
      water_signal_right <- water_signal_position + range_water_signal
    }
  }

  # Remove water signal
  for (i in water_signal_right:water_signal_left) {
    spectrum_y[i] <- 0.00000001
  }

  if (debug) {
    logf("Removed water signal")
    debuglist$ws_rm <- list(
      spectrum_length = spectrum_length,
      spectrum_x = spectrum_x,
      spectrum_x_ppm = spectrum_x_ppm,
      spectrum_y = spectrum_y,
      signal_free_region_left = signal_free_region_left,
      signal_free_region_right = signal_free_region_right,
      water_signal_position = water_signal_position,
      water_signal_position_ppm = water_signal_position_ppm,
      range_water_signal = range_water_signal,
      water_signal_left = water_signal_left,
      water_signal_right = water_signal_right
    )
  }

  # Remove negative values of spectrum by Saving the absolut values
  for (i in 1:length(spectrum_y)) {
    spectrum_y[i] <- abs(spectrum_y[i])
  }
  if (debug) {
    logf("Removed negative values")
    debuglist$neg_rm <- list(spectrum_y = spectrum_y)
  }

  # Variable Mean Filter
  smoothing_iteration <- smoothing_param[1]
  smoothing_pairs <- smoothing_param[2]

  # Check if number of smoothing pairs is uneven
  while (smoothing_pairs %% 2 == "0") {
    smoothing_pairs <- as.numeric(readline(prompt = "Number of smoothing pairs is even. Please choose uneven number: "))
  }

  for (j in 1:smoothing_iteration) {
    smoothed_spectrum_y <- c()
    for (i in 1:(length(spectrum_x))) {
      # Calculate borders
      left_border <- i - floor(smoothing_pairs / 2)
      right_border <- i + floor(smoothing_pairs / 2)

      # Calculate smoothed spectrum for borders
      if (left_border <= 0) {
        left_border <- 1
        smoothed_spectrum_y[i] <- (1 / right_border) * sum(spectrum_y[left_border:right_border])
      } else if (right_border >= length(spectrum_x)) {
        right_border <- length(spectrum_x)
        smoothed_spectrum_y[i] <- (1 / (right_border - left_border + 1)) * sum(spectrum_y[left_border:right_border])
      } else {
        # Calculate smoothed spectrum
        smoothed_spectrum_y[i] <- (1 / smoothing_pairs) * sum(spectrum_y[left_border:right_border])
      }
    }
    # Save smoothed spectrum
    spectrum_y <- smoothed_spectrum_y
  }

  if (debug) {
    logf("Smoothed spectrum")
    debuglist$smoothed <- list(spectrum_y = spectrum_y)
  }

  # Peak selection procedure

  # Calculate second derivative of spectrum
  second_derivative <- matrix(nrow = 2, ncol = length(spectrum_x) - 2)
  for (i in 2:length(spectrum_x) - 1) {
    second_derivative[1, i - 1] <- spectrum_x[i]
    second_derivative[2, i - 1] <- spectrum_y[i - 1] + spectrum_y[i + 1] - 2 * spectrum_y[i]
  }

  # Find all local minima of second derivative
  peaks_x <- c()
  peaks_index <- c()
  second_derivative_border <- ncol(second_derivative) - 1
  for (i in 2:second_derivative_border) {
    if (second_derivative[2, i] < 0) {
      if ((second_derivative[2, i] <= second_derivative[2, i - 1]) & (second_derivative[2, i] < second_derivative[2, i + 1])) {
        # if(((spectrum_y[i+1] >= spectrum_y[i]) & (spectrum_y[i+1] > spectrum_y[i+2])) | ((spectrum_y[i+1] > spectrum_y[i]) & (spectrum_y[i+1] >= spectrum_y[i+2]))){
        # Add local minima to peak list
        peaks_x <- c(peaks_x, second_derivative[1, i])
        peaks_index <- c(peaks_index, i)
      }
    }
  }

  # Find all left positions of all local minima of second derivative
  left_position <- matrix(nrow = 1, ncol = length(peaks_x))
  for (i in 1:length(peaks_x)) {
    # Save next left position of current local minima
    next_left <- peaks_index[i] + 1
    while ((peaks_index[i] < next_left) & (next_left < ncol(second_derivative))) {
      if (second_derivative[2, next_left - 1] < second_derivative[2, next_left]) {
        if (((second_derivative[2, next_left - 1] < second_derivative[2, next_left]) & (second_derivative[2, next_left + 1] <= second_derivative[2, next_left])) | ((second_derivative[2, next_left] < 0) & (second_derivative[2, next_left + 1] >= 0))) {
          left_position[i] <- next_left
          break
        } else {
          next_left <- next_left + 1
        }
      } else {
        next_left <- next_left + 1
      }
    }
  }

  # Find all right positions of all local minima of second derivative
  right_position <- matrix(nrow = 1, ncol = length(peaks_x))
  for (i in 1:length(peaks_x)) {
    # Save next right position of current local minima
    next_right <- peaks_index[i] - 1
    while ((next_right < peaks_index[i]) & (next_right >= 2)) {
      if (second_derivative[2, next_right + 1] < second_derivative[2, next_right]) {
        if (((second_derivative[2, next_right + 1] < second_derivative[2, next_right]) & (second_derivative[2, next_right - 1] <= second_derivative[2, next_right])) | ((second_derivative[2, next_right] < 0) & (second_derivative[2, next_right - 1] >= 0))) {
          right_position[i] <- next_right
          break
        } else {
          next_right <- next_right - 1
        }
      } else {
        next_right <- next_right - 1
      }
    }
  }

  if (debug) {
    logf("Selected peaks")
    debuglist$peaks_sel <- list(
      spectrum_length = spectrum_length,
      second_derivative = second_derivative,
      peaks_index = peaks_index,
      peaks_x = peaks_index,
      left_position = left_position,
      right_position = right_position
    )
  }

  # Check borders of peak triplets
  # If NA values are available, remove corresponding peak triplet
  for (i in length(left_position):1) {
    if (is.na(left_position[i]) | (is.na(right_position[i]))) {
      peaks_x <- peaks_x[-i]
      peaks_index <- peaks_index[-i]
      left_position <- left_position[-i]
      right_position <- right_position[-i]
    }
  }
  if (debug) {
    logf("Removed peaks without border")
    debuglist$peaks_wob_rm <- list( # peaks without borders removed
      peaks_x = peaks_x,
      peaks_index = peaks_index,
      left_position = left_position,
      right_position = right_position
    )
  }

  # Calculate peak triplet score to distinguish between signal and noise
  scores <- matrix(nrow = 1, ncol = length(peaks_x))
  scores_left <- matrix(nrow = 1, ncol = length(peaks_x))
  scores_right <- matrix(nrow = 1, ncol = length(peaks_x))
  for (i in 1:length(peaks_x)) {
    # Calculate left score
    left_score <- 0
    for (j in peaks_index[i]:left_position[i]) {
      left_score <- sum(left_score, abs(second_derivative[2, j]))
    }
    scores_left[i] <- left_score
    # Calculate right score
    right_score <- 0
    for (k in right_position[i]:peaks_index[i]) {
      right_score <- sum(right_score, abs(second_derivative[2, k]))
    }
    scores_right[i] <- right_score
    # Save minimum score
    scores[i] <- min(left_score, right_score)
  }

  # Calculate mean of the score and standard deviation of the score of the signal free region R
  index_left <- which(spectrum_x[peaks_index + 1] >= signal_free_region_left)
  index_right <- which(spectrum_x[peaks_index + 1] <= signal_free_region_right)

  mean_score <- mean(c(scores[index_left], scores[index_right]))
  sd_score <- stats::sd(c(scores[index_left], scores[index_right]))

  # Filter peak triplets
  filtered_peaks <- c()
  filtered_left_position <- c()
  filtered_right_position <- c()
  save_scores <- c()
  for (i in 1:length(peaks_x)) {
    if (scores[i] >= mean_score + delta * sd_score) {
      # Save peak position
      filtered_peaks <- c(filtered_peaks, peaks_index[i])
      # Save left position
      filtered_left_position <- c(filtered_left_position, left_position[i])
      # Save right position
      filtered_right_position <- c(filtered_right_position, right_position[i])
      # Save value of scores of filtered peaks
      save_scores <- c(save_scores, scores[i])
    }
  }
  if (debug) {
    logf("Filtered peaks by score")
    debuglist$peak_scores_calc <- list(
      scores = scores,
      scores_left = scores_left,
      scores_right = scores_right,
      index_left = index_left,
      index_right = index_right,
      mean_score = mean_score,
      sd_score = sd_score,
      filtered_peaks = filtered_peaks,
      filtered_left_position = filtered_left_position,
      filtered_right_position = filtered_right_position,
      save_scores = save_scores
    )
  }

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

  # Calculate parameters w, lambda and A for the initial lorentz curves
  for (i in 1:length(filtered_peaks)) {
    # Calculate position of peak triplets
    w_1 <- c(w_1, spectrum_x[filtered_left_position[i] + 1])
    w_2 <- c(w_2, spectrum_x[filtered_peaks[i] + 1])
    w_3 <- c(w_3, spectrum_x[filtered_right_position[i] + 1])

    # Calculate intensity of peak triplets
    y_1 <- c(y_1, spectrum_y[filtered_left_position[i] + 1])
    y_2 <- c(y_2, spectrum_y[filtered_peaks[i] + 1])
    y_3 <- c(y_3, spectrum_y[filtered_right_position[i] + 1])

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
    w_result <- (w_1[i]^2 * y_1[i] * y_2_3[i] + w_3[i]^2 * y_3[i] * y_1_2[i] + w_2[i]^2 * y_2[i] * (-y_1_3[i])) / (2 * w_1_2[i] * y_1[i] * y_2[i] - 2 * (w_1_3[i] * y_1[i] + (-w_2_3[i]) * y_2[i]) * y_3[i])
    w_result <- w_result + w_delta[i]
    w <- c(w, w_result)
    # Wenn y Werte nach der H?henanpassung 0 werden, so ist w_new[i] NaN
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
  lorentz_curves_initial <- matrix(nrow = length(filtered_peaks), ncol = length(spectrum_x))
  for (i in 1:length(filtered_peaks)) {
    # If A = 0, then the lorentz curve is a zero line
    if (A[i] == 0) {
      lorentz_curves_initial[i, ] <- 0
    } else {
      lorentz_curves_initial[i, ] <- abs(A[i] * (lambda[i] / (lambda[i]^2 + (spectrum_x - w[i])^2)))
    }
  }

  if (debug) {
    logf("Initialized lorentz curve parameters")
    debuglist$params_init <- list(
      A = A,
      filtered_left_position = filtered_left_position,
      filtered_peaks = filtered_peaks,
      filtered_right_position = filtered_right_position,
      lambda = lambda,
      number_iterations = number_iterations,
      spectrum_x = spectrum_x,
      spectrum_y = spectrum_y,
      w = w,
      w_delta = w_delta
    )
  }

  # Approximation of lorentz curves
  for (b in 1:number_iterations) {
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
      w_1_new <- c(w_1_new, spectrum_x[filtered_left_position[i] + 1])
      w_2_new <- c(w_2_new, spectrum_x[filtered_peaks[i] + 1])
      w_3_new <- c(w_3_new, spectrum_x[filtered_right_position[i] + 1])

      # Calculate the sum of all lorentz curves for each data point
      sum_left[i] <- sum(lorentz_curves_initial[1:length(filtered_left_position), filtered_left_position[i] + 1])
      sum_peaks[i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), filtered_peaks[i] + 1])
      sum_right[i] <- sum(lorentz_curves_initial[1:length(filtered_right_position), filtered_right_position[i] + 1])

      # Calculate the proprotion between original spectrum an the sum of the lorentz curves for each peak triplets position
      proportion_left[i] <- spectrum_y[filtered_left_position[i] + 1] / sum_left[i]
      proportion_peaks[i] <- spectrum_y[filtered_peaks[i] + 1] / sum_peaks[i]
      proportion_right[i] <- spectrum_y[filtered_right_position[i] + 1] / sum_right[i]

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
        lorentz_curves_initial[i, ] <- abs(A_new[i] * (lambda_new[i] / (lambda_new[i]^2 + (spectrum_x - w_new[i])^2)))
      }
    }

    # Calculate sum of lorentz curves
    spectrum_approx <- matrix(nrow = 1, ncol = length(spectrum_x))
    for (i in 1:length(spectrum_x)) {
      spectrum_approx[1, i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), i])
    }

    # # Calculate the difference between original spectrum and height corrections
    # difference <- c()
    # for(i in 1:length(spectrum_x)){
    #   difference[i] <- (spectrum_y[i] - spectrum_approx[i])^2
    # }
    # mse <- (1/length(difference))*sum(difference)
    # message(paste("\nMSE value of iteration", b, "is: "))
    # #print(mse)
    # #paste("MSE value of iteration", b, "is:", mse)


    # Standardization of spectra so that total area equals 1
    spectrum_y_normed <- c()
    spectrum_approx_normed <- c()

    # Standardize the spectra
    spectrum_y_normed <- spectrum_y / sum(spectrum_y)
    spectrum_approx_normed <- spectrum_approx / sum(spectrum_approx)

    # Calculate the difference between normed original spectrum and normed approximated spectrum
    difference_normed <- c()
    for (i in 1:length(spectrum_x)) {
      difference_normed[i] <- (spectrum_y_normed[i] - spectrum_approx_normed[i])^2
    }
    mse_normed <- (1 / length(difference_normed)) * sum(difference_normed)
    message(paste("\nNormed MSE value of iteration", b, "is: "))
    print(mse_normed)
  }

  # Calculate the integrals for each lorentz curve
  integrals <- matrix(nrow = 1, ncol = length(lambda_new))
  for (i in 1:length(lambda_new)) {
    integrals[1, i] <- A_new[i] * (atan((-w_new[i] + (spectrum_length / factor_x)) / lambda_new[i]) - atan((-w_new[i]) / lambda_new[i]))
  }

  if (debug) {
    logf("Approximated lorentz curve parameters")
    debuglist$params_approx <- list(
      A_new = A_new,
      lambda_new = lambda_new,
      w_new = w_new,
      integrals = integrals,
      spectrum_y_normed = spectrum_y_normed,
      spectrum_approx_normed = spectrum_approx_normed,
      difference_normed = difference_normed,
      mse_normed = mse_normed
    )
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
    peak_triplets_middle[i] <- spectrum_x_ppm[index_peak_triplets_middle[i]]
    peak_triplets_left[i] <- spectrum_x_ppm[index_peak_triplets_left[i]]
    peak_triplets_right[i] <- spectrum_x_ppm[index_peak_triplets_right[i]]
  }

  # Save values A_new, lambda_new, w_new and noise_threshold to txt document
  noise_threshold <- replicate(length(w_new), 0)
  noise_threshold[1] <- mean_score + delta * sd_score
  spectrum_info <- data.frame(rbind(w_new, lambda_new, A_new, noise_threshold))
  spectrum_output <- data.frame(spectrum_approx)
  name_info_txt <- paste(name, "parameters.txt")
  name_output_txt <- paste(name, "approximated_spectrum.txt")

  if (is.null(store_results)) {
    store_results <- get_yn_input(sprintf("Save results as text documents at %s?", getwd()))
  }
  if (store_results) {
    message(sprintf("Writing %s", norm_path(name_info_txt)))
    utils::write.table(spectrum_info, name_info_txt, sep=",", col.names=FALSE, append=FALSE)
    message(sprintf("Writing %s", norm_path(name_output_txt)))
    utils::write.table(spectrum_output, name_output_txt, sep=",", col.names=FALSE, append=FALSE)
  } else {
    message("Skipping saving of results.")
  }

  if (debug) {
    logf("Reached point where parameters were saved to txt documents previously")
    debuglist$params_saved <- list(
      index_peak_triplets_middle = index_peak_triplets_middle,
      index_peak_triplets_left = index_peak_triplets_left,
      index_peak_triplets_right = index_peak_triplets_right,
      peak_triplets_middle = peak_triplets_middle,
      peak_triplets_left = peak_triplets_left,
      peak_triplets_right = peak_triplets_right,
      noise_threshold = noise_threshold,
      spectrum_info = spectrum_info,
      spectrum_output = spectrum_output,
      name_info_txt = name_info_txt,
      name_output_txt = name_output_txt
    )
  }

  return_list <- list(
    "filename" = name,
    "spectrum_x" = spectrum_x,
    "spectrum_x_ppm" = spectrum_x_ppm,
    "spectrum_y" = spectrum_y,
    "lorentz_curves" = lorentz_curves_initial,
    "mse_normed" = mse_normed,
    "spectrum_approx" = spectrum_approx,
    "index_peak_triplets_middle" = index_peak_triplets_middle,
    "index_peak_triplets_left" = index_peak_triplets_left,
    "index_peak_triplets_right" = index_peak_triplets_right,
    "peak_triplets_middle" = peak_triplets_middle,
    "peak_triplets_left" = peak_triplets_left,
    "peak_triplets_right" = peak_triplets_right,
    "integrals" = integrals,
    "signal_free_region" = signal_free_region,
    "range_water_signal_ppm" = range_water_signal_ppm,
    "A" = A_new,
    "lambda" = lambda_new,
    "w" = w_new,
    "store_results" = store_results
  )

  if (debug) {
    return_list$debuglist <- debuglist
    logf("Finished deconvolution")
  }

  return(return_list)
}

# Plots #####

#' @export
#' @title Plot peak triplets for variable range
#' @description Plots the peak triplets for each peak detected by [MetaboDecon1D()] and stores the plots as png at `outdir`.
#' @author Martina Haeckl, Tobias Schmidt
#' @param deconv_result Saved result of the MetaboDecon1D() function
#' @param x_range Row vector with two entries consisting of the ppm start and the ppm end value to scale the range of the x-axis (optional)
#' @param y_range Row vector with two entries consisting of the ppm start and the ppm end value to scale the range of the y-axis (optional)
#' @param out_dir Directory to save the png files (optional)
#' @param ask Logical value to ask the user if the png files should be saved in the specified directory (optional)
#' @return No return value, called for side effect of plotting.
#' @seealso [MetaboDecon1D()], [calculate_lorentz_curves()], [plot_lorentz_curves_save_as_png()], [plot_spectrum_superposition_save_as_png()]
#' @examples
#' sim <- metabodecon_file("bruker/sim_subset")
#' sim_decon <- generate_lorentz_curves(sim, sfr = c(3.58, 3.42), wshw = 0, ask = FALSE, verbose = FALSE)
#' png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
#' plot_triplets(sim_decon, out_dir = png_dir, ask = FALSE)
#' dir(png_dir, full.names = TRUE)
plot_triplets <- function(deconv_result, x_range = c(), y_range = c(), out_dir = ".", ask = TRUE) {

  out_dir <- normPath(out_dir)
  if (ask) {
    continue <- get_yn_input(sprintf("Continue creating pngs in %s?", out_dir))
    if (!continue) {
      invisible(NULL)
    }
  }
  owd <- setwd(out_dir)
  on.exit(setwd(owd), add = TRUE)

  number_in_folder <- 0
  if ("number_of_files" %in% names(deconv_result)) {
    number_of_files <- deconv_result$number_of_files
  } else {
    if ("number_of_files" %in% names(deconv_result[[1]])) {
      number_of_files <- deconv_result[[1]]$number_of_files
      if (number_of_files == 1) number_in_folder <- 1
    }
  }

  if (number_of_files > 1 | number_in_folder == 1) {
    # Check if user input comprise a $ sign
    if (grepl("[$]", deparse(substitute(deconv_result)))) {
      # User would like to plot triplets of only one spectrum
      name <- deconv_result$filename
      spectrum_x_ppm <- deconv_result$x_values_ppm
      spectrum_y <- deconv_result$y_values
      index_peaks_triplets <- deconv_result$index_peak_triplets_middle
      index_left_position <- deconv_result$index_peak_triplets_left
      index_right_position <- deconv_result$index_peak_triplets_right

      message(paste("Plot triplets of ", name))

      # Check if x_range is adjusted or not
      if (is.null(x_range)) {
        filename <- paste(name, "_peak_triplets.png", sep = "")
        # Save plot as png
        grDevices::png(file = filename, width = 825, height = 525)
        plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = rev(range(spectrum_x_ppm)), ylim = y_range)
        graphics::points(spectrum_x_ppm[index_peaks_triplets], spectrum_y[index_peaks_triplets], col = "red", cex = 1.3)
        graphics::points(spectrum_x_ppm[index_right_position], spectrum_y[index_right_position], col = "green", cex = 1.3)
        graphics::points(spectrum_x_ppm[index_left_position], spectrum_y[index_left_position], col = "blue", cex = 1.3)
        graphics::legend("topright", legend = c("x_left", "x_middle", "x_right"), col = c("green", "red", "blue"), pch = 1, bty = "n", cex = 1.3)
        grDevices::dev.off()
      } else {
        filename <- paste(name, "_peak_triplets.png", sep = "")
        # Save plot as png
        grDevices::png(file = filename, width = 825, height = 525)
        plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = x_range, ylim = y_range)
        graphics::points(spectrum_x_ppm[index_peaks_triplets], spectrum_y[index_peaks_triplets], col = "red", cex = 1.3)
        graphics::points(spectrum_x_ppm[index_right_position], spectrum_y[index_right_position], col = "green", cex = 1.3)
        graphics::points(spectrum_x_ppm[index_left_position], spectrum_y[index_left_position], col = "blue", cex = 1.3)
        graphics::legend("topright", legend = c("x_left", "x_middle", "x_right"), col = c("green", "red", "blue"), pch = 1, bty = "n", cex = 1.3)
        grDevices::dev.off()
      }
    } else {
      for (l in 1:number_of_files) {
        # User would like to plot triplets of all investigated spectra
        name <- deconv_result[[l]]$filename
        spectrum_x_ppm <- deconv_result[[l]]$x_values_ppm
        spectrum_y <- deconv_result[[l]]$y_values
        index_peaks_triplets <- deconv_result[[l]]$index_peak_triplets_middle
        index_left_position <- deconv_result[[l]]$index_peak_triplets_left
        index_right_position <- deconv_result[[l]]$index_peak_triplets_right

        message(paste("Plot triplets of ", name))

        # Check if x_range is adjusted or not
        if (is.null(x_range)) {
          filename <- paste(name, "_peak_triplets.png", sep = "")
          # Save plot as png
          grDevices::png(file = filename, width = 825, height = 525)
          plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = rev(range(spectrum_x_ppm)), ylim = y_range)
          graphics::points(spectrum_x_ppm[index_peaks_triplets], spectrum_y[index_peaks_triplets], col = "red", cex = 1.3)
          graphics::points(spectrum_x_ppm[index_right_position], spectrum_y[index_right_position], col = "green", cex = 1.3)
          graphics::points(spectrum_x_ppm[index_left_position], spectrum_y[index_left_position], col = "blue", cex = 1.3)
          graphics::legend("topright", legend = c("x_left", "x_middle", "x_right"), col = c("green", "red", "blue"), pch = 1, bty = "n", cex = 1.3)
          grDevices::dev.off()
        } else {
          filename <- paste(name, "_peak_triplets.png", sep = "")
          # Save plot as png
          grDevices::png(file = filename, width = 825, height = 525)
          plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = x_range, ylim = y_range)
          graphics::points(spectrum_x_ppm[index_peaks_triplets], spectrum_y[index_peaks_triplets], col = "red", cex = 1.3)
          graphics::points(spectrum_x_ppm[index_right_position], spectrum_y[index_right_position], col = "green", cex = 1.3)
          graphics::points(spectrum_x_ppm[index_left_position], spectrum_y[index_left_position], col = "blue", cex = 1.3)
          graphics::legend("topright", legend = c("x_left", "x_middle", "x_right"), col = c("green", "red", "blue"), pch = 1, bty = "n", cex = 1.3)
          grDevices::dev.off()
        }
      }
    }
  } else {
    name <- deconv_result$filename
    spectrum_x_ppm <- deconv_result$x_values_ppm
    spectrum_y <- deconv_result$y_values
    index_peaks_triplets <- deconv_result$index_peak_triplets_middle
    index_left_position <- deconv_result$index_peak_triplets_left
    index_right_position <- deconv_result$index_peak_triplets_right

    message(paste("Plot triplets of ", name))

    # Check if x_range is adjusted or not
    if (is.null(x_range)) {
      filename <- paste(name, "_peak_triplets.png", sep = "")
      # Save plot as png
      grDevices::png(file = filename, width = 825, height = 525)
      plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = rev(range(spectrum_x_ppm)), ylim = y_range)
      graphics::points(spectrum_x_ppm[index_peaks_triplets], spectrum_y[index_peaks_triplets], col = "red", cex = 1.3)
      graphics::points(spectrum_x_ppm[index_right_position], spectrum_y[index_right_position], col = "green", cex = 1.3)
      graphics::points(spectrum_x_ppm[index_left_position], spectrum_y[index_left_position], col = "blue", cex = 1.3)
      graphics::legend("topright", legend = c("x_left", "x_middle", "x_right"), col = c("green", "red", "blue"), pch = 1, bty = "n", cex = 1.3)
      grDevices::dev.off()
    } else {
      filename <- paste(name, "_peak_triplets.png", sep = "")
      # Save plot as png
      grDevices::png(file = filename, width = 825, height = 525)
      plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = x_range, ylim = y_range)
      graphics::points(spectrum_x_ppm[index_peaks_triplets], spectrum_y[index_peaks_triplets], col = "red", cex = 1.3)
      graphics::points(spectrum_x_ppm[index_right_position], spectrum_y[index_right_position], col = "green", cex = 1.3)
      graphics::points(spectrum_x_ppm[index_left_position], spectrum_y[index_left_position], col = "blue", cex = 1.3)
      graphics::legend("topright", legend = c("x_left", "x_middle", "x_right"), col = c("green", "red", "blue"), pch = 1, bty = "n", cex = 1.3)
      grDevices::dev.off()
    }
  }
}

#' @export
#' @title Plot lorentz curves for variable range
#' @description Plots the original spectrum and all generated Lorentz curves and save the result as png under the filepath.
#' @param deconv_result Saved result of the MetaboDecon1D() function
#' @param x_range Row vector with two entries consisting of the ppm start and the ppm end value to scale the range of the x-axis (optional)
#' @param y_range Row vector with two entries consisting of the ppm start and the ppm end value to scale the range of the y-axis (optional)
#' @param out_dir Path to the directory where the png files should be saved. Default is the current working directory.
#' @param ask Logical value. Whether to ask for confirmation from the user before writing files to disk. Default is TRUE.
#' @return NULL, called for side effects.
#' @seealso [MetaboDecon1D()], [plot_triplets()], [plot_spectrum_superposition_save_as_png()]
#' @examples
#' sim <- metabodecon_file("bruker/sim_subset")
#' sim_decon <- generate_lorentz_curves(sim, sfr = c(3.58, 3.42), wshw = 0, ask = FALSE, verbose = FALSE)
#' png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
#' plot_lorentz_curves_save_as_png(sim_decon, out_dir = png_dir, ask = FALSE)
#' dir(png_dir, full.names = TRUE)
plot_lorentz_curves_save_as_png <- function(deconv_result, x_range = c(), y_range = c(), out_dir = ".", ask = TRUE) {

  out_dir <- normPath(out_dir)
  if (ask) {
    continue <- get_yn_input(sprintf("Continue creating pngs in %s?", out_dir))
    if (!continue) {
      invisible(NULL)
    }
  }
  owd <- setwd(out_dir)
  on.exit(setwd(owd), add = TRUE)

  number_in_folder <- 0
  # Get number of analyzed spectra
  if ("number_of_files" %in% names(deconv_result)) {
    number_of_files <- deconv_result$number_of_files
  } else {
    if ("number_of_files" %in% names(deconv_result[[1]])) {
      number_of_files <- deconv_result[[1]]$number_of_files
      # Check if only one spectrum is inside whole folder
      if (number_of_files == 1) {
        number_in_folder <- 1
      }
    }
  }

  # Check how many spectra are investigated
  if (number_of_files > 1 | number_in_folder == 1) {
    # Check if user input comprise a $ sign
    if (grepl("[$]", deparse(substitute(deconv_result)))) {
      # Get parameters
      name <- deconv_result$filename
      spectrum_x_ppm <- deconv_result$x_values_ppm
      spectrum_y <- deconv_result$y_values

      # Save additional parameter to know that $ sign was used by the user
      number_of_files <- 1

      # Calculate Lorentz curves
      lorentz_curves_initial <- calculate_lorentz_curves(deconv_result, number_of_files)

      message(paste("Plot Lorentz curves of", name))

      # Check if x_range is adjusted or not
      if (is.null(x_range)) {
        filename <- paste(name, "_lorentz_curves.png", sep = "")
        # Save plot as png
        grDevices::png(file = filename, width = 825, height = 525)
        plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = rev(range(spectrum_x_ppm)), ylim = y_range)
        for (i in 1:dim(lorentz_curves_initial)[1]) {
          graphics::lines(spectrum_x_ppm, lorentz_curves_initial[i, ], col = "red")
        }
        graphics::legend("topright", legend = c("Original spectrum", "Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
        grDevices::dev.off()
      } else {
        filename <- paste(name, "_lorentz_curves.png", sep = "")
        # Save plot as png
        grDevices::png(file = filename, width = 825, height = 525)
        plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = x_range, ylim = y_range)
        for (i in 1:dim(lorentz_curves_initial)[1]) {
          graphics::lines(spectrum_x_ppm, lorentz_curves_initial[i, ], col = "red")
        }
        graphics::legend("topright", legend = c("Original spectrum", "Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
        grDevices::dev.off()
      }
    } else {
      # Calculate Lorentz curves
      lorentz_curves_matrix <- calculate_lorentz_curves(deconv_result)

      for (l in 1:number_of_files) {
        # Get necessary parameters
        name <- deconv_result[[l]]$filename
        spectrum_x_ppm <- deconv_result[[l]]$x_values_ppm
        spectrum_y <- deconv_result[[l]]$y_values
        lorentz_curves_initial <- lorentz_curves_matrix[[l]]
        filtered_peaks <- deconv_result[[l]]$peak_triplets_middle

        message(paste("Plot Lorentz curves of", name))

        # Check if x_range is adjusted or not
        if (is.null(x_range)) {
          filename <- paste(name, "_lorentz_curves.png", sep = "")
          # Save plot as png
          grDevices::png(file = filename, width = 825, height = 525)
          plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = rev(range(spectrum_x_ppm)), ylim = y_range)
          for (i in 1:length(filtered_peaks)) {
            graphics::lines(spectrum_x_ppm, lorentz_curves_initial[i, ], col = "red")
          }
          graphics::legend("topright", legend = c("Original spectrum", "Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
          grDevices::dev.off()
        } else {
          filename <- paste(name, "_lorentz_curves.png", sep = "")
          # Save plot as png
          grDevices::png(file = filename, width = 825, height = 525)
          plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = x_range, ylim = y_range)
          for (i in 1:length(filtered_peaks)) {
            graphics::lines(spectrum_x_ppm, lorentz_curves_initial[i, ], col = "red")
          }
          graphics::legend("topright", legend = c("Original spectrum", "Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
          grDevices::dev.off()
        }
      }
    }
  } else {
    # Calculate Lorentz curves
    lorentz_curves_initial <- calculate_lorentz_curves(deconv_result)

    # Get parameters
    name <- deconv_result$filename
    spectrum_x_ppm <- deconv_result$x_values_ppm
    spectrum_y <- deconv_result$y_values
    filtered_peaks <- deconv_result$peak_triplets_middle

    message(paste("Plot Lorentz curves of", name))

    # Check if x_range is adjusted or not
    if (is.null(x_range)) {
      filename <- paste(name, "_lorentz_curves.png", sep = "")
      # Save plot as png
      grDevices::png(file = filename, width = 825, height = 525)
      plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = rev(range(spectrum_x_ppm)), ylim = y_range)
      for (i in 1:length(filtered_peaks)) {
        graphics::lines(spectrum_x_ppm, lorentz_curves_initial[i, ], col = "red")
      }
      graphics::legend("topright", legend = c("Original spectrum", "Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
      grDevices::dev.off()
    } else {
      filename <- paste(name, "_lorentz_curves.png", sep = "")
      # Save plot as png
      grDevices::png(file = filename, width = 825, height = 525)
      plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = x_range, ylim = y_range)
      for (i in 1:length(filtered_peaks)) {
        graphics::lines(spectrum_x_ppm, lorentz_curves_initial[i, ], col = "red")
      }
      graphics::legend("topright", legend = c("Original spectrum", "Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
      grDevices::dev.off()
    }
  }
}

#' @export
#' @title Plot spectrum approx for variable range
#' @description Plots the original spectrum and the superposition of all generated Lorentz curves and saves the result as png under the specified filepath.
#' @author Martina Haeckl
#' @param deconv_result Saved result of the MetaboDecon1D() function
#' @param x_range Row vector with two entries consisting of the ppm start and the ppm end value to scale the range of the x-axis (optional)
#' @param y_range Row vector with two entries consisting of the ppm start and the ppm end value to scale the range of the y-axis (optional)
#' @param out_dir Path to the directory where the png files should be saved. Default is the current working directory.
#' @param ask Logical value. Whether to ask for confirmation from the user before writing files to disk.
#' @return NULL, called for side effects.
#' @seealso [MetaboDecon1D()], [calculate_lorentz_curves()], [plot_triplets()], [plot_lorentz_curves_save_as_png()]
#' @examples
#' sim <- metabodecon_file("bruker/sim_subset")
#' sim_decon <- generate_lorentz_curves(sim, sfr = c(3.58, 3.42), wshw = 0, ask = FALSE, verbose = FALSE)
#' png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
#' plot_spectrum_superposition_save_as_png(sim_decon, out_dir = png_dir, ask = FALSE)
#' dir(png_dir, full.names = TRUE)
plot_spectrum_superposition_save_as_png <- function(deconv_result, x_range = c(), y_range = c(), out_dir = ".", ask = TRUE) {

  out_dir <- normPath(out_dir)
  if (ask) {
    continue <- get_yn_input(sprintf("Continue creating pngs in %s?", out_dir))
    if (!continue) {
      invisible(NULL)
    }
  }
  owd <- setwd(out_dir)
  on.exit(setwd(owd), add = TRUE)

  number_in_folder <- 0
  # Get number of analyzed spectra
  if ("number_of_files" %in% names(deconv_result)) {
    number_of_files <- deconv_result$number_of_files
  } else {
    if ("number_of_files" %in% names(deconv_result[[1]])) {
      number_of_files <- deconv_result[[1]]$number_of_files
      # Check if only one spectrum is inside whole folder
      if (number_of_files == 1) {
        number_in_folder <- 1
      }
    }
  }

  # Check how many spectra are investigated
  if (number_of_files > 1 | number_in_folder == 1) {
    # Check if user input comprise a $ sign
    if (grepl("[$]", deparse(substitute(deconv_result)))) {
      name <- deconv_result$filename
      spectrum_x_ppm <- deconv_result$x_values_ppm
      spectrum_y <- deconv_result$y_values
      spectrum_approx <- deconv_result$spectrum_superposition
      mse <- deconv_result$mse_normed
      filtered_peaks <- deconv_result$peak_triplets_middle

      message(paste("Plot superposition of", name))

      # Check if x_range is adjusted or not
      if (is.null(x_range)) {
        filename <- paste(name, "_sum_lorentz_curves.png", sep = "")
        # Save plot as png
        grDevices::png(file = filename, width = 825, height = 525)
        plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = rev(range(spectrum_x_ppm)), ylim = y_range)
        graphics::lines(spectrum_x_ppm, spectrum_approx, col = "red")
        graphics::legend("topright", legend = c("Original spectrum", "Sum of Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
        text <- paste("MSE_Normed = ", mse, sep = "")
        graphics::mtext(text, side = 3)
        grDevices::dev.off()
      } else {
        filename <- paste(name, "_sum_lorentz_curves.png", sep = "")
        # Save plot as png
        grDevices::png(file = filename, width = 825, height = 525)
        plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = x_range, ylim = y_range)
        graphics::lines(spectrum_x_ppm, spectrum_approx, col = "red")
        graphics::legend("topright", legend = c("Original spectrum", "Sum of Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
        text <- paste("MSE_Normed = ", mse, sep = "")
        graphics::mtext(text, side = 3)
        grDevices::dev.off()
      }
    } else {
      for (l in 1:number_of_files) {
        # Get necessary parameters
        name <- deconv_result[[l]]$filename
        spectrum_x_ppm <- deconv_result[[l]]$x_values_ppm
        spectrum_y <- deconv_result[[l]]$y_values
        spectrum_approx <- deconv_result[[l]]$spectrum_superposition
        mse <- deconv_result[[l]]$mse_normed
        filtered_peaks <- deconv_result[[l]]$peak_triplets_middle

        message(paste("Plot superposition of", name))

        # Check if x_range is adjusted or not
        if (is.null(x_range)) {
          filename <- paste(name, "_sum_lorentz_curves.png", sep = "")
          # Save plot as png
          grDevices::png(file = filename, width = 825, height = 525)
          plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = rev(range(spectrum_x_ppm)), ylim = y_range)
          graphics::lines(spectrum_x_ppm, spectrum_approx, col = "red")
          graphics::legend("topright", legend = c("Original spectrum", "Sum of Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
          text <- paste("MSE_Normed = ", mse, sep = "")
          graphics::mtext(text, side = 3)
          grDevices::dev.off()
        } else {
          filename <- paste(name, "_sum_lorentz_curves.png", sep = "")
          # Save plot as png
          grDevices::png(file = filename, width = 825, height = 525)
          plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = x_range, ylim = y_range)
          graphics::lines(spectrum_x_ppm, spectrum_approx, col = "red")
          graphics::legend("topright", legend = c("Original spectrum", "Sum of Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
          text <- paste("MSE_Normed = ", mse, sep = "")
          graphics::mtext(text, side = 3)
          grDevices::dev.off()
        }
      }
    }
  } else {
    name <- deconv_result$filename
    spectrum_x_ppm <- deconv_result$x_values_ppm
    spectrum_y <- deconv_result$y_values
    spectrum_approx <- deconv_result$spectrum_superposition
    mse <- deconv_result$mse_normed
    filtered_peaks <- deconv_result$peak_triplets_middle

    message(paste("Plot superposition of", name))

    # Check if x_range is adjusted or not
    if (is.null(x_range)) {
      filename <- paste(name, "_sum_lorentz_curves.png", sep = "")
      # Save plot as png
      grDevices::png(file = filename, width = 825, height = 525)
      plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = rev(range(spectrum_x_ppm)), ylim = y_range)
      graphics::lines(spectrum_x_ppm, spectrum_approx, col = "red")
      graphics::legend("topright", legend = c("Original spectrum", "Sum of Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
      text <- paste("MSE_Normed = ", mse, sep = "")
      graphics::mtext(text, side = 3)
      grDevices::dev.off()
    } else {
      filename <- paste(name, "_sum_lorentz_curves.png", sep = "")
      # Save plot as png
      grDevices::png(file = filename, width = 825, height = 525)
      plot(spectrum_x_ppm, spectrum_y, type = "l", main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]", cex = 1.3, xlim = x_range, ylim = y_range)
      graphics::lines(spectrum_x_ppm, spectrum_approx, col = "red")
      graphics::legend("topright", legend = c("Original spectrum", "Sum of Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
      text <- paste("MSE_Normed = ", mse, sep = "")
      graphics::mtext(text, side = 3)
      grDevices::dev.off()
    }
  }
}

# Helpers #####

#' @export
#' @title Calculate lorentz curves for each analyzed spectrum
#' @description Calculates the lorentz curves of each investigated spectrum.
#' @author Martina Haeckl, Tobias Schmidt
#' @param deconv_result A list as returned by [generate_lorentz_curves()] or [MetaboDecon1D].
#' @param number_of_files Number of spectra to analyze
#' @return If `deconv_result` holds the result of a single deconvolution, a matrix containing the generated Lorentz curves is returned, where each row depicts one Lorentz curve. If `deconv_result` is a list of deconvoluted spectra, a list of such matrices is returned.
#' @examples
#'
#' # Deconvolute the spectra in folder "bruker/sim_subset" into a list of
#' # Lorentz Curves (specified via the parameters A, lambda and x_0).
#' sim <- metabodecon_file("bruker/sim_subset")
#' decons <- generate_lorentz_curves(sim, sfr = c(3.58, 3.42), wshw = 0, ask = FALSE, verbose = FALSE)
#' decon1 <- decons[[1]]
#'
#' # Calculate the corresponding y values at each ppm value for each Lorentz
#' # Curve. I.e. you get a matrix of dimension n x m for each deconvolution,
#' # where n is the number of Lorentz Curves and m is the number of ppm values.
#' Xs <- calculate_lorentz_curves(decons)
#' X1 <- Xs[[1]]
#' dim(X1)
#'
#' # Visualize the 5th, 9th and 11th Lorentz curve.
#' nrs <- c(5, 9, 11)
#' col <- c("red", "blue", "orange")
#' desc <- paste("Lorentz curve", nrs)
#' plot(decon1$x_values_ppm, decon1$y_values, type = "l", lty = 2)
#' for (i in 1:3) lines(decon1$x_values_ppm, X1[nrs[i], ], col = col[i])
#' legend("topright", legend = desc, col = col, lty = 1)
#' @seealso [MetaboDecon1D()], [plot_triplets()], [plot_lorentz_curves_save_as_png()], [plot_spectrum_superposition_save_as_png()]
calculate_lorentz_curves <- function(deconv_result, number_of_files = NA) {
  number_in_folder <- 0
  if (is.na(number_of_files)) {
    # Get number of analyzed spectra
    if ("number_of_files" %in% names(deconv_result)) {
      number_of_files <- deconv_result$number_of_files
    } else {
      if ("number_of_files" %in% names(deconv_result[[1]])) {
        number_of_files <- deconv_result[[1]]$number_of_files
        # Check if only one spectrum is inside whole folder
        if (number_of_files == 1) {
          number_in_folder <- 1
        }
      }
    }
  }

  # Check if more than one spectra was analyzed
  if (number_of_files > 1 | number_in_folder == 1) {
    # Check if user input comprise a $ sign
    if (grepl("[$]", deparse(substitute(deconv_result)))) {
      spectrum_x <- deconv_result$x_values
      A_new <- deconv_result$A
      lambda_new <- deconv_result$lambda
      w_new <- deconv_result$x_0

      # Calculate lorentz curves
      lorentz_curves_initial <- matrix(nrow = length(A_new), ncol = length(spectrum_x))
      for (i in 1:length(A_new)) {
        if ((w_new[i] == 0) | (lambda_new[i] == 0) | (A_new[i] == 0)) {
          lorentz_curves_initial[i, ] <- 0
        } else {
          lorentz_curves_initial[i, ] <- abs(A_new[i] * (lambda_new[i] / (lambda_new[i]^2 + (spectrum_x - w_new[i])^2)))
        }
      }
      # Return matrix with each row contains one lorentz curve
      return(lorentz_curves_initial)
    } else {
      lorentz_curves_list <- list()
      for (l in 1:number_of_files) {
        name <- deconv_result[[l]]$filename
        spectrum_x <- deconv_result[[l]]$x_values
        A_new <- deconv_result[[l]]$A
        lambda_new <- deconv_result[[l]]$lambda
        w_new <- deconv_result[[l]]$x_0

        # Calculate lorentz curves
        lorentz_curves_initial <- matrix(nrow = length(A_new), ncol = length(spectrum_x))
        for (i in 1:length(A_new)) {
          if ((w_new[i] == 0) | (lambda_new[i] == 0) | (A_new[i] == 0)) {
            lorentz_curves_initial[i, ] <- 0
          } else {
            lorentz_curves_initial[i, ] <- abs(A_new[i] * (lambda_new[i] / (lambda_new[i]^2 + (spectrum_x - w_new[i])^2)))
          }
        }
        lorentz_curves_list[[paste0(name)]] <- lorentz_curves_initial
      }
      return(lorentz_curves_list)
    }
  } else {
    spectrum_x <- deconv_result$x_values
    A_new <- deconv_result$A
    lambda_new <- deconv_result$lambda
    w_new <- deconv_result$x_0

    # Calculate lorentz curves
    lorentz_curves_initial <- matrix(nrow = length(A_new), ncol = length(spectrum_x))
    for (i in 1:length(A_new)) {
      if ((w_new[i] == 0) | (lambda_new[i] == 0) | (A_new[i] == 0)) {
        lorentz_curves_initial[i, ] <- 0
      } else {
        lorentz_curves_initial[i, ] <- abs(A_new[i] * (lambda_new[i] / (lambda_new[i]^2 + (spectrum_x - w_new[i])^2)))
      }
    }
    # Return matrix with each row contains one lorentz curve
    return(lorentz_curves_initial)
  }
}

#' @noRd
#' @author Martina Haeckl
#' @title Show terms and conditions of license
#' @description Show terms and conditions of license
#' @param ... Not used
show_license <- function(...) {
  message("TERMS AND CONDITIONS
           0. Definitions.

           'This License' refers to version 3 of the GNU General Public License.

           'Copyright' also means copyright-like laws that apply to other kinds of works, such as semiconductor masks.

           'The Program' refers to any copyrightable work licensed under this License. Each licensee is addressed as 'you'. 'Licensees' and 'recipients' may be individuals or organizations.

           To 'modify' a work means to copy from or adapt all or part of the work in a fashion requiring copyright permission, other than the making of an exact copy. The resulting work is called a 'modified version' of the earlier work or a work 'based on' the earlier work.

           A 'covered work' means either the unmodified Program or a work based on the Program.

           To 'propagate' a work means to do anything with it that, without permission, would make you directly or secondarily liable for infringement under applicable copyright law, except executing it on a computer or modifying a private copy. Propagation includes copying, distribution (with or without modification), making available to the public, and in some countries other activities as well.

           To 'convey' a work means any kind of propagation that enables other parties to make or receive copies. Mere interaction with a user through a computer network, with no transfer of a copy, is not conveying.

           An interactive user interface displays 'Appropriate Legal Notices' to the extent that it includes a convenient and prominently visible feature that (1) displays an appropriate copyright notice, and (2) tells the user that there is no warranty for the work (except to the extent that warranties are provided), that licensees may convey the work under this License, and how to view a copy of this License. If the interface presents a list of user commands or options, such as a menu, a prominent item in the list meets this criterion.
           1. Source Code.

           The 'source code' for a work means the preferred form of the work for making modifications to it. 'Object code' means any non-source form of a work.

           A 'Standard Interface' means an interface that either is an official standard defined by a recognized standards body, or, in the case of interfaces specified for a particular programming language, one that is widely used among developers working in that language.

           The 'System Libraries' of an executable work include anything, other than the work as a whole, that (a) is included in the normal form of packaging a Major Component, but which is not part of that Major Component, and (b) serves only to enable use of the work with that Major Component, or to implement a Standard Interface for which an implementation is available to the public in source code form. A 'Major Component', in this context, means a major essential component (kernel, window system, and so on) of the specific operating system (if any) on which the executable work runs, or a compiler used to produce the work, or an object code interpreter used to run it.

           The 'Corresponding Source' for a work in object code form means all the source code needed to generate, install, and (for an executable work) run the object code and to modify the work, including scripts to control those activities. However, it does not include the work's System Libraries, or general-purpose tools or generally available free programs which are used unmodified in performing those activities but which are not part of the work. For example, Corresponding Source includes interface definition files associated with source files for the work, and the source code for shared libraries and dynamically linked subprograms that the work is specifically designed to require, such as by intimate data communication or control flow between those subprograms and other parts of the work.

The Corresponding Source need not include anything that users can regenerate automatically from other parts of the Corresponding Source.

The Corresponding Source for a work in source code form is that same work.
2. Basic Permissions.

All rights granted under this License are granted for the term of copyright on the Program, and are irrevocable provided the stated conditions are met. This License explicitly affirms your unlimited permission to run the unmodified Program. The output from running a covered work is covered by this License only if the output, given its content, constitutes a covered work. This License acknowledges your rights of fair use or other equivalent, as provided by copyright law.

You may make, run and propagate covered works that you do not convey, without conditions so long as your license otherwise remains in force. You may convey covered works to others for the sole purpose of having them make modifications exclusively for you, or provide you with facilities for running those works, provided that you comply with the terms of this License in conveying all material for which you do not control copyright. Those thus making or running the covered works for you must do so exclusively on your behalf, under your direction and control, on terms that prohibit them from making any copies of your copyrighted material outside their relationship with you.

Conveying under any other circumstances is permitted solely under the conditions stated below. Sublicensing is not allowed; section 10 makes it unnecessary.
3. Protecting Users' Legal Rights From Anti-Circumvention Law.

           No covered work shall be deemed part of an effective technological measure under any applicable law fulfilling obligations under article 11 of the WIPO copyright treaty adopted on 20 December 1996, or similar laws prohibiting or restricting circumvention of such measures.

           When you convey a covered work, you waive any legal power to forbid circumvention of technological measures to the extent such circumvention is effected by exercising rights under this License with respect to the covered work, and you disclaim any intention to limit operation or modification of the work as a means of enforcing, against the work's users, your or third parties' legal rights to forbid circumvention of technological measures.
           4. Conveying Verbatim Copies.

           You may convey verbatim copies of the Program's source code as you receive it, in any medium, provided that you conspicuously and appropriately publish on each copy an appropriate copyright notice; keep intact all notices stating that this License and any non-permissive terms added in accord with section 7 apply to the code; keep intact all notices of the absence of any warranty; and give all recipients a copy of this License along with the Program.

You may charge any price or no price for each copy that you convey, and you may offer support or warranty protection for a fee.
5. Conveying Modified Source Versions.

You may convey a work based on the Program, or the modifications to produce it from the Program, in the form of source code under the terms of section 4, provided that you also meet all of these conditions:

    a) The work must carry prominent notices stating that you modified it, and giving a relevant date.
    b) The work must carry prominent notices stating that it is released under this License and any conditions added under section 7. This requirement modifies the requirement in section 4 to 'keep intact all notices'.
    c) You must license the entire work, as a whole, under this License to anyone who comes into possession of a copy. This License will therefore apply, along with any applicable section 7 additional terms, to the whole of the work, and all its parts, regardless of how they are packaged. This License gives no permission to license the work in any other way, but it does not invalidate such permission if you have separately received it.
    d) If the work has interactive user interfaces, each must display Appropriate Legal Notices; however, if the Program has interactive interfaces that do not display Appropriate Legal Notices, your work need not make them do so.

A compilation of a covered work with other separate and independent works, which are not by their nature extensions of the covered work, and which are not combined with it such as to form a larger program, in or on a volume of a storage or distribution medium, is called an 'aggregate' if the compilation and its resulting copyright are not used to limit the access or legal rights of the compilation's users beyond what the individual works permit. Inclusion of a covered work in an aggregate does not cause this License to apply to the other parts of the aggregate.
           6. Conveying Non-Source Forms.

           You may convey a covered work in object code form under the terms of sections 4 and 5, provided that you also convey the machine-readable Corresponding Source under the terms of this License, in one of these ways:

             a) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by the Corresponding Source fixed on a durable physical medium customarily used for software interchange.
  b) Convey the object code in, or embodied in, a physical product (including a physical distribution medium), accompanied by a written offer, valid for at least three years and valid for as long as you offer spare parts or customer support for that product model, to give anyone who possesses the object code either (1) a copy of the Corresponding Source for all the software in the product that is covered by this License, on a durable physical medium customarily used for software interchange, for a price no more than your reasonable cost of physically performing this conveying of source, or (2) access to copy the Corresponding Source from a network server at no charge.
c) Convey individual copies of the object code with a copy of the written offer to provide the Corresponding Source. This alternative is allowed only occasionally and noncommercially, and only if you received the object code with such an offer, in accord with subsection 6b.
d) Convey the object code by offering access from a designated place (gratis or for a charge), and offer equivalent access to the Corresponding Source in the same way through the same place at no further charge. You need not require recipients to copy the Corresponding Source along with the object code. If the place to copy the object code is a network server, the Corresponding Source may be on a different server (operated by you or a third party) that supports equivalent copying facilities, provided you maintain clear directions next to the object code saying where to find the Corresponding Source. Regardless of what server hosts the Corresponding Source, you remain obligated to ensure that it is available for as long as needed to satisfy these requirements.
e) Convey the object code using peer-to-peer transmission, provided you inform other peers where the object code and Corresponding Source of the work are being offered to the general public at no charge under subsection 6d.

A separable portion of the object code, whose source code is excluded from the Corresponding Source as a System Library, need not be included in conveying the object code work.

A 'User Product' is either (1) a 'consumer product', which means any tangible personal property which is normally used for personal, family, or household purposes, or (2) anything designed or sold for incorporation into a dwelling. In determining whether a product is a consumer product, doubtful cases shall be resolved in favor of coverage. For a particular product received by a particular user, 'normally used' refers to a typical or common use of that class of product, regardless of the status of the particular user or of the way in which the particular user actually uses, or expects or is expected to use, the product. A product is a consumer product regardless of whether the product has substantial commercial, industrial or non-consumer uses, unless such uses represent the only significant mode of use of the product.

'Installation Information' for a User Product means any methods, procedures, authorization keys, or other information required to install and execute modified versions of a covered work in that User Product from a modified version of its Corresponding Source. The information must suffice to ensure that the continued functioning of the modified object code is in no case prevented or interfered with solely because modification has been made.

If you convey an object code work under this section in, or with, or specifically for use in, a User Product, and the conveying occurs as part of a transaction in which the right of possession and use of the User Product is transferred to the recipient in perpetuity or for a fixed term (regardless of how the transaction is characterized), the Corresponding Source conveyed under this section must be accompanied by the Installation Information. But this requirement does not apply if neither you nor any third party retains the ability to install modified object code on the User Product (for example, the work has been installed in ROM).

The requirement to provide Installation Information does not include a requirement to continue to provide support service, warranty, or updates for a work that has been modified or installed by the recipient, or for the User Product in which it has been modified or installed. Access to a network may be denied when the modification itself materially and adversely affects the operation of the network or violates the rules and protocols for communication across the network.

Corresponding Source conveyed, and Installation Information provided, in accord with this section must be in a format that is publicly documented (and with an implementation available to the public in source code form), and must require no special password or key for unpacking, reading or copying.
7. Additional Terms.

'Additional permissions' are terms that supplement the terms of this License by making exceptions from one or more of its conditions. Additional permissions that are applicable to the entire Program shall be treated as though they were included in this License, to the extent that they are valid under applicable law. If additional permissions apply only to part of the Program, that part may be used separately under those permissions, but the entire Program remains governed by this License without regard to the additional permissions.

When you convey a copy of a covered work, you may at your option remove any additional permissions from that copy, or from any part of it. (Additional permissions may be written to require their own removal in certain cases when you modify the work.) You may place additional permissions on material, added by you to a covered work, for which you have or can give appropriate copyright permission.

Notwithstanding any other provision of this License, for material you add to a covered work, you may (if authorized by the copyright holders of that material) supplement the terms of this License with terms:

  a) Disclaiming warranty or limiting liability differently from the terms of sections 15 and 16 of this License; or
b) Requiring preservation of specified reasonable legal notices or author attributions in that material or in the Appropriate Legal Notices displayed by works containing it; or
c) Prohibiting misrepresentation of the origin of that material, or requiring that modified versions of such material be marked in reasonable ways as different from the original version; or
d) Limiting the use for publicity purposes of names of licensors or authors of the material; or
e) Declining to grant rights under trademark law for use of some trade names, trademarks, or service marks; or
f) Requiring indemnification of licensors and authors of that material by anyone who conveys the material (or modified versions of it) with contractual assumptions of liability to the recipient, for any liability that these contractual assumptions directly impose on those licensors and authors.

All other non-permissive additional terms are considered 'further restrictions' within the meaning of section 10. If the Program as you received it, or any part of it, contains a notice stating that it is governed by this License along with a term that is a further restriction, you may remove that term. If a license document contains a further restriction but permits relicensing or conveying under this License, you may add to a covered work material governed by the terms of that license document, provided that the further restriction does not survive such relicensing or conveying.

If you add terms to a covered work in accord with this section, you must place, in the relevant source files, a statement of the additional terms that apply to those files, or a notice indicating where to find the applicable terms.

Additional terms, permissive or non-permissive, may be stated in the form of a separately written license, or stated as exceptions; the above requirements apply either way.
8. Termination.

You may not propagate or modify a covered work except as expressly provided under this License. Any attempt otherwise to propagate or modify it is void, and will automatically terminate your rights under this License (including any patent licenses granted under the third paragraph of section 11).

However, if you cease all violation of this License, then your license from a particular copyright holder is reinstated (a) provisionally, unless and until the copyright holder explicitly and finally terminates your license, and (b) permanently, if the copyright holder fails to notify you of the violation by some reasonable means prior to 60 days after the cessation.

Moreover, your license from a particular copyright holder is reinstated permanently if the copyright holder notifies you of the violation by some reasonable means, this is the first time you have received notice of violation of this License (for any work) from that copyright holder, and you cure the violation prior to 30 days after your receipt of the notice.

Termination of your rights under this section does not terminate the licenses of parties who have received copies or rights from you under this License. If your rights have been terminated and not permanently reinstated, you do not qualify to receive new licenses for the same material under section 10.
9. Acceptance Not Required for Having Copies.

You are not required to accept this License in order to receive or run a copy of the Program. Ancillary propagation of a covered work occurring solely as a consequence of using peer-to-peer transmission to receive a copy likewise does not require acceptance. However, nothing other than this License grants you permission to propagate or modify any covered work. These actions infringe copyright if you do not accept this License. Therefore, by modifying or propagating a covered work, you indicate your acceptance of this License to do so.
10. Automatic Licensing of Downstream Recipients.

Each time you convey a covered work, the recipient automatically receives a license from the original licensors, to run, modify and propagate that work, subject to this License. You are not responsible for enforcing compliance by third parties with this License.

An 'entity transaction' is a transaction transferring control of an organization, or substantially all assets of one, or subdividing an organization, or merging organizations. If propagation of a covered work results from an entity transaction, each party to that transaction who receives a copy of the work also receives whatever licenses to the work the party's predecessor in interest had or could give under the previous paragraph, plus a right to possession of the Corresponding Source of the work from the predecessor in interest, if the predecessor has it or can get it with reasonable efforts.

You may not impose any further restrictions on the exercise of the rights granted or affirmed under this License. For example, you may not impose a license fee, royalty, or other charge for exercise of rights granted under this License, and you may not initiate litigation (including a cross-claim or counterclaim in a lawsuit) alleging that any patent claim is infringed by making, using, selling, offering for sale, or importing the Program or any portion of it.
11. Patents.

A 'contributor' is a copyright holder who authorizes use under this License of the Program or a work on which the Program is based. The work thus licensed is called the contributor's 'contributor version'.

A contributor's 'essential patent claims' are all patent claims owned or controlled by the contributor, whether already acquired or hereafter acquired, that would be infringed by some manner, permitted by this License, of making, using, or selling its contributor version, but do not include claims that would be infringed only as a consequence of further modification of the contributor version. For purposes of this definition, 'control' includes the right to grant patent sublicenses in a manner consistent with the requirements of this License.

Each contributor grants you a non-exclusive, worldwide, royalty-free patent license under the contributor's essential patent claims, to make, use, sell, offer for sale, import and otherwise run, modify and propagate the contents of its contributor version.

In the following three paragraphs, a 'patent license' is any express agreement or commitment, however denominated, not to enforce a patent (such as an express permission to practice a patent or covenant not to sue for patent infringement). To 'grant' such a patent license to a party means to make such an agreement or commitment not to enforce a patent against the party.

If you convey a covered work, knowingly relying on a patent license, and the Corresponding Source of the work is not available for anyone to copy, free of charge and under the terms of this License, through a publicly available network server or other readily accessible means, then you must either (1) cause the Corresponding Source to be so available, or (2) arrange to deprive yourself of the benefit of the patent license for this particular work, or (3) arrange, in a manner consistent with the requirements of this License, to extend the patent license to downstream recipients. 'Knowingly relying' means you have actual knowledge that, but for the patent license, your conveying the covered work in a country, or your recipient's use of the covered work in a country, would infringe one or more identifiable patents in that country that you have reason to believe are valid.

If, pursuant to or in connection with a single transaction or arrangement, you convey, or propagate by procuring conveyance of, a covered work, and grant a patent license to some of the parties receiving the covered work authorizing them to use, propagate, modify or convey a specific copy of the covered work, then the patent license you grant is automatically extended to all recipients of the covered work and works based on it.

A patent license is 'discriminatory' if it does not include within the scope of its coverage, prohibits the exercise of, or is conditioned on the non-exercise of one or more of the rights that are specifically granted under this License. You may not convey a covered work if you are a party to an arrangement with a third party that is in the business of distributing software, under which you make payment to the third party based on the extent of your activity of conveying the work, and under which the third party grants, to any of the parties who would receive the covered work from you, a discriminatory patent license (a) in connection with copies of the covered work conveyed by you (or copies made from those copies), or (b) primarily for and in connection with specific products or compilations that contain the covered work, unless you entered into that arrangement, or that patent license was granted, prior to 28 March 2007.

Nothing in this License shall be construed as excluding or limiting any implied license or other defenses to infringement that may otherwise be available to you under applicable patent law.
12. No Surrender of Others' Freedom.

If conditions are imposed on you (whether by court order, agreement or otherwise) that contradict the conditions of this License, they do not excuse you from the conditions of this License. If you cannot convey a covered work so as to satisfy simultaneously your obligations under this License and any other pertinent obligations, then as a consequence you may not convey it at all. For example, if you agree to terms that obligate you to collect a royalty for further conveying from those to whom you convey the Program, the only way you could satisfy both those terms and this License would be to refrain entirely from conveying the Program.
13. Use with the GNU Affero General Public License.

Notwithstanding any other provision of this License, you have permission to link or combine any covered work with a work licensed under version 3 of the GNU Affero General Public License into a single combined work, and to convey the resulting work. The terms of this License will continue to apply to the part which is the covered work, but the special requirements of the GNU Affero General Public License, section 13, concerning interaction through a network will apply to the combination as such.
14. Revised Versions of this License.

The Free Software Foundation may publish revised and/or new versions of the GNU General Public License from time to time. Such new versions will be similar in spirit to the present version, but may differ in detail to address new problems or concerns.

Each version is given a distinguishing version number. If the Program specifies that a certain numbered version of the GNU General Public License 'or any later version' applies to it, you have the option of following the terms and conditions either of that numbered version or of any later version published by the Free Software Foundation. If the Program does not specify a version number of the GNU General Public License, you may choose any version ever published by the Free Software Foundation.

If the Program specifies that a proxy can decide which future versions of the GNU General Public License can be used, that proxy's public statement of acceptance of a version permanently authorizes you to choose that version for the Program.

Later license versions may give you additional or different permissions. However, no additional obligations are imposed on any author or copyright holder as a result of your choosing to follow a later version.
15. Disclaimer of Warranty.

THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM 'AS IS' WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
16. Limitation of Liability.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
17. Interpretation of Sections 15 and 16.

If the disclaimer of warranty and limitation of liability provided above cannot be given local legal effect according to their terms, reviewing courts shall apply local law that most closely approximates an absolute waiver of all civil liability in connection with the Program, unless a warranty or assumption of liability accompanies a copy of the Program in return for a fee.")
}
