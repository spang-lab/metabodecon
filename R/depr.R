# MetaboDecon1D API (Public) #####

#' @export
#'
#' @title Deconvolute 1D NMR spectrum
#'
#' @description
#' Automatic deconvolution of a 1D NMR spectrum into several Lorentz curves and
#' the integration of them. The NMR file needs to be in Bruker format or
#' jcamp-dx format.
#'
#' This function has been deprecated with metabodecon version v1.2.0 and will be
#' removed with version 2.0.0. Please use [deconvolute()] instead.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param filepath
#' Complete path of the file folder (Notice for Bruker format: filepath needs to
#' be the spectrum folder containing one or more different spectra
#' (e.g."C:/Users/Username/Desktop/spectra_from_bruker"))
#'
#' @param filename
#' Name of the NMR file. (Notice for Bruker format: filename need to be the name
#' of your spectrum which is also the name of the folder) (Default: filename =
#' NA to analyze more spectra at once)
#'
#' @param file_format
#' Format (bruker or jcampdx) of the NMR file. (Default: file_format = "bruker")
#'
#' @param number_iterations
#' Number of iterations for the approximation of the parameters for the Lorentz
#' curves (Default: number_iterations=10)
#'
#' @param range_water_signal_ppm
#' Half width of the water artefact in ppm (Default:
#' range_water_signal=0.1527692 (e.g. for urine NMR spectra))
#'
#' @param signal_free_region
#' Row vector with two entries consisting of the ppm positions for the left and
#' right border of the signal free region of the spectrum. (Default:
#' signal_free_region=c(11.44494, -1.8828))
#'
#' @param smoothing_param
#' Row vector with two entries consisting of the number of smoothing repeats for
#' the whole spectrum and the number of data points (uneven) for the mean
#' calculation (Default: smoothing_param=c(2,5))
#'
#' @param delta
#' Defines the threshold value to distinguish between signal and noise (Default:
#' delta=6.4)
#'
#' @param scale_factor
#' Row vector with two entries consisting of the factor to scale the x-axis and
#' the factor to scale the y-axis (Default: scale_factor=c(1000,1000000))
#'
#' @param debug
#' Logical value to activate the debug mode (Default: debug=FALSE)
#'
#' @param store_results
#' Specifies whether the lorentz curve parameters `A`, `lambda` and `x_0` and
#' the approximated spectrum should be stored on disk (in addition to returning
#' them). If `store_results` is `NULL` (default), the user is asked
#' interactively where the files should be stored. If FALSE, the results are not
#' stored. If TRUE, the results are stored in a subdirectory of R's per-session
#' temporary directory.
#'
#' @return
#' A `decon0` object as described in [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
#' @seealso
#' [calculate_lorentz_curves()],
#' [plot_triplets()],
#' [plot_lorentz_curves_save_as_png()],
#' [plot_spectrum_superposition_save_as_png()]
#'
#' @references
#' Haeckl, M.; Tauber, P.; Schweda, F.; Zacharias, H.U.; Altenbuchinger, M.;
#' Oefner, P.J.; Gronwald, W. An R-Package for the Deconvolution and Integration
#' of 1D NMR Data: MetaboDecon1D. Metabolites 2021, 11, 452.
#' https://doi.org/10.3390/metabo11070452
#'
#' @author
#' 2020-2021 Martina Haeckl: initial version.\cr
#' 2024-2025 Tobias Schmidt: Minor updates to pass CRAN checks. Parameters
#' `debug` and `store_results` added.
#'
#' @examples
#'
#' ## ATTENTION: using MetaboDecon1D() for deconvolution is deprecated. Please use
#' ## deconvolute() instead.
#'
#' ## The following example shows how a subset of the Sim dataset, consisting
#' ## of two spectrum objects, can be deconvoluted using `MetaboDecon1D()`. The
#' ## whole example code is wrapped into `evalwith()` to simulate user input.
#' ## When using the function interactively, you should type in the answers to
#' ## the questions manually.
#' expected_answers <- c(
#'      "10",   # Subfolder of your filepath, i.e. the experiment number?
#'      "10",   # Subsubsubfolder of filepath, i.e. the processing number?
#'      "y",    # Use same parameters for all spectra?
#'      "1",    # File to adjust all parameters.
#'      "n",    # Signal free region borders correct selected?
#'      "3.55", # Left border.
#'      "3.35", # Right border.
#'      "y",    # Signal free region borders correct selected?
#'      "n",    # Water artefact fully inside red vertical lines?
#'      "0",    # Half width range (in ppm) for the water artefact.
#'      "y",    # Water artefact fully inside red vertical lines?
#'      "n"     # Save results as text documents?
#' )
#' sim <- metabodecon_file("bruker/sim_subset")
#' evalwith(answers = expected_answers, {
#'      sim_decon <- MetaboDecon1D(sim)
#' })
#'
#' ## Deconvolute only the first spectrum of the folder "bruker/sim_subset" into
#' evalwith(answers = expected_answers[-(3:4)], {
#'      sim_decon <- MetaboDecon1D(sim, filename = "sim_01")
#' })
#'
MetaboDecon1D <- function(filepath,
                          filename = NA,
                          file_format = "bruker",
                          number_iterations = 10,
                          range_water_signal_ppm = 0.1527692,
                          signal_free_region = c(11.44494, -1.8828),
                          smoothing_param = c(2, 5),
                          delta = 6.4,
                          scale_factor = c(1000, 1000000),
                          debug = FALSE,
                          store_results = NULL) {

    stopifnot(
        is_str(filepath) && dir.exists(filepath),
        is.na(filename) || (is_str(filename) && file.exists(file.path(filepath, filename))),
        file_format %in% c("bruker", "jcampdx"),
        is_int(number_iterations, n = 1),
        is_num(range_water_signal_ppm, n = 1),
        is_num(signal_free_region, n = 2),
        is_num(smoothing_param, n = 2),
        is_num(delta, n = 1),
        is_num(scale_factor, n = 2),
        is_bool(debug, 1),
        is_bool(store_results, 1) || is.null(store_results)
    )

    if (debug) logf("Starting MetaboDecon1D")

    example <- FALSE

    local_dir(filepath)

    # If filename is, as the default value, NA, then the all spectra of the
    # filepath of the folder are analyzed
    if (is.na(filename)) {
        # Check which file_format is present
        if (file_format == "jcampdx") {
            files_list <- list.files(filepath)
            files <- c() # Save only files with .dx format
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

            # Get some values from the user
            spectroscopy_value <- readline(prompt = "What is the name of the subfolder of your filepath: \n[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10] ")
            processing_value <- readline(prompt = "What is the name of the subsubsubfolder of your filepath: \n[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10/pdata/10] ")

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
            # Get some values from the user
            spectroscopy_value <- readline(prompt = "What is the name of the subfolder of your filepath: \n[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10] ")
            processing_value <- readline(prompt = "What is the name of the subsubsubfolder of your filepath: \n[e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10/pdata/10] ")

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

#' @export
#'
#' @title Calculate lorentz curves for each analyzed spectrum
#'
#' @description
#' Helper function of [plot_lorentz_curves_save_as_png()].
#' Should not be called directly by the user.
#'
#' Calculates the lorentz curves of each investigated spectrum.
#'
#' This function has been deprecated with metabodecon version v1.4.3 and will be
#' removed with version 2.0.0.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param deconv_result
#' A list as returned by [generate_lorentz_curves()] or [MetaboDecon1D].
#'
#' @param number_of_files Number of spectra to analyze
#'
#' @return
#' If `deconv_result` holds the result of a single deconvolution, a matrix
#' containing the generated Lorentz curves is returned, where each row depicts
#' one Lorentz curve. If `deconv_result` is a list of deconvoluted spectra, a
#' list of such matrices is returned.
#'
#' @seealso
#' [MetaboDecon1D()],
#' [plot_triplets()],
#' [plot_lorentz_curves_save_as_png()],
#' [plot_spectrum_superposition_save_as_png()]
#'
#' @author
#' 2020-2021 Martina Haeckl: initial version.\cr
#' 2024-2025 Tobias Schmidt: Minor updates to pass CRAN checks
#'
#' @examples
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' ## Deconvolute the spectra in folder "bruker/sim_subset" into a list of
#' ## Lorentz Curves (specified via the parameters A, lambda and x_0).
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' decons <- generate_lorentz_curves_sim(sim[1:2])
#' decon0 <- decons[[1]]
#'
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' ## Calculate the corresponding y values at each ppm value for each Lorentz
#' ## Curve. I.e. you get a matrix of dimension n x m for each deconvolution,
#' ## where n is the number of Lorentz Curves and m is the number of ppm values.
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' yy <- calculate_lorentz_curves(decons)
#' y1 <- yy[[1]]
#' dim(y1)
#'
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' ## Visualize the 5th, 9th and 11th Lorentz curve.
#' ## -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#' nrs <- c(5, 9, 11)
#' col <- c("red", "blue", "orange")
#' desc <- paste("Lorentz curve", nrs)
#' plot(decon0$x_values_ppm, decon0$y_values, type = "l", lty = 2)
#' for (i in 1:3) lines(decon0$x_values_ppm, y1[nrs[i], ], col = col[i])
#' legend("topright", legend = desc, col = col, lty = 1)
#'
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

#' @export
#' @title Plot peak triplets for variable range
#' @description Plots the peak triplets for each peak detected by [MetaboDecon1D()] and stores the plots as png at `outdir`.
#'
#' Superseded by [plot_spectrum()] since metabodecon v1.2.0. Will be replaced with v2.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param deconv_result
#' Saved result of the MetaboDecon1D() function
#'
#' @param x_range
#' Row vector with two entries consisting of the ppm start and the ppm end value
#' to scale the range of the x-axis (optional)
#'
#' @param y_range
#' Row vector with two entries consisting of the ppm start and the ppm end value
#' to scale the range of the y-axis (optional)
#'
#' @param out_dir
#' Directory to save the png files (optional)
#'
#' @param ask
#' Logical value to ask the user if the png files should be saved in the
#' specified directory (optional)
#'
#' @return No return value, called for side effect of plotting.
#'
#' @seealso
#' [MetaboDecon1D()], [calculate_lorentz_curves()], [plot_lorentz_curves_save_as_png()],
#' [plot_spectrum_superposition_save_as_png()]
#'
#' @author
#' 2020-2021 Martina Haeckl: initial version.\cr
#' 2024-2025 Tobias Schmidt: Minor updates to pass CRAN checks.
#'
#' @examples
#' sim <- metabodecon_file("bruker/sim_subset")
#' sim_decon <- generate_lorentz_curves_sim(sim)
#' png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
#' plot_triplets(sim_decon, out_dir = png_dir, ask = FALSE)
#' dir(png_dir, full.names = TRUE)
plot_triplets <- function(deconv_result, x_range = c(), y_range = c(), out_dir = ".", ask = TRUE) {
    out_dir <- norm_path(out_dir)
    if (ask) {
        continue <- get_yn_input(sprintf("Continue creating pngs in %s?", out_dir))
        if (!continue) {
            invisible(NULL)
        }
    }
    local_dir(out_dir)

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

            message(paste("Plot triplets of", name))

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

                message(paste("Plot triplets of", name))

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

        message(paste("Plot triplets of", name))

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
#'
#' @title Plot lorentz curves for variable range
#'
#' @description
#' Plots the original spectrum and all generated Lorentz curves and save the
#' result as png under the filepath.
#'
#' Superseded by [plot_spectrum()] since metabodecon v1.2.0. Will be replaced
#' with metabodecon v2.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @param deconv_result
#' Saved result of the MetaboDecon1D() function
#'
#' @param x_range
#' Row vector with two entries consisting of the ppm start and the ppm end value
#' to scale the range of the x-axis (optional)
#'
#' @param y_range
#' Row vector with two entries consisting of the ppm start and the ppm end value
#' to scale the range of the y-axis (optional)
#'
#' @param out_dir
#' Path to the directory where the png files should be saved. Default is the
#' current working directory.
#'
#' @param ask
#' Logical value. Whether to ask for confirmation from the user before writing
#' files to disk. Default is TRUE.
#'
#' @return
#' NULL, called for side effects.
#'
#' @seealso
#' [MetaboDecon1D()], [plot_triplets()], [plot_spectrum_superposition_save_as_png()]
#'
#' @author
#' 2020-2021 Martina Haeckl: initial version.\cr
#' 2024-2025 Tobias Schmidt: Minor updates to pass CRAN checks
#'
#' @examples
#' sim <- metabodecon_file("bruker/sim_subset")
#' sim_decon <- generate_lorentz_curves_sim(sim)
#' png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
#' plot_lorentz_curves_save_as_png(sim_decon, out_dir = png_dir, ask = FALSE)
#' dir(png_dir, full.names = TRUE)
plot_lorentz_curves_save_as_png <- function(deconv_result, x_range = c(), y_range = c(), out_dir = ".", ask = TRUE) {
    out_dir <- norm_path(out_dir)
    if (ask) {
        continue <- get_yn_input(sprintf("Continue creating pngs in %s?", out_dir))
        if (!continue) {
            invisible(NULL)
        }
    }
    local_dir(out_dir)

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
#'
#' @title Plot spectrum approx for variable range
#'
#' @description
#' Plots the original spectrum and the superposition of all generated Lorentz
#' curves and saves the result as png under the specified filepath.
#'
#' Superseded by [plot_spectrum()] since metabodecon v1.2.0. Will be replaced with v2.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @author 2020-2021 Martina Haeckl: initial version.
#'
#' @param deconv_result
#' Saved result of the MetaboDecon1D() function
#'
#' @param x_range
#' Row vector with two entries consisting of the ppm start and the ppm end value
#' to scale the range of the x-axis (optional)
#'
#' @param y_range
#' Row vector with two entries consisting of the ppm start and the ppm end value
#' to scale the range of the y-axis (optional)
#'
#' @param out_dir
#' Path to the directory where the png files should be saved. Default is the
#' current working directory.
#'
#' @param ask
#' Logical value. Whether to ask for confirmation from the user before writing
#' files to disk.
#'
#' @return NULL, called for side effects.
#'
#' @seealso
#' [MetaboDecon1D()],
#' [calculate_lorentz_curves()],
#' [plot_triplets()],
#' [plot_lorentz_curves_save_as_png()]
#'
#' @examples
#' sim <- metabodecon_file("bruker/sim_subset")
#' sim_decon <- generate_lorentz_curves_sim(sim)
#' png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
#' plot_spectrum_superposition_save_as_png(sim_decon, out_dir = png_dir, ask = FALSE)
#' dir(png_dir, full.names = TRUE)
plot_spectrum_superposition_save_as_png <- function(deconv_result,
                                                    x_range = c(),
                                                    y_range = c(),
                                                    out_dir = ".",
                                                    ask = TRUE) {
    out_dir <- norm_path(out_dir)
    if (ask) {
        prompt <- sprintf("Continue creating pngs in %s?", out_dir)
        continue <- get_yn_input(prompt)
        if (!continue) return(invisible(NULL)) # styler: off
    }
    local_dir(out_dir)

    number_in_folder <- 0 # Number of analyzed spectra
    if ("number_of_files" %in% names(deconv_result)) {
        number_of_files <- deconv_result$number_of_files
    } else {
        if ("number_of_files" %in% names(deconv_result[[1]])) {
            number_of_files <- deconv_result[[1]]$number_of_files
            if (number_of_files == 1) number_in_folder <- 1
        }
    }

    # Check how many spectra are investigated
    if (number_of_files > 1 || number_in_folder == 1) {
        # Check if user input comprises a $ sign
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
        x_range <- if (is.null(x_range)) rev(range(spectrum_x_ppm)) else x_range
        message(paste("Plot superposition of", name))
        filename <- paste(name, "_sum_lorentz_curves.png", sep = "")
        grDevices::png(file = filename, width = 825, height = 525)
        plot(
            spectrum_x_ppm, spectrum_y,
            main = name, xlab = "[ppm]", ylab = "Intensity [a.u.]",
            type = "l", cex = 1.3,
            xlim = x_range,
            ylim = y_range
        )
        graphics::lines(spectrum_x_ppm, spectrum_approx, col = "red")
        graphics::legend("topright", legend = c("Original spectrum", "Sum of Lorentz curves"), col = c("black", "red"), lty = 1, bty = "n", cex = 1.3)
        text <- paste("MSE_Normed = ", mse, sep = "")
        graphics::mtext(text, side = 3)
        grDevices::dev.off()
    }
}

# MetaboDecon1D Helpers (Private) #####

#' @noRd
#' @author 2020-2021 Martina Haeckl: initial version.
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
            debuglist$data <- list(
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
            debuglist$data <- list(
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
        debuglist$wsrm <- list(
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
        debuglist$nvrm <- list(spectrum_y = spectrum_y)
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
        debuglist$smooth <- list(spectrum_y = spectrum_y)
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
        debuglist$peaksel <- list(
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
        debuglist$peakfilter <- list( # peaks without borders removed
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
        debuglist$peakscore <- list(
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
        debuglist$parinit <- list(
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
        debuglist$parapprox <- list(
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
        utils::write.table(spectrum_info, name_info_txt, sep = ",", col.names = FALSE, append = FALSE)
        message(sprintf("Writing %s", norm_path(name_output_txt)))
        utils::write.table(spectrum_output, name_output_txt, sep = ",", col.names = FALSE, append = FALSE)
    } else {
        message("Skipping saving of results.")
    }

    if (debug) {
        logf("Reached point where parameters were saved to txt documents previously")
        debuglist$parsave <- list(
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

# Deconvolution (Public) #####

#' @export
#'
#' @title Deconvolute one or more NMR spectra
#'
#' @description
#' Deconvolutes NMR spectra by modeling each detected signal within a spectrum
#' as Lorentz Curve.
#'
#' This function has been deprecated with metabodecon version v1.4.3 and will be
#' removed with version 2.0.0. Please use [deconvolute()] instead.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' @inheritParams read_spectrum
#' @inheritParams deconvolute
#'
#' @param make_rds Logical or character. If TRUE, stores results as an RDS file
#' on disk. If a character string, saves the RDS file with the specified name.
#' Should be set to TRUE if many spectra are evaluated to decrease computation
#' time.
#'
#' @return
#' A 'decon1' object if a single spectrum was provided. A 'decons1' object if
#' multiple spectra were provided. See [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html) for
#' details.
#'
#' @details
#'
#' First, an automated curvature based signal selection is performed. Each
#' signal is represented by 3 data points to allow the determination of initial
#' Lorentz curves. These Lorentz curves are then iteratively adjusted to
#' optimally approximate the measured spectrum.
#'
#' [generate_lorentz_curves_sim()] is identical to [generate_lorentz_curves()]
#' except for the defaults, which are optimized for deconvoluting the 'sim'
#' dataset, shipped with 'metabodecon'. The 'sim' dataset is a simulated
#' dataset, which is much smaller than a real NMR spectra and lacks a water
#' signal. This makes it ideal for use in examples. However, the default values
#' for `sfr`, `wshw`, and `delta` in the "normal" [generate_lorentz_curves()]
#' function are not optimal for this dataset. To avoid having to define the
#' optimal parameters repeatedly in examples, this function is provided to
#' deconvolute the "Sim" dataset with suitable parameters.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#'
#' ## Define the paths to the example datasets we want to deconvolute:
#' ## `sim_dir`: directory containing 16 simulated spectra
#' ## `sim_01`: path to the first spectrum in the `sim` directory
#' ## `sim_01_spec`: the first spectrum in the `sim` directory as a dataframe
#'
#' sim_dir        <- metabodecon_file("sim_subset")
#' sim_1_dir      <- file.path(sim_dir, "sim_01")
#' sim_2_dir      <- file.path(sim_dir, "sim_02")
#' sim_1_spectrum <- read_spectrum(sim_1_dir)
#' sim_2_spectrum <- read_spectrum(sim_2_dir)
#' sim_spectra    <- read_spectra(sim_dir)
#'
#'
#' ## Show that `generate_lorentz_curves()` and `generate_lorentz_curves_sim()`
#' ## produce the same results:
#'
#' sim_1_decon0 <- generate_lorentz_curves(
#'     data_path = sim_1_dir, # Path to directory containing spectra
#'     sfr = c(3.55, 3.35),   # Borders of signal free region (SFR) in ppm
#'     wshw = 0,              # Half width of water signal (WS) in ppm
#'     ask = FALSE,           # Don't ask for user input
#'     verbose = FALSE        # Suppress status messages
#' )
#' sim_1_decon1 <- generate_lorentz_curves_sim(sim_1_dir)
#' stopifnot(all.equal(sim_1_decon0, sim_1_decon1))
#'
#'
#' ## Show that passing a spectrum produces the same results as passing the
#' ## the corresponding directory:
#'
#' decon_from_spectrum_dir <- generate_lorentz_curves_sim(sim_1_dir)
#' decon_from_spectrum_obj <- generate_lorentz_curves_sim(sim_1_spectrum)
#' decons_from_spectra_obj <- generate_lorentz_curves_sim(sim_spectra)
#' decons_from_spectra_dir <- generate_lorentz_curves_sim(sim_dir)
#'
#' most.equal <- function(x1, x2) {
#'     ignore <- which(names(x1) %in% c("number_of_files", "filename"))
#'     equal <- all.equal(x1[-ignore], x2[-ignore])
#'     invisible(stopifnot(isTRUE(equal)))
#' }
#'
#' all.equal(  decon_from_spectrum_dir, decon_from_spectrum_obj     )
#' all.equal(  decons_from_spectra_dir, decons_from_spectra_obj     )
#' most.equal( decon_from_spectrum_dir, decons_from_spectra_obj[[1]])
#' most.equal( decon_from_spectrum_dir, decons_from_spectra_dir[[1]])
generate_lorentz_curves <- function(data_path,
    file_format="bruker", make_rds=FALSE, expno=10, procno=10, raw=TRUE,
    nfit=10, smopts=c(2,5), delta=6.4, sfr=c(11.44494,-1.8828), wshw=0.1527692,
    ask=TRUE, force=FALSE, verbose=TRUE, nworkers=1
) {
    # Check inputs
    stopifnot(
        is_existing_path(data_path) ||
        is_spectrum(data_path) || is_spectra(data_path),
        is_char(file_format,1,"(bruker|jcampdx)"),
        is_bool(make_rds,1) || is_char(make_rds,1),
        is_int(expno,1),    is_int(procno,1),   is_int(nfit,1),
        is_int(smopts,2),   is_num(delta,1),    is_num(sfr,2),
        is_num(wshw,1),     is_bool(ask,1),     is_bool(force,1),
        is_bool(verbose,1), is_int(nworkers,1)
    )

    # Read spectra
    spectra <- as_spectra(
        data_path, file_format, expno, procno, raw,
        silent = !verbose, force = force
    )

    # Deconvolute
    decons1 <- deconvolute_spectra(spectra,
        nfit, smopts, delta, sfr, wshw,
        ask, force, verbose, bwc=1,
        use_rust=FALSE, nw=nworkers, igr=list(), rtyp="decon1"
    )

    # Store and return
    store_as_rds(decons1, make_rds, data_path)
    if (length(decons1) == 1) decons1[[1]] else decons1
}

#' @export
#' @rdname generate_lorentz_curves
generate_lorentz_curves_sim <- function(data_path,
    file_format="bruker", make_rds=FALSE, expno=10, procno=10, raw=TRUE,
    nfit=10, smopts=c(2,5), delta=6.4, sfr=c(3.55,3.35), wshw=0,
    ask=FALSE, force=FALSE, verbose=FALSE, nworkers=1
) {
    generate_lorentz_curves(
        data_path, file_format, make_rds, expno, procno, raw,
        nfit, smopts, delta, sfr, wshw,
        ask, force, verbose, nworkers
    )
}

# Plotting (Private) #####

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
plot_si_mat <- function(X = generate_lorentz_curves_sim("bruker/sim"),
                        lgdcex = "auto",
                        main = NULL,
                        mar = par("mar")) {
    n <- nrow(X) # Number of Spectra
    p <- ncol(X) # Number of Datapoints
    cols <- rainbow(n) # Colors
    dpis <- seq_len(p) # Datapoint Indices
    spis <- seq_len(n) # Spectrum Indices
    ltxt <- paste("Spectrum", spis)
    xlab <- "Datapoint Number"
    ylab <- "Signal Intensity"
    xlim <- c(1, p)
    ylim <- c(0, max(X))
    args <- named(x = NA, type = "n", xlim, ylim, xlab, ylab)
    local_par(mar = mar)
    do.call(plot, args)
    for (i in spis) lines(x = dpis, y = X[i, ], col = cols[i])
    if (lgdcex == "auto") lgdcex <- 1 / max(1, log(n, base = 8))
    legend(x = "topright", legend = ltxt, col = cols, lty = 1, cex = lgdcex)
    if (!is.null(main)) title(main = main)
}

#' @noRd
#' @title Plots a spectrum simulated with [get_sim_params()]
#' @param simspec A simulated spectrum as returned by [get_sim_params()].
#' @author 2024-2025 Tobias Schmidt: initial version.
#' @examples
#' simspec <- get_sim_params()
#' plot_sim_spec(simspec)
plot_sim_spec <- function(simspec = get_sim_params()) {
    X <- simspec$X
    top_ticks <- seq(from = min(X$cs), to = max(X$cs), length.out = 5)
    top_labels <- round(seq(from = min(X$fq), to = max(X$fq), length.out = 5))
    line1_text <- sprintf("Name: %s", gsub("blood", "sim", simspec$filename))
    line2_text <- sprintf("Base: %s ", simspec$filename)
    line3_text <- sprintf("Range: 3.6-3.4 ppm")
    plot(
        x = X$cs, y = X$si_raw * 1e-6, type = "l",
        xlab = "Chemical Shift [PPM]", ylab = "Signal Intensity [AU]",
        xlim = c(3.6, 3.4), col = "black"
    )
    lines(x = X$cs, y = X$si_smooth, col = "blue")
    lines(x = X$cs, y = X$si_sim, col = "red")
    legend(
        x = "topleft", lty = 1,
        col = c("black", "blue", "red"),
        legend = c("Original Raw / 1e6", "Original Smoothed", "Simulated")
    )
    axis(3, at = top_ticks, labels = top_labels, cex.axis = 0.75)
    mtext(sprintf("Frequency [Hz]"), side = 3, line = 2)
    mtext(line1_text, side = 3, line = -1.1, col = "red", cex = 1, adj = 0.99)
    mtext(line2_text, side = 3, line = -2.1, col = "red", cex = 1, adj = 0.99)
    mtext(line3_text, side = 3, line = -3.1, col = "red", cex = 1, adj = 0.99)
}

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
plot_noise_methods <- function(siRND, siSFR, n = 300, start = 5000) {
    ymin <- min(c(min(siRND), min(siSFR)))
    ymax <- max(c(max(siRND), max(siSFR)))
    ylim <- c(ymin, ymax)
    local_par(mfrow = c(5, 1), mar = c(3, 4, 0, 1))
    for (i in 1:5) {
        redT <- rgb(1, 0, 0, 0.1)
        bluT <- rgb(0, 0, 1, 0.1)
        idx <- ((i - 1) * n + 1):(i * n) + start
        ysiRND <- siRND[idx]
        ysiSFR <- siSFR[idx]
        plot(1:n, ysiRND, type = "n", ylim = ylim, ylab = "", xlab = "", xaxt = "n")
        axis(1, at = seq(1, n, by = 50), labels = idx[seq(1, n, by = 50)])
        points(1:n, ysiRND, col = "red", pch = 20)
        points(1:n, ysiSFR, col = "blue", pch = 20)
        lines(1:n, ysiRND, col = "red")
        lines(1:n, ysiSFR, col = "blue")
        lines(1:n, rep(0, n), col = "black", lty = 2)
        polygon(c(1:n, n:1), c(ysiRND, rep(0, n)), col = redT, border = NA)
        polygon(c(1:n, n:1), c(ysiSFR, rep(0, n)), col = bluT, border = NA)
        legend("topleft", NULL, c("RND", "SFR"), col = c("red", "blue"), lty = 1)
    }
}

# Alignment (Private) #####

#' @noRd
#' @author
#' 2021-2024 Wolfram Gronwald: wrote initial version as part of (gen_feat_mat[]).\cr
#' 2024-2025 Tobias Schmidt: refactored into a separate function.
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
#' @param feat
#' Output of `gen_feat_mat()`.
#'
#' @param maxShift
#' Maximum number of points along the "ppm-axis" which a value can be moved by
#' speaq package e.g. 50. 50 is a suitable starting value for plasma spectra
#' with a digital resolution of 128K. Note that this parameter has to be
#' individually optimized depending on the type of analyzed spectra and the
#' digital resolution. For urine which is more prone to chemical shift
#' variations this value most probably has to be increased.
#'
#' @param spectrum_data
#' Output of `generate_lorentz_curves()`.
#'
#' @param si_size_real_spectrum
#' Number of real data points in your original spectra.
#'
#' @return
#' A matrix containing the aligned integral values of all spectra. Each row
#' contains the data of each spectrum and each column corresponds to one data
#' point. Each entry corresponds to the integral of a deconvoluted signal with
#' the signal center at this specific position after alignment by speaq.
#'
#' @author 2021-2024 Wolfram Gronwald: initial version.
#'
#' @examples
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

# Data Management (Private) #####

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
simulate_from_decon <- function(x,
                                # Simulation Params
                                cs_min = 3.4,
                                cs_max = 3.6,
                                pk_min = 3.45,
                                pk_max = 3.55,
                                noise_method = "SFR",
                                # Output params
                                show = TRUE,
                                save = FALSE,
                                verbose = TRUE,
                                ...) {
    if (!isFALSE(verbose)) logf("Checking function inputs")
    stopifnot(is_num(cs_min, 1))
    stopifnot(is_num(cs_max, 1))
    stopifnot(is_num(pk_min, 1))
    stopifnot(is_num(pk_max, 1))
    stopifnot(is_char(noise_method, 1, "(RND|SFR)"))
    stopifnot(is_bool(show, 1))
    stopifnot(is_bool(save, 1))
    stopifnot(is_bool(verbose, 1))
    d <- as_decon1(x)
    logv <- if (verbose) logf else function(...) NULL

    logv("Simulating spectrum from '%s' (noise method = '%s')", d$filename, noise_method)
    s <- as_spectrum(d)

    logv("Throwing away datapoints outside of %.1f to %.1f", cs_min, cs_max)
    ix <- which(s$cs >= cs_min & s$cs <= cs_max)
    s$cs <- s$cs[ix]
    s$si <- s$si[ix]
    if (!is.null(s$meta$fq)) s$meta$fq <- s$meta$fq[ix]

    logv("Keeping only within %.1f to %.1f", pk_min, pk_max)
    ip <- which(d$x_0_ppm <= pk_max & d$x_0_ppm >= pk_min)
    s$simpar <- list(
        A      = -d$A_ppm[ip],
        x_0    = d$x_0_ppm[ip],
        lambda = -d$lambda_ppm[ip]
    )

    logv("Calculating simulated signal intensities (si_sim) as superposition of lorentz curves")
    rownames(s) <- NULL
    si <- lorentz_sup(x = s$cs, lcpar = s$simpar)

    logv("Adding %s noise to simulated data", noise_method)
    if (noise_method == "RND") {
        # ToSc, 2024-08-02: noise_method 'RND' turned out to be a bad idea,
        # because the true noise is not normally distributed, but has long
        # stretches of continuous increase or decrease that would be incredibly
        # unlikely to occur with a normal distribution. This can be seen by
        # running `analyze_noise_methods()`.
        sfr_sd <- sd(d$y_values_raw[c(1:10000, 121073:131072)] * 1e-6) # (1)
        # (1) The true SFR covers approx. the first and last 20k datapoints.
        # However, to be on the safe side, only use the first and last 10k
        # datapoints for the calculation.
        noise <- rnorm(n = length(si), mean = 0, sd = sfr_sd)
        si <- si + noise
    } else {
        idx <- ix - min(ix) + 5000 # Use SI of datapoints 5000:6308 for noise
        noise <- d$y_values_raw[idx] * 1e-6
        si <- si + noise
    }

    logv("Discretize signal intensities to allow efficient storage as integers")
    s$si <- as.integer(si * 1e6) / 1e6 # (2)
    # (2) Convert to integers and back. We do this to not lose precision later
    # on when the data gets written to disc as integers.

    if (show) plot_sim_spec(s)
    if (store) save_spectrum(s, verbose = verbose, ...)

    logv("Finished spectrum simulation")
    invisible(s)
}

#' @noRd
#' @title Count Stretches of Increases and Decreases
#'
#' @description
#' Counts the lengths of consecutive increases and decreases in a numeric
#' vector.
#'
#' @param x A numeric vector.
#'
#' @return
#' A numeric vector containing the lengths of stretches of increases and
#' decreases.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' #
#' # Example Data (x)
#' # | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
#' # |-------------------------------|
#' # |   |   |   |###|   |   |   |###|
#' # |   |   |###|###|   |   |###|###|
#' # |   |###|###|###|###|   |###|###|
#' # |###|###|###|###|###|###|###|###|
#' # |-------------------------------|
#' # |   +   +   +   -   -   +   +   | + is increase
#' # |-------------------------------| - is decrease
#' #
#' x <- c(1.0, 2.2, 3.0, 4.4, 2.0, 1.0, 3.0, 4.0)
#' count_stretches(x) # Returns c(3, 2, 2) because we have + + + - - + +
count_stretches <- function(x) {
    if (length(x) < 2) {
        return(integer(0))
    }
    ss <- numeric(length(x))
    s <- 1
    inc <- x[2] > x[1]
    for (i in 3:length(x)) {
        if ((inc && x[i] > x[i - 1]) || (!inc && x[i] < x[i - 1])) {
            s <- s + 1
        } else {
            ss[i - 1] <- s
            s <- 1
            inc <- x[i] > x[i - 1]
        }
    }
    ss[i] <- s
    return(ss[ss != 0])
}

#' @noRd
#' @description
#' Used during development of `simulate_spectra()` to find a realistic method
#' for noise generation.
#' @author 2024-2025 Tobias Schmidt: initial version.
analyze_noise_methods <- function(ask = TRUE) {
    download_example_datasets()
    blood_1 <- datadir("example_datasets/bruker/blood/blood_01")
    deconv <- generate_lorentz_curves(blood_1, ask = FALSE)
    si_raw <- deconv$y_values_raw
    sd_sfr <- sd(si_raw[c(1:10000, 121073:131072)] * 1e-6)
    siRND <- rnorm(n = 10000, mean = 0, sd = sd_sfr)
    siSFR <- si_raw[1:10000] * 1e-6

    logf("Visualizing raw SIs for noise methods RND and SFR")
    plot_noise_methods(siRND, siSFR)
    if (!get_yn_input("Continue?")) {
        return()
    }

    logf("Visualizing smoothed SIs for noise methods RND and SFR")
    siRND_sm <- smooth_signals(siRND)$y_smooth
    siSFR_sm <- smooth_signals(siSFR)$y_smooth
    plot_noise_methods(siRND_sm, siSFR_sm)
    if (!get_yn_input("Continue?")) {
        return()
    }

    logf("Visualizing lengths of intervals of continuous increase and/or decrease")
    slRND <- count_stretches(siRND) # stretch lengths of siRND
    slSFR <- count_stretches(siSFR) # stretch lengths of SFR
    table(slRND)
    table(slSFR)
    local_par(mfrow = c(2, 1))
    hist(slRND, breaks = 0:20 + 0.5, xlim = c(0, 20))
    hist(slSFR, breaks = 0:20 + 0.5, xlim = c(0, 20))
}

