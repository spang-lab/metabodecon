# Main #####

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
#' @noRd
generate_lorentz_curves_v2 <- function(data_path,
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
