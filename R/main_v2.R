#' @export
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
#' @examples \dontrun{
#' spectrum_data <- generate_lorentz_curves(data_path = system.file(package = "metabodecon"), file_format = "jcampdx", filename = "urine.dx")
#' }
#' @details First, an automated curvature based signal selection is performed. Each signal is represented by 3 data points to allow the determination of initial Lorentz curves. These Lorentz curves are then iteratively adjusted to optimally approximate the measured spectrum. For each spectrum two text files will be created in the parent folder i.e. the folder given in data path. The spectrum approximated from all deconvoluted signals and a parameter file that contains all numerical values of the deconvolution. Furthermore, the numerical values of the deconvolution are also stored in a data_frame.
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

  # Checks arguments
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
    spectroscopy_value <- readline(prompt="What is the name of the subfolder of your filepath: (e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10) ")
    processing_value <- readline(prompt="What is the name of the subsubsubfolder of your filepath: (e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10/pdata/10) ")
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
  return_list <- list()
  for (i in 1:length(files)) {
    # Bruker Path Example
    #   fullpath:  C:\Users\tobi\.local\share\R\metabodecon\urine\urine_2\10\pdata\10\procs
    #   filepath:  C:\Users\tobi\.local\share\R\metabodecon\urine\urine_2\10
    #   data_path  C:\Users\tobi\.local\share\R\metabodecon\urine
    #   name:                                                     urine_2
    #   spectroscopy_value:                                               10
    #   constant:                                                            pdata    procs
    #   processing_value:                                                          10
    # Jcampdx Path Example
    #   fullpath: C:\Users\tobi\.local\share\R\metabodecon\urine\urine_2.dx
    #   filepath: C:\Users\tobi\.local\share\R\metabodecon\urine
    #   name:                                                    urine_2
    name <- files[i] # bruker: `urine_2`, jcampdx: `urine_2.dx`
    filepath <- switch(file_format, "bruker" = paste(data_path, name, spectroscopy_value, sep="/"), "jcampdx" = data_path, stop("Invalid file format"))
    return_list[[files[i]]] <- deconvolute_spectrum(filepath, name, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, current_filenumber = i, number_of_files = length(files))
    range_water_signal_ppm <- list_file$range_water_signal_ppm # Save `range_water_signal` and `signal_free_region` for next loop passage as those might have been adjusted interactively by the user
    signal_free_region <- list_file$signal_free_region
  }

  # Save results
  if (make_rds) {
    saveRDS(object = spectrum_data, file = file.path(data_path, "spectrum_data.rds"))
  }
  return(spectrum_data)
}


deconvolute_spectrum <- function(filepath, name, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, current_filenumber, number_of_files) {
  message(paste("Start deconvolution of ", name, ":", sep = ""))
  x <- deconvolution(filepath, name, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, current_filenumber)
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
    "A" = x$A,
    "lambda" = x$lambda,
    "x_0" = x$w # [17] rename
  )
  return(y)
}

######### General Utils #######################################################

#' @title Get numeric input from user
#' @description Prompts the user for input and checks if the input is a number between a minimum and maximum value. If the input is not valid, it keeps asking the user for input until they provide a valid response.
#' @param prompt The prompt to display to the user.
#' @param min The minimum valid value. Default is -Inf.
#' @param max The maximum valid value. Default is Inf.
#' @param int Whether the input should be an integer. Default is FALSE.
#' @return The user's input as a numeric value.
#' @examples \dontrun{
#' get_num_input("Enter a number between 1 and 10: ", min = 1, max = 10)
#' }
get_num_input <- function(prompt, min = -Inf, max = Inf, int = FALSE) {
  pat <- if (int) "^[0-9]+$" else "^[0-9]*\\.[0-9]+$"
  typ <- if (int) "number" else "value"
  x <- trimws(readline(prompt = prompt))
  while (!(grepl(pat, x) && as.numeric(x) >= min && as.numeric(x) <= max)) {
    message("Error. Please enter a ", typ, " between ", min, " and ", max, ".")
    x <- trimws(readline(prompt = prompt))
  }
  x <- as.numeric(x)
  return(x)
}

#' @title Get string input from user
#' @description Prompts the user for input and checks if the input is in a list of valid responses. If the input is not valid, it keeps asking the user for input until they provide a valid response.
#' @param prompt The prompt to display to the user.
#' @param valid A vector of valid responses.
#' @return The user's input.
#' @examples \dontrun{
#' get_str_input("Enter a, b or c: ", c("a", "b", "c"))
#' }
get_str_input <- function(prompt, valid) {
  x <- readline(prompt = prompt)
  n <- length(valid)
  valid_str <- if (n == 1) valid else paste(collapse(valid[-n]), "or", valid[n])
  while (!(x %in% valid)) {
    message("Error. Please type either ", valid_str, ".")
    x <- readline(prompt = prompt)
  }
  return(x)
}

#' @title Get yes/no input from user
#' @description Prompts the user for input until they enter either 'y' or no 'n'. Returns TRUE if the user entered 'y' and FALSE if they entered 'n'.
#' @param prompt The prompt to display to the user.
#' @return TRUE if the user entered 'y' and FALSE if they entered 'n'.
#' @examples \dontrun{
#' show_dir <- get_yn_input("List dir content? (y/n) ")
#' if (show_dir) print(dir())
#' }
get_yn_input <- function(prompt) {
  if (!grepl("(y/n)", prompt, fixed = TRUE)) {
    prompt <- paste0(prompt, " (y/n) ")
  }
  x <- get_str_input(prompt, c("y", "n"))
  y <- if (x == "y") TRUE else FALSE
  return(y)
}

#' @title Collapse a vector into a string
#' @description Collapses a vector into a single string, with elements separated by a specified separator. Essentially a shorthand for `paste(x, collapse = sep)`.
#' @param x A vector to collapse.
#' @param sep A string to use as the separator between elements. Default is ", ".
#' @return A single string with elements of x separated by sep.
#' @examples
#' collapse(c("a", "b", "c")) # "a, b, c"
#' collapse(1:5, sep = "-")   # "1-2-3-4-5"
collapse <- function(x, sep = ", ") {
  paste(x, collapse = ", ")
}

#' @title Check if a value is NULL or NA
#' @description Checks if a given value is NULL or NA.
#' @param x The value to check.
#' @return TRUE if x is NULL or x is NA, else FALSE.
is.none <- function(x) {
  return(is.null(x) || is.na(x))
}
