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
