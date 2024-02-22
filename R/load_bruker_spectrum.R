#' @title Load single Bruker Spectrum
#' @description Loads a single Bruker spectrum file and returns the spectrum data in ppm.
#' @param filepath The path of the file to be loaded.
#' @param processing_value The processing value for the file.
#' @param scale_factor A vector of two elements to scale the x and y axis values. Default is c(1, 1).
#' @return A list containing the spectrum data.
#' @examples
#' \dontrun{
#'   spectrum_data <- load_bruker_spectrum("spectrum.bruker", 1, c(2, 2))
#'   print(spectrum_data)
#' }
#' @noRd
load_bruker_spectrum <- function(filepath, processing_value, scale_factor = c(1, 1)) {
  # Get file paths of all necessary documents
  acqus_file <- file.path(filepath, "acqus")
  spec_file <- file.path(filepath, paste("pdata", processing_value, "1r", sep = "/"))
  procs_file <- file.path(filepath, paste("pdata", processing_value, "procs", sep = "/"))

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
  factor_x <- scale_factor[1]
  factor_y <- scale_factor[2]

  # Calculate length of spectrum
  spectrum_length <- length(spectrum_y)

  # Calculate x axis
  spectrum_x_ppm <- seq(ppm_highest_value, ppm_highest_value - ppm_range, by = -(ppm_range / (spectrum_length - 1)))

  # Calculate spectrum_x and recalculate spectum_y
  spectrum_x <- seq((spectrum_length - 1) / factor_x, 0, -0.001)
  spectrum_y <- spectrum_y / factor_y

  # Return the spectrum data
  return(list(spectrum_x = spectrum_x,
              spectrum_y = spectrum_y,
              spectrum_x_ppm = spectrum_x_ppm,
              spectrum_length = spectrum_length,
              ppm_range = ppm_range,
              ppm_highest_value = ppm_highest_value,
              ppm_lowest_value = ppm_lowest_value))
}
