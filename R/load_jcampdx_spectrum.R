#' @title Load single JCAMP-DX Spectrum
#' @description Loads a single JCAMP-DX spectrum file and returns the spectrum data in ppm.
#' @param name The path of the file to be loaded.
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
load_jcampdx_spectrum <- function(name, scale_factor = c(1, 1)) {

    # Import data using readJDX package
    data <- readJDX::readJDX(file = name, SOFC = TRUE, debug = 0) # reading urine_1.dx (~1MB) takes ~30s on machine r31

    # Scale factors for x and y axis
    factor_x <- scale_factor[1]
    factor_y <- scale_factor[2]

    # Get the length of the spectrum
    spectrum_length <- length(data[[4]]$x) - 1

    # Scale the x and y axis values
    spectrum_x <- seq((spectrum_length / factor_x), 0, -1 / factor_x)
    spectrum_y <- (data[[4]]$y) / factor_y

    # Extract spectral width in ppm
    ppm_range_index <- which(startsWith(data[[2]], "##$SW="))
    ppm_range <- as.numeric(gsub("\\D+", "", data[[2]][ppm_range_index]))

    # Extract highest and lowest ppm values
    ppm_highest_index <- which(startsWith(data[[2]], "##$OFFSET="))
    ppm_highest_value <- as.numeric(gsub("\\D+", "", data[[2]][ppm_highest_index]))
    ppm_lowest_value <- ppm_highest_value - ppm_range

    # Generate ppm sequence for x axis
    spectrum_x_ppm <- seq(ppm_highest_value, ppm_lowest_value, (-ppm_range / spectrum_length))

    # Return the spectrum data
    return(list(spectrum_x = spectrum_x,
                spectrum_y = spectrum_y,
                spectrum_x_ppm = spectrum_x_ppm,
                spectrum_length = spectrum_length,
                ppm_range = ppm_range,
                ppm_highest_value = ppm_highest_value,
                ppm_lowest_value = ppm_lowest_value))
}