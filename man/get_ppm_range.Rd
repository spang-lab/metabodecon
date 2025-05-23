% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align.R
\name{get_ppm_range}
\alias{get_ppm_range}
\title{Get PPM Range covered by Spectra}
\usage{
get_ppm_range(spectrum_data, full_range = FALSE)
}
\arguments{
\item{spectrum_data}{A list of deconvoluted spectra as returned by \code{\link[=generate_lorentz_curves]{generate_lorentz_curves()}}.}

\item{full_range}{If TRUE, the full range of the spectra is returned. If FALSE, only the range
from the lowest to the highest peak center is returned.}
}
\value{
A vector containing the lowest and highest ppm value over all peaks of the
provided deconvoluted spectra.
}
\description{
Helper function of \code{align()}. Should not be called directly by the user.

Returns the ppm range across all peaks of the provided deconvoluted spectra.

Direct usage of this function has been deprecated with metabodecon version
1.4.3 and will be removed with metabodecon version 2.0.0.

\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
}
\examples{
spectrum_data <- generate_lorentz_curves(
    data_path = sim[1:2],
    nfit = 3,
    sfr = c(3.55, 3.35),
    wshw = 0,
    ask = FALSE,
    verbose = FALSE
)
ppm_rng <- get_ppm_range(spectrum_data)
print(ppm_rng)
}
\author{
2021-2024 Wolfram Gronwald: initial version.\cr
2024-2025 Tobias Schmidt: refactored initial version.
}
