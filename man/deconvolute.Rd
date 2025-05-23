% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decon.R
\name{deconvolute}
\alias{deconvolute}
\title{Deconvolute one or more NMR spectra}
\usage{
deconvolute(
  x,
  nfit = 3,
  smopts = c(2, 5),
  delta = 6.4,
  sfr = NULL,
  wshw = 0,
  ask = FALSE,
  force = FALSE,
  verbose = TRUE,
  nworkers = 1,
  use_rust = FALSE
)
}
\arguments{
\item{x}{A \code{spectrum} or \code{spectra} object as described in \href{https://spang-lab.github.io/metabodecon/articles/Classes.html}{Metabodecon Classes}.}

\item{nfit}{Integer. Number of iterations for approximating the parameters
for the Lorentz curves. See 'Details'.}

\item{smopts}{Numeric vector with two entries: the number of smoothing
iterations and the number of data points to use for smoothing (must be odd).
See 'Details'.}

\item{delta}{Threshold for peak filtering. Higher values result in more peaks
being filtered out. A peak is filtered if its score is below \eqn{\mu +
\sigma \cdot \delta}{mu + s * delta}, where \eqn{\mu}{mu} is the average
peak score in the signal-free region (SFR), and \eqn{\sigma}{s} is the
standard deviation of peak scores in the SFR. See 'Details'.}

\item{sfr}{Numeric vector with two entries: the ppm positions for the left
and right border of the signal-free region of the spectrum. See 'Details'.}

\item{wshw}{Half-width of the water artifact in ppm.  See 'Details'.}

\item{ask}{Logical. Whether to ask for user input during the deconvolution
process. If FALSE, the provided default values will be used.}

\item{force}{If FALSE, the function stops with an error message if no peaks
are found in the signal free region (SFR), as these peaks are required as a
reference for peak filtering. If TRUE, the function instead proceeds without
peak filtering, potentially increasing runtime and memory usage
significantly.}

\item{verbose}{Logical. Whether to print log messages during the
deconvolution process.}

\item{nworkers}{Number of workers to use for parallel processing. If
\code{"auto"}, the number of workers will be determined automatically. If a number
greater than 1, it will be limited to the number of spectra.}

\item{use_rust}{Logical. Whether to use the Rust backend for deconvolution.
Requires the \href{https://github.com/spang-lab/mdrb}{mdrb} package. If TRUE and
mdrb is missing, an error is thrown. If FALSE, the R implementation is used.
If NULL, the Rust backend is used if available, otherwise the R implementation
is used.}
}
\value{
A 'decon2' object as described in \href{https://spang-lab.github.io/metabodecon/articles/Classes.html}{Metabodecon Classes}.
}
\description{
Deconvolutes NMR spectra by modeling each detected signal within
a spectrum as Lorentz Curve.
}
\details{
First, an automated curvature based signal selection is performed. Each
signal is represented by 3 data points to allow the determination of initial
Lorentz curves. These Lorentz curves are then iteratively adjusted to
optimally approximate the measured spectrum.
}
\examples{
## Deconvolute a single spectrum
spectrum <- sim[1]
decon <- deconvolute(spectrum)

## Read multiple spectra from disk and deconvolute at once
spectra_dir <- metabodecon_file("sim_subset")
spectra <- read_spectra(spectra_dir)
decons <- deconvolute(spectra, sfr = c(3.55, 3.35))
}
\author{
2024-2025 Tobias Schmidt: initial version.
}
