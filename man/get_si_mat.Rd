% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align.R
\name{get_si_mat}
\alias{get_si_mat}
\title{Extract Matrix of aligned Signal Intensities}
\usage{
get_si_mat(x)
}
\arguments{
\item{x}{An object of type \code{aligns}.}
}
\value{
A matrix of aligned signal intensities.
}
\description{
Takes an object of type \code{aligns}, i.e., a list of deconvoluted and aligned
spectra, extracts the vector of aligned signal integrals for each spectrum
and returns them as a matrix with datapoints in rows and spectra in columns.
}
\examples{
decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
aligns <- align(decons)
si_mat <- get_si_mat(aligns) # 2048 x 2 matrix (2048 datapoints, 2 spectra)
}
\author{
2024-2025 Tobias Schmidt: initial version.
}
