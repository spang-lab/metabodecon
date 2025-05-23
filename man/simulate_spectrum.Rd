% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spectrum.R
\name{simulate_spectrum}
\alias{simulate_spectrum}
\title{Simulate a 1D NMR Spectrum}
\usage{
simulate_spectrum(
  name = "sim_00",
  seed = sum(utf8ToInt(name)),
  ndp = 2048,
  npk = 10,
  csres = 0.00015,
  cs = seq(from = 3.6, length.out = ndp, by = -csres),
  pkr = quantile(cs, c(0.25, 0.75)),
  fqref = 600252806.95,
  x0 = sort(runif(npk, pkr[1], pkr[2])),
  A = runif(npk, 2.5, 20) * 1000,
  lambda = runif(npk, 0.9, 1.3)/1000,
  noise = rnorm(length(cs), sd = 1200)
)
}
\arguments{
\item{name}{The name of the spectrum.}

\item{seed}{The seed for the random number generator.}

\item{ndp}{The number of data points in the spectrum.}

\item{npk}{The number of peaks in the spectrum.}

\item{csres}{The chemical shift resolution in PPM.}

\item{cs}{The vector of chemical shifts in PPM.}

\item{pkr}{The start and stop of the peak region in PPM.}

\item{fqref}{The reference frequency in Hz.}

\item{x0}{The peak center positions in PPM.}

\item{A}{The peak area parameter.}

\item{lambda}{The peak width parameter.}

\item{noise}{The noise to add to the spectrum.}
}
\value{
A \code{spectrum} object as described in \href{https://spang-lab.github.io/metabodecon/articles/Classes.html}{Metabodecon Classes}.
}
\description{
Simulates a 1D NMR spectrum based on the provided parameters.

\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
}
\examples{
simA <- simulate_spectrum("simA")
simA_copy <- simulate_spectrum("simA")
simB <- simulate_spectrum("simB")
simC <- simulate_spectrum("simC", npk = 20)
plot_spectrum(simC)
if (!identical(simA, simA_copy)) stop()
if ( identical(simA, simB     )) stop()
}
\author{
2024-2025 Tobias Schmidt: initial version.
}
