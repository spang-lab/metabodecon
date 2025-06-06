% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depr.R
\name{plot_lorentz_curves_save_as_png}
\alias{plot_lorentz_curves_save_as_png}
\title{Plot lorentz curves for variable range}
\usage{
plot_lorentz_curves_save_as_png(
  deconv_result,
  x_range = c(),
  y_range = c(),
  out_dir = ".",
  ask = TRUE
)
}
\arguments{
\item{deconv_result}{Saved result of the MetaboDecon1D() function}

\item{x_range}{Row vector with two entries consisting of the ppm start and the ppm end value
to scale the range of the x-axis (optional)}

\item{y_range}{Row vector with two entries consisting of the ppm start and the ppm end value
to scale the range of the y-axis (optional)}

\item{out_dir}{Path to the directory where the png files should be saved. Default is the
current working directory.}

\item{ask}{Logical value. Whether to ask for confirmation from the user before writing
files to disk. Default is TRUE.}
}
\value{
NULL, called for side effects.
}
\description{
Plots the original spectrum and all generated Lorentz curves and save the
result as png under the filepath.

Superseded by \code{\link[=plot_spectrum]{plot_spectrum()}} since metabodecon v1.2.0. Will be replaced
with metabodecon v2.

\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
}
\examples{
sim <- metabodecon_file("bruker/sim_subset")
sim_decon <- generate_lorentz_curves_sim(sim)
png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
plot_lorentz_curves_save_as_png(sim_decon, out_dir = png_dir, ask = FALSE)
dir(png_dir, full.names = TRUE)
}
\seealso{
\code{\link[=MetaboDecon1D]{MetaboDecon1D()}}, \code{\link[=plot_triplets]{plot_triplets()}}, \code{\link[=plot_spectrum_superposition_save_as_png]{plot_spectrum_superposition_save_as_png()}}
}
\author{
2020-2021 Martina Haeckl: initial version.\cr
2024-2025 Tobias Schmidt: Minor updates to pass CRAN checks
}
