% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decon.R, R/depr.R
\name{deconvolute}
\alias{deconvolute}
\alias{generate_lorentz_curves_sim}
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

generate_lorentz_curves_sim(
  data_path,
  file_format = "bruker",
  make_rds = FALSE,
  expno = 10,
  procno = 10,
  raw = TRUE,
  nfit = 10,
  smopts = c(2, 5),
  delta = 6.4,
  sfr = c(3.55, 3.35),
  wshw = 0,
  ask = FALSE,
  force = FALSE,
  verbose = FALSE,
  nworkers = 1
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

\item{data_path}{The path of the file/folder containing the spectrum data. E.g.
\code{"example_datasets/jcampdx/urine/urine_1.dx"} or
\code{"example_datasets/bruker/urine/urine"}.}

\item{file_format}{The file_format of the spectrum file. E.g. \code{"bruker"} or \code{"jcampdx"}.}

\item{make_rds}{Logical or character. If TRUE, stores results as an RDS file
on disk. If a character string, saves the RDS file with the specified name.
Should be set to TRUE if many spectra are evaluated to decrease computation
time.}

\item{expno, procno}{The experiment/processing number for the file. E.g. \code{"10"}. Only relevant if
\code{file_format} equals \code{"bruker"}. For details see section \href{https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure}{File Structure} in
the metabodecon FAQ.}

\item{raw}{If \code{FALSE}, scales the returned signal intensities based on information
available in the spectrum metadata, in particular \code{NC_proc}. For details see
\code{processing-reference.pdf}, available at \url{https://www.bruker.com/en.html} at
section 'Services & Support > Documentation & Manuals > Magnetic Resonance >
Acquisition & Processing > TopSpin Processing Commands and Parameters'
(requires login).}
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

\code{\link[=generate_lorentz_curves_sim]{generate_lorentz_curves_sim()}} is identical to \code{\link[=generate_lorentz_curves]{generate_lorentz_curves()}}
except for the defaults, which are optimized for deconvoluting the 'Sim'
dataset, shipped with 'metabodecon'. The 'Sim' dataset is a simulated
dataset, which is much smaller than a real NMR spectra and lacks a water
signal. This makes it ideal for use in examples or as a default value for
functions. However, the default values for \code{sfr}, \code{wshw}, and \code{delta} in the
"normal" \code{\link[=generate_lorentz_curves]{generate_lorentz_curves()}} function are not optimal for this
dataset. To avoid having to define the optimal parameters repeatedly in
examples, this function is provided to deconvolute the "Sim" dataset with
suitable parameters.

In \code{\link[=generate_lorentz_curves]{generate_lorentz_curves()}} the parameters \code{nfit}, \code{smopts}, \code{delta},
\code{sfr} and \code{wshw} must be fully specified. In \code{\link[=deconvolute]{deconvolute()}}, these
parameters can be set to \code{NULL} (the default value). In this case, the
function will try to determine the optimal values for these parameters
automatically. The values chosen are stored in field \code{args} of the returned
\code{decon2} object.
}
\examples{

## Define the paths to the example datasets we want to deconvolute:
## `sim_dir`: directory containing 16 simulated spectra
## `sim_01`: path to the first spectrum in the `sim` directory
## `sim_01_spec`: the first spectrum in the `sim` directory as a dataframe

sim_dir <- metabodecon_file("sim_subset")
sim_1_dir <- file.path(sim_dir, "sim_01")
sim_2_dir <- file.path(sim_dir, "sim_02")
sim_1_spectrum <- read_spectrum(sim_1_dir)
sim_2_spectrum <- read_spectrum(sim_2_dir)
sim_spectra <- read_spectra(sim_dir)


## Show that `generate_lorentz_curves()` and `generate_lorentz_curves_sim()`
## produce the same results:

sim_1_decon0 <- generate_lorentz_curves(
    data_path = sim_1_dir, # Path to directory containing spectra
    sfr = c(3.55, 3.35), # Borders of signal free region (SFR) in ppm
    wshw = 0, # Half width of water signal (WS) in ppm
    ask = FALSE, # Don't ask for user input
    verbose = FALSE # Suppress status messages
)
sim_1_decon1 <- generate_lorentz_curves_sim(sim_1_dir)
stopifnot(all.equal(sim_1_decon0, sim_1_decon1))


## Show that passing a spectrum produces the same results as passing the
## the corresponding directory:

decon_from_spectrum_dir <- generate_lorentz_curves_sim(sim_1_dir)
decon_from_spectrum_obj <- generate_lorentz_curves_sim(sim_1_spectrum)
decons_from_spectra_obj <- generate_lorentz_curves_sim(sim_spectra)
decons_from_spectra_dir <- generate_lorentz_curves_sim(sim_dir)

most.equal <- function(x1, x2) {
    ignore <- which(names(x1) \%in\% c("number_of_files", "filename"))
    equal <- all.equal(x1[-ignore], x2[-ignore])
    invisible(stopifnot(isTRUE(equal)))
}

all.equal(  decon_from_spectrum_dir, decon_from_spectrum_obj     )
all.equal(  decons_from_spectra_dir, decons_from_spectra_obj     )
most.equal( decon_from_spectrum_dir, decons_from_spectra_obj[[1]])
most.equal( decon_from_spectrum_dir, decons_from_spectra_dir[[1]])
}
\author{
2024-2025 Tobias Schmidt: initial version.
}
