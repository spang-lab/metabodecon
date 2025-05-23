% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{metabodecon_file}
\alias{metabodecon_file}
\title{Return Path to File or Directory in metabodecon Package}
\usage{
metabodecon_file(name = "sim_01")
}
\arguments{
\item{name}{The name to search for.}
}
\value{
The file or directory path.
}
\description{
Recursively searches for files or directories within the
'metabodecon' package that match the given name.
}
\examples{
# Unambiguous paths
metabodecon_file("urine_1")
metabodecon_file("urine_1.dx")
metabodecon_file("sim/sim_01")

# Ambiguous paths (i.e. multiple matches)
metabodecon_file("sim")
metabodecon_file("urine")

# Non-existing paths (i.e. a character vector of length zero gets returned)
metabodecon_file("asdfasdf")
}
\author{
2024-2025 Tobias Schmidt: initial version.
}
