% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{get_data_dir}
\alias{get_data_dir}
\title{Retrieve directory path of an example dataset}
\usage{
get_data_dir(dataset_name = c("", "blood", "test", "urine"), warn = TRUE)
}
\arguments{
\item{dataset_name}{Either \code{""}, \code{"test"}, \code{"blood"} or \code{"urine"}.}

\item{warn}{Whether to print a warning message when the example folders do
not yet exist, i.e. \code{\link[=download_example_datasets]{download_example_datasets()}} has not been called yet.}
}
\description{
Returns the path to the directory storing the example files
shipped with metabodecon.

Deprecated. Please use \code{\link[=datadir]{datadir()}} instead. See examples below for usage.
``
}
\examples{
x <- get_data_dir("urine")                     # Deprecated
y <- datadir("example_datasets/bruker/urine")  # Preferred
}
\seealso{
\code{\link[=download_example_datasets]{download_example_datasets()}}
}
