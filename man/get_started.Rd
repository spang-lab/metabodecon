% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{get_started}
\alias{get_started}
\alias{aaa_get_started}
\title{Get URL of Metabodecon "Get Started" Page}
\usage{
get_started(open_browser = interactive())

aaa_get_started(open_browser = interactive())
}
\arguments{
\item{open_browser}{If TRUE, the "Get Stated" page is opened in the default
browser.}
}
\value{
A character string containing the URL of the "Get Started" page.
}
\description{
\code{get_started} and \code{aaa_get_started} both return (and optionally open) the URL
of the "Get Started" page of the metabodecon documentation. The
\code{aaa_get_started} version exists, because functions are listed alphabetically
in the reference manual and we want \code{get_started} to be shown at the top of
the list.
}
\examples{
get_started(open_browser = FALSE)
get_started()
}
