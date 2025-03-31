# Imports ######################################################################

# Standard Lib
#' @import grDevices
#' @import stats
#' @import utils
#' @import parallel
#' @import stats
#' @import graphics
#' @import grDevices
#' @import mathjaxr

# 3rd Party
#' @import withr
#' @import toscutil

# Docs ########################################################################

#' @keywords internal
"_PACKAGE"

#' @export
#' @title Get URL of Metabodecon "Get Started" Page
#' @description
#' `get_started` and `aaa_get_started` both return (and optionally open) the URL
#' of the "Get Started" page of the metabodecon documentation. The
#' `aaa_get_started` version exists, because functions are listed alphabetically
#' in the reference manual and we want `get_started` to be shown at the top of
#' the list.
#' @param open_browser If TRUE, the "Get Stated" page is opened in the default
#' browser.
#' @return A character string containing the URL of the "Get Started" page.
#' @examples
#' get_started(open_browser = FALSE)
#' get_started()
get_started <- function(open_browser = interactive()) {
    url <- "https://spang-lab.github.io/metabodecon/articles/Get_Started.html"
    if (open_browser) utils::browseURL(url)
    url
}

#' @export
#' @rdname get_started
`aaa_get_started` <- get_started
