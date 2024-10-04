#' @export
`[.spectra` <- function(x, i, ...) {
    result <- NextMethod("[")
    class(result) <- class(x)
    result
}
