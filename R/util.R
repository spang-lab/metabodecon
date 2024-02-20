`%||%` <- function (x, y)
{
    if (is.null(x)) y else x
}

msg <- function(..., sep = " ", appendLF = TRUE)
{
    message(paste(..., sep = sep), appendLF = appendLF)
}

msgf <- function(fmt, ..., appendLF = TRUE)
{
    message(sprintf(fmt, ...), appendLF = appendLF)
}
