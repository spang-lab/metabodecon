#' @noRd
#' @title Combine Peaks
#' @description Even after calling `speaq_align()`, the alignment of individual signals is not always perfect, as 'speaq' performs a segment-wise alignment i.e. groups of signals are aligned. For further improvements, partly filled neighboring columns are merged.
#' @param shifted_mat The matrix returned by `speaq_align()`.
#' @param range Amount of adjacent columns which are permitted to be used for improving the alignment.
#' @param lower_bound Minimum amount of non-zero elements per column to trigger the alignment improvement.
#' @param spectrum_data The list of deconvoluted spectra as returned by `generate_lorentz_curves()` that was used to generate `shifted_mat`. No longer required since version 1.2 of Metabodecon.
#' @param data_path If not NULL, the returned dataframes `long` and `short` are written to `data_path` as "aligned_res_long.csv" and "aligned_res_short.csv".
#' @return A list containing two data frames `long` and `short`. The first data frame contains one one column for each data point in the original spectrum. The second data frame contains only columns where at least one entry is non-zero.
#' @details
#' Example of what the function does:
#'
#' ```txt
#' |            | 3.56 | 3.54 | 3.51 | 3.51 | 3.50 |
#' |----------- |------|------|------|------|------|
#' | Spectrum 1 | 0.13 | 0    | 0.11 | 0    | 0    |
#' | Spectrum 2 | 0.13 | 0    | 0.12 | 0    | 0    |
#' | Spectrum 3 | 0.07 | 0    | 0    | 0    | 0    |
#' | Spectrum 4 | 0.08 | 0    | 0    | 0.07 | 0    |
#' | Spectrum 5 | 0.04 | 0    | 0.04 | 0    | 0    |
#'
#' becomes
#'
#' |            | 3.56 | 3.54 | 3.51 | 3.50 |
#' |----------- |------|------|------|------|
#' | Spectrum 1 | 0.13 | 0    | 0.11 | 0    |
#' | Spectrum 2 | 0.13 | 0    | 0.12 | 0    |
#' | Spectrum 3 | 0.07 | 0    | 0    | 0    |
#' | Spectrum 4 | 0.08 | 0    | 0.07 | 0    |
#' | Spectrum 5 | 0.04 | 0    | 0.04 | 0    |
#'
#' I.e. column 3 and 4 get merged, because they are in `range` of each other and have no common non-zero entries.
#' ```
#'
#' @author Initial version from Wolfram Gronwald. Refactored by Tobias Schmidt in 2024.
#' @examples
#' spectrum_data <- generate_lorentz_curves_sim("bruker/sim")
#' shifted_mat <- speaq_align(spectrum_data = spectrum_data, verbose = FALSE)
#' range <- 5
#' lower_bound <- 1
#' obj <- combine_peaks(shifted_mat, range, lower_bound)
#' str(obj)
combine_peaks <- function(shifted_mat = speaq_align(),
                          range = 5,
                          lower_bound = 1,
                          spectrum_data = NULL,
                          data_path = NULL) {
    M <- replace(shifted_mat, is.na(shifted_mat), 0)
    U <- M != 0 # Unequal zero matrix. U[i,j] is TRUE if M[i,j] is nonzero, else FALSE.
    uu <- colSums(U) # Unequal zero vector. u[j] gives the amount of nonzero elements in M[,j].
    for (i in (nrow(M) - 1):lower_bound) {
        for (j in which(uu == i)) {
            if (uu[j] == 0) next
            nn <- c((j - range):(j - 1), (j + 1):(j + range)) # Neighbors of j.
            cc <- combine_scores(U, uu, j, nn)
            while (!all(cc == 0)) {
                n <- nn[which.max(cc)]
                M[, j] <- M[, j] + M[, n]
                U[, j] <- U[, j] | U[, n]
                uu[j] <- uu[j] + uu[n]
                M[, n] <- 0
                U[, n] <- FALSE
                uu[n] <- 0
                cc <- combine_scores(U, uu, j, nn)
            }
        }
    }
    if (!is.null(spectrum_data)) colnames(M) <- spectrum_data[[1]]$x_values_ppm
    S <- M[, uu > 0]
    if (!is.null(data_path)) {
        long.csv <- file.path(data_path, "aligned_res_long.csv")
        short.csv <- file.path(data_path, "aligned_res_short.csv")
        utils::write.csv2(M, file = long.csv)
        utils::write.csv2(S, file = short.csv)
    }
    list("short" = S, "long" = M)
}

#' @noRd
#' @description Calculates a "combine score" for a set of columns `nn` of a Matrix `M`.
#' The score describes how "beneficial" it is to merge column `n` (from `nn`) into column `j`.
#' A column `n` is considered "combinable" if there is no row where columns `n` and `j` both have a non-zero element.
#' If a column is not combinable, its combine score is 0.
#' If a column is combinable, its combine score is the amount of non-zero elements in the column.
#' This function calculates the combine score for all columns `n` in `nn` and then returns the calculated combine scores as vector.
#' @param U Matrix describing nonzero entries of `M.` `U[i,j]` should be `TRUE` if `M[i,j]` is nonzero, else `FALSE.` I.e. for any matrix `M` you can generate `U` as `U <- M != 0`.
#' @param uu Vector describing the amount of nonzero elements in each column of `M`. I.e. for any matrix `M` you can generate `uu` as `uu <- colSums(M != 0)`. Called `uu` because it denotes the amount of elements that are *unequal zero*.
#' @param j Index of the column where the neighbor columns `nn` should be combined into.
#' @param nn Indices of the columns to calculate the combine score for. Called `nn` because this should be a set of *neighbor columns*.
#' @details Since we only need to know whether an element is nonzero in order to calculate the combine scores, the function takes the matrix `U` and the vector `uu` as input instead of their common ancestor `M`. This is much more efficient if the function is called multiple times, as the conversion must not be done multiple times.
#' @examples
#' M <- rbind(
#'     c(2, 0, 2, 2, 0),
#'     c(2, 1, 0, 2, 0),
#'     c(0, 1, 0, 2, 0),
#'     c(0, 0, 3, 0, 1)
#' )
#' U <- M != 0
#' uu <- colSums(U) # 2 2 2 3 1
#'
#' cc <- combine_scores(U, uu, j = 2, nn = c(1, 3))
#' cc[1] == 0  # M[,1] and M[, 2] are not combinable
#' cc[2] == 2  # M[,3] and M[, 2] are combinable and M[, 3] has two nonzero elements
#'
#' cc <- combine_scores(U, uu, j = 1, nn = 2:5)
#' cc[1] == 0  # M[, 2] and M[, 1] are not combinable
#' cc[2] == 0  # M[, 3] and M[, 1] are not combinable
#' cc[3] == 0  # M[, 4] and M[, 1] are not combinable
#' cc[4] == 1  # M[, 5] and M[, 1] are combinable and M[, 5] has one nonzero element
combine_scores <- function(U, uu, j, nn) {
    sapply(seq_along(nn), function(k) {
        n <- nn[k] # Index of neighboring column
        combinable <- sum(U[, n] & U[, j]) == 0
        if (combinable) uu[n] else 0
    })
}
