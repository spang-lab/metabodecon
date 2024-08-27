#' @export
#' @title Combine Peaks
#' @description Even after calling `speaq_align()`, the alignment of individual signals is not always perfect, as 'speaq' performs a segment-wise alignment i.e. groups of signals are aligned. For further improvements, partly filled neighboring columns are merged.
#' @param shifted_mat The matrix returned by `speaq_align()`.
#' @param range Amount of adjacent columns which are permitted to be used for improving the alignment.
#' @param lower_bound Minimum amount of non-zero elements per column to trigger the alignment improvement.
#' @param spectrum_data The list of deconvoluted spectra as returned by `generate_lorentz_curves()`.
#' @param data_path If not NULL, the returned dataframes `long` and `short` are written to `data_path` as "aligned_res_long.csv" and "aligned_res_short.csv".
#' @return A list containing two data frames `long` and `short`. The first one contains one for one column for each data point in the original spectrum. The second one contains only columns where at least one entry is non-zero. The returned data.frames are additionally written to disk as .csv files `aligned_res_long.csv` and `aligned_res_short.csv`.
#' @examples
#' shifted_mat <- speaq_align(spectrum_data = glc_sim("bruker/sim"))
#' system.time(M1 <- combine_peaks_v1(shifted_mat, spectrum_data = glc_sim("bruker/sim")))
#' system.time(M2 <- combine_peaks(shifted_mat, spectrum_data = glc_sim("bruker/sim")))
#' all.equal(M1, M2)
combine_peaks <- function(shifted_mat = speaq_align(spectrum_data = spectrum_data, show = TRUE),
                          range = 5,
                          lower_bound = 1,
                          spectrum_data = glc_sim("bruker/sim"),
                          data_path = NULL) {
    M <- replace(shifted_mat, is.na(shifted_mat), 0) # Shifted matrix with NAs replaced by 0
    s <- nrow(M) # Number of spectra
    nz <- apply(M, 2, function(x) sum(x != 0)) # Number of Non-zero Entries (NZs) per column
    nzIDs <- which(nz >= lower_bound) # Columns IDs (CIDs) with at least `lower_bound` NZs.
    nzIDs <- nzIDs[order(nz[nzIDs], decreasing = TRUE)] # Sort column IDs with most NZs first
    for (j in nzIDs) {
        if (nz[j] == 0) next # Can happen if column j has already been absorbed by another column
        while (TRUE) {
            nbIDs <- c((j - range):(j - 1), (j + 1):(j + range)) # Indices of Neighbour columns
            nbIDs <- nbIDs[nbIDs > 0 & nbIDs <= ncol(M)] # Remove out-of-bounds indices
            cvs <- nz[nbIDs] # Number of Combinable Values (CVs) per column
            if (all(cvs %in% c(0, s))) break # Stop early if neighbour cols are empty or completely filled already
            for (i in seq_along(cvs)[cvs > 0 & cvs < s]) {
                n <- nbIDs[i] # Index of neighbour column
                if (any(M[, n] != 0 & M[, j] != 0)) cvs[n] <- 0
                # Example:
                # c(0, 0, 3, 4, 0) -> M[, n]
                # c(2, 3, 0, 2, 1) -> M[, j]
                # c(F, F, F, T, F) -> M[, n] != 0 & M[, j] != 0
                # Not combinable because there is a common non-zero element at index 4
            }
            if (all(cvs %in% c(0, nrow(M)))) break
            x <- nbIDs[which.max(cvs)] # Index of column to combine with current column
            if (any(M[, x] == 0.03336328)) browser()
            M[, j] <- M[, j] + M[, x]
            M[, x] <- 0
            nz[x] <- 0
        }
    }
    ppm_rounded <- round(spectrum_data[[1]]$x_values_ppm, 4)
    colnames(M) <- as.character(ppm_rounded)
    nz <- apply(M, 2, function(x) sum(x != 0))
    nzIDs <- which(nz >= lower_bound)
    S <- M[, nzIDs] # Shifted matrix without zero-only-columns
    # utils::write.csv2(M2, file = file.path(data_path, "aligned_res_long.csv"))
    # utils::write.csv2(M, file = file.path(data_path, "aligned_res_short.csv"))
    list("short" = S, "long" = M)
}