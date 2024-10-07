# Structures #####

str_urine_1_deconvoluted <- function(cf = 1, nf = 1, dx = FALSE, nested = TRUE, ni = 10) {
    fn <- if (dx) "urine_1.dx" else "urine_1"
    elemstr <- c(
        sprintf("number_of_files           : int %d", nf),
        sprintf('filename                  : chr "%s"', fn),
        sprintf("x_values                  : num [1:131072] 131 131 131 131 131 ..."),
        sprintf("x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ..."),
        sprintf("y_values                  : num [1:131072] 0.000831 0.000783 0.000743 0.000717 0.00065 ..."),
        sprintf("spectrum_superposition    : num [1, 1:131072] 3.51e-05 3.51e-05 3.51e-05 3.51e-05 3.52e-05 ..."),
        sprintf("mse_normed                : num 3.92e-11"),
        sprintf("index_peak_triplets_middle: num [1:1227] 36159 37149 37419 37435 38943 ..."),
        sprintf("index_peak_triplets_left  : num [1:1227] 36161 37160 37423 37438 38949 ..."),
        sprintf("index_peak_triplets_right : num [1:1227] 36156 37140 37415 37432 38938 ..."),
        sprintf("peak_triplets_middle      : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ..."),
        sprintf("peak_triplets_left        : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ..."),
        sprintf("peak_triplets_right       : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ..."),
        sprintf("integrals                 : num [1, 1:1227] 0.000501 0.026496 0.000402 0.000375 0.008274 ..."),
        # sprintf("signal_free_region        : num [1:2] 109.1 21.9"),
        # sprintf("range_water_signal_ppm    : num 0.153"),
        sprintf("A                         : num [1:1227] -0.00016 -0.008436 -0.000128 -0.000119 -0.002634 ..."),
        sprintf("lambda                    : num [1:1227] -0.00775 -0.02188 -0.00675 -0.00562 -0.01343 ..."),
        sprintf("x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ...")
    )
    ne <- length(elemstr)
    plainstr <- c(
        paste0("List of ", ne),
        paste0(" $ ", elemstr)
    )
    nestedstr <- c(
        paste0("List of ", nf),
        paste0(" $ ", fn, ":List of ", ne),
        paste0("  ..$ ", elemstr)
    )
    if (nested) nestedstr else plainstr
}

str_urine_2_deconvoluted <- function(cf = 2, nf = 2, dx = FALSE, nested = TRUE, ni = 10) {
    fn <- if (dx) "urine_2.dx" else "urine_2"
    elemstr <- c(
        sprintf("number_of_files           : int %d", nf),
        sprintf('filename                  : chr "%s"', fn),
        sprintf("x_values                  : num [1:131072] 131 131 131 131 131 ..."),
        sprintf("x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ..."),
        sprintf("y_values                  : num [1:131072] 0.00586 0.00578 0.00569 0.00557 0.00548 ..."),
        sprintf("spectrum_superposition    : num [1, 1:131072] 4.21e-05 4.21e-05 4.21e-05 4.21e-05 4.21e-05 ..."),
        sprintf("mse_normed                : num 2.86e-11"),
        sprintf("index_peak_triplets_middle: num [1:1393] 36290 37241 38346 38826 39025 ..."),
        sprintf("index_peak_triplets_left  : num [1:1393] 36297 37244 38349 38835 39028 ..."),
        sprintf("index_peak_triplets_right : num [1:1393] 36285 37234 38343 38823 39019 ..."),
        sprintf("peak_triplets_middle      : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ..."),
        sprintf("peak_triplets_left        : num [1:1393] 9.26 9.12 8.95 8.87 8.85 ..."),
        sprintf("peak_triplets_right       : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ..."),
        sprintf("integrals                 : num [1, 1:1393] 0.00683 0.00504 0.00322 0.0174 0.00274 ..."),
        # sprintf("signal_free_region        : num [1:2] 109.1 21.9"),
        # sprintf("range_water_signal_ppm    : num 0.153"),
        sprintf("A                         : num [1:1393] -0.002176 -0.001604 -0.001025 -0.005541 -0.000872 ..."),
        sprintf("lambda                    : num [1:1393] -0.0189 -0.0168 -0.014 -0.0291 -0.0146 ..."),
        sprintf("x_0                       : num [1:1393] 94.8 93.8 92.7 92.2 92 ...")
    )
    ne <- length(elemstr)
    plainstr <- c(
        paste0("List of ", ne),
        paste0(" $ ", elemstr)
    )
    nestedstr <- c(
        paste0("List of ", nf),
        paste0(" $ ", fn, ":List of ", ne),
        paste0("  ..$ ", elemstr)
    )
    if (nested) nestedstr else plainstr
}

str_urine_deconvoluted <- function(nf = 2, dx = FALSE, nested = TRUE, ni = 10) {
    u1 <- str_urine_1_deconvoluted(nf = 2, dx = dx, nested = TRUE, ni = ni)
    u2 <- str_urine_2_deconvoluted(nf = 2, dx = dx, nested = TRUE, ni = ni)
    c("List of 2", u1[2:length(u1)], u2[2:length(u2)])
}

# Compare #####

#' @noRd
#' @title Check the quality of a deconvolution by comparing with the true
#' parameters
#'
#' @description
#' Checks the quality of a deconvolution by comparing the deconvolution results
#' with the true parameters (which are only known for simulated spectra). See
#' 'Details' for more information on how the quality is assessed.
#'
#' @param decon A list containing the deconvolution results, as returned by
#' [generate_lorentz_curves()].
#'
#' @param lcpar A data frame containing the true parameters of the peaks, as
#' returned by `get_sim_params()`.
#'
#' @return A quality score for the deconvolution. In addition, a plot is created
#' to visualize the deconvolution results.
#'
#' @details
#' The function first plots the deconvolution results for visual inspection and
#' then returns a quality score for the deconvolution.
#'
#' The plotting is done as follows:
#'
#' 1. Plot the deconvoluted spectrum using `plot_spectrum()`.
#' 2. Draw green circles around found peaks [^1].
#' 3. Draw red circles around missed peaks [^1].
#' 4. Draw red rectangles around falsely detected peaks. [^2]
#'
#' [^1]: we consider a peak as 'found' if there is at least one detected peak
#'       center within 0.001 ppm of the true peak position. If this is not the
#'       case, the peak is considered as 'missed'.
#' [^2]: we consider a peak as 'falsely detected' if there is no true peak
#'       center within 0.001 ppm of the detected peak position.
#'
#' In addition, a quality score is calculated as follows:
#'
#' quality    = peak_ratio * area_ratio
#' peak_ratio = min(peaks_true, peaks_found) / max(peaks_true, peaks_found)
#' area_ratio = min(area_true,  area_found)  / max(area_true,  area_found)
#'
#' I.e., the score is close to 1 if the number of peaks and the area of the
#' peaks are similar in the true and found spectra and the score is close to 0
#' if the number of peaks and/or the area of the peaks are very different.
check_decon_quality <- function(decon, lcpar) {
    true <- lcpar # True Params
    found <- data.frame( # Found Params
        A      = decon$A_ppm,
        x_0    = decon$x_0_ppm,
        lambda = decon$lambda_ppm
    )

    # Calculate y values at found x_0 position
    ppm <- decon$x_values_ppm
    y <- decon$y_values
    i_0 <- convert_pos(found$x_0, ppm, 1:length(ppm))
    i_0_floor <- floor(i_0)
    i_0_ceil <- ceiling(i_0)
    i_0_frac <- i_0 - i_0_floor
    y_floor <- y[i_0_floor]
    y_ceil <- y[i_0_ceil]
    found$y_0 <- y_floor + (y_ceil - y_floor) * i_0_frac

    # Check which peaks are found correctly and which were missed
    for (i in seq_along(true$x_0)) {
        true[i, "closest"] <- j <- which.min(abs(found$x_0 - true$x_0[i]))
        true[i, "dist"] <- d <- abs(found$x_0[j] - true$x_0[i])
        true[i, "correct"] <- d < 0.001
    }
    for (i in seq_along(found$x_0)) {
        found[i, "closest"] <- j <- which.min(abs(true$x_0 - found$x_0[i]))
        found[i, "dist"] <- d <- abs(true$x_0[j] - found$x_0[i])
        found[i, "correct"] <- d < 0.001
    }

    plot_spectrum(decon, foc_rgn = c(0.2, 0.8), sub_show = FALSE)
    points(
        x = found$x_0,
        y = found$y_0,
        pch = 21,
        cex = 3,
        bg = "transparent",
        col = ifelse(found$correct, "green", "red")
    )
    abline(v = found$x_0)

    symbols(
        x = found$x_0,
        y = found$y_0,
        circles = rep(0.01, length(found$x_0)), # Radius is half the width
        fg = "green"
    )

    for (i in seq_len(nrow(true))) {
        rect(
            xleft = true$x_0[i] - true$lambda[i] / 2,
            ybottom = 0,
            xright = true$x_0[i] + true$lambda[i] / 2,
            ytop = par("usr")[4], # Full y height
            col = rgb(0.5, 0.5, 0.5, 0.1), # Transparent grey
            border = NA
        )
    }

    found$x_0_closest <- sapply(found$x_0, function(x_0) which.min(abs(true$x_0 - x_0)))
    true$x_0_assign <- sapply(seq_len(nrow(true)), function(i) {
        collapse(which(found$x_0_closest == i))
    })
    plot_spectrum(found, foc_rgn = c(0, 1), sub_show = FALSE)

    # Draw vertical lines at T$x_0
    opar <- par(xpd = TRUE)
    on.exit(par(opar), add = TRUE)
    y_range <- range(par("usr")[3:4])
    for (i in seq_along(T$x_0)) {
        x0 <- T$x_0[i]
        y0 <- y_range[1]
        x1 <- T$x_0[i]
        y1 <- y_range[2]
        segments(x0, y0, x1, y1, lty = 2)
        # text(T$x_0[i], y_range[1], labels = i, pos = 3, col = "blue")
    }
}

pairwise_identical <- function(x) {
    lapply(seq_len(length(x) - 1), function(i) identical(x[[i]], x[[i + 1]]))
}

#' @noRd
#' @title Compare two vectors
#' @description Checks if `x` and `y` are identical, all.equal or different
#' vectors. If differences are greated than expected, the differing elements are
#' printed.
#' @param x First vector.
#' @param y Second vector.
#' @param xpct Expected result. 0==identical, 1==all.equal, 2==different,
#' 3==error.
#' @param silent Logical indicating whether to print the results.
#' @return 0==identical, 1==all.equal, 2==different, 3==error
#' @details The function compares the vectors `x` and `y` and prints the
#' results. If the vectors are different, the differing elements are printed as
#' follows:
#'
#' ```R
#'              x       y  i  a                b                     z
#' class  integer numeric  .  .                .                     .
#' length      10      10  2  2                2                     2
#' head        1L       1 5L 5L                5 -8.88178419700125e-16
#' tail       10L      10 9L 9L 8.99999999999999  1.06581410364015e-14
#' ```
#'
#' Where
#'
#' - `i` == the indices of the differing elements
#' - `a` == the differing elements from `x`
#' - `b` == the differing elements from `y`
#' - `z` == the difference between `a` and `b`
#'
#' @examples
#' x <- 1:10
#' y1 <- 1:10
#' y2 <- c(1:4, 5.000000000000001, 6:8, 8.99999999999999, 10)
#' y3 <- c(1:4, 5.1, 6:8, 8.9, 10)
#' y4 <- list("a", "b", "c")
#' vcomp(x, y1)
#' vcomp(x, y2)
#' vcomp(x, y3)
#' vcomp(x, y4)
vcomp <- function(x, y, xpct = 0, silent = FALSE) {
    isvec <- function(x) is.vector(x) && !is.list(x)
    callstr <- paste(deparse(sys.call()), collapse = "")
    callstr <- gsub("\\s+", " ", callstr)
    status <- if (identical(x, y)) 0 else if (isTRUE(all.equal(x, y))) 1 else if (isvec(x) && isvec(y)) 2 else 3
    if (!silent) {
        msg <- c("identical", "all.equal", "different", "error")[status + 1]
        if (status == 3) msg <- paste0(msg, ": class(x)==", class(x), ", class(y)==", class(y))
        col <- c(GREEN, YELLOW, RED, RED)[status + 1]
        cat2(callstr, " ", col, msg, RESET, sep = "")
        if (status %in% 1:2 && status != xpct) {
            line <- "----------------------------------------"
            cat2(col, line, RESET, sep = "")
            i <- which(x != y)
            a <- x[i]
            b <- y[i]
            z <- a - b # nolint: object_usage_linter
            rows <- c("class", "length", "head", "tail")
            cols <- c("x", "y", "i", "a", "b", "z")
            objs <- list(x = x, y = y, i = i, a = a, b = b, z = z)
            df <- as.data.frame(matrix(".", length(rows), length(cols), dimnames = list(rows, cols)))
            df["class", c("x", "y")] <- c(class(x), class(y))
            df["length", c("x", "y", "i")] <- c(length(x), length(y), length(i))
            df["length", c("a", "b", "z")] <- c(".", ".", ".")
            df["head", ] <- sapply(objs, function(obj) dput2(head(obj, 1)))
            df["tail", ] <- sapply(objs, function(obj) dput2(tail(obj, 1)))
            print(df)
            cat2(col, line, RESET, sep = "")
        }
    }
    invisible(status)
}

#' @noRd
#' @description Compares a spectrum deconvoluted with
#' [generate_lorentz_curves_v12()] with a spectrum deconvoluted with
#' [MetaboDecon1D()].
#' @param x Result of [generate_lorentz_curves_v12()].
#' @param y Result of [MetaboDecon1D()].
#' @examples \donttest{
#' sim_01 <- metabodecon_file("sim_01")[1]
#' new <- generate_lorentz_curves_sim(sim_01)
#' old <- MetaboDecon1D(sim_01)
#' compare_spectra(new, old)
#' }
compare_spectra <- function(new, old, silent = FALSE) {
    # Define comparison function
    comp <- vcomp
    arglist <- formals(comp)
    arglist$silent <- silent
    formals(comp) <- arglist
    r <- numeric()

    # Extract old data
    o1 <- old$debuglist$args # nolint: object_usage_linter.
    o2 <- old$debuglist$data
    o3 <- old$debuglist$wsrm
    o4 <- old$debuglist$nvrm
    o5 <- old$debuglist$smooth
    o6 <- old$debuglist$peaksel
    o7 <- old$debuglist$peakfilter
    o8 <- old$debuglist$peakscore
    o9 <- old$debuglist$parinit
    o10 <- old$debuglist$parapprox
    o11 <- old$debuglist$parsave # nolint: object_usage_linter.

    # [1-5] spectra <- read_spectra(data_path, file_format, expno, procno, ask, sf, bwc = TRUE)
    # [6-7] spectra <- get_sfr(spectra, sfr, ask, adjno)
    # [8-9] spectra <- get_wshw(spectra, wshw, ask, adjno)
    # [10] spec <- rm_water_signal(spec)
    # [11] spec <- rm_negative_signals(spec)
    # [12] spec <- smooth_signals(spec, reps = smopts[1], k = smopts[2])
    r[1] <- comp(new$y_raw, o2$spectrum_y_raw)
    r[2] <- comp(new$y_scaled, o2$spectrum_y)
    r[3] <- comp(new$n, o3$spectrum_length)
    r[4] <- comp(new$sdp, o3$spectrum_x)
    r[5] <- comp(new$ppm, o3$spectrum_x_ppm)
    r[6] <- comp(new$sfr$left_sdp, o3$signal_free_region_left)
    r[7] <- comp(new$sfr$right_sdp, o3$signal_free_region_right)
    r[8] <- comp(new$wsr$left_dp, o3$water_signal_left)
    r[9] <- comp(new$wsr$right_dp, o3$water_signal_right)
    r[10] <- comp(new$y_nows, o3$spectrum_y)
    r[11] <- comp(new$y_pos, o4$spectrum_y)
    r[12] <- comp(new$y_smooth, o5$spectrum_y)

    # spec <- find_peaks(spec)
    r[13] <- comp(new$n, o6$spectrum_length)
    r[14] <- comp(new$d, c(NA, o6$second_derivative[2, ], NA)) # (1)
    r[16] <- comp(new$peak$center, as.integer(o6$peaks_index + 1)) # (1)
    r[17] <- comp(new$peak$center, as.integer(o6$peaks_x + 1)) # (1)
    r[18] <- comp(new$peak$right, as.integer(o6$left_position[1, ]) + 1) # (1)
    r[19] <- comp(new$peak$left, as.integer(o6$right_position[1, ]) + 1) # (1)
    # (1) MetaboDecon1D did not store NAs at the border, which is bad, because
    # you need to shift every index by one when you switch from
    # `second_derivative` to any other vector like `x_ppm` or `y_au`.

    # spec <- filter_peaks(spec, delta)
    border_is_na <- which(is.na(new$peak$left) | is.na(new$peak$right)) # (A)
    # (A) The original MetaboDecon1D implementation throws away NAs, so for
    # comparsion we need to do the same
    new_peak2 <- if (length(border_is_na) > 0) new$peak[-border_is_na, ] else new$peak
    r[20] <- comp(new_peak2$center, as.integer(o7$peaks_index + 1)) # (1) see above
    r[21] <- comp(new_peak2$right, as.integer(o7$left_position) + 1) # (1)
    r[22] <- comp(new_peak2$left, as.integer(o7$right_position) + 1) # (1)
    r[23] <- comp(new$sdp[new_peak2$center], o7$peaks_x)
    r[24] <- comp(mean(new_peak2$score[new_peak2$region %in% c("sfrl", "sfrr")]), o8$mean_score)
    r[25] <- comp(sd(new_peak2$score[new_peak2$region %in% c("sfrl", "sfrr")]), o8$sd_score)
    r[26] <- comp(new_peak2$score, o8$scores[1, ])
    r[27] <- comp(which(new_peak2$region == "sfrl"), o8$index_left)
    r[28] <- comp(which(new_peak2$region == "sfrr"), o8$index_right)
    r[29] <- comp(new$peak$center[new$peak$high], as.integer(o8$filtered_peaks + 1)) # (1)
    r[30] <- comp(new$peak$right[new$peak$high], o8$filtered_left_position + 1) # (1)
    r[31] <- comp(new$peak$left[new$peak$high], o8$filtered_right_position + 1) # (1)
    r[32] <- comp(new$peak$score[new$peak$high], o8$save_scores)

    # spec <- init_lorentz_curves_v12(spec)
    r[33] <- comp(new$sdp, o9$spectrum_x, xpct = 0)
    r[34] <- comp(new$y_smooth, o9$spectrum_y, xpct = 0)
    r[35] <- comp(new$peak$center[new$peak$high], as.integer(o9$filtered_peaks + 1), xpct = 0) # (1)
    r[36] <- comp(new$peak$right[new$peak$high], o9$filtered_left_position + 1, xpct = 0) # (1)
    r[37] <- comp(new$peak$left[new$peak$high], o9$filtered_right_position + 1, xpct = 0) # (1)
    r[38] <- comp(new$lc$A, o9$A, xpct = 1)
    r[39] <- comp(new$lc$lambda, o9$lambda, xpct = 1)
    r[40] <- comp(new$lc$w, o9$w, xpct = 0)

    # spec <- refine_lorentz_curves_v12(spec, nfit)
    r[42] <- comp(new$lcr$w_new, o10$w_new, xpct = 0)
    r[43] <- comp(new$lcr$lambda_new, o10$lambda_new, xpct = 1)
    r[44] <- comp(new$lcr$A_new, o10$A_new, xpct = 1)
    r[45] <- comp(new$lcr$integrals[1, ], o10$integrals[1, ], xpct = 1)

    # spec <- add_return_list_v12(spec, n, nam, debug)
    r[46] <- comp(new$ret$number_of_files, as.integer(old$number_of_files)) # int 1
    r[47] <- comp(new$ret$filename, old$filename) # chr[1] "urine_1"
    r[48] <- comp(new$ret$x_values, old$x_values) # num[131072] 131.071 131.070 131.069 ...
    r[49] <- comp(new$ret$x_values_ppm, old$x_values_ppm) # num[131072] 14.80254 14.80239 14.80223 ...
    r[50] <- comp(new$ret$y_values, old$y_values) # num[131072] 0.000831 0.000783 0.000743 ...
    r[51] <- comp(new$ret$spectrum_superposition, old$spectrum_superposition[1, ], 1) # num [1:131072] 3.29e-05 3.29e-05 3.29e-05 3.29e-05 3.29e-05 ...
    r[52] <- comp(new$ret$mse_normed, old$mse_normed) # num 4.1e-11
    r[53] <- comp(new$ret$index_peak_triplets_middle, as.integer(old$index_peak_triplets_middle)) # int[1227] 36158 37148 37418 37434 38942 38959 39030 39046 39091 39109 ...
    r[54] <- comp(new$ret$index_peak_triplets_left, old$index_peak_triplets_left) # num[1227] 36160 37159 37422 37437 38948 ...
    r[55] <- comp(new$ret$index_peak_triplets_right, old$index_peak_triplets_right) # num[1227] 36155 37139 37414 37431 38937 ...
    r[59] <- comp(new$ret$integrals, old$integrals[1, ], 1) # num[1227] 0.000502 0.026498 0.000396 0.000379 0.00732 ...
    if (!is.null(old$signal_free_region)) r[60] <- comp(new$ret$signal_free_region, old$signal_free_region) # num[2] 109.1 21.9
    if (!is.null(old$range_water_signal_ppm)) r[61] <- comp(new$ret$range_water_signal_ppm, old$range_water_signal_ppm) # num 0.153
    r[62] <- comp(new$ret$A, old$A, 1) # num[1227] -0.00016 -0.008437 -0.000126 -0.000121 -0.00233 ...
    r[63] <- comp(new$ret$lambda, old$lambda, 1) # num[1227] -0.00775 -0.02188 -0.00672 -0.00566 -0.01252 ...
    r[64] <- comp(new$ret$x_0, old$x_0, 0) # num[1227] 94.9 93.9 93.7 93.6 92.1 ...

    r[56] <- comp(new$ret$peak_triplets_middle, old$peak_triplets_middle) # NULL
    r[57] <- comp(new$ret$peak_triplets_left, old$peak_triplets_left) # NULL
    r[58] <- comp(new$ret$peak_triplets_right, old$peak_triplets_right) # NULL

    r[is.na(r)] <- 4
    msg <- "Identical: %s, Equal: %s, Different: %s, Error: %s, Empty: %s"
    if (!silent) {
        logf(msg, sum(r == 0), sum(r == 1), sum(r == 2), sum(r == 3), sum(r == 4))
    }

    invisible(r)
}

#' @noRd
#' @description Compares a spectrum deconvoluted with
#' [generate_lorentz_curves_v12()] with a spectrum deconvoluted with
#' [MetaboDecon1D()].
#' @param x Result of [generate_lorentz_curves_v12()].
#' @param y Result of [MetaboDecon1D()].
#' @examples
#' sim_01 <- metabodecon_file("sim_01")[1]
#' new <- generate_lorentz_curves_sim(sim_01)
#' old <- MetaboDecon1D(sim_01)
#' r <- compare_spectra_v13(new, old, silent = FALSE)
compare_spectra_v13 <- function(new, old, silent = FALSE) {

    # styler: off
    # Define comparison functions
    ident <- update_defaults(vcomp, xpct = 0, silent = silent)
    equal <- update_defaults(vcomp, xpct = 1, silent = silent)

    # Read and smooth spectrum
    o2 <- old$debuglist$data
    o3 <- old$debuglist$wsrm
    o4 <- old$debuglist$nvrm
    o5 <- old$debuglist$smooth
    r    <-  ident(new$y_raw,         o2$spectrum_y_raw)
    r[2] <-  ident(new$y_scaled,      o2$spectrum_y)
    r[3] <-  ident(new$n,             o3$spectrum_length)
    r[4] <-  ident(new$sdp,           o3$spectrum_x)
    r[5] <-  ident(new$ppm,           o3$spectrum_x_ppm)
    r[6] <-  ident(new$sfr$left_sdp,  o3$signal_free_region_left)
    r[7] <-  ident(new$sfr$right_sdp, o3$signal_free_region_right)
    r[8] <-  ident(new$wsr$left_dp,   o3$water_signal_left)
    r[9] <-  ident(new$wsr$right_dp,  o3$water_signal_right)
    r[10] <- ident(new$y_nows,        o3$spectrum_y)
    r[11] <- ident(new$y_pos,         o4$spectrum_y)
    r[12] <- ident(new$y_smooth,      o5$spectrum_y)

    # Find peaks
    o6 <- old$debuglist$peaksel
    r[13] <- ident(new$n,           o6$spectrum_length)
    r[16] <- equal(new$peak$center, o6$peaks_x + 1) # (1)
    r[15] <- equal(new$peak$center, o6$peaks_index + 1)
    r[17] <- equal(new$peak$right,  o6$left_position[1, ] + 1)
    r[18] <- equal(new$peak$left,   o6$right_position[1, ] + 1)
    r[14] <- equal(new$d,           c(NA, o6$second_derivative[2, ], NA))

    # Filter peaks
    o7 <- old$debuglist$peakfilter
    o8 <- old$debuglist$peakscore
    border_is_na <- which(is.na(new$peak$left) | is.na(new$peak$right)) # (2)
    new_peak2    <- if (length(border_is_na) > 0) new$peak[-border_is_na, ] else new$peak
    left_pos     <- if (is.matrix(o7$left_position)) o7$left_position[1, ] else o7$left_position
    right_pos    <- if (is.matrix(o7$right_position)) o7$right_position[1, ] else o7$right_position
    new_idx_left   <- which(new_peak2$region == "sfrl")
    new_idx_right  <- which(new_peak2$region == "sfrr")
    new_idx_sfr    <- c(new_idx_left, new_idx_right)
    new_mean_score <- mean(c(new_peak2$score[new_idx_sfr]))
    new_mean_sd    <- sd(new_peak2$score[new_idx_sfr])
    r[19] <- equal(new_peak2$center,               o7$peaks_index + 1) # (1)
    r[20] <- equal(new_peak2$right,                left_pos + 1)
    r[21] <- equal(new_peak2$left,                 right_pos + 1)
    r[22] <- equal(new$sdp[new_peak2$center],      o7$peaks_x)
    r[25] <- ident(new_peak2$score,                o8$scores[1, ])
    r[26] <- ident(new_idx_left,                   o8$index_left)
    r[27] <- ident(new_idx_right,                  o8$index_right)
    r[23] <- ident(new_mean_score,                 o8$mean_score)
    r[24] <- ident(new_mean_sd,                    o8$sd_score)
    r[28] <- equal(new$peak$center[new$peak$high], o8$filtered_peaks + 1) # (1)
    r[29] <- ident(new$peak$right[new$peak$high],  o8$filtered_left_position + 1)
    r[30] <- ident(new$peak$left[new$peak$high],   o8$filtered_right_position + 1)
    r[31] <- ident(new$peak$score[new$peak$high],  o8$save_scores)

    # Lorent Curve init values
    o9    <- old$debuglist$parinit
    r[32] <- equal(new$lci$P$ic,   o9$filtered_peaks + 1) # (1)
    r[33] <- ident(new$lci$P$ir,   o9$filtered_left_position + 1)
    r[34] <- ident(new$lci$P$il,   o9$filtered_right_position + 1)
    r[35] <- equal(new$lci$A,      o9$A)
    r[36] <- equal(new$lci$lambda, o9$lambda)
    r[37] <- equal(new$lci$w,      o9$w)

    # Lorent Curve refined values
    o10 <- old$debuglist$parapprox
    r[38] <- ident(new$lcr$w,         o10$w_new         )
    r[39] <- equal(new$lcr$lambda,    o10$lambda_new    )
    r[40] <- equal(new$lcr$A,         o10$A_new         )
    r[41] <- equal(new$lcr$integrals, o10$integrals[1, ])

    # Return List
    old_sfr <- old$signal_free_region %||% new$ret$signal_free_region
    old_rws <- old$range_water_signal_ppm %||% new$ret$range_water_signal_ppm
    r[42] <- equal(new$ret$number_of_files,            old$number_of_files)
    r[43] <- ident(new$ret$filename,                   old$filename)
    r[44] <- ident(new$ret$x_values,                   old$x_values)
    r[45] <- ident(new$ret$x_values_ppm,               old$x_values_ppm)
    r[46] <- ident(new$ret$y_values,                   old$y_values)
    r[47] <- equal(new$ret$spectrum_superposition,     old$spectrum_superposition[1, ])
    r[48] <- ident(new$ret$mse_normed,                 old$mse_normed)
    r[49] <- equal(new$ret$index_peak_triplets_middle, old$index_peak_triplets_middle)
    r[50] <- equal(new$ret$index_peak_triplets_left,   old$index_peak_triplets_left)
    r[51] <- equal(new$ret$index_peak_triplets_right,  old$index_peak_triplets_right)
    r[52] <- ident(new$ret$peak_triplets_middle,       old$peak_triplets_middle)
    r[53] <- ident(new$ret$peak_triplets_left,         old$peak_triplets_left)
    r[54] <- ident(new$ret$peak_triplets_right,        old$peak_triplets_right)
    r[55] <- equal(new$ret$integrals,                  old$integrals[1, ])
    r[56] <- equal(new$ret$A,                          old$A)
    r[57] <- equal(new$ret$lambda,                     old$lambda)
    r[58] <- ident(new$ret$x_0,                        old$x_0)
    r[59] <- equal(new$ret$signal_free_region,         old_sfr)
    r[60] <- equal(new$ret$range_water_signal_ppm,     old_rws)

    # Print summary
    if (!silent) {
        msg <- "Identical: %s, Equal: %s, Different: %s, Error: %s, Empty: %s"
        logf(msg, sum(r == 0), sum(r == 1), sum(r == 2), sum(r == 3), sum(r == 4))
    }

    # styler: on
    # (1) MetaboDecon1D did not store NAs at the border, which is bad, because
    #     you need to shift every index by one when you switch from
    #     `second_derivative` to any other vector like `x_ppm` or `y_au`.
    # (2) The original MetaboDecon1D implementation throws away NAs, so for
    #     comparsion we need to do the same

    # Return results
    r[is.na(r)] <- 4
    invisible(r)
}

#' @noRd
#' @description Helper of [compare_spectra_v13()].
update_defaults <- function(func, ...) {
    kwargs <- list(...)
    defaults <- formals(func)
    for (name in names(kwargs)) {
        defaults[[name]] <- kwargs[[name]]
    }
    formals(func) <- defaults
    func
}
