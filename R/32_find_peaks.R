# Private Main #####

find_peaks_v12 <- function(spec) {
    logf("Starting peak selection")
    d <- spec$d <- calc_second_derivative_v12(y = spec$y_smooth)
    a <- abs(d)
    m <- length(d)
    dl <- c(NA, d[-m]) # dl[i] == d[i-1]
    dr <- c(d[-1], NA) # dr[i] == d[i+1]
    center <- which(d < 0 & d <= dl & d < dr)
    spec$peak <- data.frame(left = NA, center = center, right = NA, score = NA)
    for (i in seq_along(center)) {
        j <- center[i]
        l <- spec$peak$left[i] <- get_left_border_v12(j, d)
        r <- spec$peak$right[i] <- get_right_border_v12(j, d, m)
        spec$peak$score[i] <- get_peak_score_v12(j, l, r, a)
    }
    logf("Detected %d peaks", length(center))
    return(spec)
}

# Private Helpers #####

calc_second_derivative_v12 <- function(y) {
    n <- length(y)
    x <- c(NA, y[-n]) # x[i] == y[i-1]
    z <- c(y[-1], NA) # z[i] == y[i+1]
    d <- x + z - 2 * y
    d
}

get_right_border_v12 <- function(j, d, m) {
    r <- j + 1
    while (r < m) { # use r<m instead of r<=m because c4 requires d[r+1]
        c1 <- d[r] > d[r - 1]
        c2 <- d[r] >= d[r + 1]
        c3 <- d[r] < 0
        c4 <- d[r + 1] >= 0
        is_right_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_right_border)) {
            return(r)
        }
        r <- r + 1
    }
    return(NA)
}

get_left_border_v12 <- function(j, d) {
    l <- j - 1
    while (l > 1) { # use l>1 instead of l>=1 because c4 requires d[l-1]
        c1 <- d[l] > d[l + 1]
        c2 <- d[l] >= d[l - 1]
        c3 <- d[l] < 0
        c4 <- d[l - 1] >= 0
        is_left_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_left_border)) {
            return(l)
        }
        l <- l - 1
    }
    return(NA)
}

#' @noRd
#' @title Get Peak Score
#' @description Calculate the score of a peak based on the sum of absolute second derivative values of its datapoints.
#' @param j <- Index of the peak center
#' @param l <- Index of the left border
#' @param r <- Index of the right border
#' @param a <- Absolute values of the second derivative for all data points
#' @return The score of the peak.
#' @examples
#' y <- c( 0, 1, 2, 3, 4, 3, 3, 2, 1, 0, 0, 1, 2, 3, 1, 1, 0  )
#' #      ____________________________________________________
#' #     |____2________5___________9_______12____14____16_____|
#' #     |             x                                      |
#' #     |          x  x  x  x                    x           |
#' #     |       x  x  x  x  x  x              x  x           |
#' #     |_.__x__x__x__x__x__x__x__x__.__.__x__x__x__x__x__.__|
#' a <- c(NA, 0, 0, 0, 2, 1, 1, 0, 0, 1, 1, 0, 0, 3, 2, 1, NA )
#' #          |----2---|-----4-----|        |--3--|--6--|
#' all.equal(a, abs(calc_second_derivative_v12(y)))
#'
#' s1 <- get_peak_score_v12( 5, 2,   9, a)
#' s2 <- get_peak_score_v12(14, 12, 16, a)
#' stopifnot(s1 == min(sum(a[2:5]), sum(a[5:9])))
#' stopifnot(s2 == min(sum(a[12:14]), sum(a[14:16])))
get_peak_score_v12 <- function(j, l, r, a) {
    if (any(is.na(a[c(l, j, r)]))) {
        0
    } else {
        min(sum(a[l:j]), sum(a[j:r]))
    }
}
