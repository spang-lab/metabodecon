# Private Main #####

find_peaks_v12 <- function(spec) {
    cat3("Starting peak selection")
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
    cat3("Detected", length(center), "peaks")
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

get_peak_score_v12 <- function(j, l, r, a) {
    if (any(is.na(a[c(l, j, r)]))) {
        0
    } else {
        min(sum(a[l:j]), sum(a[j:r]))
    }
}
