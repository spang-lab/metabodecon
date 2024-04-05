refine_lorentz_curves_v13_example <- function(spec, nfit = 3) {
    x <- spec$sdp
    y <- spec$y_smooth
    pc <- as.integer(spec$peak$center[spec$peak$high] - 1)
    pl <- spec$peak$right[spec$peak$high] - 1
    pr <- spec$peak$left[spec$peak$high] - 1
    A <- spec$lc$A
    lambda <- spec$lc$lambda
    w <- spec$lc$w
    refine_lorentz_curves_v13(x, y, pc, pl, pr, A, lambda, w, nfit = 3)
}

refine_lorentz_curves_v13 <- function(x, y, pc, pl, pr, A, lambda, w, nfit = 3) {
    cat3("Refining Lorentz curves")

    p <- length(pc)
    n <- length(x)
    lc <- matrix(nrow = n, ncol = p) # 131072 x 1227 for urine_1
    for (j in seq_along(pc)) {
        lc[, j] <- if (A[j] == 0) 0 else lc[, j] <- abs(A[j] * (lambda[j] / (lambda[j]^2 + (x - w[j])^2)))
    }

    # Approximation of lorentz curves
    
    for (b in 1:nfit) {
        # Calculate new heights of peak triplets
        xl <- c()
        xc <- c()
        xr <- c()
        yl <- c()
        yc <- c()
        yr <- c()
        xlc <- c()
        xlr <- c()
        xcr <- c()
        ylc <- c()
        ylr <- c()
        ycr <- c()
        wd <- c()
        w <- c()
        l <- c()
        A_new <- c()
        sum_left <- c()
        sum_peaks <- c()
        sum_right <- c()
        proportion_left <- c()
        proportion_peaks <- c()
        proportion_right <- c()

        for (i in seq_along(pc)) {
            # Calculate the position of the peak triplets
            xl <- c(xl, x[pl[i] + 1])
            xc <- c(xc, x[pc[i] + 1])
            xr <- c(xr, x[pr[i] + 1])

            # Calculate the sum of all lorentz curves for each data point
            sum_left[i] <- sum(lc[seq_along(pl), pl[i] + 1])
            sum_peaks[i] <- sum(lc[seq_along(pc), pc[i] + 1])
            sum_right[i] <- sum(lc[seq_along(pr), pr[i] + 1])

            # Calculate the proprotion between original spectrum an the sum of the lorentz curves for each peak triplets position
            proportion_left[i] <- y[pl[i] + 1] / sum_left[i]
            proportion_peaks[i] <- y[pc[i] + 1] / sum_peaks[i]
            proportion_right[i] <- y[pr[i] + 1] / sum_right[i]

            # Calculate the new heights of the peak triplets
            yl[i] <- lc[i, pl[i] + 1] * proportion_left[i]
            yc[i] <- lc[i, pc[i] + 1] * proportion_peaks[i]
            yr[i] <- lc[i, pr[i] + 1] * proportion_right[i]

            # Calculate mirrored points if necesccary
            # For ascending shoulders
            if ((yl[i] < yc[i]) && (yc[i] < yr[i])) {
                xr[i] <- 2 * xc[i] - xl[i]
                yr[i] <- yl[i]
            }
            # For descending shoulders
            if ((yl[i] > yc[i]) && (yc[i] > yr[i])) {
                xl[i] <- 2 * xc[i] - xr[i]
                yl[i] <- yr[i]
            }

            # Move triplet to zero position
            wd[i] <- xl[i]
            xl[i] <- xl[i] - wd[i]
            xc[i] <- xc[i] - wd[i]
            xr[i] <- xr[i] - wd[i]

            # Calculate difference of peak triplet positions
            xlc <- c(xlc, xl[i] - xc[i])
            xlr <- c(xlr, xl[i] - xr[i])
            xcr <- c(xcr, xc[i] - xr[i])

            # Calculate difference of new intensity values of peak triplets
            ylc <- c(ylc, yl[i] - yc[i])
            ylr <- c(ylr, yl[i] - yr[i])
            ycr <- c(ycr, yc[i] - yr[i])

            # Calculate w for each peak triplet
            w_result <- (xl[i]^2 * yl[i] * ycr[i] + xr[i]^2 * yr[i] * ylc[i] + xc[i]^2 * yc[i] * (-ylr[i])) / (2 * xlc[i] * yl[i] * yc[i] - 2 * (xlr[i] * yl[i] + (-xcr[i]) * yc[i]) * yr[i])
            w_result <- w_result + wd[i]
            w <- c(w, w_result)

            # If y values are getting 0 after height adjustment, then w[i]=NaN
            if (is.nan(w[i])) {
                w[i] <- 0
            }

            # Calculate lambda for each peak triplet
            lambda_result <- -((sqrt(abs(((-xc[i]^4 * yc[i]^2 * ylr[i]^2 - xl[i]^4 * yl[i]^2 * ycr[i]^2 - xr[i]^4 * ylc[i]^2 * yr[i]^2 + 4 * xc[i] * xr[i]^3 * yc[i] * ((-yl[i]) + yc[i]) * yr[i]^2 + 4 * xc[i]^3 * xr[i] * yc[i]^2 * yr[i] * ((-yl[i]) + yr[i]) + 4 * xl[i]^3 * yl[i]^2 * ycr[i] * (xc[i] * yc[i] - xr[i] * yr[i]) + 4 * xl[i] * yl[i] * (xc[i]^3 * yc[i]^2 * ylr[i] - xc[i] * xr[i]^2 * yc[i] * (yl[i] + yc[i] - 2 * yr[i]) * yr[i] + xr[i]^3 * ylc[i] * yr[i]^2 - xc[i]^2 * xr[i] * yc[i] * yr[i] * (yl[i] - 2 * yc[i] + yr[i])) + 2 * xc[i]^2 * xr[i]^2 * yc[i] * yr[i] * (yl[i]^2 - 3 * yc[i] * yr[i] + yl[i] * (yc[i] + yr[i])) + 2 * xl[i]^2 * yl[i] * (-2 * xc[i] * xr[i] * yc[i] * yr[i] * (-2 * yl[i] + yc[i] + yr[i]) + xr[i]^2 * yr[i] * (yl[i] * (yc[i] - 3 * yr[i]) + yc[i] * (yc[i] + yr[i])) + xc[i]^2 * yc[i] * (yl[i] * (-3 * yc[i] + yr[i]) + yr[i] * (yc[i] + yr[i]))))))))) / (2 * sqrt((xl[i] * yl[i] * ycr[i] + xr[i] * ylc[i] * yr[i] + xc[i] * yc[i] * ((-yl[i]) + yr[i]))^2))

            # If y and w are 0, then 0/0=NaN
            if (is.nan(lambda_result)) {
                lambda_result <- 0
            }
            l <- c(l, lambda_result)

            # Calculate scaling factor A for each peak triplet
            A_result <- (-4 * xlc[i] * xlr[i] * xcr[i] * yl[i] * yc[i] * yr[i] * (xl[i] * yl[i] * ycr[i] + xr[i] * yr[i] * ylc[i] + xc[i] * yc[i] * (-ylr[i])) * l[i]) / (xlc[i]^4 * yl[i]^2 * yc[i]^2 - 2 * xlc[i]^2 * yl[i] * yc[i] * (xlr[i]^2 * yl[i] + xcr[i]^2 * yc[i]) * yr[i] + (xlr[i]^2 * yl[i] - xcr[i]^2 * yc[i])^2 * yr[i]^2)

            # If y and w are 0, then 0/0=NaN
            if (is.nan(A_result)) {
                A_result <- 0
            }
            A_new <- c(A_new, A_result)

            # Calculate new lorentz curves
            # If y values are zero, then lorentz curves should also be zero
            if ((w[i] == 0) || (l[i] == 0) || (A_new[i] == 0)) {
                lc[i, ] <- 0
            } else {
                lc[i, ] <- abs(A_new[i] * (l[i] / (l[i]^2 + (x - w[i])^2)))
            }
        }

        # Calculate sum of lorentz curves
        spectrum_approx <- matrix(nrow = 1, ncol = n)
        for (i in seq_along(x)) {
            spectrum_approx[1, i] <- sum(lc[seq_along(pc), i])
        }

        # Standardize the spectra
        y_normed <- y / sum(y)
        spectrum_approx_normed <- spectrum_approx / sum(spectrum_approx)

        # Calculate the difference between normed original spectrum and normed approximated spectrum
        difference_normed <- c()
        for (i in seq_along(x)) {
            difference_normed[i] <- (y_normed[i] - spectrum_approx_normed[i])^2
        }
        mse_normed <- (1 / length(difference_normed)) * sum(difference_normed)
        cat3(sprintf("Normed MSE after iteration %d: %.22f", b, mse_normed))
    }

    # Calculate the integrals for each lorentz curve
    integrals <- matrix(nrow = 1, ncol = length(l))
    for (i in seq_along(l)) {
        integrals[1, i] <- A_new[i] * (atan((-w[i] + (spec$n / spec$sfx)) / l[i]) - atan((-w[i]) / l[i]))
    }

    spec$lcr <- list(
        A_new = A_new,
        l = l,
        w = w,
        spectrum_approx = spectrum_approx,
        y_normed = y_normed,
        spectrum_approx_normed = spectrum_approx_normed,
        difference_normed = difference_normed,
        mse_normed = mse_normed,
        integrals = integrals
    )
    spec
}

# Helpers #####

benchmark_calc_lc <- function(spec, nfit = 3) {
    x <- spec$sdp
    A <- spec$lc$A
    lambda <- spec$lc$lambda
    w <- spec$lc$w
    fnams <- c("calc_lc_pxn", "calc_lc_nxp", "calc_lc_repx", "calc_lc_repx2", "calc_lc_cpp1", "calc_lc_cpp2", "calc_lc_cpp3")
    colors <- c("red", "blue", "green", "orange", "purple", "brown", "black")
    fobjs <- sapply(fnams, function(fnam) get(fnam))
    runtimes <- matrix(nrow = length(fnams), ncol = 5)
    for (i in seq_along(fnams)) {
        fobj <- fobjs[[i]]
        cat2("Benchmarking ", fnams[i])
        for (j in seq_len(5)) {
            cat3("Run ", j)
            runtimes[i, j] <- system.time(fobj(A, lambda, w, x))[3]
        }
    }
    runtimes_df <- as.data.frame(t(runtimes))
    names(runtimes_df) <- fnams

    # Create the boxplot with a log2 scale on the y-axis
    pdf("benchmark_calc_lc.pdf", width = 14, height = 7)
    boxplot(runtimes_df, main = "Runtime of different functions", xlab = "Function", ylab = "Runtime (s)", col = colors, log = "y")
    legend("bottomright", legend = fnams, col = colors, lty = 1)
    dev.off()
}

calc_lc_pxn <- function(A, lambda, w, x) {
    p <- length(A)
    n <- length(x)
    lc <- matrix(nrow = p, ncol = n) # 1227 x 131072 for urine_1
    for (i in seq_along(pc)) {
        if (A[i] == 0) { # If A = 0, then the lorentz curve is a zero line
            lc[i, ] <- 0
        } else {
            lc[i, ] <- abs(A[i] * (lambda[i] / (lambda[i]^2 + (x - w[i])^2)))
        }
    }
}

calc_lc_nxp <- function(A, lambda, w, x) {
    p <- length(A)
    n <- length(x)
    lc <- matrix(nrow = n, ncol = p) # 131072 x 1227 for urine_1
    for (j in seq_along(pc)) {
        if (A[j] == 0) { # If A = 0, then the lorentz curve is a zero line
            lc[, j] <- 0
        } else {
            lc[, j] <- abs(A[j] * (lambda[j] / (lambda[j]^2 + (x - w[j])^2)))
        }
    }
}

calc_lc_repx <- function(A, lambda, w, x) {
    p <- length(A)
    n <- length(x)
    xx <- rep(x, each = p)
    lcvec <- abs(A * (lambda / (lambda^2 + (xx - w)^2)))
    matrix(lcvec, nrow = p, ncol = n, byrow = FALSE)
}

calc_lc_repx2 <- function(A, lambda, w, x) {
    p <- length(A)
    n <- length(x)
    matrix(abs(A * (lambda / (lambda^2 + (rep(x, each = p) - w)^2))), nrow = p, ncol = n, byrow = FALSE)
}
