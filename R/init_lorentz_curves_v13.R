# x <- spec$sdp
# y <- spec$y_smooth
# pc <- spec$peak$center[spec$peak$high]
# pl <- spec$peak$right[spec$peak$high]
# pr <- spec$peak$left[spec$peak$high]
init_lorentz_curves_v13 <- function(x, y, pc, pl, pr) {

    msg("Initializing Lorentz curves")

    xl <- x[pl] # position of peak triplets
    xc <- x[pc]
    xr <- x[pr]

    yl <- y[pl] # intensity of peak triplets
    yc <- y[pc]
    yr <- y[pr]

    # Calculate mirrored points for ascending/descending shoulders
    i <- which((yl < yc) & (yc < yr)) # ascending shoulders
    j <- which((yl > yc) & (yc > yr)) # descending shoulders
    xr[i] <- 2 * xc[i] - xl[i]
    xl[j] <- 2 * xc[j] - xr[j]
    yr[i] <- yl[i]
    yl[j] <- yr[j]

    # Move triplet to zero position
    xd <- xl
    xl <- xl - xd
    xc <- xc - xd
    xr <- xr - xd

    # Calculate difference of position of peak triplets
    xlc <- xl - xc
    xlr <- xl - xr
    xcr <- xc - xr

    # Calculate difference of intensity values of peak triplets
    ylc <- yl - yc
    ylr <- yl - yr
    ycr <- yc - yr

    t1 <- xl^2 * yl * ycr
    t2 <- xr^2 * yr * ylc
    t3 <- xc^2 * yc * ylr
    t4 <- 2 * xlc * yl * yc
    t5 <- 2 * xcr * yc * yr
    t6 <- 2 * xlr * yl * yr
    w <- (t1 + t2 - t3) / (t4 + t5 - t6) + xd
    w[is.nan(w)] <- 0 # If (t4 + t5 - t6) is 0, then w is NaN. In this case we set w to 0.

    lambda <- -((sqrt(abs((-xc^4 * yc^2 * ylr^2 - xl^4 * yl^2 * ycr^2 - xr^4 * ylc^2 * yr^2 + 4 * xc * xr^3 * yc * ((-yl) + yc) * yr^2 + 4 * xc^3 * xr * yc^2 * yr * ((-yl) + yr) + 4 * xl^3 * yl^2 * ycr * (xc * yc - xr * yr) + 4 * xl * yl * (xc^3 * yc^2 * ylr - xc * xr^2 * yc * (yl + yc - 2 * yr) * yr + xr^3 * ylc * yr^2 - xc^2 * xr * yc * yr * (yl - 2 * yc + yr)) + 2 * xc^2 * xr^2 * yc * yr * (yl^2 - 3 * yc * yr + yl * (yc + yr)) + 2 * xl^2 * yl * (-2 * xc * xr * yc * yr * (-2 * yl + yc + yr) + xr^2 * yr * (yl * (yc - 3 * yr) + yc * (yc + yr)) + xc^2 * yc * (yl * (-3 * yc + yr) + yr * (yc + yr)))))))) / (2 * sqrt((xl * yl * ycr + xr * ylc * yr + xc * yc * ((-yl) + yr))^2))
    lambda[is.nan(lambda)] <- 0 

    A <- (-4 * xlc * xlr * xcr * yl * yc * yr * (xl * yl * ycr + xr * yr * ylc + xc * yc * (-ylr)) * lambda) / (xlc^4 * yl^2 * yc^2 - 2 * xlc^2 * yl * yc * (xlr^2 * yl + xcr^2 * yc) * yr + (xlr^2 * yl - xcr^2 * yc)^2 * yr^2)
    A[is.nan(A)] <- 0

    # TODO: break up calculation of w, lambda, and A and add explanation

    lc <- list(A = A, lambda = lambda, w = w, w_delta = xd)
    lc
}
