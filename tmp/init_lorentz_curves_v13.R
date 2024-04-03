# x <- spec$sdp
# y <- spec$y_smooth
# pc <- spec$peak$center[spec$peak$high]
# pl <- spec$peak$right[spec$peak$high]
# pr <- spec$peak$left[spec$peak$high]
init_lorentz_curves_v13 <- function(spec) {

    x <- spec$sdp; y <- spec$y_smooth; p <- spec$peak
    pl <- p$right[p$high]; pc <- p$center[p$high]; pr <- p$left[p$high];

    cat3("Initializing Lorentz curves")

    p <- length(pc); n <- length(x)
    wl <- x[pl]; wc <- x[pc]; wr <- x[pr] # position of peak triplets
    yl <- y[pl]; yc <- y[pc]; yr <- y[pr] # intensity of peak triplets
    P1 <- data.frame(pl = pl, pc = pc, pr = pr, wl = wl, wc = wc, wr = wr, yl = yl, yc = yc, yr = yr)

    # Replace ascending/descending shoulders (i/j) with mirrored points
    i <- (yl < yc) & (yc < yr)
    j <- (yl > yc) & (yc > yr)
    wr[i] <- 2 * wc[i] - wl[i]
    wl[j] <- 2 * wc[j] - wr[j]
    yr[i] <- yl[i];
    yl[j] <- yr[j]
    P2 <- data.frame(pl = pl, pc = pc, pr = pr, wl = wl, wc = wc, wr = wr, yl = yl, yc = yc, yr = yr, i = i, j = j)

    # Express wr/wc/wl as "distance to right border" and calculate x and y distances between triplet positions
    wd <- wl; wl <- wl - wd; wc <- wc - wd; wr <- wr - wd
    xlc <- wl - wc; xlr <- wl - wr; xcr <- wc - wr
    ylc <- yl - yc; ylr <- yl - yr; ycr <- yc - yr
    P3 <- data.frame(pl = pl, pc = pc, pr = pr, wl = wl, wc = wc, wr = wr, yl = yl, yc = yc, yr = yr, i = i, j = j, xlc = xlc, xlr = xlr, xcr = xcr, ylc = ylc, ylr = ylr, ycr = ycr, wd = wd)

    # Calculate parameters of Lorentz Curves from peak triplets
    w <- calc_w(wl, wc, wr, yl, yc, yr, xlc, xlr, xcr, ylc, ylr, ycr, wd)
    lambda <- calc_lambda(wl, wc, wr, yl, yc, yr, xlc, xlr, xcr, ylc, ylr, ycr, wd)
    A <- calc_A(wl, wc, wr, yl, yc, yr, xlc, xlr, xcr, ylc, ylr, ycr, lambda, w, wd)
    Y <- matrix(nrow = n, ncol = p) # 131072 x 1227 for urine_1
    for (j in seq_along(pc)) Y[, j] <- if (A[j] == 0) 0 else abs(A[j] * (lambda[j] / (lambda[j]^2 + (x - w[j])^2)))

    cat3("Done")

    list(P1 = P1, P2 = P2, P3 = P3, A = A, lambda = lambda, w = w, w_delta = wd, Y = Y)
}

calc_w <- function(wl, wc, wr, yl, yc, yr, xlc, xlr, xcr, ylc, ylr, ycr, wd) {
    t1 <- wl^2 * yl * ycr
    t2 <- wr^2 * yr * ylc
    t3 <- wc^2 * yc * ylr
    t4 <- 2 * xlc * yl * yc
    t5 <- 2 * xcr * yc * yr
    t6 <- 2 * xlr * yl * yr
    w <- (t1 + t2 - t3) / (t4 + t5 - t6) + wd
    w[is.nan(w)] <- 0 # If (t4 + t5 - t6) is 0, then w is NaN. In this case we set w to 0.
    w
}

calc_lambda <- function(wl, wc, wr, yl, yc, yr, xlc, xlr, xcr, ylc, ylr, ycr, wd) {
    lambda <- -((sqrt(abs((-wc^4 * yc^2 * ylr^2 - wl^4 * yl^2 * ycr^2 - wr^4 * ylc^2 * yr^2 + 4 * wc * wr^3 * yc * ((-yl) + yc) * yr^2 + 4 * wc^3 * wr * yc^2 * yr * ((-yl) + yr) + 4 * wl^3 * yl^2 * ycr * (wc * yc - wr * yr) + 4 * wl * yl * (wc^3 * yc^2 * ylr - wc * wr^2 * yc * (yl + yc - 2 * yr) * yr + wr^3 * ylc * yr^2 - wc^2 * wr * yc * yr * (yl - 2 * yc + yr)) + 2 * wc^2 * wr^2 * yc * yr * (yl^2 - 3 * yc * yr + yl * (yc + yr)) + 2 * wl^2 * yl * (-2 * wc * wr * yc * yr * (-2 * yl + yc + yr) + wr^2 * yr * (yl * (yc - 3 * yr) + yc * (yc + yr)) + wc^2 * yc * (yl * (-3 * yc + yr) + yr * (yc + yr)))))))) / (2 * sqrt((wl * yl * ycr + wr * ylc * yr + wc * yc * ((-yl) + yr))^2))
    lambda[is.nan(lambda)] <- 0
    lambda
}

calc_A <- function(wl, wc, wr, yl, yc, yr, xlc, xlr, xcr, ylc, ylr, ycr, lambda, w, wd) {
    A <- (-4 * xlc * xlr * xcr * yl * yc * yr * (wl * yl * ycr + wr * yr * ylc + wc * yc * (-ylr)) * lambda) / (xlc^4 * yl^2 * yc^2 - 2 * xlc^2 * yl * yc * (xlr^2 * yl + xcr^2 * yc) * yr + (xlr^2 * yl - xcr^2 * yc)^2 * yr^2)
    A[is.nan(A)] <- 0
    A
}