library(testthat)

skip("Manual execution only")

simx_pick_idx <- function(npk, lo, hi, mingap) {
    repeat {
        idx <- sort(sample(lo:hi, npk))
        if (all(diff(idx) >= mingap)) {
            return(idx)
        }
    }
}

simx_make_one <- function(name, grp, cs, x0, A, lambda, didx, csres) {
    shift <- sample(-50:50, 1) * csres
    x0i <- x0 + shift + rnorm(length(x0), sd = 2 * csres)
    Ai <- A * exp(rnorm(length(A), sd = 0.08))
    li <- lambda * exp(rnorm(length(lambda), sd = 0.05))
    if (grp == "AKI") {
        x0i[didx] <- x0i[didx] + c(-4, 5, -3, 4) * csres
        Ai[didx] <- Ai[didx] * c(2.6, 0.35, 1.9, 0.45)
        li[didx] <- li[didx] * c(1.20, 0.85, 1.15, 0.90)
    }
    simulate_spectrum(
        name = name,
        cs = cs,
        x0 = x0i,
        A = Ai,
        lambda = li,
        noise = rnorm(length(cs), sd = 700)
    )
}

make_simx <- function(seed = 1, n_per_grp = 10) {
    set.seed(seed)
    csres <- 0.00015
    cs <- seq(from = 3.59, length.out = 2048, by = -csres)
    idx <- simx_pick_idx(npk = 12, lo = 260, hi = 1780, mingap = 85)
    x0 <- cs[idx]
    A <- runif(length(x0), 5, 18) * 1e3
    lambda <- runif(length(x0), 0.9, 1.3) / 1e3
    didx <- c(2, 5, 8, 11)

    ctrl <- lapply(seq_len(n_per_grp), function(i) {
        simx_make_one(
            name = sprintf("simx_ctrl_%02d", i),
            grp = "Control",
            cs = cs,
            x0 = x0,
            A = A,
            lambda = lambda,
            didx = didx,
            csres = csres
        )
    })
    aki <- lapply(seq_len(n_per_grp), function(i) {
        simx_make_one(
            name = sprintf("simx_aki_%02d", i),
            grp = "AKI",
            cs = cs,
            x0 = x0,
            A = A,
            lambda = lambda,
            didx = didx,
            csres = csres
        )
    })
    sp <- as_spectra(c(ctrl, aki))
    y <- factor(rep(c("Control", "AKI"), each = n_per_grp))
    list(
        spectra = sp,
        y = y,
        sfr = c(3.55, 3.35),
        x0_case = x0[didx],
        csres = csres
    )
}

get_nonzero_features <- function(m) {
    cf <- as.matrix(coef(m))
    keep <- cf[, 1] != 0
    nm <- rownames(cf)[keep]
    nm <- setdiff(nm, "(Intercept)")
    as.numeric(nm)
}

count_feature_hits <- function(obs, exp, tol) {
    hits <- vapply(exp, function(x0) {
        any(abs(obs - x0) <= tol)
    }, logical(1))
    sum(hits)
}

test_that("cv_mdlm and benchmark_mdlm work on simx", {
    simx <- make_simx(seed = 11, n_per_grp = 12)
    grid <- mdlm_grid(
        nfit = 5,
        smit = 2,
        smws = 5,
        delta = 4,
        maxShift = 50,
        maxCombine = 20
    )

    cv <- cv_mdlm(
        spectra = simx$spectra,
        y = simx$y,
        sfr = simx$sfr,
        use_rust = TRUE,
        grid = grid,
        nworkers = 2,
        verbose = FALSE,
        nfolds = 4,
        seed = 3
    )

    nz <- get_nonzero_features(cv)
    tol <- 20 * simx$csres
    nhit <- count_feature_hits(nz, simx$x0_case, tol)

    expect_s3_class(cv, "mdlm")
    expect_true(length(nz) > 0)
    expect_gte(nhit, 3)

    bm <- benchmark_mdlm(
        spectra = simx$spectra,
        y = simx$y,
        sfr = simx$sfr,
        nwo = 2,
        nwi = 2,
        verbose = FALSE,
        use_rust = TRUE,
        grid = grid,
        nfo = 4,
        nfi = 4,
        seed = 3
    )

    expect_length(bm$models, 4)
    expect_equal(nrow(bm$predictions), length(simx$y))
    expect_false(anyNA(bm$predictions$prob))
    expect_true(all(levels(bm$predictions$pred) == levels(simx$y)))
})
