test_that("get_si_mat returns a matrix of the correct dimensions", {

    withr::local_output_sink(nullfile())

    # 'speaq' requires 'MassSpecWavelet' and 'impute' to be installed
    deps <- c("MassSpecWavelet", "impute")
    inst <- sapply(deps, requireNamespace, quietly = TRUE)
    if (!all(inst)) skip(paste("Missing deps:", collapse(deps[!inst])))

    decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
    aligns <- align(decons)
    si_mat <- get_si_mat(aligns)
    expect_equal(dim(si_mat), c(2048, 2))
    expect_equal(as.numeric(rownames(si_mat)), aligns[[1]]$cs)
    expect_equal(colnames(si_mat), get_names(aligns))
})

test_that("get_si_mat drop_zero removes all-zero rows", {

    withr::local_output_sink(nullfile())

    deps <- c("MassSpecWavelet", "impute")
    inst <- sapply(deps, requireNamespace, quietly = TRUE)
    if (!all(inst)) skip(paste("Missing deps:", collapse(deps[!inst])))

    decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
    aligns <- align(decons)
    full <- get_si_mat(aligns)
    compact <- get_si_mat(aligns, drop_zero = TRUE)
    expect_lt(nrow(compact), nrow(full))
    expect_true(all(rowSums(compact != 0) > 0))
    expect_equal(compact, full[rowSums(full != 0) > 0, , drop = FALSE])
})

test_that("get_si_mat with maxCombine returns reduced matrix", {

    withr::local_output_sink(nullfile())

    deps <- c("MassSpecWavelet", "impute")
    inst <- sapply(deps, requireNamespace, quietly = TRUE)
    if (!all(inst)) skip(paste("Missing deps:", collapse(deps[!inst])))

    decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
    al <- align(decons, verbose = FALSE)

    mat_raw   <- get_si_mat(al)
    mat_20dp  <- get_si_mat(al, maxCombine = 20, drop_zero = TRUE)
    mat_40dp  <- get_si_mat(al, maxCombine = 40, drop_zero = TRUE)

    # After combining and dropping all-zero rows, matrix should be smaller
    expect_lt(nrow(mat_20dp), nrow(mat_raw))
    expect_lt(nrow(mat_40dp), nrow(mat_raw))
    expect_equal(ncol(mat_20dp), ncol(mat_raw))
    expect_equal(ncol(mat_40dp), ncol(mat_raw))
})

make_aligns_for_get_si_mat_test <- function(cs, peaks_list, areas_list) {
    objs <- lapply(seq_along(peaks_list), function(i) {
        x0 <- cs[peaks_list[[i]]]
        al <- numeric(length(cs))
        al[peaks_list[[i]]] <- areas_list[[i]] * pi
        structure(list(
            cs = cs,
            lcpar = list(x0 = x0, x0al = x0, A = areas_list[[i]]),
            sit = list(al = al),
            meta = list(name = sprintf("spec_%d", i))
        ), class = c("align", "decon2", "spectrum"))
    })
    structure(objs, class = c("aligns", "decons2", "spectra"))
}

test_that("si_mat with peakPos reproduces the docs example", {
    # Reproduce the example from ?si_mat verbatim.
    # cs = 9:1 so colnames are decreasing and column index == position label.
    cs <- 9:1
    aligns <- make_aligns_for_get_si_mat_test(
        cs = cs,
        peaks_list = list(c(4, 8), c(3, 8), c(3, 4, 7), c(5, 9), c(5, 9)),
        areas_list = list(c(2, 4), c(3, 4), c(2, 4, 5), c(3, 3), c(2, 3))
    )

    mat <- si_mat(aligns, maxCombine = 1, peakPos = c(3, 4, 9))

    expected <- rbind(
        c(0, 0, 0, 2, 0, 0, 0, 0, 4),
        c(0, 0, 3, 0, 0, 0, 0, 0, 4),
        c(0, 0, 2, 4, 0, 0, 5, 0, 0),
        c(0, 0, 0, 3, 0, 0, 0, 0, 3),
        c(0, 0, 0, 2, 0, 0, 0, 0, 3)
    ) * pi
    rownames(expected) <- sprintf("spec_%d", 1:5)
    colnames(expected) <- as.character(cs)

    expect_equal(mat, expected)
})

test_that("si_mat with peakPos shifts only the closest peak on ties", {
    # Peaks at columns 7 and 8, peakPos=9, maxCombine=2: both have closest
    # peakPos = 9, but only the peak at 8 (dist 1) shifts; peak at 7 stays.
    cs <- 9:1
    aligns <- make_aligns_for_get_si_mat_test(
        cs = cs,
        peaks_list = list(c(7, 8)),
        areas_list = list(c(10, 20))
    )

    mat <- si_mat(aligns, maxCombine = 2, peakPos = 9)

    expected <- matrix(0, nrow = 1, ncol = 9)
    expected[1, 7] <- 10 * pi   # stays
    expected[1, 9] <- 20 * pi   # shifted from 8
    rownames(expected) <- "spec_1"
    colnames(expected) <- as.character(cs)

    expect_equal(mat, expected)
})

test_that("si_mat with peakPos preserves the full cs grid", {
    cs <- 9:1
    aligns <- make_aligns_for_get_si_mat_test(
        cs = cs,
        peaks_list = list(c(4, 8)),
        areas_list = list(c(2, 4))
    )

    mat <- si_mat(aligns, maxCombine = 1, peakPos = c(3, 4, 9))

    expect_equal(ncol(mat), length(cs))
    expect_equal(colnames(mat), as.character(cs))
})

test_that("si_mat round-trips peakPos: training grid matches prediction grid", {
    # Mimic the fit_mdm() / predict.mdm() pipeline:
    #   train: M_train <- si_mat(al, maxCombine=k); pp <- which(colSums(M_train != 0) > 0)
    #   test:  M_test  <- si_mat(al_new, maxCombine=k, peakPos=pp)
    # The training feature columns must equal those addressable in test output.
    cs <- 20:1
    train <- make_aligns_for_get_si_mat_test(
        cs = cs,
        peaks_list = list(c(3, 7, 13), c(4, 7, 14)),
        areas_list = list(c(1, 2, 3), c(4, 5, 6))
    )
    test <- make_aligns_for_get_si_mat_test(
        cs = cs,
        peaks_list = list(c(3, 8, 13)),
        areas_list = list(c(7, 8, 9))
    )

    M_train <- si_mat(train, maxCombine = 1)
    pp <- which(colSums(M_train != 0) > 0)
    M_test <- si_mat(test, maxCombine = 1, peakPos = pp)

    # Both matrices share the cs grid, so reducing to peakPos columns must
    # yield matrices with identical column labels.
    expect_equal(colnames(M_train[, pp, drop = FALSE]),
                 colnames(M_test[, pp, drop = FALSE]))
    # The test peak at column 8 is within maxCombine=1 of pp (column 7),
    # so it must snap there.
    expect_equal(M_test[1, 7], (8) * pi)
    expect_equal(M_test[1, 3], 7 * pi)
    expect_equal(M_test[1, 13], 9 * pi)
})
