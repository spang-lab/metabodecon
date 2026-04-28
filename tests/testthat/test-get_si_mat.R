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
    al <- align(decons, maxCombine = 0, verbose = FALSE)

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
            lcpar = list(x0 = x0, x0_al = x0, A = areas_list[[i]]),
            sit = list(al = al),
            meta = list(name = sprintf("spec_%d", i))
        ), class = "align")
    })
    structure(objs, class = "aligns")
}

test_that("get_si_mat with maxCombine sums areas within snapped intervals", {
    cs <- 10:1
    aligns <- make_aligns_for_get_si_mat_test(
        cs = cs,
        peaks_list = list(c(2, 5, 8), c(1, 3, 4, 6, 8, 10), c(2, 5, 7)),
        areas_list = list(c(10, 20, 30), c(1, 2, 3, 4, 5, 6), c(7, 8, 9))
    )

    mat <- get_si_mat(aligns, maxCombine = 1, ref = aligns[[1]])

    expected <- cbind(
        c(10, 20, 30),
        c(1 + 2, 3 + 4, 5),
        c(7, 8, 9)
    ) * pi
    rownames(expected) <- cs[c(2, 5, 8)]
    colnames(expected) <- c("spec_1", "spec_2", "spec_3")

    expect_equal(mat, expected)
})

test_that("get_si_mat with maxCombine assigns midpoint ties to the left", {
    cs <- 6:1
    aligns <- make_aligns_for_get_si_mat_test(
        cs = cs,
        peaks_list = list(c(2, 4), c(3)),
        areas_list = list(c(10, 20), c(5))
    )

    mat <- get_si_mat(aligns, maxCombine = 1, ref = aligns[[1]])

    expected <- cbind(c(10, 20), c(5, 0)) * pi
    rownames(expected) <- cs[c(2, 4)]
    colnames(expected) <- c("spec_1", "spec_2")

    expect_equal(mat, expected)
})
