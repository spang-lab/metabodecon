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

test_that("get_si_mat with maxSnap returns reduced matrix", {

    withr::local_output_sink(nullfile())

    deps <- c("MassSpecWavelet", "impute")
    inst <- sapply(deps, requireNamespace, quietly = TRUE)
    if (!all(inst)) skip(paste("Missing deps:", collapse(deps[!inst])))

    decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
    al <- align(decons, maxCombine = 0, verbose = FALSE)

    mat_raw <- get_si_mat(al)
    mat_1hw <- get_si_mat(al, maxSnap = 1)
    mat_2hw <- get_si_mat(al, maxSnap = 2)

    # Snapped matrices should have fewer rows (one per ref peak)
    expect_lt(nrow(mat_1hw), nrow(mat_raw))
    expect_lt(nrow(mat_2hw), nrow(mat_raw))
    expect_equal(ncol(mat_1hw), ncol(mat_raw))
    expect_equal(ncol(mat_2hw), ncol(mat_raw))
})
