test_that("si_mat returns a matrix of the correct dimensions", {

    withr::local_output_sink(nullfile())

    # 'speaq' requires 'MassSpecWavelet' and 'impute' to be installed
    deps <- c("MassSpecWavelet", "impute")
    inst <- sapply(deps, requireNamespace, quietly = TRUE)
    if (!all(inst)) skip(paste("Missing deps:", collapse(deps[!inst])))

    decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
    aligns <- align(decons)
    mat <- si_mat(aligns)
    # si_mat() returns spectra in rows, chemical shifts in columns.
    expect_equal(dim(mat), c(2, 2048))
    expect_equal(as.numeric(colnames(mat)), aligns[[1]]$cs)
    expect_equal(rownames(mat), get_names(aligns))
})

test_that("si_mat drop_zero removes all-zero columns", {

    withr::local_output_sink(nullfile())

    deps <- c("MassSpecWavelet", "impute")
    inst <- sapply(deps, requireNamespace, quietly = TRUE)
    if (!all(inst)) skip(paste("Missing deps:", collapse(deps[!inst])))

    decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
    aligns <- align(decons)
    full <- si_mat(aligns)
    compact <- si_mat(aligns, drop_zero = TRUE)
    expect_lt(ncol(compact), ncol(full))
    expect_true(all(colSums(compact != 0) > 0))
    expect_equal(compact, full[, colSums(full != 0) > 0, drop = FALSE])
})

make_aligns_for_si_mat_test <- function(cs, peaks_list, areas_list) {
    objs <- lapply(seq_along(peaks_list), function(i) {
        pcial <- as.integer(peaks_list[[i]])
        x0 <- cs[pcial]
        structure(list(
            cs = cs,
            lcpar = data.frame(x0 = x0, x0al = x0, pcial = pcial,
                               A = areas_list[[i]]),
            sit = list(),
            meta = list(name = sprintf("spec_%d", i))
        ), class = c("align", "decon2", "spectrum"))
    })
    names(objs) <- vapply(objs, function(o) o$meta$name, character(1))
    structure(objs, class = c("aligns", "decons2", "spectra"))
}

test_that("si_mat places each peak at its column and sums collisions", {
    cs <- 9:1
    aligns <- make_aligns_for_si_mat_test(
        cs = cs,
        peaks_list = list(c(4, 8), c(3, 8), c(3, 4, 7), c(5, 9), c(5, 9)),
        areas_list = list(c(2, 4), c(3, 4), c(2, 4, 5), c(3, 3), c(2, 3))
    )

    mat <- si_mat(aligns)

    expected <- rbind(
        c(0, 0, 0, 2, 0, 0, 0, 4, 0),
        c(0, 0, 3, 0, 0, 0, 0, 4, 0),
        c(0, 0, 2, 4, 0, 0, 5, 0, 0),
        c(0, 0, 0, 0, 3, 0, 0, 0, 3),
        c(0, 0, 0, 0, 2, 0, 0, 0, 3)
    ) * pi
    rownames(expected) <- sprintf("spec_%d", 1:5)
    colnames(expected) <- as.character(cs)

    expect_equal(mat, expected)
})

test_that("si_mat preserves the full cs grid", {
    cs <- 9:1
    aligns <- make_aligns_for_si_mat_test(
        cs = cs,
        peaks_list = list(c(4, 8)),
        areas_list = list(c(2, 4))
    )

    mat <- si_mat(aligns)

    expect_equal(ncol(mat), length(cs))
    expect_equal(colnames(mat), as.character(cs))
})
