library(testthat)

# Helper: build a minimal aligns-class object from cs + per-spectrum
# (pcial, A, lambda) triples. Bypasses deconvolution so unit tests can
# drive snap_to_ref() on known peak lists.
make_aligns <- function(cs, peaks_list, areas_list, lambdas_list=NULL) {
    if (is.null(lambdas_list)) {
        lambdas_list <- lapply(peaks_list, function(p) rep(0.5, length(p)))
    }
    objs <- lapply(seq_along(peaks_list), function(i) {
        pcial <- as.integer(peaks_list[[i]])
        x0 <- cs[pcial]
        lcpar <- data.frame(
            x0=x0, A=areas_list[[i]], lambda=lambdas_list[[i]],
            pcide=pcial, x0al=x0, pcial=pcial
        )
        structure(list(
            cs=cs,
            lcpar=lcpar,
            sit=list(),
            meta=list(name=sprintf("spec_%d", i))
        ), class=c("align", "decon2", "spectrum"))
    })
    names(objs) <- vapply(objs, function(o) o$meta$name, character(1))
    structure(objs, class=c("aligns", "decons2", "spectra"))
}

test_that("snap_to_ref adds pcisn/x0sn and preserves all other fields", {
    cs <- 9:1
    a <- make_aligns(
        cs=cs,
        peaks_list=list(c(2, 5), c(3, 6)),
        areas_list=list(c(2, 3), c(4, 5))
    )
    out <- snap_to_ref(a, ref=a[[1]], maxCombine=1)
    # Spectrum 1 (= ref): both peaks snap to themselves.
    expect_equal(out[[1]]$lcpar$pcisn, c(2L, 5L))
    expect_equal(out[[1]]$lcpar$x0sn,  cs[c(2L, 5L)])
    # Spectrum 2: peak at 3 -> nearest ref 2, peak at 6 -> nearest 5.
    expect_equal(out[[2]]$lcpar$pcisn, c(2L, 5L))
    expect_equal(out[[2]]$lcpar$x0sn,  cs[c(2L, 5L)])
    # x0, A, lambda, pcide, pcial preserved on every row.
    for (i in seq_along(out)) {
        expect_equal(out[[i]]$lcpar$x0,     a[[i]]$lcpar$x0)
        expect_equal(out[[i]]$lcpar$A,      a[[i]]$lcpar$A)
        expect_equal(out[[i]]$lcpar$lambda, a[[i]]$lcpar$lambda)
        expect_equal(out[[i]]$lcpar$pcide,  a[[i]]$lcpar$pcide)
        expect_equal(out[[i]]$lcpar$pcial,  a[[i]]$lcpar$pcial)
    }
})

test_that("snap_to_ref marks out-of-range peaks with NA pcisn/x0sn", {
    cs <- 9:1
    a <- make_aligns(
        cs=cs,
        peaks_list=list(c(2), c(5)),
        areas_list=list(c(10), c(20))
    )
    out <- snap_to_ref(a, ref=a[[1]], maxCombine=1)
    # Spectrum 2's peak at 5 is distance 3 from the only ref peak (2):
    # row is kept but pcisn / x0sn are NA.
    expect_equal(nrow(out[[2]]$lcpar), 1L)
    expect_equal(out[[2]]$lcpar$pcisn, NA_integer_)
    expect_equal(out[[2]]$lcpar$x0sn,  NA_real_)
    # x0, A, lambda preserved.
    expect_equal(out[[2]]$lcpar$x0,     a[[2]]$lcpar$x0)
    expect_equal(out[[2]]$lcpar$A,      a[[2]]$lcpar$A)
    expect_equal(out[[2]]$lcpar$lambda, a[[2]]$lcpar$lambda)
})

test_that("snap_to_ref does not sum amplitudes on collision", {
    cs <- 9:1
    # Two peaks (cols 2 and 4) both nearest to ref col 3, both within
    # maxCombine=2. Both rows must survive with their original A; only
    # si_mat() should sum them.
    a <- make_aligns(
        cs=cs,
        peaks_list=list(c(3), c(2, 4)),
        areas_list=list(c(1), c(2, 3))
    )
    out <- snap_to_ref(a, ref=a[[1]], maxCombine=2)
    expect_equal(out[[2]]$lcpar$pcisn, c(3L, 3L))
    expect_equal(out[[2]]$lcpar$A,     c(2, 3))  # NOT summed
    # si_mat does the summing at rasterisation time.
    m <- si_mat(out)
    # Column for cs == 7 (i.e. position 3) of spectrum 2 should be (2+3)*pi.
    expect_equal(m["spec_2", "7"], (2 + 3) * pi)
})

test_that("snap_to_ref with maxCombine=0 is a no-op", {
    cs <- 9:1
    a <- make_aligns(
        cs=cs,
        peaks_list=list(c(2, 5), c(3, 6)),
        areas_list=list(c(2, 3), c(4, 5))
    )
    out <- snap_to_ref(a, ref=a[[1]], maxCombine=0)
    expect_equal(out, a)
})

test_that("snap_to_ref picks reference automatically when ref=NULL", {
    skip_if_speaq_deps_missing()
    decons <- deconvolute(sim[1:3], sfr=c(3.55, 3.35), verbose=FALSE)
    aligned <- align(decons, maxShift=20, maxCombine=0, verbose=FALSE)
    snapped <- snap_to_ref(aligned, maxCombine=10)
    # Every non-NA pcisn must sit on the chosen reference's pcial grid.
    ref <- find_ref(aligned)
    pp <- ref$lcpar$pcial
    expect_s3_class(snapped, "aligns")
    for (s in snapped) {
        pcisn <- s$lcpar$pcisn
        ok <- !is.na(pcisn)
        expect_true(all(pcisn[ok] %in% pp))
    }
})
