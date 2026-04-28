library(metabodecon)

s <- metabodecon:::read_aki_data()$spectra[[1]]
profvis::profvis({
    for (i in 1:10) {
        cat(sprintf("%s %d/10\n", Sys.time(), i))
        metabodecon::deconvolute(s, nfit = 10, use_rust = 0.5, verbose = FALSE)
    }
})
