library(metabodecon)

out <- "tmp/profile.out"
dir.create("tmp", showWarnings = FALSE)

s <- metabodecon:::read_aki_data()$spectra[[1]]
Rprof(out, interval = 0.01, line.profiling = FALSE)
for (i in 1:20) {
    a <- Sys.time()
    cat(sprintf("%s %d/10", format(a, "%Y-%m-%d %H:%M:%S"), i))
    metabodecon::deconvolute(s, nfit=10, use_rust=0, verbose=FALSE)
    d <- Sys.time() - a
    cat(sprintf(" (%.2f seconds)\n", as.numeric(d, units="secs")))
}
Rprof(NULL)

sm <- summaryRprof(out)
cat("\n=== Self time (top 20) ===\n")
print(head(sm$by.self,  20))
cat("\n=== Total time (top 20) ===\n")
print(head(sm$by.total, 20))

write.csv(sm$by.self,  "tmp/profile_self.csv")
write.csv(sm$by.total, "tmp/profile_total.csv")
