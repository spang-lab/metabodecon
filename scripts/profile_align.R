library(metabodecon)

out <- "tmp/profile.out"
dir.create("tmp", showWarnings = FALSE)

s <- metabodecon:::read_aki_data()$spectra[[1]]
dd <- deconvolute(s[1:4], nfit=10, use_rust=0, verbose=FALSE)

Rprof(out, interval = 0.005, line.profiling = FALSE)
for (i in 1:10) {
    cat(sprintf("%s %d/10\n", Sys.time(), i))
    aa <- metabodecon::align(dd)
}
Rprof(NULL)

sm <- summaryRprof(out)
cat("\n=== Self time (top 20) ===\n")
print(head(sm$by.self,  20))
cat("\n=== Total time (top 20) ===\n")
print(head(sm$by.total, 20))

write.csv(sm$by.self,  "tmp/align_profile_self.csv")
write.csv(sm$by.total, "tmp/align_profile_total.csv")
