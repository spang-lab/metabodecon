---
title: "Urine Analysis"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Urine Analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<style>
.sourceCode.txt {
  background-color: #f1f3f5;
}
figure {
  border: 1px solid #dee2e6;
  border-radius: .375rem;
  padding: 0;
  background-color: #f1f3f5;
}
figcaption {
  text-align: center;
  font-style: italic;
  border-top: 1px solid #dee2e6;
  padding: 0.75rem;
  background-color: #ffffff;
}
</style>

```R
urine_dir <- metabodecon_file("bruker/urine")
urine_01_dir <- file.path(urine_dir, "urine_1")
urine_01_deconv <- generate_lorentz_curves(
    urine_01_dir, ask = FALSE
)
plot_spectrum(urine_01_deconv, foc_rgn = c(-1.9, -1.8), foc_unit = "ppm", trp_pch = c(17, 4, 4, 124))
plot_spectrum(urine_01_deconv, foc_rgn = c(0.25, 0.50))

sim_dir <- metabodecon::metabodecon_file("bruker/sim")
sim_01_dir <- file.path(sim_dir, "sim_01")
sim_01_deconv <- metabodecon::generate_lorentz_curves(
    sim_01_dir, sfr = c(3.42, 3.58), wshw = 0, smopts = c(1, 1), delta = 0.1, ask = FALSE
)
metabodecon::plot_spectrum(sim_01_deconv, foc_rgn = c(0, 1), sub_show = FALSE)
metabodecon::plot_spectrum(sim_01_deconv, foc_rgn = c(0.25, 0.50))
metabodecon::plot_spectrum(sim_01_deconv, foc_rgn = c(0.00, 0.25))
```
