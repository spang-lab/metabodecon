---
title: "PRARPX vs MSE"
output: rmarkdown::html_vignette
css: styles.css
vignette: >
  %\VignetteIndexEntry{PRARPX vs MSE}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r set-defaults, echo=FALSE, results=FALSE, message=FALSE}
knitr::opts_chunk$set(
  fig.dim = c(5, 5),
  fig.show = "hold",
  out.height = "auto",
  eval = FALSE,
  echo = FALSE,
  results = FALSE,
  message = FALSE
)
```

This vignette illustrates the relationship between PRARPX and MSE for a
simulated spectrum. PRARPX incorporates both peak recovery and residual error,
so it better reflects deconvolution quality than MSE alone. In practice, MSE
may favor overly complex solutions, whereas PRARPX penalizes overfitting by
accounting for peak correctness. For the Sim spectrum, the global optimum under
MSE coincides with the PRARPX optimum, which is reassuring when using MSE as a
proxy during parameter search.

# Parameter grid

The following code evaluates PRARPX and MSE over the grid used by
`find_best_params()`. We use `sim[[2]]`, which is representative of the
simulated spectra and yields a matching optimum for PRARPX and MSE on this
grid.

```{r prarpx-grid, eval=FALSE}
smopts_grid <- list(
  c(1, 3), c(1, 5), c(2, 3), c(2, 5),
  c(3, 3), c(3, 5), c(3, 7)
)
delta_grid <- unique(pmax(0.1, 6.4 + c(-1, -0.5, 0, 0.5, 1)))
nfit_grid <- c(2, 3, 4)

spec <- metabodecon::sim[[2]]
df <- metabodecon:::prarpx_grid(
    spectrum = spec,
    smopts_grid = smopts_grid,
    delta_grid = delta_grid,
    nfit_grid = nfit_grid,
    verbose = TRUE
)

best <- metabodecon:::prarpx_best_rows(df)
print(best$best_mse)
print(best$best_prarpx)

metabodecon:::prarpx_plot(df, "figs/prarpx_mse_sim.png")
```

```{r prarpx-plot, echo=FALSE}
knitr::include_graphics("figs/prarpx_mse_sim.png")
```

In this grid, the best PRARPX and best MSE occur at the same parameter
combination, supporting the use of MSE for the coarse search while retaining
PRARPX for final validation.
