---
title: "Deconvolution Details"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{metabodecon}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r knitr-setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)
```

# The Problem

Consider a spectrum of a mixture of water, ethanol, acetate, alanine and lactate. The spectrum could look like this:

PICTURE

As you can see, the peaks of ethanol and lactate overlap, making it hard to identify the individual absorption patterns of these two metabolites. To solve this problem, one can try to model the spectrum as a sum of peaks, where each peak has the shape of a Lorentzian function. For an explanation of why lorentz curves are the best choice for this task, see TODO from TODO et al.

E.g. the spectrum above could be written as:

$$
TODO
$$

where $A_i$ is the amplitude of the $i$-th peak, $f_i$ is the frequency of the $i$-th peak, $w_i$ is the width of the $i$-th peak and $c$ is a constant offset.

That means, the goal is to find the list of peak positions $w_0$ and their corresponding amplitudes $A_0$ that best fit the spectrum. Often times, this problem is divided into two subproblems:

1. Finding the peak positions $w_0$
2. Finding the corresponding amplitudes $A$ and $\lambda$ values that optimize the MSE

Subproblem 1 could be solved naively by searching for local maxima in the spectrum and then directly using these as inputs for subproblem 2. This is what's done in the LEVENBERG-MARQUARDT algorithm, which is often found in commercial NMR software. However, this approach has them disadvantages, which are explained in section [Levenberg-Marqardt-Peak-selection](#Levenberg-Marqardt-Peak-selection). The alternative approach, which is taken by metabodecon is explained in section [Koh-Peak-Selection] and is based on searching for local minima of the second derivative, i.e. the curvature, of the spectrum. Subproblem 2 is then solved using the iterative procedure described in section [Koh-Parameter-Approximation] and was first published together with the Koh-Peak-Selection algorithm in 2009 by Koh et al.

# Levenberg-Marqardt-Peak-selection

TODO

# Koh-et-al-Peak-Selection

# Koh-et-al-Parameter-Approximation
