---
title: "NMR Basics"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{NMR Basics}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r knitr-setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)
```

# NMR Basics

This article attempts to give a brief introduction to 1D NMR spectroscopy, highlight the difficulties involved and explain how metabodecon can be used to overcome these difficulties.

## Whats the general idea of NMR

1.  Metabolites are molecules.
2.  Molecules are constructed from atoms. For instance, water (H2O) is composed of two hydrogen atoms and one oxygen atom.
3.  Atoms possess physical properties such as mass, charge, and nuclear spin. A hydrogen atom (H1), for instance, has a nuclear spin of 1/2, while an oxygen atom (O16) has a nuclear spin of 0.
4.  Nuclei with a spin not equal to 0 behave like tiny magnets when placed in a magnetic field of strength B, i.e., aligning in opposite direction is energetically worse then aligning in the same direction as B. Therefore, the proportion of nuclei aligned in the same direction as the magnetic field is slightly higher than the proportion of nuclei aligned in the opposite direction. The amount of nuclei occupying the higher energy state depends on the strength of the magnetic field and the temperature of the system. The higher the magnetic field strength and the lower the temperature, the more nuclei will align in the low energy state.
5.  If we direct photons, e.g. in the form of radiowaves, at these nuclei, a nucleus may transition into a "high-energy" state by absorbing a photon and flipping its orientation. This transition requires the photon to have the exact amount of energy needed for the switch.
6.  This change in orientation induces a shift in the magnetic field, which in turn induces a small current in the probe that can be measured.
7.  By directing photons of varying energies at a mixture of nuclei, we can identify which photons are absorbed by measuring this current. This results in a frequency/absorption spectrum.
8.  If we place pure water in an NMR device operating at a field strength B = 2MT, we can observe the unique absorption pattern of water, which consists of two peaks in the NMR spectrum at frequencies of 2.3 and 4.7 MHz. The reason for this exact pattern is, that only the hydrogen nuclei will absorb photons, as the oxygen nuclei have a spin of 0 and therefore are unaffected by the magnetic field.
9.  If we place pure ethanol in an NMR device operating at a field strength B = 2MT, we will observe the unique absorption pattern of ethanol instead, which consists of three peaks at frequencies of 1.2, 3.6, and 6.0 MHz. This is because ethanol is composed of two different types of hydrogen nuclei, one attached to the carbon atom and one attached to the oxygen atom.
10. If we place a mixture of water and ethanol in an NMR device operating at a field strength B = 2MT, we will observe a combination of the above absorption patterns. The amplitude of the peaks reflects the concentration of the respective metabolites in the mixture, allowing us to determine the concentration of metabolites from a mixture.

## What are the difficulties

### 1. Overlapping peaks

1.  As soon as we add more and more different molecules in a our mixture, the resulting spectrum gets more and more crowded. At some point, certain peaks will overlap, making it hard to distinguish the individual absorption pattern of the metabolites from another.

2.  Consider the example from section [Whats the general idea of NMR](#whats-the-general-idea-of-nmr), where we have a mixture of water and ethanol.

    PICTURE

3.  Now let's add three more metabolites, e.g. acetate, alanine and lactate, to the mixture.

    PICTURE

4.  As you can see, the peaks of ethanol and lactate overlap, making it hard to identify the individual absorption patterns of these two metabolites.

    PICTURE

5.  The untangling of overlapping peaks into individual peaks is called deconvolution and is one of the problems metabodecon tries to solve. For the exact algorithm, that metabodecon uses to deconvolute spectra, see section [Deconvolution](how-#how-does-metabodecon-work) [How does metabodecon work](#how-does-metabodecon-work).


