# 0. Script inits
# ===============

devtools::load_all()
get_started_dir <- pkg_file("vignettes/Get_Started")
sim_01_param_checking.svg = file.path(get_started_dir, "sim_01_param_checking.svg")
sim_01_spectrum.svg = file.path(get_started_dir, "sim_01_spectrum.svg")
sim_01_peak_triplets.svg = file.path(get_started_dir, "sim_01_peak_triplets.svg")
sim_01_lorentz_curves.svg = file.path(get_started_dir, "sim_01_lorentz_curves.svg")
sim_01_superposition.svg = file.path(get_started_dir, "sim_01_superposition.svg")

# 1. Deconvolute spectra
# ======================

sim_dir <- metabodecon_file("bruker/sim")
ewobj <- evalwith(
    pars = list(mfrow = c(1, 2)),
    answers = c("y", "1", "y", "y"),
    generate_lorentz_curves(data_path = sim_dir, sfr = c(3.42, 3.58), wshw = 0, smopts = c(1, 5), delta = 2, verbose = FALSE),
    plot = svg(sim_01_param_checking.svg, width = 10, height = 5),
    output = "captured",
    message = "captured"
)
deconvs <- ewobj$rv

# 2. Visualize deconvoluted spectra
# =================================

# After devonvoluting the spectra, we might want to visualize the deconvoluted
# signals, to get an impression of the quality of the deconvolution. Metabodecon
# provides four functions for this purpose:

# 1. `plot_spectrum()` for plotting the the spectrum before the deconvolution
# 2. `plot_triplets()` for plotting the peak centers and borders detected by the
#    deconvolution algorithm
# 3. `plot_lorentz_curves_save_as_png()` for plotting the individual
#    deconvoluted signals
# 4. `plot_spectrum_superposition_save_as_png()` for plotting the superposition
#    of all deconvoluted signals

# The following code snippet shows how these functions can be used for the
# deconvoluted Sim dataset. The resulting plots for the first spectrum of the
# Sim dataset, are shown in [Figure 2](#fig-sim1-plots).

decon <- deconvs[[1]] # Generated in section 'Deconvolute spectra'
svg(sim_01_spectrum.svg, width = 5, height = 5)
plot_spectrum(decon, foc_rgn = c(0.4, 0.6))
dev.off()

plot_triplets(deconvs[1])

plot_lorentz_curves_save_as_png()
plot_spectrum_superposition_save_as_png()
