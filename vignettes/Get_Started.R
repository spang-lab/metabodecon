# Deconvolute spectra

devtools::load_all()

sim_dir <- metabodecon_file("bruker/sim")
sim_01_decon.svg = file.path(pkg_file("vignettes/Get_Started"), "sim_01_decon.svg")
ewobj <- evalwith(
    pars = list(mfrow = c(1, 2)),
    answers = c("y", "1", "y", "y"),
    generate_lorentz_curves(data_path = sim_dir, sfr = c(3.42, 3.58), wshw = 0, smopts = c(1, 5), delta = 2, verbose = FALSE),
    plot = svg(sim_01_decon.svg, width = 10, height = 5),
    output = "captured",
    message = "captured"
)
deconvs <- ewobj$rv

# Visualize deconvoluted spectra

devtools::load_all()

decon <- deconvs[[1]] # Generated in section 'Deconvolute spectra'
sim_01_spectrum.svg = file.path(pkg_file("vignettes/Get_Started"), "sim_01_spectrum.svg")
svg(sim_01_spectrum.svg, width = 5, height = 5)
plot_spectrum(decon, foc_rgn = c(0.4, 0.6))
dev.off()
