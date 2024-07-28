# Deconvolute spectra

devtools::load_all()

sim_dir <- metabodecon_file("bruker/sim")
svg_path = file.path(pkg_file("vignettes/Get_Started"), "sim_01_decon.svg")
ewobj <- evalwith(
    pars = list(mfrow = c(1, 2)),
    answers = c("y", "1", "y", "y"),
    generate_lorentz_curves(data_path = sim_dir, sfr = c(3.42, 3.58), wshw = 0, smopts = c(1, 5), delta = 2, verbose = FALSE),
    plot = svg(svg_path, width = 10, height = 5),
    output = "captured",
    message = "captured"
)
deconvs <- ewobj$rv

# Visualize deconvoluted spectra

devtools::load_all()

decon <- deconvs[[1]] # Generated in section 'Deconvolute spectra'
plot_spectrum(decon, foc_rgn = c(0.4, 0.6))
