# Script inits #####

devtools::load_all()
img_dir <- pkg_file("vignettes/Get_Started")
sim_01_param_checking.svg <- file.path(img_dir, "sim_01_param_checking.svg")
sim_01_spectrum.svg <- file.path(img_dir, "sim_01_spectrum.svg")
sim_02_spectrum.svg <- file.path(img_dir, "sim_02_spectrum.svg")


# Deconvolute spectra #####

sim_dir <- metabodecon_file("bruker/sim_subset")
ewobj <- evalwith(
    # output = "captured", message = "captured",
    plot = svg(sim_01_param_checking.svg, width = 10, height = 5),
    pars = list(mfrow = c(1, 2)),
    answers = c("y", "1", "y", "y"),
    expr = generate_lorentz_curves(
        data_path = sim_dir,
        sfr = c(3.42, 3.58), wshw = 0,
        smopts = c(2, 5), delta = 0.3,
        verbose = TRUE, nworkers = 1
    )
)
deconvs <- ewobj$rv

# Visualize deconvoluted spectra #####

svg(sim_01_spectrum.svg, width = 7, height = 7)
plot_spectrum(deconvs[[1]], foc_rgn = c(0.5, 0.2))
dev.off()

svg(sim_02_spectrum.svg, width = 7, height = 7)
plot_spectrum(deconvs[[2]], foc_rgn = c(0.5, 0.2))
dev.off()
