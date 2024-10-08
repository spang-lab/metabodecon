# Script inits #####

img_dir <- pkg_file("vignettes/Get_Started")
img_path <- function(file) file.path(img_dir, file)

# Deconvolute spectra #####

sim_dir <- metabodecon_file("bruker/sim")
ewobj <- evalwith(
    # output = "captured", message = "captured",
    plot = svg(img_path("sim_01_param_checking.svg"), width = 10, height = 5),
    pars = list(mfrow = c(1, 2)),
    answers = c("y", "1", "y", "y"),
    expr = generate_lorentz_curves(
        data_path = sim_dir,
        sfr = c(3.35, 3.55), wshw = 0,
        smopts = c(2, 5), delta = 0.1,
        verbose = FALSE, nworkers = 1
    )
)
deconvs <- ewobj$rv

# Visualize deconvoluted spectra #####

svg(img_path("sim_01_spectrum.svg"), width = 7, height = 7)
plot_spectrum(deconvs[[1]], foc_rgn = c(0.5, 0.2))
dev.off()

svg(img_path("sim_02_spectrum.svg"), width = 7, height = 7)
plot_spectrum(deconvs[[2]], foc_rgn = c(0.5, 0.2))
dev.off()

# Align deconvoluted spectra

## Look up spectrum range

svg(img_path("sim_spectra_overlayed.svg"), width = 10, height = 5)
ppm_range <- get_ppm_range(spectrum_data = deconvs, show = TRUE)
dev.off()
