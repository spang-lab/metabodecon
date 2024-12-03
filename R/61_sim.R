# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Sim #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

docs <- "see below"
#' @title The Sim Dataset
#'
#' @description
#' A simulated dataset generated from the [Blood](
#' https://spang-lab.github.io/metabodecon/articles/Datasets.html#Blood)
#' dataset.
#'
#' @format A [spectra](metabodecon_classes) object consisting of 16 [spectrum]
#' objects, where each spectrum object contains 2048 datapoints ranging from
#' 3.60 to 3.29 ppm
"sim"

update_sim <- function(overwrite = FALSE, verbose = FALSE) {
    sim <- make_sim()
    if (isTRUE(overwrite)) {
        path <- file.path(pkg_file("example_datasets/bruker"), "sim")
        if (verbose) logf("Overwrite is TRUE. Updating %s." , path)
        save_spectra(sim, path, force = TRUE, verbose = verbose)
        usethis::use_data(sim, overwrite = TRUE)
    } else {
        path <- file.path(tmpdir(subdir = TRUE, create = TRUE), "sim")
        if (verbose) logf("Overwrite is FALSE. Writing spectra to %s." , path)
        save_spectra(sim, path, force = FALSE, verbose = verbose)
        save(sim, file = file.path(path, "sim.rda"), compress = "bzip2")
    }
    path
}

# "C:/Users/tobi/AppData/Local/Temp/RtmpE3Ma4o/metabodecon/4b0857617065/sim/sim.rda"

make_sim <- function() {
    decons <- deconvolute_blood()
    simpars <- lapply(decons, get_sim_params, pkr = c(3.52, 3.37))
    simpars[[2]] <- get_sim_params(decons[[2]], pkr = c(3.51, 3.36)) # (1)
    # (1) Blood_02 is shifted approx. 0.01 ppm to the right, so we also need to
    # shift the interval from which we pick our peaks by 0.01 to end up with
    # signals from the same metabolites. This was determined by visual
    # inspection after deconvoluting the blood dataset.
    sim <- lapply(simpars, function(simpar) {
        simulate_spectrum(
            name = gsub("blood", "sim", simpar$name),
            ndp = 2048,
            x0 = simpar$x0,
            A = simpar$A,
            lambda = simpar$lambda,
            noise = rnorm(2048, sd = simpar$noiseSD)
        )
    })
    names(sim) <- gsub("blood", "sim", names(sim))
    class(sim) <- "spectra"
    sim
}

#' @noRd
#'
#' @title Deconvolute the Blood Dataset
#'
#' @description
#' Downloads and deconvolutes the Blood Dataset.
#'
#' @param read_cache If TRUE the function will try to read the results from a
#' cache file instead of downloading and deconvoluting the dataset again.
#'
#' @param write_cache If TRUE, the results will be written to the cache file
#' path. If the cache file already exists and the results from the deconvolution
#' and the stored object differ, a warning will be issued and the new results
#' will only be stored if force is TRUE.
#'
#' @param force If TRUE, the results will be written to the cache file path even
#' if the results from the deconvolution and the stored object differ.
#'
#' @return
#' A `decons2` object containing the deconvolution results. For details
#' see [deconvolute()].
#'
#' @examples
#' deconvolute_blood() # Use cache file if it exists
#' deconvolute_blood(read_cache = FALSE) # Deconvolute from scratch
#' deconvolute_blood(force = TRUE) # Force update of cache file
deconvolute_blood <- function(read_cache = TRUE,
                              write_cache = TRUE,
                              force = FALSE,
                              verbose = TRUE,
                              nworkers = 1) {
    rds <- file.path(cachedir(), "deconvolute_blood.rds")
    new <- old <- if (file.exists(rds)) readRDS(rds) else NULL # 62.6 MB ~= 0.6s
    if (!read_cache || is.null(new)) {
        download_example_datasets()
        path <- datadir("example_datasets/bruker/blood")
        new <- generate_lorentz_curves(
            data_path = path,
            nworkers = nworkers,
            verbose = verbose,
            ask = FALSE
        )
    }
    if (!(identical(new, old) || is.null(old))) {
        warning("Cache and deconvolution differ.", immediate. = TRUE)
    }
    if (write_cache && (is.null(old) || force)) {
        logf("Writing results to cache file %s", rds)
        saveRDS(new, rds)
    }
    return(new)
}

#' @noRd
#' @examples
#' xdir <- download_example_datasets()
#' path <- file.path(xdir, "bruker/blood/blood_01")
#' spec <- read_spectrum(path)
#' x <- deconvolute_ispec(spec)
#' sim_params <- get_sim_params(x)
get_sim_params <- function(x, pkr = c(3.4, 3.5)) {
    d <- as_decon1(x)
    p <- which(d$x_0_ppm <= max(pkr) & d$x_0_ppm >= min(pkr))
    A <- -d$A_ppm[p] * 1e6
    x0 <- d$x_0_ppm[p]
    lambda <- -d$lambda_ppm[p]
    noiseSD <- sd(d$y_values_raw[c(1:10000, 121073:131072)])
    name <- d$filename
    named(A, x0, lambda, noiseSD, name)
}

get_sim1_decon1 <- function() {
    generate_lorentz_curves(
        data_path = sim[[1]],
        sfr = c(3.55, 3.35),
        wshw = 0,
        ask = FALSE,
        verbose = FALSE
    )
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Deprecated #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

simulate_from_decon <- function(x,
                                # Simulation Params
                                cs_min = 3.4,
                                cs_max = 3.6,
                                pk_min = 3.45,
                                pk_max = 3.55,
                                noise_method = "SFR",
                                # Output params
                                show = TRUE,
                                save = FALSE,
                                verbose = TRUE,
                                ...) {
    if (!isFALSE(verbose)) logf("Checking function inputs")
    stopifnot(is_num(cs_min, 1))
    stopifnot(is_num(cs_max, 1))
    stopifnot(is_num(pk_min, 1))
    stopifnot(is_num(pk_max, 1))
    stopifnot(is_char(noise_method, 1, "(RND|SFR)"))
    stopifnot(is_bool(show, 1))
    stopifnot(is_bool(save, 1))
    stopifnot(is_bool(verbose, 1))
    d <- as_decon1(x)
    logv <- if (verbose) logf else function(...) NULL

    logv("Simulating spectrum from '%s' (noise method = '%s')", d$filename, noise_method)
    s <- as_spectrum(d)

    logv("Throwing away datapoints outside of %.1f to %.1f", cs_min, cs_max)
    ix <- which(s$cs >= cs_min & s$cs <= cs_max)
    s$cs <- s$cs[ix]
    s$si <- s$si[ix]
    if (!is.null(s$meta$fq)) s$meta$fq <- s$meta$fq[ix]

    logv("Keeping only within %.1f to %.1f", pk_min, pk_max)
    ip <- which(d$x_0_ppm <= pk_max & d$x_0_ppm >= pk_min)
    s$simpar <- list(
        A      = -d$A_ppm[ip],
        x_0    = d$x_0_ppm[ip],
        lambda = -d$lambda_ppm[ip]
    )

    logv("Calculating simulated signal intensities (si_sim) as superposition of lorentz curves")
    rownames(s) <- NULL
    si <- lorentz_sup(x = s$cs, lcpar = s$simpar)

    logv("Adding %s noise to simulated data", noise_method)
    if (noise_method == "RND") {
        # ToSc, 2024-08-02: noise_method 'RND' turned out to be a bad idea,
        # because the true noise is not normally distributed, but has long
        # stretches of continuous increase or decrease that would be incredibly
        # unlikely to occur with a normal distribution. This can be seen by
        # running `analyze_noise_methods()`.
        sfr_sd <- sd(d$y_values_raw[c(1:10000, 121073:131072)] * 1e-6) # (1)
        # (1) The true SFR covers approx. the first and last 20k datapoints.
        # However, to be on the safe side, only use the first and last 10k
        # datapoints for the calculation.
        noise <- rnorm(n = length(si), mean = 0, sd = sfr_sd)
        si <- si + noise
    } else {
        idx <- ix - min(ix) + 5000 # Use SI of datapoints 5000:6308 for noise
        noise <- d$y_values_raw[idx] * 1e-6
        si <- si + noise
    }

    logv("Discretize signal intensities to allow efficient storage as integers")
    s$si <- as.integer(si * 1e6) / 1e6 # (2)
    # (2) Convert to integers and back. We do this to not lose precision later
    # on when the data gets written to disc as integers.

    if (show) plot_sim_spec(s)
    if (store) save_spectrum(s, verbose = verbose, ...)

    logv("Finished spectrum simulation")
    invisible(s)
}

#' @noRd
#' @title Count Stretches of Increases and Decreases
#' @description Counts the lengths of consecutive increases and decreases in a
#' numeric vector.
#' @param x A numeric vector.
#' @return A numeric vector containing the lengths of stretches of increases and
#' decreases.
#' @examples
#' #
#' # Example Data (x)
#' # | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
#' # |-------------------------------|
#' # |   |   |   |###|   |   |   |###|
#' # |   |   |###|###|   |   |###|###|
#' # |   |###|###|###|###|   |###|###|
#' # |###|###|###|###|###|###|###|###|
#' # |-------------------------------|
#' # |   +   +   +   -   -   +   +   | + is increase
#' # |-------------------------------| - is decrease
#' #
#' x <- c(1.0, 2.2, 3.0, 4.4, 2.0, 1.0, 3.0, 4.0)
#' count_stretches(x) # Returns c(3, 2, 2) because we have + + + - - + +
count_stretches <- function(x) {
    if (length(x) < 2) {
        return(integer(0))
    }
    ss <- numeric(length(x))
    s <- 1
    inc <- x[2] > x[1]
    for (i in 3:length(x)) {
        if ((inc && x[i] > x[i - 1]) || (!inc && x[i] < x[i - 1])) {
            s <- s + 1
        } else {
            ss[i - 1] <- s
            s <- 1
            inc <- x[i] > x[i - 1]
        }
    }
    ss[i] <- s
    return(ss[ss != 0])
}

#' @noRd
#' @description Used during development of `simulate_spectra()` to find a
#' realistic method for noise generation.
analyze_noise_methods <- function(ask = TRUE) {
    download_example_datasets()
    blood_1 <- datadir("example_datasets/bruker/blood/blood_01")
    deconv <- generate_lorentz_curves(blood_1, ask = FALSE)
    si_raw <- deconv$y_values_raw
    sd_sfr <- sd(si_raw[c(1:10000, 121073:131072)] * 1e-6)
    siRND <- rnorm(n = 10000, mean = 0, sd = sd_sfr)
    siSFR <- si_raw[1:10000] * 1e-6

    logf("Visualizing raw SIs for noise methods RND and SFR")
    plot_noise_methods(siRND, siSFR)
    if (!get_yn_input("Continue?")) {
        return()
    }

    logf("Visualizing smoothed SIs for noise methods RND and SFR")
    siRND_sm <- smooth_signals(list(y_pos = siRND))$y_smooth
    siSFR_sm <- smooth_signals(list(y_pos = siSFR))$y_smooth
    plot_noise_methods(siRND_sm, siSFR_sm)
    if (!get_yn_input("Continue?")) {
        return()
    }

    logf("Visualizing lengths of intervals of continuous increase and/or decrease")
    slRND <- count_stretches(siRND) # stretch lengths of siRND
    slSFR <- count_stretches(siSFR) # stretch lengths of SFR
    table(slRND)
    table(slSFR)
    opar <- par(mfrow = c(2, 1))
    on.exit(par(opar), add = TRUE)
    hist(slRND, breaks = 0:20 + 0.5, xlim = c(0, 20))
    hist(slSFR, breaks = 0:20 + 0.5, xlim = c(0, 20))
}
