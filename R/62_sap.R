
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Module Description #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

# Functions for creating "Simple as Possible" (SAP) spectrum objects and
# corresponding deconvolution objects. The SAP spectra should be so simple, that
# the expected deconvolution results can be either:
#
# 1. Calculated by hand, i.e. without calling a deconvolution function, or
# 2. Checked for correctness by inspecting the deconvolution object. Inspection
#    can be done by looking at the object structure or by plotting the object.
#
# This way we get a a set of objects that we can use to test wether our
# deconvolution and conversion functions produce the expected results.


# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# SAP2 #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @noRd
#' @title Create the Simple-As-Possible-Spectrum-v2
#' @details
#'
#' ```
#' -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#' |      SFR      |               w               |     SFR      |
#' |               |  x           www       p      |              |
#  |               | xxxa        wwwww     ppp     |              |
#  |               |xxxaaa      wwwwwww  ppppppp   |              |
#' |~-~-~-~-~-~-~-~|~-|-|-~-~-~-~-~|~-~-~-~-|-~-~-~-~-~-~-~-~-~-~-~
#' |               |  | |          |        |      |
#' 6.4             |  | 2.24       0.047    -2.22  -3.2
#'                 |  2.61
#'                 3.2
#' ```
#'
get_sap2_spectrum <- function() {
    spectrum <- simulate_spectrum(
        name   = "sap2",
        cs     = (64:-63) / 10,
        x0     = c(2.81, 2.10, 0.047, -2.22),
        A      = c(2000,  500,  5000,   1500),
        lambda = c( 0.4,  0.3,   0.6,   0.5),
        noise  = round(rnorm(128, sd = 50)),
        fqref  = 6e8,
        pkr    = c(3.2, -3.2)
    )
}

get_sap2_ispec <- function() {
    as_ispec(get_sap2_spectrum())
}

get_sap2_idecon <- function(rmws = FALSE, bwc = 2) {
    if (rmws %!=% FALSE|| bwc %!=% 2) {
        call <- deparse(match.call())
        stopf("%s is not implemented yet. Use rmws=FALSE and bwc=2.", call, call. = FALSE)
    }
    sap_decon <- deconvolute_ispec(
        ispec = as_ispec(get_sap2_spectrum()),
        sfr = c(3.2, -3.2),
        smopts = c(0, 3),
        delta = 1,
        verbose = FALSE
    )
}

get_optimal_values <- function(spectrum) {
    # Meta Data
    name   <- spectrum$sap1
    sf     <- c(1e3, 1e6)
    fq_ref <- convert_pos(0, spectrum$cs, spectrum$meta$fq)
    # Lorentz params
    x0_ppm     <- spectrum$meta$simpar$x0
    A_raw_ppm  <- spectrum$meta$simpar$A
    lambda_ppm <- spectrum$meta$simpar$lambda
    # Axis values
    cs             <- spectrum$cs
    noise          <- spectrum$meta$simpar$noise
    si_sup         <- lorentz_sup(cs, x0_ppm, A_raw_ppm, lambda_ppm)
    si_sup_norm    <- si_sup / sum(si_sup)
    si_sup_sc      <- si_sup / sf[2]
    si_sup_sc_norm <- si_sup_sc / sum(si_sup_sc)
    si             <- si_sup + noise
    si_norm        <- si / sum(si)
    si_sc          <- si / sf[2]
    si_sc_norm     <- si_sc / sum(si_sc)
    n              <- length(si)
    dpi            <- 1:n
    dpn            <- (n-1):0
    sdp            <- dpn / sf[1]
    fq             <- fq_ref - (fq_ref / 1e6) * cs
    cs_step        <- 0.1
    dpn_step       <- 1
    fq_step        <- cs_step * (fq_ref / 1e6)
    sdp_step       <- dpn_step / sf[1]
    # Peak values
    peak_ppm <- data.frame(
        left = x0_ppm + lambda_ppm,
        center = x0_ppm,
        right = x0_ppm - lambda_ppm
    )
    peak_idx <- data.frame(
        left = round(convert_pos(peak_ppm$left, cs, dpi)),
        center = round(convert_pos(peak_ppm$center, cs, dpi)),
        right = round(convert_pos(peak_ppm$right, cs, dpi))
    )
    peak_sdp <- data.frame(
        left = convert_pos(peak_ppm$left, cs, sdp),
        center = convert_pos(peak_ppm$center, cs, sdp),
        right = convert_pos(peak_ppm$right, cs, sdp)
    )
    # Lorentz params converted
    x0_dp  <- convert_pos(x0_ppm, cs, dpn) # x0_dp  = 7
    x0_sdp <- convert_pos(x0_ppm, cs, sdp) # x0_sdp = 0.007
    x0_hz  <- convert_pos(x0_ppm, cs, fq)  # x0_hz  = 599999880
    A_raw_dp  <- A_raw_ppm * (dpn_step  / cs_step) # A_raw_dp  = 20000
    A_raw_sdp <- A_raw_ppm * (sdp_step / cs_step)  # A_raw_sdp = 20
    A_raw_hz  <- A_raw_ppm * (fq_step  / cs_step)  # A_raw_hz  = 1200000
    A_sc_ppm  <- A_raw_ppm / sf[2] # A_sc_ppm  = 0.002
    A_sc_dp   <- A_raw_dp  / sf[2] # A_sc_dp   = 0.02
    A_sc_sdp  <- A_raw_sdp / sf[2] # A_sc_sdp  = 0.00002
    A_sc_hz   <- A_raw_hz  / sf[2] # A_sc_hz   = 1.2
    lambda_dp  <- convert_width(lambda_ppm, cs, dpn)     # lambda_dp  = 0.8
    lambda_sdp <- convert_width(lambda_ppm, cs, sdp)     # lambda_sdp = 0.0008
    lambda_hz  <- abs(convert_width(lambda_ppm, cs, fq)) # lambda_hz  = 48
    # Deconv params
    nfit <- 3
    smopts <- c(0, 3)
    delta <- 6.4
    wshw <- 0
    sfr_ppm <- c(1, -1) # outside of spectrum range
    sfr_dp  <- convert_pos(sfr_ppm, cs, dpn) # c(15, -5)
    sfr_sdp <- convert_pos(sfr_ppm, cs, sdp) # c(0.015, -0.005)
    sfr_sdp_0 <- sfr_in_sdp_bwc(sfr_ppm, cs, sf) # c(15, -5)
    sfr_hz  <- convert_pos(sfr_ppm, cs, fq)  # c(599999400, 600000600)
    args <- named(nfit, smopts, delta, sfr = sfr_ppm, wshw)
    # Integrals
    # (use n-step intervals instead of n-1-steps to reproduce (gy) results from v0.x)
    cs_range_0  <- c(min(cs ), max(cs ) + cs_step ) # cs_range_0  = c(-0.5, 0.6)
    dpn_range_0 <- c(min(dpn), max(dpn) + dpn_step) # dpn_range_0 = c(0, 11)
    sdp_range_0 <- c(min(sdp), max(sdp) + sdp_step) # sdp_range_0 = c(0, 0.011)
    fq_range_0  <- c(min(fq ), max(fq ) + fq_step ) # fq_range_0  = c(599999700, 600000360)
    integrals_rawsi_ppm   <- lorentz_int(x0_ppm, A_raw_ppm, lambda_ppm, limits = range(cs ))   # integrals_rawsi_ppm   = 5534.39650939749
    integrals_rawsi_dp    <- lorentz_int(x0_dp,  A_raw_dp,  lambda_dp,  limits = range(dpn))   # integrals_rawsi_dp    = 55343.9650939749
    integrals_rawsi_sdp   <- lorentz_int(x0_sdp, A_raw_sdp, lambda_sdp, limits = range(sdp))   # integrals_rawsi_sdp   = 55.3439650939749
    integrals_rawsi_hz    <- lorentz_int(x0_hz,  A_raw_hz,  lambda_hz,  limits = range(fq ))   # integrals_rawsi_hz    = 3320637.90563849
    integrals_rawsi_ppm_0 <- lorentz_int(x0_ppm, A_raw_ppm, lambda_ppm, limits = cs_range_0)   # integrals_rawsi_ppm_0 = 5660.81017319241
    integrals_rawsi_dp_0  <- lorentz_int(x0_dp,  A_raw_dp,  lambda_dp,  limits = dpn_range_0)  # integrals_rawsi_dp_0  = 56608.1017319241
    integrals_rawsi_sdp_0 <- lorentz_int(x0_sdp, A_raw_sdp, lambda_sdp, limits = sdp_range_0)  # integrals_rawsi_sdp_0 = 56.6081017319241
    integrals_rawsi_hz_0  <- lorentz_int(x0_hz,  A_raw_hz,  lambda_hz,  limits = fq_range_0)   # integrals_rawsi_hz_0  = 3337585.93122155
    integrals_scsi_ppm   <- integrals_rawsi_ppm   / sf[2] # integrals_sc_ppm   = 0.00553439650939
    integrals_scsi_dp    <- integrals_rawsi_dp    / sf[2] # integrals_sc_dp    = 0.05534396509397
    integrals_scsi_sdp   <- integrals_rawsi_sdp   / sf[2] # integrals_sc_sdp   = 0.00005534396509
    integrals_scsi_hz    <- integrals_rawsi_hz    / sf[2] # integrals_sc_hz    = 3.32063790563849
    integrals_scsi_ppm_0 <- integrals_rawsi_ppm_0 / sf[2] # integrals_sc_ppm_0 = 0.00566081017319
    integrals_scsi_dp_0  <- integrals_rawsi_dp_0  / sf[2] # integrals_sc_dp_0  = 0.05660810173192
    integrals_scsi_sdp_0 <- integrals_rawsi_sdp_0 / sf[2] # integrals_sc_sdp_0 = 0.00005660810173
    integrals_scsi_hz_0  <- integrals_rawsi_hz_0  / sf[2] # integrals_sc_hz_0  = 3.33758593122155
    # MSEs
    vae_raw_si  <- abs(si - si_sup)           # Vector of absolute errors for raw SIs
    vae_norm_si <- abs(si_norm - si_sup_norm) # Vector of absolute errors for normalized SIs
    vse_raw_si  <- vae_raw_si^2               # Vector of squared errors for raw SIs
    vse_norm_si <- vae_norm_si^2              # Vector of squared errors for normalized SIs
    mae_raw_si  <- mean(vae_raw_si)           # Mean absolute error for raw SIs
    mae_norm_si <- mean(vae_norm_si)          # Mean absolute error for normalized SIs
    mse_raw_si  <- mean(vse_raw_si)           # Mean squared  error for raw SIs
    mse_norm_si <- mean(vse_norm_si)          # Mean squared  error for normalized SIs
    # Compound Objects
    simpar <- named(name = "sap1", cs, pkr = NULL, x0 = x0_ppm, A = A_raw_ppm, lambda = lambda_ppm, noise)
    lcpar  <- named(A = A_sc_ppm, lambda = lambda_ppm, x0 = x0_ppm)
    args   <- named(nfit, smopts, delta, sfr = sfr_ppm, wshw)
    meta   <- named(name, fq, simpar)
    sit    <- named(wsrm = si, nvrm = si, sm = si, sup = si_sup, al = NULL)
    mse    <- named(raw = mse_raw_si, norm = mse_norm_si, sm = mse_raw_si, smnorm = mse_norm_si)
    peak   <- peak_idx
    # Return Objects
    locals()
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# SAP1 #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @noRd
#'
#' @title Create, update or get SAP1
#'
#' @description
#' SAP1 is a list containing a "simple as possible" (SAP) spectrum object in
#' `spectrum` and `ispec` format as well as the corresponding deconvolution
#' object in `decon0`, `decon1`, `decon2` and `idecon` format.
#'
#' Use
#' `make_sap1()` to create the SAP1 object from scratch,
#' `update_sap1()` to save it to the package and
#' `get_sap1()` to load it from the package.
#'
#'
#' @details
#' The SAP spectrum contains 11 datapoints ranging from -0.5 to 0.5 ppm with one
#' peak at 0.2 ppm, i.e. at the 4th datapoint. At the far right a huge noise
#' peak is added to achieve a non-zero MSE. The spectrum has no water signal and
#' no signal free region. It's most important properties are shown below, with
#' the frequency (fq) is given as difference to the reference frequency (6e8).
#' I.e. fq* = +60 means fq = 600000060 Hz.
#'
#' ```
#'                  left  centr  rigt
#' si   1660  3448  9756  25000  9756  3448  1660   962   624   437   322
#' cs    0.5   0.4   0.3    0.2   0.1     0  -0.1  -0.2  -0.3  -0.4  -0.5
#' dpi     1     2     3      4     5     6     7     8     9    10    11
#' dpn    10     9     8      7     6     5     4     3     2     1     0
#' sdp  .010  .009  .008   .007  .006  .005  .004  .003  .002  .001  .000
#' fq*  -300  -240  -180   -120   -60     0   +60  +120  +180  +240  +300
#' ```
#'
#' On a Windows Machine with Samsung 980 Pro SSD and AMD Ryzen 7 5700U CPU the
#' following runtimes were observed:
#'
#' 1. Creation from scratch via [make_sap1()]: approx. 0.11 seconds.
#' 2. Reading from disk via [get_sap1()]: approx. 0.01 seconds
#'
make_sap1 <- function() {
    # Meta Data
    name   <- "sap1"
    sf     <- c(1e3, 1e6)
    fq_ref <- 6e8
    # Lorentz params
    x0_ppm     <- 0.2
    A_raw_ppm  <- 2000
    lambda_ppm <- 0.08
    # Axis values
    cs             <- (5:-5) / 10
    noise          <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10000)
    si_sup         <- lorentz_sup(cs, x0_ppm, A_raw_ppm, lambda_ppm)
    si_sup_norm    <- si_sup / sum(si_sup)
    si_sup_sc      <- si_sup / sf[2]
    si_sup_sc_norm <- si_sup_sc / sum(si_sup_sc)
    si             <- si_sup + noise
    si_norm        <- si / sum(si)
    si_sc          <- si / sf[2]
    si_sc_norm     <- si_sc / sum(si_sc)
    n              <- length(si)
    dpi            <- 1:n
    dpn            <- (n-1):0
    sdp            <- dpn / sf[1]
    fq             <- fq_ref - (fq_ref / 1e6) * cs
    cs_step        <- 0.1
    dpn_step       <- 1
    fq_step        <- cs_step * (fq_ref / 1e6)
    sdp_step       <- dpn_step / sf[1]
    # Peak values
    peak_idx <- data.frame(left = 3,     center = 4,     right = 5)
    peak_ppm <- data.frame(left = 0.3,   center = 0.2,   right = 0.1)
    peak_sdp <- data.frame(left = 0.008, center = 0.007, right = 0.006)
    # Lorentz params converted
    x0_dp  <- convert_pos(x0_ppm, cs, dpn) # x0_dp  = 7
    x0_sdp <- convert_pos(x0_ppm, cs, sdp) # x0_sdp = 0.007
    x0_hz  <- convert_pos(x0_ppm, cs, fq)  # x0_hz  = 599999880
    A_raw_dp  <- A_raw_ppm * (dpn_step  / cs_step) # A_raw_dp  = 20000
    A_raw_sdp <- A_raw_ppm * (sdp_step / cs_step)  # A_raw_sdp = 20
    A_raw_hz  <- A_raw_ppm * (fq_step  / cs_step)  # A_raw_hz  = 1200000
    A_sc_ppm  <- A_raw_ppm / sf[2] # A_sc_ppm  = 0.002
    A_sc_dp   <- A_raw_dp  / sf[2] # A_sc_dp   = 0.02
    A_sc_sdp  <- A_raw_sdp / sf[2] # A_sc_sdp  = 0.00002
    A_sc_hz   <- A_raw_hz  / sf[2] # A_sc_hz   = 1.2
    lambda_dp  <- convert_width(lambda_ppm, cs, dpn)     # lambda_dp  = 0.8
    lambda_sdp <- convert_width(lambda_ppm, cs, sdp)     # lambda_sdp = 0.0008
    lambda_hz  <- abs(convert_width(lambda_ppm, cs, fq)) # lambda_hz  = 48
    # Deconv params
    nfit <- 3
    smopts <- c(0, 3)
    delta <- 6.4
    wshw <- 0
    sfr_ppm <- c(1, -1) # outside of spectrum range
    sfr_dp  <- convert_pos(sfr_ppm, cs, dpn) # c(15, -5)
    sfr_sdp <- convert_pos(sfr_ppm, cs, sdp) # c(0.015, -0.005)
    sfr_sdp_0 <- sfr_in_sdp_bwc(sfr_ppm, cs, sf) # c(15, -5)
    sfr_hz  <- convert_pos(sfr_ppm, cs, fq)  # c(599999400, 600000600)
    args <- named(nfit, smopts, delta, sfr = sfr_ppm, wshw)
    # Integrals
    # (use n-step intervals instead of n-1-steps to reproduce (gy) results from v0.x)
    cs_range_0  <- c(min(cs ), max(cs ) + cs_step ) # cs_range_0  = c(-0.5, 0.6)
    dpn_range_0 <- c(min(dpn), max(dpn) + dpn_step) # dpn_range_0 = c(0, 11)
    sdp_range_0 <- c(min(sdp), max(sdp) + sdp_step) # sdp_range_0 = c(0, 0.011)
    fq_range_0  <- c(min(fq ), max(fq ) + fq_step ) # fq_range_0  = c(599999700, 600000360)
    integrals_rawsi_ppm   <- lorentz_int(x0_ppm, A_raw_ppm, lambda_ppm, limits = range(cs ))   # integrals_rawsi_ppm   = 5534.39650939749
    integrals_rawsi_dp    <- lorentz_int(x0_dp,  A_raw_dp,  lambda_dp,  limits = range(dpn))   # integrals_rawsi_dp    = 55343.9650939749
    integrals_rawsi_sdp   <- lorentz_int(x0_sdp, A_raw_sdp, lambda_sdp, limits = range(sdp))   # integrals_rawsi_sdp   = 55.3439650939749
    integrals_rawsi_hz    <- lorentz_int(x0_hz,  A_raw_hz,  lambda_hz,  limits = range(fq ))   # integrals_rawsi_hz    = 3320637.90563849
    integrals_rawsi_ppm_0 <- lorentz_int(x0_ppm, A_raw_ppm, lambda_ppm, limits = cs_range_0)   # integrals_rawsi_ppm_0 = 5660.81017319241
    integrals_rawsi_dp_0  <- lorentz_int(x0_dp,  A_raw_dp,  lambda_dp,  limits = dpn_range_0)  # integrals_rawsi_dp_0  = 56608.1017319241
    integrals_rawsi_sdp_0 <- lorentz_int(x0_sdp, A_raw_sdp, lambda_sdp, limits = sdp_range_0)  # integrals_rawsi_sdp_0 = 56.6081017319241
    integrals_rawsi_hz_0  <- lorentz_int(x0_hz,  A_raw_hz,  lambda_hz,  limits = fq_range_0)   # integrals_rawsi_hz_0  = 3337585.93122155
    integrals_scsi_ppm   <- integrals_rawsi_ppm   / sf[2] # integrals_sc_ppm   = 0.00553439650939
    integrals_scsi_dp    <- integrals_rawsi_dp    / sf[2] # integrals_sc_dp    = 0.05534396509397
    integrals_scsi_sdp   <- integrals_rawsi_sdp   / sf[2] # integrals_sc_sdp   = 0.00005534396509
    integrals_scsi_hz    <- integrals_rawsi_hz    / sf[2] # integrals_sc_hz    = 3.32063790563849
    integrals_scsi_ppm_0 <- integrals_rawsi_ppm_0 / sf[2] # integrals_sc_ppm_0 = 0.00566081017319
    integrals_scsi_dp_0  <- integrals_rawsi_dp_0  / sf[2] # integrals_sc_dp_0  = 0.05660810173192
    integrals_scsi_sdp_0 <- integrals_rawsi_sdp_0 / sf[2] # integrals_sc_sdp_0 = 0.00005660810173
    integrals_scsi_hz_0  <- integrals_rawsi_hz_0  / sf[2] # integrals_sc_hz_0  = 3.33758593122155
    # MSEs
    vae_raw_si  <- abs(si - si_sup)           # Vector of absolute errors for raw SIs
    vae_norm_si <- abs(si_norm - si_sup_norm) # Vector of absolute errors for normalized SIs
    vse_raw_si  <- vae_raw_si^2               # Vector of squared errors for raw SIs
    vse_norm_si <- vae_norm_si^2              # Vector of squared errors for normalized SIs
    mae_raw_si  <- mean(vae_raw_si)           # Mean absolute error for raw SIs
    mae_norm_si <- mean(vae_norm_si)          # Mean absolute error for normalized SIs
    mse_raw_si  <- mean(vse_raw_si)           # Mean squared  error for raw SIs
    mse_norm_si <- mean(vse_norm_si)          # Mean squared  error for normalized SIs
    # Compound Objects
    simpar <- named(name = "sap1", cs, pkr = NULL, x0 = x0_ppm, A = A_raw_ppm, lambda = lambda_ppm, noise)
    lcpar  <- named(A = A_sc_ppm, lambda = lambda_ppm, x0 = x0_ppm)
    args   <- named(nfit, smopts, delta, sfr = sfr_ppm, wshw)
    meta   <- named(name, fq, simpar)
    sit    <- named(wsrm = si, nvrm = si, sm = si, sup = si_sup, al = NULL)
    mse    <- named(raw = mse_raw_si, norm = mse_norm_si, sm = mse_raw_si, smnorm = mse_norm_si)
    peak   <- peak_idx
    # Return Objects
    spectrum <- structure(class = "spectrum", named(si, cs, meta))
    decon0 <- list(
        number_of_files            = 1L,
        filename                   = "sap1",
        x_values                   = sdp,
        x_values_ppm               = cs,
        y_values                   = si_sc,
        spectrum_superposition     = structure(si_sup_sc, dim = c(1L, 11L)),
        mse_normed                 = mse_norm_si,
        index_peak_triplets_middle = peak_idx[1, 2],
        index_peak_triplets_left   = peak_idx[1, 3],  # (1)
        index_peak_triplets_right  = peak_idx[1, 1],  # (1)
        peak_triplets_middle       = peak_ppm[1, 2],
        peak_triplets_left         = peak_ppm[1, 3],  # (1)
        peak_triplets_right        = peak_ppm[1, 1],  # (1)
        integrals                  = structure(integrals_scsi_sdp_0, dim = c(1, 1)), # (2)
        signal_free_region         = sfr_sdp_0,    # (3)
        range_water_signal_ppm     = 0,
        A                          = - A_sc_sdp,   # (4)
        lambda                     = - lambda_sdp, # (4)
        x_0                        = x0_sdp
    )
    decon1 <- structure(class = "decon1", c(decon0, list(
        y_values_raw           = si,
        x_values_hz            = fq,
        mse_normed_raw         = mse_norm_si,
        signal_free_region_ppm = sfr_ppm,
        x_0_hz                 = x0_hz,
        x_0_dp                 = x0_dp,
        x_0_ppm                = x0_ppm,
        A_hz                   = A_sc_hz,
        A_dp                   = - A_sc_dp,    # (4)
        A_ppm                  = - A_sc_ppm,   # (4)
        lambda_hz              = lambda_hz,
        lambda_dp              = - lambda_dp,  # (4)
        lambda_ppm             = - lambda_ppm  # (4)
    )))
    decon2 <- structure(class = "decon2", named(cs, si, meta, args, sit, peak, lcpar, mse))
    ispec <- structure(class = "ispec", list(
        y_raw     = si,
        y_scaled  = si_sc,
        n         = n,
        dpn       = dpn,
        sdp       = sdp,
        sf        = sf,
        ppm       = cs,
        hz        = fq,
        ppm_range = diff(range(cs)),
        ppm_max   = max(cs),
        ppm_min   = min(cs),
        ppm_step  = cs_step,
        ppm_nstep = diff(range(cs)) / length(cs),
        name      = "sap1",
        meta      = meta
    ))
    idecon <- structure(class = "idecon", c(ispec, list(
        TODO     = "TODO: Finish idecon implementation",
        sf       = sf,
        y_nows   = si,
        y_pos    = si,
        Z        = NULL,
        y_smooth = si,
        d        = calc_second_derivative(si),
        peak     = NULL,
        lci      = NULL,
        lca      = NULL,
        lcr      = NULL
    )))
    # (1): right and left are switched in decon0/decon1
    # (2): integrals are calculated over [0, n] instead of [0, n-1] in decon0/decon1
    # (3): conversion of ppm to sdp is slightly off in decon0/decon1
    # (4): x_0 and lambda are returned as negative values in decon0/decon1
    named(spectrum, decon0, decon1, decon2, ispec, idecon)
}

update_sap1 <- function(overwrite = FALSE, verbose = TRUE) {
    sap1 <- make_sap1()
    if (isTRUE(overwrite)) {
        xds_path  <- pkg_file("example_datasets")
        rds_path  <- file.path(xds_path, "rds") # might not exist
        sap1_path <- file.path(rds_path, "sap1.rds")
        if (!dir.exists(rds_path)) dir.create(rds_path)
        if (verbose) logf("Overwrite is TRUE. Updating %s." , sap1_path)
        saveRDS(sap1, sap1_path)
    } else {
        sap1_path <- file.path(tmpdir(subdir = TRUE, create = TRUE), "sap1.rds")
        if (verbose) logf("Overwrite is FALSE. Writing spectra to %s." , sap1_path)
        saveRDS(sap1, sap1_path)
    }
    sap1_path
}

get_sap1 <- function() {
    xds_path  <- pkg_file("example_datasets")
    rds_path  <- file.path(xds_path, "rds")
    sap1_path <- file.path(rds_path, "sap1.rds")
    sap1 <- readRDS(sap1_path)
    sap1
}


