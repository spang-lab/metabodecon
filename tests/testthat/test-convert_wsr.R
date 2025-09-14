test_that("enrich_wshw works", {
    ppm <- c(4.7, 3.4, 2.1, 0.8, -0.5, -1.8, -3.1, -4.4)
    n <- length(ppm)
    ppm_nstep <- (max(ppm) - min(ppm)) / (n)
    ispec <- structure(list(ppm = ppm, n = n, ppm_nstep = ppm_nstep), class = "ispec")
    wsr <- enrich_wshw(wshw = 0.2, ispec)
    # Below results are incorrect but expected nonetheless to maintain backwards
    # compatibility. For details see 'CHECK-3: water signal calculation' from
    # TODOS.md (no longer part of the repository). To retrieve it, fetch a
    # commit older than 2025-09-14.
    expect_equal(rev(wsr), list(
        hwidth_dp = 0.175824175824176, # hwidth_ppm / spec$ppm_nstep
        center_dp = 4, # spec$n / 2
        right_dp = 4.17582417582418, # center_dp + hwidth_dp
        left_dp = 3.82417582417582, # center_dp - hwidth_dp
        hwidth_ppm = 0.2, # wshw
        center_ppm = 0.8, # spec$ppm[center_dp]
        right_ppm = 0.8, # spec$ppm[right_dp]
        left_ppm = 2.1 # spec$ppm[left_dp]
    ))
})
