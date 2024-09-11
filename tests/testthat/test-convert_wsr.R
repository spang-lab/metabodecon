test_that("enrich_wshw works", {
    ppm <- c(4.7, 3.4, 2.1, 0.8, -0.5, -1.8, -3.1, -4.4)
    n <- length(ppm)
    ppm_nstep <- (max(ppm) - min(ppm)) / (n)
    spec <- list(ppm = ppm, n = n, ppm_nstep = ppm_nstep)
    wsr <- enrich_wshw(spec, wshw = 0.2)
    # Below results are incorrect but expected nonetheless to maintain backwards compatibility. For details see [CHECK-3: water signal calculation](TODOS.md).
    expect_equal(rev(wsr), list(
        hwidth_ppm = 0.2, # wshw
        hwidth_dp = 0.175824175824176, # hwidth_ppm / spec$ppm_nstep
        center_dp = 4, # spec$n / 2
        right_dp = 4.17582417582418, # center_dp + hwidth_dp
        left_dp = 3.82417582417582, # center_dp - hwidth_dp
        center_ppm = 0.8, # spec$ppm[center_dp]
        right_ppm = 0.8, # spec$ppm[right_dp]
        left_ppm = 2.1 # spec$ppm[left_dp]
    ))
})
