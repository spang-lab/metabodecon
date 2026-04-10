test_that("human_readable formats difftime >= 1 hour as h + min", {
    x <- as.difftime(63, units = "mins")
    expect_equal(metabodecon:::human_readable(x), "1h 3min")
})

test_that("human_readable formats difftime >= 1 min as min + s", {
    x <- as.difftime(14 * 60 + 2, units = "secs")
    expect_equal(metabodecon:::human_readable(x), "14min 2s")
})

test_that("human_readable formats difftime < 1 min with 2 decimals", {
    x <- as.difftime(34.7, units = "secs")
    expect_equal(metabodecon:::human_readable(x), "34.70s")
})

test_that("human_readable formats negative difftime values", {
    x <- as.difftime(-75, units = "secs")
    expect_equal(metabodecon:::human_readable(x), "-1min 15s")
})
