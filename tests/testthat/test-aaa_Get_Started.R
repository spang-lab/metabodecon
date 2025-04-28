test_that("aaa_Get_Started works", {
    obj <- evalwith(
        output = "captured",
        opts = list(browser = cat2),
        expr = x <- aaa_Get_Started(open_browser = TRUE)
    )
    url <- "https://spang-lab.github.io/metabodecon/articles/Get_Started.html"
    expect_equal(x, url)
    expect_equal(obj$output, url)
})
