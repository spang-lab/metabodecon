test_that("get_names() works", {

    xx <- list( a = list(name="c", meta=list(name="e")),
                b = list(name="d", meta=list(name="f")) )
    expect_equal(get_names(xx), c("e", "f"))

    xx <- list( a = list(name="foo"),
                b = list(name="bar"))
    expect_equal(get_names(xx), c("foo", "bar"))

    xx <- list(a = list(), b = list())
    expect_equal(get_names(xx), c("a", "b"))

    xx <- list(a = 1, b = list(name = "foo"))
    expect_equal(get_names(xx), c("a", "foo"))

    xx <- list(a = 1, b = 2)
    expect_equal(get_names(xx), c("a", "b"))

    xx <- c(a = 1, b = 2)
    expect_equal(get_names(xx), c("a", "b"))

    xx <- c(a = 1, 2)
    expect_equal(get_names(xx), c("a", "spectrum_2"))

})

