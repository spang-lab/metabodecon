test_that("redirect works", {

    # Setup test env
    wd <- file.path(testdir(), "redirect/1")
    mkdirs(wd)
    owd <- setwd(wd)
    on.exit(setwd(owd), add = TRUE)

    # Run test
    redirects <- redirect(output = "captured", message = "captured", plots = "tmp.pdf")
    cat2("Hello")
    cat2("from cat")
    message("Goodbye")
    message("from message")
    plot(1:20)
    restore()

    # Check results
    expect_equal(redirects$output$text, c("Hello", "from cat"))
    expect_equal(redirects$message$text, c("Goodbye", "from message"))
    expect_equal(file.size("tmp.pdf"), 5413)
})
