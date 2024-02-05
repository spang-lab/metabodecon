test_that("with()", {

    before_wd <- getwd()

    x <- with(
        testdir = "with/1",
        answers = c("y", "n"),
        output = "captured", message = "captured", plots = "plots.pdf",
        inputs = c(urine.dx = "jcampdx/urine/urine.dx"),
        expr = {
            cat2("Helloworld!") # output is captured
            readline("Continue?") # readline is mocked
            readline("Continue?") # readline is mocked
            message("Roar") # messages are captured
            warning("Blub") # warnings are transformed to messages and captured as well
            test_wd <- getwd() # working dir is set to '{testdir()}/with/1'
            z <- 3 # vars can be assigned
            2 # return value of expression is captured in rv
        }
    )

    after_wd <- getwd()

    expect_true(file.exists(file.path(x$testdir, "plots.pdf")))
    expect_true(file.exists(file.path(x$testdir, "urine.dx")))
    expect_equal(x$rv, 2)
    expect_equal(z, 3)
    expect_true(x$runtime <= 1)
    expect_equal(x$output$text, "Helloworld!")
    expect_equal(x$message$text, c("Continue?y", "Continue?n", "Roar", "Warning: Blub"))
    expect_equal(test_wd, file.path(testdir(), "with/1"))
    expect_equal(after_wd, before_wd)
})
