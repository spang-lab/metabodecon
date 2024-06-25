test_that("evalwith works", {
    before_wd <- getwd()
    x <- evalwith(
        testdir = "with/1",
        answers = c("y", "n"),
        output = "captured", message = "captured", plot = "plots.pdf",
        inputs = "jcampdx/urine/urine_1.dx",
        expr = {
            cat2("Helloworld!") # output is captured
            readline("Continue?") # readline is mocked
            readline("Continue?") # readline is mocked
            message("Roar") # messages are captured
            warning("Blub") # warnings are transformed to messages and captured as well
            test_wd <- getwd() # working dir is set to '{testdir()}/with/1'
            y <- 2
            z <- 3 # vars can be assigned
            list(y = y, z = z) # return value of expression is captured in rv
        }
    )
    after_wd <- getwd()
    expect_true(file.exists(file.path(test_wd, "plots.pdf")))
    expect_true(file.exists(file.path(test_wd, "urine_1.dx")))
    expect_equal(x$rv, list(y = 2, z = 3))
    expect_equal(z, 3)
    expect_true(x$runtime <= 1)
    expect_equal(x$output, "Helloworld!")
    expect_equal(x$message, c("Continue?y", "Continue?n", "Roar", "Warning: Blub"))
    expect_equal(test_wd, file.path(testdir(), "with/1"))
    expect_equal(after_wd, before_wd)
})
