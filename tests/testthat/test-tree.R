test_that("tree show.counts counts direct children", {
    root <- tempfile("tree-test-")
    dir.create(root)
    dir.create(file.path(root, "sample"), recursive = TRUE)
    dir.create(file.path(root, "sample", "10"), recursive = TRUE)
    file.create(file.path(root, "meta.txt"))

    out <- capture.output(tree(root, max.level = 1, show.counts = TRUE))
    txt <- paste(out, collapse = "\n")

    expect_true(grepl("\\[2\\]$", out[1]))
    expect_match(txt, "sample/ \\[1\\]")
})

test_that("tree can list files before directories", {
    root <- tempfile("tree-test-")
    dir.create(root)
    dir.create(file.path(root, "a_dir"), recursive = TRUE)
    file.create(file.path(root, "z_file.txt"))

    out <- capture.output(tree(root, max.level = 1, files.first = TRUE))
    file_line <- grep("z_file\\.txt", out)
    dir_line <- grep("a_dir/", out)

    expect_true(length(file_line) == 1)
    expect_true(length(dir_line) == 1)
    expect_true(file_line < dir_line)
})

test_that("tree shows omitted entry count", {
    root <- tempfile("tree-test-")
    dir.create(root)
    file.create(file.path(root, c("a.txt", "b.txt", "c.txt", "d.txt")))

    out <- capture.output(tree(root, max.level = 1, max.entries = 2))
    txt <- paste(out, collapse = "\n")

    expect_true(grepl("... [2]", txt, fixed = TRUE))
})

test_that("tree truncation with odd max.entries uses top-biased split", {
    root <- tempfile("tree-test-")
    dir.create(root)
    file.create(file.path(root, c("a.txt", "b.txt", "c.txt", "d.txt", "e.txt", "f.txt")))

    out <- capture.output(tree(root, max.level = 1, max.entries = 3, files.first = TRUE))
    shown <- out[-1]

    expect_true(any(grepl("a.txt$", shown)))
    expect_true(any(grepl("b.txt$", shown)))
    expect_true(any(grepl("f.txt$", shown)))
    expect_true(any(grepl("... [3]", shown, fixed = TRUE)))
    expect_false(any(grepl("c.txt$", shown)))
    expect_false(any(grepl("d.txt$", shown)))
    expect_false(any(grepl("e.txt$", shown)))
})

test_that("tree truncation with even max.entries splits equally", {
    root <- tempfile("tree-test-")
    dir.create(root)
    file.create(file.path(root, c("a.txt", "b.txt", "c.txt", "d.txt", "e.txt", "f.txt")))

    out <- capture.output(tree(root, max.level = 1, max.entries = 4, files.first = TRUE))
    shown <- out[-1]

    expect_true(any(grepl("a.txt$", shown)))
    expect_true(any(grepl("b.txt$", shown)))
    expect_true(any(grepl("e.txt$", shown)))
    expect_true(any(grepl("f.txt$", shown)))
    expect_true(any(grepl("... [2]", shown, fixed = TRUE)))
    expect_false(any(grepl("c.txt$", shown)))
    expect_false(any(grepl("d.txt$", shown)))
})

test_that("tree_preview uses compact preview defaults", {
    root <- tempfile("tree-test-")
    dir.create(root)
    file.create(file.path(root, c("a.txt", "b.txt", "c.txt", "d.txt")))
    dir.create(file.path(root, "z_dir"), recursive = TRUE)
    dir.create(file.path(root, "y_dir"), recursive = TRUE)

    out <- capture.output(tree_preview(root))
    txt <- paste(out, collapse = "\n")

    expect_true(grepl("\\[6\\]$", out[1]))
    expect_true(length(grep("a\\.txt", out)) == 1)
    expect_true(length(grep("y_dir/", out)) == 1)
    expect_true(length(grep("z_dir/", out)) == 1)
    expect_true(length(grep("d\\.txt", out)) == 1)
    expect_false(grepl("...", txt, fixed = TRUE))
})
