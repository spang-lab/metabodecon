### Create new version

1. Clone the repo
2. Create a new branch
3. Make your desired changes
4. Open an R terminal inside the repository root folder
5. Use command `devtools::load_all()` to load the new package from source
6. Try the changed functions
7. Increase the version in [DESCRIPTION](DESCRIPTION)
8. Use `devtools::document()` to update the documentation files based on roxygen comments describing each function
9. Use `devtools::check()` to run all tests required by CRAN
10. Push you changes to gitlab
11. Create a merge request on gitlab
12. If all automatic checks have passed, accept the merge request

### Submit to CRAN

According to <https://r-pkgs.org/release.html> the following steps are necessary

```R
devtools::document() # Update documentation
rcmdcheck::rcmdcheck(
    # Run `R CMD check` for this package
    args=c("--no-manual", "--as-cran"),
    build_args=c("--no-manual"),
    check_dir="check"
)
if (FALSE) {
    # package is currently broken
    revdepcheck::revdep_check(num_workers = detectCores(logical=FALSE)) #1
    #1 Run `R CMD check` for all dependencies
    #1 For first time setup use `usethis::use_revdep()`
} else {
    # use function from core team instead
    if (!dir.exists("dist")) {
        dir.create("dist")
    }
    devtools::build(path="dist/")
    tools::check_packages_in_dir(
        dir = "dist/",
        check_args = "--as-cran",
        reverse = list(
            repos = getOption("repos")["CRAN"],
            which = "all",
            recursive = TRUE
        ),
        clean = FALSE
    )
    tools::summarize_check_packages_in_dir_results("dist/")
    tools::summarize_check_packages_in_dir_timings("dist/")
    tools::check_packages_in_dir_details("dist/")
}
# Update cran-comments.md
devtools::spell_check() # Check spelling of package
devtools::release() # Builds, tests and submits the package to CRAN.
# Manual submission can be done at: https://cran.r-project.org/submit.html
```
