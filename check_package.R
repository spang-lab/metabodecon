# Used by .gitlab-ci.yml
devtools::install_deps(upgrade = FALSE)
options(repos =  c(CRAN = "https://cloud.r-project.org"))
devtools::check()
