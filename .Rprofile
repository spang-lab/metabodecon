cat("Setting up languageserver.formatting_style\n")
options(
    languageserver.formatting_style = function(options) {
        transformers <- styler::tidyverse_style(
            indent_by = options$tabSize
        )
        unindent_fun_dec <- transformers$indention$unindent_fun_dec
        formals(unindent_fun_dec)$indent_by <- options$tabSize
        transformers$indention$unindent_fun_dec <- unindent_fun_dec
        transformers$space$set_no_space_around_eq_sub <- function(pd) {
            is_eq_sub <- pd$token == "EQ_SUB"
            is_before_eq_sub <- c(is_eq_sub[-1], FALSE)
            pd$spaces[is_eq_sub | is_before_eq_sub] <- 0L
            pd
        }
        transformers
    }
)
print(options("languageserver.formatting_style"))
