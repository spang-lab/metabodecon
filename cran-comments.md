
## R CMD check results

0 errors | 0 warnings | 1 note

## Resubmission

This is a resubmission.

In this version all issues identified by CRAN Maintainer Konstanze Lauseker concerning MetaboDecon1D v1.2.3 (submitted on Thu 16/01/2025) have been addressed by the following changes:

1. Explained acronyms in the description text (in particular 'NMR').
2. Added `\value` descriptions to `is_metabodecon_class.Rd` and `print_methods.Rd`
3. Improved examples for following functions:
   - `deconvolute()`: removed the `\dontrun{}` part of the example as it did not add value to the previous code snippets. (It only used a different dataset as input and set `ask=TRUE` instead of `ask=FALSE`, which does not really justify a separate code snippet.)
   - `MetaboDecon1D()`: Replaced the previous `\dontrun{}` example with a new one that is more informative and does not require user input. (The previous version required the user to answer 10 questions, which is why it was wrapped in `\dontrun{}`.)
   - `read_spectrum()`: removed the `\dontrun{}` part of the example as it did not add value to the previous code snippets. (It only showed how to read in a different file format, which required setting `file_format = "jcampdx"` instead of `file_format = "bruker"`, which does not really justify a separate code snippet.)
