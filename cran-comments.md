
## R CMD check results

0 errors | 0 warnings | 1 note

## Resubmission

This is a resubmission.

In this version all issues identified by CRAN Maintainer Konstanze Lauseker concerning MetaboDecon1D v1.2.3 (submitted on Thu 16/01/2025) have been addressed by the following changes:

1. Explained acronyms in the description text (in particular 'NMR').
2. Added `\value` descriptions to `is_metabodecon_class.Rd` and `print_methods.Rd`
3. Improved examples for following functions:
   - `read_spectrum()`: replaced `\dontrun{}` with `\donttest{}`. (Unfortunatelly the example is, although indeed executable, very slow: 30s on my machine. A faster example can not easily be provided, as it is supposed to showcase the reading of JDX files which are not trivial to generate from scratch. The original example file we provide is output by Bruker TopSpin software, which is proprietary and requires a subscription.)
   - `deconvolute()`: removed the `\dontrun{}` part of the example as it did not add value to the previous code snippets. (It essentially used a different dataset as input and set `ask=TRUE` instead of `ask=FALSE`, which does not really justify a separate code snippet.)
   - `MetaboDecon1D()`: Replaced the previous `\dontrun{}` example with a new one that is more informative and does not require user input. (The previous version required the user to answer 10 questions, which is why it was wrapped in `\dontrun{}`.)
