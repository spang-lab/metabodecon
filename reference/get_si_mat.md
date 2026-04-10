# Extract Matrix of aligned Signal Intensities

Extracts a peak-area matrix from aligned spectra. Rows are
chemical-shift positions, columns are spectra. When `maxCombine = 0` the
raw aligned integral vectors are returned. When `maxCombine > 0`, peaks
are combined using one of two methods controlled by `combineMethod`.

## Usage

``` r
get_si_mat(x, ref = NULL, drop_zero = FALSE, maxCombine = 0, combineMethod = 2)
```

## Arguments

- x:

  An object of type `aligns`.

- ref:

  A single `align` or `decon2` object whose peaks define the rows of the
  output matrix. If `NULL` (default), the reference is auto-detected
  from `x` via `find_ref()`. Used by `combineMethod = 2` only.

- drop_zero:

  If `TRUE`, rows where all values are zero are removed.

- maxCombine:

  Controls peak combining in datapoints. `0` (default): off. Any
  positive number enables peak combining. For `combineMethod = 1` this
  is the `range` parameter of
  [`combine_peaks()`](https://spang-lab.github.io/metabodecon/reference/combine_peaks.md).
  For `combineMethod = 2` this is the maximum distance to the nearest
  reference peak. See 'Details'.

- combineMethod:

  Peak-combining backend. `1`: merge adjacent columns via
  [`combine_peaks()`](https://spang-lab.github.io/metabodecon/reference/combine_peaks.md).
  `2` (default): nearest-neighbour snapping to reference peaks.

## Value

A numeric matrix with chemical shifts as rownames and spectrum names as
colnames.

## Details

**Method 1** (`combineMethod = 1`) passes the transposed integral matrix
to
[`combine_peaks()`](https://spang-lab.github.io/metabodecon/reference/combine_peaks.md)
with `range = maxCombine`. Partly-filled neighbouring columns are merged
when they share no common non-zero row.

**Method 2** (`combineMethod = 2`) maps every peak to the nearest
reference peak and keeps it only if the distance is at most `maxCombine`
datapoints. Areas of peaks that map to the same reference peak are
summed.

Example for method 2 with ref peaks at indices 5, 9, 20 and
`maxCombine = 4`:

    Step 1 – Build intervals [ref ± maxCombine]:

      ref:     5              9                    20
               |              |                     |
      int:  [1 ····· 9]   [5 ···· 13]        [16 ···· 24]
                overlap!

    Step 2 – Shrink overlapping neighbours to midpoint.
             Refs 1 & 2 overlap → mid = floor((5+9)/2) = 7.
             Refs 2 & 3 don't  → keep maxCombine boundary.

      ref:     5         9                         20
               |         |                          |
      int:  [1 ··· 7] [8 ·· 13]   gap        [16 ···· 24]

    Step 3 – Assign peaks to nearest ref; keep if ≤ maxCombine:

      | Peak | Nearest | Dist | ≤ 4? | Action |
      |------|---------|------|------|--------|
      |    3 | ref 1   |    2 | yes  | keep   |
      |    6 | ref 1   |    1 | yes  | keep   |
      |    8 | ref 2   |    1 | yes  | keep   |
      |   11 | ref 2   |    2 | yes  | keep   |
      |   15 | ref 3   |    5 | no   | DROP   |
      |   21 | ref 3   |    1 | yes  | keep   |

    Step 4 – Sum areas (A × pi) per reference peak:

      ref 1 ← A(3) + A(6)    (peaks at 3, 6)
      ref 2 ← A(8) + A(11)   (peaks at 8, 11)
      ref 3 ← A(21)          (peak 15 dropped)

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
if (interactive()) {
    decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
    aligns <- align(decons, maxCombine = 0)
    si_mat_0 <- get_si_mat(aligns)
    si_mat_1 <- get_si_mat(aligns, maxCombine = 20)
    si_mat_2 <- get_si_mat(aligns, maxCombine = 20, combineMethod = 1)
}
```
