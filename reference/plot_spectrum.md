# Plot Spectrum

Plot a spectrum and zoom in on a specific region.

**\[experimental\]**

## Usage

``` r
plot_spectrum(
  x,
  ...,
  obj = as_v12_singlet(x),
  foc_frac = get_foc_frac(obj),
  foc_rgn = get_foc_rgn(obj, foc_frac),
  sub1 = TRUE,
  sub2 = FALSE,
  sub3 = width(foc_rgn) < width(obj$cs),
  mar = NULL,
  frame = FALSE,
  con_lines = TRUE
)
```

## Arguments

- x:

  An object of type `spectrum`, `decon0`, `decon1`, `decon2` or `align`.
  For details see [Metabodecon
  Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).

- ...:

  Additional arguments passed to
  [`draw_spectrum()`](https://spang-lab.github.io/metabodecon/reference/draw_spectrum.md)
  for **every** sub figure. See 'Details'.

- obj:

  An object of type `spectrum` or `decon2`. Usually auto generated from
  `x`, but can be set manually in case the default conversion is not
  sufficient.

- foc_frac:

  A numeric vector specifying the start and end of the focus region as
  fraction of the full spectrum width. Only used if `foc_rgn` is set to
  NULL.

- foc_rgn:

  A numeric vector specifying the start and end of the focus region in
  ppm. If set to NULL, `foc_frac` is used to determine the focus region.
  If both `foc_rgn` and are set to NULL, a suitable focus region is
  chosen automatically. Takes precedence over `foc_frac`.

- sub1, sub2, sub3:

  List of arguments passed to
  [`draw_spectrum()`](https://spang-lab.github.io/metabodecon/reference/draw_spectrum.md)
  when drawing sub figure 1-3. See 'Details'.

- mar:

  A numeric vector of length 4 specifying the margins of the plot.
  Passed to [`par()`](https://rdrr.io/r/graphics/par.html). If set to
  `NULL`, a suitable value is chosen automatically.

- frame:

  A list of values passed to
  [`box()`](https://rdrr.io/r/graphics/box.html) when drawing the frame
  around plot region. If set to `NULL`, no frame is drawn.

- con_lines:

  A list of values passed to
  [`lines()`](https://rdrr.io/r/graphics/lines.html) when drawing the
  connecting lines between sub figure 1 and the focus rectangle in sub
  figure 3. See 'Details'. If set to `NULL`, the connecting lines are
  not drawn.

## Value

NULL. Called for side effect of plotting as sketched in 'Details'.

## Details

This function first initializes a new plotting canvas. After that it
calls
[`draw_spectrum()`](https://spang-lab.github.io/metabodecon/reference/draw_spectrum.md)
multiple times to draw the following sub figures onto the plotting
canvas:

1.  The signal intensities in the focus region

2.  The second derivative in the focus region

3.  The signal intensities over all datapoints

The argument lists for the individual calls to
[`draw_spectrum()`](https://spang-lab.github.io/metabodecon/reference/draw_spectrum.md)
are determined at runtime and depend on the arguments passed to
`plot_spectrum()` as well as the currently active graphics device. To
customize the appearance of the individual sub plots, you can overwrite
each value passed to
[`draw_spectrum()`](https://spang-lab.github.io/metabodecon/reference/draw_spectrum.md)
by providing a corresponding named element in `sub1`, `sub2` or `sub3`.

A sketch of the resulting figure is shown below.

     __________________________________________
    |        ______________1_____________      |
    |       | Sub1: Signal Intensity in  |     |
    |       | Focus Region               |     |
    |       |             /\             |     |
    |       |            /  \            |     |
    |       |           /    \  /\       |     |
    |     11|          /      \/  \      |7    |
    |       |     /\  /            \     |     |
    |       |    /  \/              \    |     |
    |       |   /                    \   |     |
    |       |__/___________0__________\__|     |
    |       | Sub2: Second Derivative    |     |
    |     11| in Focus Region            |7    |
    |       |____________________________|     |
    |                      3                   |
    |    __________________3_________________  |
    |   |  Sub3: Signal Intensity over all   | |
    |   |  Datapoints     ________________   | |
    | 5 |                | Focus Rectangle|  |1|
    |   |     /\         |       /\       |  | |
    |   |    /  \        |      /  \/\    |  | |
    |   |   /    \   /\  |   /\/      \   |  | |
    |   |__/______\_/__\_|__/__________\__|__| |
    |______________________5___________________|

Note that the figure created by `plot_spectrum()` can be part of a
multi-figure configuration as created when setting `mfrow` or `mfcol`
via [`par()`](https://rdrr.io/r/graphics/par.html). Example:

    _______________________________________
    | Plot Spectrum with   | Other Figure  |
    | sub3 = TRUE          | Other Figure  |
    |      ___________     |  ___________  |
    |     | Sub Fig 1 |    | | x      x  | |
    |     |___________|    | |      x    | |
    |     |_Sub_Fig_2_|    | |      x    | |
    |   _________________  | |   x     x | |
    |  |    Sub Fig 3    | | |      x    | |
    |  |_________________| | |___________| |
    |______________________|_______________|
    | Some other Figure    | Plot Spectrum |
    |                      | sub3 = FALSE  |
    |  _________________   |  ___________  |
    | |     ___         |  | | Sub Fig 1 | |
    | | ___/   \___     |  | |           | |
    | |/           \____|  | |___________| |
    | |                 |  | | Sub Fig 2 | |
    | |_________________|  | |___________| |
    |______________________|_______________|

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
## 1. Prepare a deconvoluted spectrum as input

spec <- read_spectrum(metabodecon_file("sim/sim_01"))
decon <- generate_lorentz_curves_sim(spec)

## 2.1. Plot the full (non-deconvoluted) spectrum
## 2.2. Remove connecting lines, and focus on a specific region specified in ppm
## 2.3. Show second derivative and focus on a specific region specified as fraction
## 2.4. Change color of focus rectangle and margins of sub figure 1
## 2.5. Hide xlab and show second derivative
## 2.6. Change the figure region for sub figure 1

plot_spectrum(spec, sub1 = FALSE)

plot_spectrum(decon, foc_rgn = c(3.49, 3.45), con_lines = FALSE)

plot_spectrum(decon, sub2 = TRUE, foc_frac = c(0.40, 0.30))

plot_spectrum(decon,
    sub1 = list(mar = c(3, 6, 3, 6), lt_axis = list(col = "violet")),
    foc_rect = list(border = "violet", col = transp("violet")),
    con_lines = list(col = "violet")
)

plot_spectrum(decon,
    sub2 = TRUE,
    sub3 = list(bt_text = list(text = "")),
    frame = TRUE,
    con_lines = FALSE
)

plot_spectrum(decon, sub1 = list(fig_rgn_npc = c(0,1,.3,1), mar = c(0,5,0,0)))

```
