Read the functions `fit_mdm()`, `cv_mdm()` and `benchmark_mdm()`. They have
underdone quite some change and are a currently a work in progress.

I would like to fix/finalize them now, by combining `fit_mdm()` and `cv_mdm()`
into a single function `fit_mdm()` (which essentially behaves like `cv_mdm()`).
I.e., the (only) supported interface for training should accept grids containing
all deconvolutions and alignment parameters at once. If a single model is to be
trained, the grid can just be one row long (as returned by
`get_pgrid('default')`). `benchmark_mdm` should be renamed to `benchmark` and
accept any model fitting function as input.

Side-note: We should probably rename `get_pgrid` to `get_mdm_pgrid` or
`get_dap_grid` or `get_dap` meaning get deconvolution+alignment-parameter-grid.
This prevents confusion with the deconvolution-parameter-grid used by
`grid_deconvolute_spectra` / `grid_deconvolute_spectrum`.

Regarding refactoring `fit_mdm()` / `cv_mdm() `/ `benchmark_mdm()`: the end goal
is the following design:

## `fit_mdm(x, y, dap)`

1. Takes as required arguments:
   - the input spectra x (which can already by deconvoluted and/or aligned)
   - the class labels y as factor and
   - the grid of deconvolutioin+alignment params `dap`
   - the remaining args like sfr, use_rust, nworkers, verbosity, seed,
     nfolds, are optional with defaults

2. Calls `x <- grid_deconvolute_spectra(x, dap)` if `dap` contain any `npmax`
   values greater than 0. This means we need to make the
   deconvolution-param-grid (maybe called `dep`?) the second argument of
   `deconvolute_spectra()` and `deconvolute_spectrum()`. It should be a subset
   of the columns of `dap` so it should be possible to just pass `dap` directly
   and `deconvolute_spectr*()` just takes the columns it needs (nfit, smws,
   smit, delta)

3. Then, using the "enriched" spectra (containing their `dep`-grid incl.
   residual-to-spectra area ratios) we can start the actual grid-search over
   `dap`. For that we iterate over each row of `dap`, deconvolute and align the
   spectra. Extract the signal-integral-matrix (`si_mat()`) as train a model on
   that matrix using `cv_glmnet`. In addition to the accuracy and AUC of each
   cv.glmnet model we need to store the best performing model plus its
   peak-positions to return at the end.

## `benchmark(x, y, k = 3, fun = fit_mdm, ...)`

1. Calls `x <- grid_deconvolute_spectra(x, dap)` if `fun == fit_mdm` and
   `dots$dap` contains any `npmax` values greater than 0 (`dots <- list(...)`).
2. Splits `x` and `y` into `k` folds, i.e., k pairs of `xte`, `yte` and their
   `xtr`, `ytr` counterparts.
3. Calls `fun(xtr, ytr, ...)` for each `xtr`/`ytr` pair to obtain a
   corresponding `model`
4. Calls `predict(model, yte)` to obtain corresponding predictions
5. Collects all predictions, calculates AUC and acc.
6. Prints a performance summary
7. Return predictions and performances.
8. In addition to `fit_mdm`, fun should also accept "binning models", as
   sketched in `misc/stash/bm.R`


Any objections or suggestions from your side before starting the implementation?


-------------------------------------------------------------------------------


> fun == fit_mdm special-case in benchmark() — fragile (breaks if user wraps
> fit_mdm, e.g. with purrr::partial or a closure that fixes dap). Cleaner
> alternatives:
>
> (preferred) Drop the branch entirely. Let users call x <-
>     grid_deconvolute_spectra(x, sfr=...) themselves before benchmark() when
>     they want grid-search caching. Keeps benchmark() model-agnostic.
>
> Or add an explicit preprocess callback argument. predict(model, yte) in your
>     spec is a typo — should be xte. Just confirming.
Mhm, I can see your point, but then we should do the same thing in
`deconvolute()`, i.e., if users want to call something like
`deconvolute(x,sfr=c(10,-2), npmax=1000)`, they need to call
`x<-grid_deconvolute_spectra(x,sfr=c(10,-2))` first.

This introduces two problems:

1. We need to provide a public interface for `grid_deconvolute_spectra()` (which
   should probably be called `grid_deconvolute()`)

2. We need to make sure that the optional deconvolution parameter (`use_rust`,
   `sfr`, `igrs`, etc.) match between the calls to `grid_deconvolute` and
   `deconvolute()`. E.g. in

   ```R
   x <- grid_deconvolute(x, sfr=c(10,-2), use_rust=FALSE)
   x <- deconvolute(x, sfr=c(11,-3), use_rust=TRUE, npmax=100)
   ```

   The grid attac hed by `grid_deconvolute` is not valid.

   Right now this cannot happen, because we attach the grid as part of
   `deconvolute()`, so it always fits.

That being said, both problems above could be solved. But it might be easier to
just keep the call to `grid_deconvolute_spectra()` as part of `benchmark()`.
Instead of accepting `fun = fit_mdm`  and/or `fun = fit_bm` we could be more
specific and accept the functions as strings, e.g.
`benchmark("fit_mdm", x, y, k=3, dap=dagrid("default"))`
or
`benchmark("fit_bm", x, y, k=3)`
Then we could easily do something
`if (fun=="fit_mdm" && max(dap$npmax)>0) x <- grid_deconvolute_spectra(x, dap)`

> dap as 2nd arg of deconvolute_spectr*() — what's the semantics if dap has more
> than one row? My assumption: deconvolute_spectr*() only ever consumes a 1-row
> dap (scalar params). The multi-row iteration stays in fit_mdm. If so, fine —
> but I'd suggest the parameter name be dap = NULL (default keeps current scalar
> API) rather than positional, so deconvolute_spectra(x, sfr=...) still works
> without a grid.
Sorry, that was explained wrong in my instructions. I meant we need to make the
`dep` (Deconvolution Parameter Grid) the second argument of
`grid_deconvolute_spectr*`. So we can easily pass the `dap`
(Deconvolution+Alignment Grid) to `grid_deconvolute_spectra()`.

Now that I've written it a few times. Should we use the terms `deg` and `dag`
instead of `dep` and `dap`?

> Re-deconvolution across dap rows when npmax == 0 — when the grid varies only
> maxShift/maxCombine, decon is identical across rows. Worth deduplicating
> decon-tuples (nfit,smit,smws,delta) and caching results across the inner loop.
> Currently relies on the somewhat-implicit metabodecon.ds.cache option; a small
> explicit cache keyed by tuple would be cleaner.
Mhm, good point. The caching in `deconvolute_spectra()` and `align_decons()` was
added for exactly this reason. To speedup multiple calls to
`deconvolute_spectra()` / `align_decons()` where the inputs haven't changed. But
I can see that it would simplify the caching logic if it wasn't split across
three files. If it doesn't make `cv_mdm` a lot more complicated, I would be
happy to do it. That would mean we need to keep track of the last deconvolution
and alignment in each iteration, right?

> benchmark(... k=3 ...) — current default is nfo=5. Suggest keeping k=5 (3
> outer folds is statistically noisy).
Ok, let's keep 5.

> Name benchmark — generic verb; collides with bench::benchmark / common usage.
> Consider nested_cv() or benchmark_model(). Not blocking, just a heads-up.
Let's keep it. Users can call `metabodecon::benchmark()` and
`bench::benchmark()` to distinguish if they want to use both.

> Outer-fold parallelism story — benchmark_bm parallelizes outer folds;
> benchmark_mdm does inner parallelism only. The unified benchmark() should pick
> a consistent default. I'd default to sequential outer + delegate parallelism
> to fun(...) (cleaner logging, simpler), with no outer parallelism for now. Add
> later if needed.
Yes, exactly. Parallelism should one be done in the inner folds. This keeps
things much easier.


-------------------------------------------------------------------------------


> Final plan (ordered)

> 1. Rename get_pgrid → get_dag. Add maxShift, maxCombine to its grids properly.
>    Sort by decon-tuple.
Should we maybe even name it `mog` or `mfg` (Model Fitting Grid)? Reasoning:

- Deconvolution Grid (deg): nfit, smit, smws, delta, npmax
- Deconvolution+Alignment Grid: adds maxShift
- Model-Fitting-Grid: adds maxCombine

We must sort by deconvolution params, followed by maxShift, followed by
maxCombine. This way we can cache the last deconvolution AND the last alignment.
If only maxCombine changes, we don't need a new deconvolution nor do we need a
new alignment.

> 2. Refactor grid_deconvolute_spectrum/spectra to take deg as 2nd positional
>    arg (was hardcoded grid). deg may be dag — extract unique(deg[deg$npmax>0,
>    c("nfit","smit","smws","delta")]) internally. Keep current default if deg
>    missing.
Good.

> 3. Drop metabodecon.ds.cache, metabodecon.ad.cache, full=FALSE, hash attribute
>    plumbing in deconvolute_spectra / align_decons. Caching moves into fit_mdm.
full=FALSE is still useful to skip costly calculations and to not store costly
members. If we cache the last deconvolutuion results plus the last alignment
results, we will have times where the spectra and two deconvolutions and two
alignments are in RAM. Then distributed among processes (which on Windows have
their own RAM requirements). Theirfore, we should try to keep the RAM used by
the deconvoluted and/or aligned spectra as low as possible and not spend time on
the calculation of fields that are never used. THis is not a high priority, but
we shouldn't remove the full argument.

> 4. New fit_mdm(x, y, dag = get_dag("default"), sfr=NULL, use_rust=0.5,
>    nworkers=1, verbosity=2, seed=1, nfolds=10, check=TRUE) — replaces both old
>    fit_mdm and cv_mdm. Logic per the sketch above. Returns mdm with $model,
>    $ref, $meta, $dag (with acc/auc columns).
Very good. Also make sure that acc and auc don't need to be pre-initalized when
we pass the dag argument to the function. That would be confusing in the docs.

> 5. Delete cv_mdm.
Good.

> 6. New benchmark(fun, x, y, k=5, seed=1, verbosity=2, ...) — fun is a string.
>    If fun == "fit_mdm" and dots$dag has any npmax>0, run
>    grid_deconvolute_spectra(x, dots$dag, sfr=dots$sfr) once up front. Split,
>    fit, predict, aggregate, return list(models, predictions, performance).
>    Replaces benchmark_mdm.
Good.

> 7. Update predict.mdm to read params from meta after the new fit_mdm shape
>    (mostly already correct).
Good. Maybe we can find a better name then `meta`?

> 8. Add fit_bm / predict.bm (move from bm.R into mdm.R or new R/bm.R) so
>    benchmark("fit_bm", ...) works.
I think we can move it into the `mdm.R` file. It's just two functions and we
only have the "binning models" (bm) so we can compare our "metabodecon models"
(mdm) against them. So it's justifyable to have them in the same file.

> 9. Fix the small bugs noticed: pc → row vars in old cv_mdm body becomes moot;
>    get_pgrid ord referencing missing cols is fixed by adding the cols.
Good.

> 10. Update doc/aliases (@rdname mdm, mdm_methods, NAMESPACE export changes),
>     grep for callers of cv_mdm/benchmark_mdm/get_pgrid and update.
Good.

> A couple of decisions to confirm before I touch the code:

> (a) OK to delete the RAM-cache options + full/hash plumbing in
>     deconvolute_spectra/align_decons? (Replaced by in-loop cache in fit_mdm.)
Yes, that's ok.

Just to verify: i.e., we will only have explicit, non-options-based caching,
right?

1. `x <- grid_deconvolute_spectra(x)` attaches the `deg` grids to its input
   spectra `x`, so that subsequent calls to `deconvolute_spectra(x)` can access
   them.

2. `fit_mdm(x)` keeps track of deconvolution and/or alignment of the last
   iteration and only re-deconvolutes/re-alignes if any relevant input has
   changed compared to the last iteration.

> (b) Should fit_bm live in mdm.R or its own R/bm.R?
In `mdm.R`.

> (c) benchmark() signature — benchmark(fun, x, y, ...) or benchmark(x, y,
>     fun="fit_mdm", ...) (matches your benchmark("fit_mdm", x, y, ...) example,
>     which has fun first)?
I think I like `benchmark(x, y, ..., fun="fit_mdm")` the most. `fit_mdm` is what
the vast majority of users will use, so they directly write
`benchmark(x, y, dag)`. What's your recommendation?


-------------------------------------------------------------------------------


> Grid name → mfg (Model Fitting Grid): Yes, that's a cleaner taxonomy. So:
> deg (Deconvolution Grid): nfit, smit, smws, delta, npmax
> dag (Decon+Align Grid): deg + maxShift
> mfg (Model Fitting Grid): dag + maxCombine
I actually like `mog` more. It's pronouncable as "mog" where "mfg" must be
pronounced like the german "em-ef-tschi".

> grid_deconvolute_spectr*() takes deg as 2nd arg (accepts dag/mfg too — extra
> cols ignored). fit_mdm() takes mfg. Sort mfg rows by (npmax, nfit, smit, smws,
> delta, maxShift, maxCombine) so caching works as you described.
> Constructor: get_mfg(conf="default").
Good, but let's use `mog`.

> full=FALSE stays: Agreed, keep it. We just remove the option-based caching
> layer (metabodecon.ds.cache/metabodecon.ad.cache) and the hash attribute
> plumbing. full is orthogonal and useful.
Good.

> Caching confirmation — yes, exactly two layers:
> - grid_deconvolute_spectra(x) attaches $grid (the per-spectrum deg performance
>   table) to each spectrum. Idempotent. Persists on x.
> - fit_mdm(x, y, mfg) keeps two locals: last_decons (keyed by decon-tuple) and
>   last_aligns (keyed by decon-tuple + maxShift). Recomputes only when the key
>   changes. Sorted mfg ensures runs of identical keys cluster.
Good. But I'm wondering, do we even need seperate variables for that? Can't we just
have something like:

```R
if (decon_pars_have_changed) d <- deconvolute_spectra(x, ...)
if (align_pars_have_changed) a <- align(d, ...)
```

This would keep RAM usage as low as possible.

> acc/auc not pre-initialized in mfg: Confirmed — get_mfg() returns only the
> input columns. fit_mdm() adds acc/auc to the returned $mfg.
Good.

> meta rename: params is the most obvious replacement. config would also work.
> I'd go with params — short, accurate (it stores the chosen scalar parameters
> of the best model). mdm$params reads well in code. OK?
Ok.

> benchmark() signature — benchmark(x, y, ..., fun="fit_mdm"): Agreed, this is
> better. Matches typical-case ergonomics (benchmark(x, y, mfg)) and fun= stays
> explicit when you need it. My recommendation matches yours.
Good.

> fit_bm arg name: It currently uses regions/binwidth (no grid). For consistency
> with the new benchmark() signature, I'll keep it as is — benchmark(x, y,
> regions=..., binwidth=..., fun="fit_bm").
Change to takes `igrs` (Ignore Regions) and `nbin` (Number of bins). This is
closer to the `deconvolute` signature, so we can easily reuse the same ignore
regions.

> Open clarifications before I start

> 1. deconvolute_spectra(x, deg=NULL, ...) 2nd-arg shift — wait, on rereading
>    you said deg should be 2nd arg of grid_deconvolute_spectrum/spectra, not
>    deconvolute_spectra. Confirming: deconvolute_spectra keeps its current
>    signature; only grid_deconvolute_spectr* gets deg as 2nd positional.
Exactly.
`deconvolute_spectra()` keeps it signature.
`grid_deconvolute_spectra()` and `grid_deconvolute_spectrum()` get `deg` as
second argument.

> 2. benchmark() performance reporting — I'll print acc and AUC overall (pooled
>    across folds). Per-fold also? Suggest: pooled in summary line, per-fold in
>    returned performance data.frame.
Also print per fold. This way users can control the verbosity.
If verbosity is 2, they see the individual `cv_mdm` logs.
If verbosity is 1, they see only the acc/auc after each fold and the final acc/auc.
If verbosity is 0, they see nothing.

> 3. Old exports to remove from NAMESPACE/Rd: cv_mdm, benchmark_mdm, get_pgrid.
>    New exports: fit_mdm (signature change), benchmark, get_mfg, fit_bm. Predict
>    methods predict.bm need S3 registration. OK?
That's ok, but use `mog` instead of `mfg`.

