1. Is it possible to store the profiling results from `profile_deconvolute.R` in
   a human readable format? E.g. as json, markdown or even plain text?

2. According to the results from `profile_deconvolute.R`, deconvolution of the
   first AKI spectrum takes oin average 7460 milliseconds, with 4110
   milliseconds spent on `lorentz_sup` calls executed from
   `deconvolute_spectrum_r()` and 2120 ms spent on lorentz_sup calls executed
   from `deconvolute_spectrum_r()`/`fit_lorentz_curves2()`. Please check whether there is
   1. Any way how we can make these `lorentz_sup()` calls faster (e.g. by improving the C code)
   2. Whether we really need the `x$sit$sup <- loretz_sup(...)` assignment from `deconvolute_spectrum_r()`. As far as I can see the `$sup` element is mostly used during the alignment of spectra as part of `align`. But maybe the alignment algorithm doesn't really need the actual superposition or can be rewritten in a way that I must not be calculated explicitely? Please verify this.

Don't implement anything. Just analyze and tell me your findings.


-------------------------------------------------------------------------------

> Option A — Skip profvis, use base R Rprof() directly. The output is a
> plain-text file, and summaryRprof() returns readable tables of self/total time
> per call:

Good. Please rewrite `profile_deconvolute.R` accordingly.

> 1. OpenMP on the outer loop. The work per i is independent and writes only
> pr[i]. Adding #pragma omp parallel for (with the standard R-package OpenMP dance
> in Makevars / Makevars.win) typically gives near-linear speed-up on the outer
> loop. On a 4-core machine this alone could turn the 4110 ms full-sup call into
> ~1.1–1.5 s.

Multiprocessing should be done over spectra from R, not within a single spectrum
in C.

> Distance-based pruning. Lorentzians decay as 1/d². If we keep x0, Al, l2 sorted
> by x0 and pre-compute a per-peak cutoff dmax[j] such that Al[j] / (l2[j] +
> dmax[j]^2) < eps * scale, then for each xi we can binary-search the relevant
> peak window and only iterate that. With ~hundreds of peaks spread over many ppm
> and individual half-widths in the 10⁻³ ppm range, the effective inner-loop
> length per xi drops from np to a small constant. This is the single biggest
> expected speed-up. Note: it changes results slightly unless eps is chosen
> tighter than current numerical noise; using e.g. eps = 1e-12 * max(Al/l2) is
> safe.

Sounds too complicated. Calculating relevant peaks and maintaining the
additional logic is not worth the gain.

> SIMD-friendliness. The inner loop is already vectorizable; making sure it's
> compiled with -O3 (and possibly -funsafe-math-optimizations if you accept
> reordering of FP adds) helps. Adding __restrict__ qualifiers on the pointers and
> a small static inline helper is a low-risk change. Auto-vectorization is more
> reliable if the inner loop body is written with a contiguous accumulator and
> constant strides — it already is.

This sounds very good. What exactly does __restrict__ do and are there other
ways we can optimize the compiled code?

> Block over peaks for cache reuse. With np small the three peak arrays fit in
> L1, so this matters less. Skip unless 1+2 don't suffice.

Sounds also interesting. How exactly would you implement this. Explain in detail.

> Algorithmic alternative for the fit-loop call. Inside fit_lorentz_curves2()
> (decon.R:429-431), lorentz_sup is called on rs_x of length 3*np, but each
> stencil point only "really" needs the contribution from peaks near it (same
> argument as #2). With pruning, the fit-loop call's cost essentially becomes
> O(np) instead of O(np²) — that's where the 2120 ms goes.

Again, not worth the effort.