#include <R.h>       /* R_xlen_t, R_NilValue, etc.  */
#include <Rinternals.h> /* SEXP, REAL, XLENGTH, ...    */
#include <string.h>  /* memset, memcpy               */

/*
 * lorentz_sup_c  —  Superposition of Lorentz curves
 * ===================================================
 *
 * Evaluates the sum
 *
 *   pr[i] = sum_{j=0}^{np-1}  Al[j] / (l2[j] + (x[i] - x0[j])^2)
 *
 * for every position x[i], where Al[j] = |A_j * lambda_j| and
 * l2[j] = lambda_j^2 are pre-computed by the caller to avoid redundant
 * multiplications inside the hot loop.
 *
 * Parameters
 * ----------
 * x   REALSXP of length nx  — evaluation positions (ppm grid)
 * x0  REALSXP of length np  — peak centres
 * Al  REALSXP of length np  — |A * lambda| per peak
 * l2  REALSXP of length np  — lambda^2 per peak
 *
 * Returns
 * -------
 * A freshly allocated REALSXP of length nx.
 *
 *
 * R API primer
 * ------------
 * SEXP      The universal R object type.  Every R value — vector, list,
 *           function, environment — is a pointer to an SEXP struct.
 *
 * REAL(v)   Returns a `double *` pointing to the data array of a numeric
 *           SEXP.  No copy is made; the pointer is valid as long as `v`
 *           is protected from garbage collection (see Rf_protect below).
 *
 * XLENGTH(v)
 *           Returns the number of elements in SEXP `v` as an R_xlen_t
 *           (a signed integer type, at least 64 bits on long-vector
 *           builds).  Equivalent to R's `length(v)`.
 *
 * Rf_allocVector(type, n)
 *           Allocates a new R vector of the given type (REALSXP = double)
 *           and length n.  The contents are uninitialised.  The returned
 *           SEXP is on the R heap and will be collected by the garbage
 *           collector unless protected.
 *
 * Rf_protect(v) / Rf_unprotect(n)
 *           R's garbage collector can run at any allocation point.  Every
 *           SEXP that was just allocated and is not yet reachable from R's
 *           own symbol table must be pushed onto the "protection stack" via
 *           Rf_protect().  Rf_unprotect(n) pops the n most recently pushed
 *           entries.  Forgetting to protect leads to use-after-free bugs;
 *           forgetting to unprotect causes a stack overflow on long runs.
 *           Here we protect exactly one object (res) and unprotect it
 *           before returning.
 *
 * memset(ptr, byte, n)
 *           C standard library function.  Fills `n` bytes starting at
 *           `ptr` with the value `byte` (interpreted as unsigned char).
 *           `memset(tmp, 0, blen * sizeof(double))` zero-initialises the
 *           accumulator block.  The bit pattern of all-zeros is defined
 *           by IEEE 754 to be +0.0 for doubles.
 *
 * memcpy(dst, src, n)
 *           Copies exactly `n` bytes from `src` to `dst`.  The regions
 *           must not overlap (use memmove otherwise).  Used here to
 *           flush the stack-allocated `tmp` block back to the heap
 *           output array after each block of positions is finished.
 *
 *
 * Optimisation techniques applied
 * --------------------------------
 *
 * 1. -O3 compilation (see src/Makevars)
 *    R packages compile C with -O2 by default.  -O3 additionally enables
 *    auto-vectorisation (-ftree-vectorize), loop unrolling, and aggressive
 *    instruction scheduling.  With -O3 the compiler has more freedom to
 *    issue multiple FP operations per clock on modern out-of-order cores.
 *
 * 2. `const` and `restrict` on every pointer
 *
 *    `const double *p`  declares that the data pointed to by `p` will
 *    not be modified through `p`.  This allows the compiler to treat the
 *    values as immutable for the lifetime of the function, enabling
 *    constant-folding and load-elimination.
 *
 *    `restrict` (C99) is a promise to the compiler that no two restrict-
 *    qualified pointers in the same scope point to overlapping memory.
 *    Without it the compiler must conservatively assume that writing
 *    `pr[i]` could overwrite `px0[j]`, `pAl[j]`, or `pl2[j]` — a false
 *    aliasing hazard that prevents register promotion and vectorisation.
 *    With `restrict` on all five pointers, the compiler knows that a
 *    store to `prb[i]` (or `tmp[i]`) cannot invalidate the cached values
 *    of x0j / alj / l2j, so it keeps them in scalar FP registers for
 *    the entire inner loop.
 *
 * 3. Loop order swap (outer = peaks, inner = positions)
 *
 *    Original order:  for i (positions)  for j (peaks)
 *    New order:       for i0 (blocks)  for j (peaks)  for i (block)
 *
 *    In the original order, for each position i the compiler must load
 *    px0[j], pAl[j], pl2[j] from memory on every j step.  Even with
 *    restrict, if np is large these arrays may not stay in L1 cache.
 *    More importantly, pr[i] is written once per i-iteration but read
 *    back as an accumulator only after iterating all j — the compiler
 *    cannot safely use a register accumulator because it would have to
 *    prove no aliasing over the entire i-loop, which is hard for the
 *    loop-carried dependency.
 *
 *    With peaks as the outer loop, the three scalars x0j / alj / l2j
 *    are loop-invariant for the entire inner (position) loop and are
 *    trivially promoted to registers.  The inner loop then reduces to
 *    a pure sequential read of pxb[i] and read-modify-write of tmp[i],
 *    both of which are friendly to vectorisation.
 *
 * 4. Cache blocking over positions (block size BLOCK)
 *
 *    nx can be ~131072 (≈1 MB for double), far exceeding L1 (≈32 KB).
 *    In the swapped-loop order, for each peak j the entire pr[] array
 *    would be streamed through cache once.  With np peaks, pr[] is
 *    written np times sequentially, but on the first pass after a gap
 *    each cache line must be fetched from L2/L3.
 *
 *    Cache blocking limits the working set:
 *      - `tmp[0..BLOCK]`       512 bytes (BLOCK=64 doubles)
 *      - `pxb[0..BLOCK]`       512 bytes (alias of px + i0)
 *
 *    Both fit in L1 alongside the three peak arrays (np doubles each;
 *    for np=1000 that is 8 KB per array — comfortably in L1+L2).
 *    After processing all np peaks for the current block, tmp[] is
 *    flushed to pr[] via memcpy and we advance to the next block.
 *    This means pr[] is written only once per block, as a sequential copy,
 *    rather than being touched once per peak-iteration.
 *
 *    BLOCK = 64 was chosen so that two 64-double arrays (tmp + pxb)
 *    occupy 1 KB — about 16 cache lines — leaving the rest of a typical
 *    32 KB L1 for the peak arrays and compiler spills.
 *
 * Author: 2026 Tobias Schmidt.
 */

#define BLOCK 64  /* number of positions processed per cache block */

SEXP lorentz_sup_c(SEXP x, SEXP x0, SEXP Al, SEXP l2) {

    R_xlen_t nx = XLENGTH(x), np = XLENGTH(x0);
    if (XLENGTH(Al) != np || XLENGTH(l2) != np) {
        Rf_error("x0, Al, and l2 must have the same length.");
    }

    /* restrict + const: compiler may keep pAl[j]/pl2[j]/px0[j] in
     * registers and auto-vectorise the inner loop (see comment above). */
    const double * restrict px  = REAL(x);
    const double * restrict px0 = REAL(x0);
    const double * restrict pAl = REAL(Al);
    const double * restrict pl2 = REAL(l2);

    /* Allocate output on the R heap and protect it from GC. */
    SEXP res = Rf_protect(Rf_allocVector(REALSXP, nx));
    double * restrict pr = REAL(res);

    /* Stack-allocated accumulator block.  Lives in L1 for the whole
     * j-sweep of each position block (see cache-blocking comment). */
    double tmp[BLOCK];

    /* Outer loop: advance through positions in blocks of BLOCK. */
    for (R_xlen_t i0 = 0; i0 < nx; i0 += BLOCK) {
        R_xlen_t blen = (i0 + BLOCK <= nx) ? BLOCK : (nx - i0);
        const double * restrict pxb = px + i0; /* start of this block */
        double       * restrict prb = pr + i0;

        /* Zero the accumulator for this block.
         * memset with byte=0 is safe for IEEE 754 doubles (+0.0). */
        memset(tmp, 0, (size_t)blen * sizeof(double));

        /* Middle loop: iterate over all peaks.  The three scalars
         * x0j / alj / l2j are hoisted into registers here because
         * they are loop-invariant for the inner i-loop and restrict
         * guarantees no aliasing. */
        for (R_xlen_t j = 0; j < np; j++) {
            double x0j = px0[j], alj = pAl[j], l2j = pl2[j];

            /* Inner loop: accumulate contribution of peak j into tmp.
             * Sequential read of pxb[i] + read-modify-write of tmp[i];
             * both arrays fit in L1 for the lifetime of this i0 block.
             * With -O3 + restrict the compiler can often vectorise this
             * loop. */
            for (R_xlen_t i = 0; i < blen; i++) {
                double d = pxb[i] - x0j;
                tmp[i] += alj / (l2j + d * d);
            }
        }

        /* Flush the finished block from the stack to the R heap output.
         * memcpy gives the compiler and C library a simple sequential copy. */
        memcpy(prb, tmp, (size_t)blen * sizeof(double));
    }

    Rf_unprotect(1);
    return res;
}
