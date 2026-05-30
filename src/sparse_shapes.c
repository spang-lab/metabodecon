#include <R.h>
#include <Rinternals.h>
#include <string.h>

/*
 * sparse_shapes.c — Fast peak-shape superpositions on the cssh grid
 * ==================================================================
 *
 * Two narrow-support alternatives to `lorentz_sup_c` used as supsh
 * input for CluPA's FFT cross-correlation:
 *
 *   triangle_sup_c   isoceles triangle, full-width = lambda (in
 *                    datapoints), peak height = A. value at offset
 *                    d from the center is A * max(0, 1 - |d| / hw)
 *                    where hw = lambda / 2 (in datapoints).
 *
 *   rect_sup_c       constant A on [c - hw, c + hw], zero outside.
 *
 * Both run in O(npeaks * hw) — independent of the cssh length — and
 * are >10x faster than the Lorentzian on a 128k-point grid. The
 * tradeoff is sharper localization (no heavy tails extending past
 * neighboring peaks), which is exactly what FFT cross-correlation
 * benefits from: a shift of k columns becomes detectable as a clean
 * peak in the correlation function instead of a broad maximum.
 *
 * Centers `c` and half-widths `hw` are 1-based integer indices on the
 * cssh column space (the caller has already rounded x0 to its nearest
 * cssh column and converted lambda from ppm to datapoints). This is
 * the round-first contract: the supsh peak position and the integer
 * peakList index passed into hclust_align() are byte-equal.
 *
 * Author: 2026 Tobias Schmidt.
 */

SEXP triangle_sup_c(SEXP c_, SEXP hw_, SEXP A_, SEXP n_) {
    R_xlen_t np = XLENGTH(c_);
    if (XLENGTH(hw_) != np || XLENGTH(A_) != np) {
        Rf_error("c, hw, and A must have the same length.");
    }
    int n = Rf_asInteger(n_);
    if (n <= 0) Rf_error("n must be > 0.");

    const int    * restrict pc  = INTEGER(c_);
    const int    * restrict phw = INTEGER(hw_);
    const double * restrict pA  = REAL(A_);

    SEXP res = Rf_protect(Rf_allocVector(REALSXP, n));
    double * restrict pr = REAL(res);
    memset(pr, 0, (size_t)n * sizeof(double));

    for (R_xlen_t i = 0; i < np; i++) {
        int    ci  = pc[i];
        int    hwi = phw[i];
        double Ai  = pA[i];
        if (hwi <= 0 || Ai == 0.0) continue;
        int lo = ci - hwi; if (lo < 1) lo = 1;
        int hi = ci + hwi; if (hi > n) hi = n;
        double inv_hw = 1.0 / (double)hwi;
        for (int j = lo; j <= hi; j++) {
            int d = j - ci;
            double ad = (d < 0) ? -(double)d : (double)d;
            pr[j - 1] += Ai * (1.0 - ad * inv_hw);
        }
    }
    Rf_unprotect(1);
    return res;
}

SEXP rect_sup_c(SEXP c_, SEXP hw_, SEXP A_, SEXP n_) {
    R_xlen_t np = XLENGTH(c_);
    if (XLENGTH(hw_) != np || XLENGTH(A_) != np) {
        Rf_error("c, hw, and A must have the same length.");
    }
    int n = Rf_asInteger(n_);
    if (n <= 0) Rf_error("n must be > 0.");

    const int    * restrict pc  = INTEGER(c_);
    const int    * restrict phw = INTEGER(hw_);
    const double * restrict pA  = REAL(A_);

    SEXP res = Rf_protect(Rf_allocVector(REALSXP, n));
    double * restrict pr = REAL(res);
    memset(pr, 0, (size_t)n * sizeof(double));

    for (R_xlen_t i = 0; i < np; i++) {
        int    ci  = pc[i];
        int    hwi = phw[i];
        double Ai  = pA[i];
        if (hwi <= 0 || Ai == 0.0) continue;
        int lo = ci - hwi; if (lo < 1) lo = 1;
        int hi = ci + hwi; if (hi > n) hi = n;
        for (int j = lo; j <= hi; j++) {
            pr[j - 1] += Ai;
        }
    }
    Rf_unprotect(1);
    return res;
}
