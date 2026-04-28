#include <R.h>
#include <Rinternals.h>
#include <string.h>

/* Superposition of Lorentz curves using pre-computed Al = |A*lambda| and
 * l2 = lambda^2. This avoids recomputing those inside the refinement loop.
 *
 * Arguments (all numeric vectors):
 *   x   - positions at which to evaluate (length nx)
 *   x0  - peak centres (length np)
 *   Al  - |A * lambda| per peak (length np)
 *   l2  - lambda^2 per peak (length np)
 *
 * Returns a REALSXP vector of length nx.
 *
 * Author: 2026 Tobias Schmidt.
 */
SEXP lorentz_sup_c(SEXP x, SEXP x0, SEXP Al, SEXP l2) {
    int nx = Rf_length(x), np = Rf_length(x0);
    double *px  = REAL(x),  *px0 = REAL(x0);
    double *pAl = REAL(Al), *pl2 = REAL(l2);
    SEXP res = Rf_protect(Rf_allocVector(REALSXP, nx));
    double *pr = REAL(res);
    memset(pr, 0, nx * sizeof(double));
    for (int i = 0; i < nx; i++) {
        double s = 0.0, xi = px[i];
        for (int j = 0; j < np; j++) {
            double d = xi - px0[j];
            s += pAl[j] / (pl2[j] + d * d);
        }
        pr[i] = s;
    }
    Rf_unprotect(1);
    return res;
}
