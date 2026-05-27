#include <R.h>
#include <Rinternals.h>
#include <stddef.h>

/*
 * Needleman-Wunsch pairwise alignment of NMR peak lists.
 *
 * Inputs (all R-side preallocated):
 *   M  -- nx x ny REAL matrix, column-major; M[i,j] = match cost(p[i], q[j])
 *   gp -- nx REAL vector;     gp[i] = cost of leaving p[i] unpaired
 *   gq -- ny REAL vector;     gq[j] = cost of leaving q[j] unpaired
 *
 * Returns: list(cost = numeric(1), alignment = integer matrix [L x 2]) where
 * each alignment row is one column of the alignment (i, j) with NA_INTEGER
 * marking a gap on that side.
 */
SEXP align_dp_c(SEXP M_, SEXP gp_, SEXP gq_) {
    if (!isReal(M_)) error("M must be a numeric matrix");
    if (!isReal(gp_) || !isReal(gq_)) error("gap costs must be numeric");
    int nx = LENGTH(gp_);
    int ny = LENGTH(gq_);
    SEXP dims = getAttrib(M_, R_DimSymbol);
    if (dims == R_NilValue || LENGTH(dims) != 2) error("M must be a 2D matrix");
    if (INTEGER(dims)[0] != nx) error("nrow(M) must equal length(gp)");
    if (INTEGER(dims)[1] != ny) error("ncol(M) must equal length(gq)");

    double *M  = REAL(M_);
    double *gp = REAL(gp_);
    double *gq = REAL(gq_);

    size_t n1 = (size_t)nx + 1;
    size_t m1 = (size_t)ny + 1;
    double *dp = (double*)R_alloc(n1 * m1, sizeof(double));
    char   *bt = (char*)  R_alloc(n1 * m1, sizeof(char));

    /* Column-major layout: cell (i, j) at offset i + n1 * j. */
    #define DP(i, j) dp[(size_t)(i) + n1 * (size_t)(j)]
    #define BT(i, j) bt[(size_t)(i) + n1 * (size_t)(j)]

    /* dp[0,0]=0; first column = all-gap-in-q, first row = all-gap-in-p. */
    DP(0, 0) = 0.0; BT(0, 0) = 0;
    for (int i = 1; i <= nx; i++) {
        DP(i, 0) = DP(i-1, 0) + gp[i-1];
        BT(i, 0) = 2; /* gap in q */
    }
    for (int j = 1; j <= ny; j++) {
        DP(0, j) = DP(0, j-1) + gq[j-1];
        BT(0, j) = 3; /* gap in p */
    }

    /* Forward fill, column-major to keep M accesses contiguous. */
    for (int j = 1; j <= ny; j++) {
        size_t Mcol = (size_t)(j-1) * (size_t)nx;
        for (int i = 1; i <= nx; i++) {
            double cm = DP(i-1, j-1) + M[(size_t)(i-1) + Mcol];
            double cg = DP(i-1, j  ) + gp[i-1];
            double cp = DP(i,   j-1) + gq[j-1];
            double best = cm; char bm = 1;
            if (cg < best) { best = cg; bm = 2; }
            if (cp < best) { best = cp; bm = 3; }
            DP(i, j) = best;
            BT(i, j) = bm;
        }
    }

    /* Backtrack into temp buffers (reverse order). */
    int max_len = nx + ny;
    int *ti = max_len ? (int*)R_alloc((size_t)max_len, sizeof(int)) : NULL;
    int *tj = max_len ? (int*)R_alloc((size_t)max_len, sizeof(int)) : NULL;
    int len = 0;
    {
        int i = nx, j = ny;
        while (i > 0 || j > 0) {
            char m = BT(i, j);
            if (m == 1)       { ti[len] = i;           tj[len] = j;           i--; j--; }
            else if (m == 2)  { ti[len] = i;           tj[len] = NA_INTEGER;  i--;      }
            else              { ti[len] = NA_INTEGER;  tj[len] = j;           j--;      }
            len++;
        }
    }

    /* Build return: list(cost, alignment [L x 2 integer matrix]). */
    SEXP cost_     = PROTECT(ScalarReal(DP(nx, ny)));
    SEXP align_mat = PROTECT(allocMatrix(INTSXP, len, 2));
    int *am = INTEGER(align_mat);
    for (int k = 0; k < len; k++) {
        am[k]       = ti[len - 1 - k]; /* column "i" */
        am[k + len] = tj[len - 1 - k]; /* column "j" */
    }
    SEXP colnames = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(colnames, 0, mkChar("i"));
    SET_STRING_ELT(colnames, 1, mkChar("j"));
    SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, R_NilValue);
    SET_VECTOR_ELT(dimnames, 1, colnames);
    setAttrib(align_mat, R_DimNamesSymbol, dimnames);

    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, cost_);
    SET_VECTOR_ELT(out, 1, align_mat);
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("cost"));
    SET_STRING_ELT(names, 1, mkChar("alignment"));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(6);
    return out;

    #undef DP
    #undef BT
}
