#include <R.h>
#include <Rinternals.h>
#include <math.h>

/* Find peaks in a signal via its second derivative.
 *
 * A point i is a candidate for each role when these conditions hold
 * (d_p = d[i-1], d_n = d[i+1]; NA comparisons are always false):
 *
 *   lbc (left-border candidate):
 *     (a) sign change from left: d[i] < 0  AND  d[i-1] >= 0, OR
 *     (b) local max:             d[i] > d[i+1]  AND  d[i] >= d[i-1]
 *
 *   pcc (peak-centre candidate):
 *     local minimum of d that is negative:
 *       d[i] < 0  AND  d[i] <= d[i-1]  AND  d[i] < d[i+1]
 *
 *   rbc (right-border candidate):
 *     (a) sign change to right:  d[i] < 0  AND  d[i+1] >= 0, OR
 *     (b) local max:             d[i-1] < d[i]  AND  d[i+1] <= d[i]
 *
 * State machine (single left-to-right pass):
 *
 *   FREE:
 *     lbc  ->  pending_lb = i,  goto HAVE_LB
 *
 *   HAVE_LB:
 *     lbc  ->  pending_lb = i   (update; still HAVE_LB)
 *     pcc  ->  pending_pc = i,  goto HAVE_PC
 *
 *   HAVE_PC:
 *     rbc  ->  emit (pending_lb, pending_pc, i)
 *              if also lbc: pending_lb = i, goto HAVE_LB
 *              else:                        goto FREE
 *
 *
 * Argument:
 *   y  - signal intensities (REALSXP, length n).
 *
 * Returns a data.frame with integer columns left/center/right (1-based,
 * NA_INTEGER where absent) and double column score.
 *
 * Author: 2026 Tobias Schmidt.
 */

/* Append one peak to the output buffers.
 * Indices lb0, pc0, rb0 are 0-based; rb0 = -1 signals a missing right border.
 * cs is the prefix-sum array of |d| (length n+1). */
static void emit_peak(
    int lb0, int pc0, int rb0, double *cs,
    int *lb_buf, int *pc_buf, int *rb_buf, double *sc_buf, int *npc
) {
    double sum_lj = cs[pc0 + 1] - cs[lb0];
    double sum_jr = (rb0 < 0) ? -1.0 : (cs[rb0 + 1] - cs[pc0]);
    double sc = (sum_jr < 0.0) ? 0.0
              : (sum_lj < sum_jr ? sum_lj : sum_jr);
    lb_buf[*npc] = lb0 + 1;                       /* convert to 1-based */
    pc_buf[*npc] = pc0 + 1;
    rb_buf[*npc] = (rb0 < 0) ? NA_INTEGER : rb0 + 1;
    sc_buf[*npc] = sc;
    (*npc)++;
}

SEXP find_peaks_c(SEXP y_in) {

    int n = Rf_length(y_in);
    double *y = REAL(y_in);

    /* ---- second derivative -------------------------------------------- */
    double *d = (double*)R_alloc(n < 1 ? 1 : n, sizeof(double));
    if (n >= 1) d[0] = NA_REAL;
    for (int i = 1; i < n - 1; i++)
        d[i] = y[i - 1] + y[i + 1] - 2.0 * y[i];
    if (n >= 2) d[n - 1] = NA_REAL;

    /* ---- prefix sums of |d| for score computation --------------------- */
    double *cs = (double*)R_alloc(n + 1, sizeof(double));
    cs[0] = 0.0;
    for (int i = 0; i < n; i++)
        cs[i + 1] = cs[i] + (ISNAN(d[i]) ? 0.0 : fabs(d[i]));

    /* ---- output buffers (at most n/3 peaks possible, n is safe bound) - */
    int    *lb_buf = (int*)   R_alloc(n < 1 ? 1 : n, sizeof(int));
    int    *pc_buf = (int*)   R_alloc(n < 1 ? 1 : n, sizeof(int));
    int    *rb_buf = (int*)   R_alloc(n < 1 ? 1 : n, sizeof(int));
    double *sc_buf = (double*)R_alloc(n < 1 ? 1 : n, sizeof(double));
    int npc = 0;

    /* ---- state machine: single left-to-right pass --------------------- */
    /* States: 0=FREE, 1=HAVE_LB, 2=HAVE_PC */
    int state = 0, pending_lb = -1, pending_pc = -1;

    for (int i = 0; i < n; i++) {

        double di  = d[i];
        double d_p = (i > 0)     ? d[i - 1] : NA_REAL;
        double d_n = (i < n - 1) ? d[i + 1] : NA_REAL;

        int is_lbc =
            (!ISNAN(di) && !ISNAN(d_p) && di < 0.0 && d_p >= 0.0) ||
            (!ISNAN(di) && !ISNAN(d_p) && !ISNAN(d_n)
             && di > d_n && di >= d_p);

        int is_pcc =
            !ISNAN(di) && !ISNAN(d_p) && !ISNAN(d_n) &&
            di < 0.0 && di <= d_p && di < d_n;

        int is_rbc =
            (!ISNAN(di) && !ISNAN(d_n) && di < 0.0 && d_n >= 0.0) ||
            (!ISNAN(di) && !ISNAN(d_p) && !ISNAN(d_n)
             && d_p < di && d_n <= di);

        if (state == 0) {
            if (is_lbc) { pending_lb = i; state = 1; }

        } else if (state == 1) {
            if (is_pcc)      { pending_pc = i; state = 2; }
            else if (is_lbc) { pending_lb = i; }

        } else {
            if (is_rbc) {
                emit_peak(pending_lb, pending_pc, i,
                          cs, lb_buf, pc_buf, rb_buf, sc_buf, &npc);
                if (is_lbc) { pending_lb = i; state = 1; }
                else        { state = 0; }
            }
        }
    }

    /* ---- build data.frame --------------------------------------------- */
    SEXP left_s   = Rf_protect(Rf_allocVector(INTSXP,  npc));
    SEXP center_s = Rf_protect(Rf_allocVector(INTSXP,  npc));
    SEXP right_s  = Rf_protect(Rf_allocVector(INTSXP,  npc));
    SEXP score_s  = Rf_protect(Rf_allocVector(REALSXP, npc));

    for (int k = 0; k < npc; k++) {
        INTEGER(left_s)[k]   = lb_buf[k];
        INTEGER(center_s)[k] = pc_buf[k];
        INTEGER(right_s)[k]  = rb_buf[k];
        REAL(score_s)[k]     = sc_buf[k];
    }

    SEXP df = Rf_protect(Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(df, 0, left_s);
    SET_VECTOR_ELT(df, 1, center_s);
    SET_VECTOR_ELT(df, 2, right_s);
    SET_VECTOR_ELT(df, 3, score_s);

    SEXP nms = Rf_protect(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, Rf_mkChar("left"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("center"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("right"));
    SET_STRING_ELT(nms, 3, Rf_mkChar("score"));
    Rf_setAttrib(df, R_NamesSymbol, nms);

    /* Compact row-name representation used by R's data.frame: c(NA, -npc) */
    SEXP rn = Rf_protect(Rf_allocVector(INTSXP, 2));
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -npc;
    Rf_setAttrib(df, R_RowNamesSymbol, rn);

    SEXP cls = Rf_protect(Rf_mkString("data.frame"));
    Rf_setAttrib(df, R_ClassSymbol, cls);

    Rf_unprotect(8); /* left, center, right, score, df, nms, rn, cls */
    return df;
}
