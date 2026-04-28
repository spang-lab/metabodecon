#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP lorentz_sup_c(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallMethods[] = {
    {"lorentz_sup_c", (DL_FUNC) &lorentz_sup_c, 4},
    {NULL, NULL, 0}
};

void R_init_metabodecon(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
