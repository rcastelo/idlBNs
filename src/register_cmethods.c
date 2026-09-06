#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP
C_iBIC_node_score(SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP
C_iBIC_score(SEXP, SEXP, SEXP, SEXP, SEXP);

/* registration of C-entry points */

static const R_CallMethodDef CallEntries[] = {
    {"C_iBIC_node_score", (DL_FUNC) &C_iBIC_node_score, 5},
    {"C_iBIC_score",      (DL_FUNC) &C_iBIC_score,      5},
    {NULL, NULL, 0}
};

void
R_init_idlBNs(DllInfo* dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
