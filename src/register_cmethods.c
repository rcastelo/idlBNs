#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP
C_iBIC_score(SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP
C_iBGe_score(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP
C_iBIC_nh_scores(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP
C_iBGe_nh_scores(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP
C_unif_index(SEXP);

/* registration of C-entry points */

static const R_CallMethodDef CallEntries[] = {
    {"C_iBIC_score",      (DL_FUNC) &C_iBIC_score,      5},
    {"C_iBGe_score",      (DL_FUNC) &C_iBGe_score,      6},
    {"C_iBIC_nh_scores",  (DL_FUNC) &C_iBIC_nh_scores,  8},
    {"C_iBGe_nh_scores",  (DL_FUNC) &C_iBGe_nh_scores,  9},
    {"C_unif_index",      (DL_FUNC) &C_unif_index,      1},
    {NULL, NULL, 0}
};

void
R_init_idlBNs(DllInfo* dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
