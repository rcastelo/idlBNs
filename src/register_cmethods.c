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

extern SEXP C_dag_new(SEXP);
extern SEXP C_dag_apply_move(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_dag_pasets(SEXP);
extern SEXP C_dag_pasets_sorted(SEXP);
extern SEXP C_dag_edgeM(SEXP);
extern SEXP C_dag_nedges(SEXP);
extern SEXP C_dag_anc(SEXP);
extern SEXP C_dag_desc(SEXP);
extern SEXP C_dag_check(SEXP);
extern SEXP C_dag_nh(SEXP, SEXP, SEXP);
extern SEXP C_dag_cedges(SEXP, SEXP);

/* installs the external-pointer tag symbol; see src/dag_R.c */
extern void idl_dag_R_init(void);

/* registration of C-entry points */

static const R_CallMethodDef CallEntries[] = {
    {"C_iBIC_score",      (DL_FUNC) &C_iBIC_score,      5},
    {"C_iBGe_score",      (DL_FUNC) &C_iBGe_score,      6},
    {"C_iBIC_nh_scores",  (DL_FUNC) &C_iBIC_nh_scores,  8},
    {"C_iBGe_nh_scores",  (DL_FUNC) &C_iBGe_nh_scores,  9},
    {"C_unif_index",      (DL_FUNC) &C_unif_index,      1},
    {"C_dag_new",           (DL_FUNC) &C_dag_new,           1},
    {"C_dag_apply_move",    (DL_FUNC) &C_dag_apply_move,    4},
    {"C_dag_pasets",        (DL_FUNC) &C_dag_pasets,        1},
    {"C_dag_pasets_sorted", (DL_FUNC) &C_dag_pasets_sorted, 1},
    {"C_dag_edgeM",         (DL_FUNC) &C_dag_edgeM,         1},
    {"C_dag_nedges",        (DL_FUNC) &C_dag_nedges,        1},
    {"C_dag_anc",           (DL_FUNC) &C_dag_anc,           1},
    {"C_dag_desc",          (DL_FUNC) &C_dag_desc,          1},
    {"C_dag_check",         (DL_FUNC) &C_dag_check,         1},
    {"C_dag_nh",            (DL_FUNC) &C_dag_nh,            3},
    {"C_dag_cedges",        (DL_FUNC) &C_dag_cedges,        2},
    {NULL, NULL, 0}
};

void
R_init_idlBNs(DllInfo* dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    idl_dag_R_init();
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
