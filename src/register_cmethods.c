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
C_iBIC_nh_argmax(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                 SEXP);

extern SEXP
C_iBGe_nh_scores(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP
C_iBGe_nh_argmax(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                 SEXP, SEXP);

extern SEXP
C_unif_index(SEXP);

extern SEXP C_dag_new(SEXP);
extern SEXP C_dag_apply_move(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_dag_pasets(SEXP);
extern SEXP C_dag_pasets_sorted(SEXP);
extern SEXP C_dag_edgeM(SEXP);
extern SEXP C_dag_nedges(SEXP);
extern SEXP C_dag_pastamp(SEXP);
extern SEXP C_dag_anc(SEXP);
extern SEXP C_dag_desc(SEXP);
extern SEXP C_dag_check(SEXP);
extern SEXP C_dag_nh(SEXP, SEXP, SEXP);
extern SEXP C_dag_cedges(SEXP, SEXP);
extern SEXP C_dag_rcar(SEXP, SEXP, SEXP);
extern SEXP C_dag_imec_sample(SEXP, SEXP, SEXP);
extern SEXP C_dag_imec_size(SEXP, SEXP, SEXP);
extern SEXP C_dag_imec_members(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_dag_set_edges(SEXP, SEXP);
extern SEXP C_cp_amo_count(SEXP);
extern SEXP C_cp_amo_sample(SEXP);
extern SEXP C_sccache_new(SEXP);
extern SEXP C_sccache_dump(SEXP);
extern SEXP C_sccache_stats(SEXP);

/* installs the external-pointer tag symbol; see src/dag_R.c */
extern void idl_dag_R_init(void);

/* installs the score cache's external-pointer tag; see src/sccache.c */
extern void idl_sc_R_init(void);

/* registration of C-entry points */

static const R_CallMethodDef CallEntries[] = {
    {"C_iBIC_score",      (DL_FUNC) &C_iBIC_score,      5},
    {"C_iBGe_score",      (DL_FUNC) &C_iBGe_score,      6},
    {"C_iBIC_nh_scores",  (DL_FUNC) &C_iBIC_nh_scores,  8},
    {"C_iBGe_nh_scores",  (DL_FUNC) &C_iBGe_nh_scores,  9},
    {"C_iBIC_nh_argmax",  (DL_FUNC) &C_iBIC_nh_argmax, 10},
    {"C_iBGe_nh_argmax",  (DL_FUNC) &C_iBGe_nh_argmax, 11},
    {"C_unif_index",      (DL_FUNC) &C_unif_index,      1},
    {"C_dag_new",           (DL_FUNC) &C_dag_new,           1},
    {"C_dag_apply_move",    (DL_FUNC) &C_dag_apply_move,    4},
    {"C_dag_pasets",        (DL_FUNC) &C_dag_pasets,        1},
    {"C_dag_pasets_sorted", (DL_FUNC) &C_dag_pasets_sorted, 1},
    {"C_dag_edgeM",         (DL_FUNC) &C_dag_edgeM,         1},
    {"C_dag_nedges",        (DL_FUNC) &C_dag_nedges,        1},
    {"C_dag_pastamp",       (DL_FUNC) &C_dag_pastamp,       1},
    {"C_dag_anc",           (DL_FUNC) &C_dag_anc,           1},
    {"C_dag_desc",          (DL_FUNC) &C_dag_desc,          1},
    {"C_dag_check",         (DL_FUNC) &C_dag_check,         1},
    {"C_dag_nh",            (DL_FUNC) &C_dag_nh,            3},
    {"C_dag_cedges",        (DL_FUNC) &C_dag_cedges,        2},
    {"C_dag_rcar",          (DL_FUNC) &C_dag_rcar,          3},
    {"C_dag_imec_sample",   (DL_FUNC) &C_dag_imec_sample,   3},
    {"C_dag_imec_size",     (DL_FUNC) &C_dag_imec_size,     3},
    {"C_dag_imec_members",  (DL_FUNC) &C_dag_imec_members,  4},
    {"C_dag_set_edges",     (DL_FUNC) &C_dag_set_edges,     2},
    {"C_cp_amo_count",      (DL_FUNC) &C_cp_amo_count,      1},
    {"C_cp_amo_sample",     (DL_FUNC) &C_cp_amo_sample,     1},
    {"C_sccache_new",       (DL_FUNC) &C_sccache_new,       1},
    {"C_sccache_dump",      (DL_FUNC) &C_sccache_dump,      1},
    {"C_sccache_stats",     (DL_FUNC) &C_sccache_stats,     1},
    {NULL, NULL, 0}
};

void
R_init_idlBNs(DllInfo* dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    idl_dag_R_init();
    idl_sc_R_init();
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
