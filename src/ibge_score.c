#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "cache_key.h"

/* prototypes */

double
iBGe_node_score(const double* TNj, int p, const int* pa, int lp, int node,
               double awpN_i, double gsP, const double* scoreconstvec_i);

/*
 * C_iBGe_node_score
 *
 * Computes the iBGe score contribution for one node, replicating the
 * inner body of the iBGe() R loop using LAPACK Cholesky + triangular
 * solve. This is an R-facing API C wrapper for the iBGe_node_score()
 * function, which does the actual computation.
 *
 * Arguments
 * ---------
 * TNj_R              REALSXP  the p x p matrix global.sufstats$TN[[i]],
 *                             stored in column-major order
 * pa_R               INTSXP   parent variable indices, 1-based (may be
 *                             length 0)
 * node_R             INTSXP   scalar: the R loop variable i (1-based,
 *                             range 1..p) -- unlike iBIC's S, TN has no
 *                             intercept row/column, so node_R is used as
 *                             a plain 1-based index, converted to 0-based
 *                             internally (no offset-cancellation trick)
 * awpN_i_R           REALSXP  scalar: global.sufstats$awpN[i]
 * gsP_R              REALSXP  scalar: global.sufstats$p (the total number
 *                             of vertices, same for every node)
 * scoreconstvec_i_R  REALSXP  global.sufstats$scoreconstvec[[i]], length p
 *
 * Returns a length-1 REALSXP containing the node score s.
 */
SEXP
C_iBGe_node_score(SEXP TNj_R, SEXP pa_R, SEXP node_R, SEXP awpN_i_R,
                  SEXP gsP_R, SEXP scoreconstvec_i_R) {

    if (TYPEOF(TNj_R) != REALSXP)
        error("C_iBGe_node_score: 'TNj' must be a numeric matrix");
    R_xlen_t TNj_len = XLENGTH(TNj_R);

    int p      = (int) sqrt((double) TNj_len); /* dim of square TN  */

    if ((R_xlen_t) p * p != TNj_len)
        error("C_iBGe_node_score: 'TNj' must be a square matrix, but length(TNj)=%lld is not a perfect square",
              (long long) TNj_len);
    if (TYPEOF(pa_R) != INTSXP)
        error("C_iBGe_node_score: 'pa' must be an integer vector");

    int lp     = LENGTH(pa_R);                 /* number of parents */

    if (TYPEOF(node_R) != INTSXP || XLENGTH(node_R) != 1)
        error("C_iBGe_node_score: 'node' must be an integer scalar");

    int node   = INTEGER(node_R)[0];           /* 1-based node idx  */

    double awpN_i = REAL(awpN_i_R)[0];
    double gsP    = REAL(gsP_R)[0];

    const double* TNj = REAL(TNj_R);
    const int*    pa  = INTEGER(pa_R);
    const double* scoreconstvec_i = REAL(scoreconstvec_i_R);

    double s = iBGe_node_score(TNj, p, pa, lp, node, awpN_i, gsP,
                               scoreconstvec_i);

    /* return as length-1 numeric vector */
    SEXP result = PROTECT(allocVector(REALSXP, 1));
    REAL(result)[0] = s;
    UNPROTECT(1);

    return result;
}

/*
 * iBGe_node_score
 *
 * Computes the iBGe score contribution for one node, replicating the
 * inner body of the iBGe() R loop using LAPACK Cholesky + triangular
 * solve.
 *
 * Arguments
 * ---------
 * TNj             double*  the p x p matrix TN[[i]], column-major,
 *                          symmetric by construction (T0 diagonal +
 *                          covariance + rank-1 outer product)
 * p               int      dimension of TNj (== global.sufstats$p)
 * pa              int*     parent variable indices, 1-based (may be
 *                          length 0)
 * lp              int      number of parents
 * node            int      1-based node index, range 1..p (converted to
 *                          0-based internally -- no offset-cancellation
 *                          trick here, unlike iBIC_node_score(), since TN
 *                          has no intercept row/column)
 * awpN_i          double   global.sufstats$awpN[i]
 * gsP             double   global.sufstats$p (same for every node)
 * scoreconstvec_i double*  global.sufstats$scoreconstvec[[i]], length p
 *
 * Returns the node score s.
 */
double
iBGe_node_score(const double* TNj, int p, const int* pa, int lp, int node,
                double awpN_i, double gsP, const double* scoreconstvec_i) {
    if (node <= 0 || node > p)
        error("iBGe_node_score: node index %d out of range [1,%d]", node, p);

    int    node0 = node - 1; /* 0-based row/col index for node */
    double A = TNj[(size_t) node0 * p + node0];
    double awpNd2 = (awpN_i - gsP + lp + 1) / 2.0;

    /* closed-form fast path for a root node (no parents): no Cholesky
       needed at all, matching the R code's if (lp == 0L) branch */
    if (lp == 0)
        return scoreconstvec_i[0] - awpNd2 * log(A);

    /* build 0-based parent index array, with the same bounds safeguards
       as iBIC_node_score() */
    int* idx = (int *) R_alloc(lp, sizeof(int));
    for (int k = 0; k < lp; k++) {
        int pk = pa[k];

        if (pk == NA_INTEGER || pk <= 0 || pk > p || pk == node)
            error("iBGe_node_score: invalid parent index (%d) for node %d and p=%d",
                  pk, node, p);
        idx[k] = pk - 1;
    }

    /* extract D = TNj[pa, pa] (lp x lp), column-major */
    double* D = (double *) R_alloc((size_t) lp * lp, sizeof(double));
    for (int c = 0; c < lp; c++)
        for (int r = 0; r < lp; r++)
            D[c * lp + r] = TNj[(size_t) idx[c] * p + idx[r]];

    /* extract B = TNj[node, pa] (lp x 1); TNj is symmetric, so this
       equals TNj[pa, node] too */
    double* B = (double *) R_alloc(lp, sizeof(double));
    for (int r = 0; r < lp; r++)
        B[r] = TNj[(size_t) idx[r] * p + node0];

    /* Cholesky factorisation of D (upper triangular in-place) */
    int info;
    F77_CALL(dpotrf)("U", &lp, D, &lp, &info FCONE);
    if (info != 0)
        error("iBGe_node_score: dpotrf failed (info=%d); "
              "D is not positive definite", info);

    /* log determinant of D via 2*sum(log(diag(R))), R = Cholesky factor
       now held in D's upper triangle (including the diagonal) */
    double logdetD = 0.0;
    for (int k = 0; k < lp; k++)
        logdetD += log(D[k * lp + k]);
    logdetD *= 2.0;

    /* triangular solve: R^T * x = B, result overwrites B */
    /* R code: backsolve(R, B, transpose=TRUE), matching dtrtrs "U","T","N" */
    int nrhs = 1;
    F77_CALL(dtrtrs)("U", "T", "N", &lp, &nrhs, D, &lp, B, &lp, &info
                     FCONE FCONE FCONE);
    if (info != 0)
        error("iBGe_node_score: dtrtrs failed (info=%d)", info);

    double ssq = 0.0;
    for (int k = 0; k < lp; k++)
        ssq += B[k] * B[k];

    double logdetpart2 = log(A - ssq);
    double s = scoreconstvec_i[lp] - awpNd2 * logdetpart2 - logdetD / 2.0;

    return s;
}

/*
 * C_iBGe_score
 *
 * Computes the iBGe score for a whole DAG, replicating the inner body of
 * the iBGe() R loop (cache lookup, cache-miss computation via
 * iBGe_node_score(), cache write-back) entirely in C, so that a single
 * .Call() replaces what used to be a per-node R for loop. The caching
 * logic and key format are unchanged from .cached_scores_key()/
 * cached.scores[[i]][[k]] -- only their implementation moved from R to C
 * (see cache_key.h, shared with C_iBIC_score()).
 *
 * Arguments
 * ---------
 * TN_R            VECSXP   global.sufstats$TN: a list of p p x p
 *                          matrices, one per vertex
 * pasets_R        VECSXP   a list of p integer vectors, one per vertex,
 *                          the 1-based parent indices of that vertex
 * awpN_R          REALSXP  global.sufstats$awpN, length p
 * gsP_R           REALSXP  scalar: global.sufstats$p
 * scoreconstvec_R VECSXP   global.sufstats$scoreconstvec: a list of p
 *                          numeric vectors, each of length p
 * cached_scores_R VECSXP   a list of p environments, or R_NilValue if no
 *                          caching is requested (cached.scores=NULL)
 *
 * Returns a length-1 REALSXP containing the total score (sum over nodes).
 */
SEXP
C_iBGe_score(SEXP TN_R, SEXP pasets_R, SEXP awpN_R, SEXP gsP_R,
             SEXP scoreconstvec_R, SEXP cached_scores_R) {
    if (TYPEOF(pasets_R) != VECSXP)
        error("C_iBGe_score: 'pasets' must be a list");

    int p = LENGTH(pasets_R);

    if (TYPEOF(TN_R) != VECSXP || LENGTH(TN_R) != p)
        error("C_iBGe_score: 'TN' must be a list of length %d", p);
    if (TYPEOF(scoreconstvec_R) != VECSXP || LENGTH(scoreconstvec_R) != p)
        error("C_iBGe_score: 'scoreconstvec' must be a list of length %d", p);
    if (TYPEOF(awpN_R) != REALSXP || LENGTH(awpN_R) != p)
        error("C_iBGe_score: 'awpN' must be a numeric vector of length %d", p);
    if (TYPEOF(gsP_R) != REALSXP || LENGTH(gsP_R) != 1)
        error("C_iBGe_score: 'gsP' must be a numeric scalar");

    double gsP = REAL(gsP_R)[0];
    int has_cache = (cached_scores_R != R_NilValue);

    if (has_cache) {
        if (TYPEOF(cached_scores_R) != VECSXP || LENGTH(cached_scores_R) != p)
            error("C_iBGe_score: 'cached_scores' must be either NULL or a list of length %d",
                  p);
    }
    double total = 0.0;

    for (int i = 0; i < p; i++) {
        void *vmax = vmaxget(); /* bound R_alloc accumulation to one node
                                   at a time, instead of the whole loop */

        SEXP pa_R = VECTOR_ELT(pasets_R, i);
        const int *pa = INTEGER(pa_R);
        int lp = LENGTH(pa_R);

        SEXP env = R_NilValue;
        SEXP sym = R_NilValue;
        double s = 0.0;
        int found = 0;

        if (has_cache) {
            env = VECTOR_ELT(cached_scores_R, i);
            char *key = build_cache_key(pa, lp);
            sym = Rf_install(key); /* symbols are GC-safe unprotected */
            if (R_existsVarInFrame(env, sym)) {
                SEXP val = R_getVar(sym, env, FALSE);
                s = REAL(val)[0];
                found = 1;
            }
        }

        if (!found) {
            SEXP TNj_R = VECTOR_ELT(TN_R, i);
            int tp = (int) sqrt((double) LENGTH(TNj_R));
            double awpN_i = REAL(awpN_R)[i];
            const double *scoreconstvec_i = REAL(VECTOR_ELT(scoreconstvec_R, i));
            s = iBGe_node_score(REAL(TNj_R), tp, pa, lp, i + 1, awpN_i, gsP,
                               scoreconstvec_i);
            if (has_cache) {
                /* Rf_install() above and Rf_ScalarReal() here are kept as
                   separate statements, and the ScalarReal() result is
                   protected before defineVar(): nesting both allocating
                   calls as sibling arguments to defineVar() would leave
                   the unprotected ScalarReal() result exposed to GC under
                   C's unspecified argument-evaluation order */
                SEXP val = PROTECT(Rf_ScalarReal(s));
                Rf_defineVar(sym, val, env);
                UNPROTECT(1);
            }
        }
        total += s;

        vmaxset(vmax); /* reclaim this node's scratch space (idx/D/B
                           inside iBGe_node_score(), plus the key buffer) */
    }

    return Rf_ScalarReal(total);
}
