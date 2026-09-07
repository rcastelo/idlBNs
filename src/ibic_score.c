#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "cache_key.h"
#include "cache_ref.h"
#include "nh_scores.h"

/* prototypes */

double
iBIC_node_score(const double* Sj, int p1, const int* pa, int lp, int node,
                double Nj, double n);

/*
 * iBIC_node_score
 *
 * Computes the iBIC score contribution for one node, replicating the inner
 * body of the iBIC() R loop using LAPACK Cholesky + triangular solve.
 *
 * Arguments
 * ---------
 * Sj      double*  the (p+1)x(p+1) sufficient-statistics matrix S[[i]],
 *                  stored in column-major order (as R matrices always are)
 * p1      int      p+1, dim of square S
 * pa      int*     parent variable indices, 1-based (may be length 0)
 * lp      int      number of parents
 * node    int      the R loop variable i (1-based, range 1..p)
 *                  NOTE: the 0-based column of the response in S is node=i,
 *                  not i-1, because R accesses S at column i+1 (1-based) which
 *                  is column i in 0-based indexing
 * Nj      double   number of non-intervened observations for this node
 * n       double   total observation count n (used to compute lambda)
 *
 * Returns the node score s.
 */
double
iBIC_node_score(const double* Sj, int p1, const int* pa, int lp, int node,
                double Nj, double n) {
    int m = lp + 1;                  /* intercept + parents */

    /* build 0-based index array: [0, pa[0], pa[1], ...] */
    int* idx = (int *) R_alloc((size_t) m, sizeof(int));
    idx[0] = 0;
    if (node <= 0 || node >= p1)
        error("iBIC_node_score: node index %d out of range [1,%d]", node, p1-1);

    for (int k = 0; k < lp; k++) {
        int pk = pa[k];

        if (pk == NA_INTEGER || pk <= 0 || pk >= p1 || pk == node)
            error("iBIC_node_score: invalid parent index (%d) for node %d and p1=%d",
                  pk, node, p1);
        idx[k + 1] = pk;   /* pa[k] is 1-based R variable index, it equals correct
                              0-based S column */
    }

    /* extract ZtZ (mxm), column-major */
    double* ZtZ = (double *) R_alloc((size_t) m * m, sizeof(double));
    for (int c = 0; c < m; c++)
        for (int r = 0; r < m; r++)
            ZtZ[c * m + r] = Sj[(size_t)idx[c] * p1 + idx[r]];

    /* extract ZtY (mx1) */
    double* ZtY = (double *) R_alloc((size_t) m, sizeof(double));
    for (int r = 0; r < m; r++)
        ZtY[r] = Sj[(size_t)node * p1 + idx[r]];

    /* extract YtY scalar */
    double YtY = Sj[(size_t)node * p1 + node];

    /* Cholesky factorisation of ZtZ (upper triangular in-place) */
    int info;
    F77_CALL(dpotrf)("U", &m, ZtZ, &m, &info FCONE);
    if (info != 0)
        error("iBIC_node_score: dpotrf failed (info=%d); "
              "ZtZ is not positive definite", info);

    /* triangular solve: R^T * cc = ZtY, result overwrites ZtY */
    /* R code: cc <- backsolve(R, ZtY, transpose=TRUE)
       backsolve with transpose=TRUE solves R^T x = b, matching dtrtrs "U","T","N" */
    int nrhs = 1;
    F77_CALL(dtrtrs)("U", "T", "N", &m, &nrhs, ZtZ, &m, ZtY, &m, &info
                     FCONE FCONE FCONE);
    if (info != 0)
        error("iBIC_node_score: dtrtrs failed (info=%d)", info);

    /* RSS = YtY - sum(cc^2) */
    double rss = YtY;
    for (int k = 0; k < m; k++)
        rss -= ZtY[k] * ZtY[k];

    /* iBIC node score */
    double lambda = 0.5 * log(n);
    double s = -0.5 * Nj * (1.0 + log(rss / Nj)) - lambda * (1.0 + lp);

    return s;
}

/*
 * C_iBIC_score
 *
 * Computes the iBIC score for a whole DAG, replicating the inner body of
 * the iBIC() R loop (cache lookup, cache-miss computation via
 * iBIC_node_score(), cache write-back) entirely in C, so that a single
 * .Call() replaces what used to be one .Call() per cache-miss node inside
 * an R for loop. The caching logic and key format are unchanged from
 * .cached_scores_key()/cached.scores[[i]][[k]] -- only their
 * implementation moved from R to C.
 *
 * Arguments
 * ---------
 * S_R             VECSXP   global.sufstats$S: a list of p (p+1)x(p+1)
 *                          sufficient-statistics matrices, one per vertex
 * pasets_R        VECSXP   a list of p integer vectors, one per vertex,
 *                          the 1-based parent indices of that vertex
 * data_count_R    REALSXP  global.sufstats$data.count, length p
 * n_R             REALSXP  scalar: global.sufstats$n
 * cached_scores_R VECSXP   a list of p environments, or R_NilValue if no
 *                          caching is requested (cached.scores=NULL)
 *
 * Returns a length-1 REALSXP containing the total score (sum over nodes).
 */
SEXP
C_iBIC_score(SEXP S_R, SEXP pasets_R, SEXP data_count_R, SEXP n_R,
            SEXP cached_scores_R) {
    int p = LENGTH(pasets_R);
    double n = REAL(n_R)[0];
    idl_cache_ref cr = idl_cache_resolve(cached_scores_R, p, "C_iBIC_score");
    double total = 0.0;

    for (int i = 0; i < p; i++) {
        void *vmax = vmaxget(); /* bound R_alloc accumulation to one node
                                   at a time, instead of the whole loop */

        SEXP pa_R = VECTOR_ELT(pasets_R, i);
        const int *pa = INTEGER(pa_R);
        int lp = LENGTH(pa_R);

        double s = 0.0;
        int found = 0;

        /* the cache is keyed on the ASCENDING parent set, while the score
           itself is computed from 'pa' in its given order -- see
           cached_node_score() in nh_scores.c for why that distinction
           matters. p sorts per call is nothing here, unlike in the
           neighbourhood driver where it would be p sorts per candidate. */
        int *key = (int *) R_alloc((size_t) (lp > 0 ? lp : 1), sizeof(int));
        if (lp > 0) {
            memcpy(key, pa, (size_t) lp * sizeof(int));
            qsort(key, (size_t) lp, sizeof(int), cmp_int);
        }
        if (cr.kind != IDL_CACHE_NONE)
            found = idl_cache_get(&cr, i, key, lp, &s);

        if (!found) {
            SEXP Sj_R = VECTOR_ELT(S_R, i);
            int p1 = (int) sqrt((double) LENGTH(Sj_R));
            double Nj = REAL(data_count_R)[i];
            s = iBIC_node_score(REAL(Sj_R), p1, pa, lp, i + 1, Nj, n);
            if (cr.kind != IDL_CACHE_NONE)
                idl_cache_put(&cr, i, key, lp, s);
        }
        total += s;

        vmaxset(vmax); /* reclaim this node's scratch space (idx/ZtZ/ZtY
                           inside iBIC_node_score(), plus the key buffer) */
    }

    return Rf_ScalarReal(total);
}

/*
 * iBIC_ctx / iBIC_node_score_thunk
 *
 * Binds the iBIC global sufficient statistics to the generic
 * node_score_fn signature that nh_scores_driver() calls back into, so
 * that the delta-scoring driver stays score-function agnostic.
 */
typedef struct {
    SEXP          S_R;        /* global.sufstats$S, a list of p matrices  */
    const double *data_count; /* global.sufstats$data.count, length p     */
    double        n;          /* global.sufstats$n                        */
} iBIC_ctx;

static double
iBIC_node_score_thunk(void *ctxv, const int *pa, int lp, int node) {
    iBIC_ctx *ctx = (iBIC_ctx *) ctxv;
    SEXP Sj_R = VECTOR_ELT(ctx->S_R, node - 1);
    int p1 = (int) sqrt((double) LENGTH(Sj_R));

    return iBIC_node_score(REAL(Sj_R), p1, pa, lp, node,
                           ctx->data_count[node - 1], ctx->n);
}

/*
 * C_iBIC_nh_scores
 *
 * Scores a whole neighborhood of candidate moves against the current DAG
 * in a single .Call(), returning one total iBIC score per candidate. See
 * nh_scores.c for how the per-candidate totals are derived from the
 * current DAG's per-vertex terms.
 *
 * Arguments
 * ---------
 * S_R             VECSXP   global.sufstats$S, a list of p (p+1)x(p+1)
 *                          sufficient-statistics matrices
 * pasets_R        VECSXP   a list of p integer vectors, the 1-based
 *                          parent indices of each vertex in the current
 *                          DAG
 * data_count_R    REALSXP  global.sufstats$data.count, length p
 * n_R             REALSXP  scalar: global.sufstats$n
 * cached_scores_R VECSXP   a list of p environments, or R_NilValue
 * op_R            INTSXP   move operation codes, length k
 * u_R             INTSXP   move tail vertices, 1-based, length k
 * v_R             INTSXP   move head vertices, 1-based, length k
 *
 * Returns a REALSXP of length k.
 */
SEXP
C_iBIC_nh_scores(SEXP S_R, SEXP pasets_R, SEXP data_count_R, SEXP n_R,
                 SEXP cached_scores_R, SEXP op_R, SEXP u_R, SEXP v_R) {
    if (TYPEOF(pasets_R) != VECSXP)
        error("C_iBIC_nh_scores: 'pasets' must be a list");

    int p = LENGTH(pasets_R);

    if (TYPEOF(S_R) != VECSXP || LENGTH(S_R) != p)
        error("C_iBIC_nh_scores: 'S' must be a list of length %d", p);
    if (TYPEOF(data_count_R) != REALSXP || LENGTH(data_count_R) != p)
        error("C_iBIC_nh_scores: 'data_count' must be a numeric vector of length %d",
              p);
    if (TYPEOF(n_R) != REALSXP || LENGTH(n_R) != 1)
        error("C_iBIC_nh_scores: 'n' must be a numeric scalar");

    iBIC_ctx ctx;
    ctx.S_R        = S_R;
    ctx.data_count = REAL(data_count_R);
    ctx.n          = REAL(n_R)[0];

    return nh_scores_driver(pasets_R, cached_scores_R, op_R, u_R, v_R,
                            iBIC_node_score_thunk, &ctx);
}

/*
 * C_iBIC_nh_argmax
 *
 * As C_iBIC_nh_scores(), but returns only the winning candidate --
 * list(index=, total=, band=, worst=) -- found through the error-bounded
 * candidate band rather than by summing all p vertex terms for every
 * candidate. See the band commentary in nh_scores.c.
 *
 * verify_R  LGLSXP  when TRUE, additionally scores every candidate exactly
 *                   and checks the band against it. O(p * k), for testing.
 */
SEXP
C_iBIC_nh_argmax(SEXP S_R, SEXP pasets_R, SEXP data_count_R, SEXP n_R,
                  SEXP cached_scores_R, SEXP op_R, SEXP u_R, SEXP v_R,
                  SEXP stamp_R, SEXP verify_R) {
    if (TYPEOF(pasets_R) != VECSXP)
        error("C_iBIC_nh_argmax: 'pasets' must be a list");

    iBIC_ctx ctx;
    ctx.S_R        = S_R;
    ctx.data_count = REAL(data_count_R);
    ctx.n          = REAL(n_R)[0];

    return nh_argmax_driver(pasets_R, cached_scores_R, op_R, u_R, v_R,
                            stamp_R, asLogical(verify_R) == TRUE, iBIC_node_score_thunk, &ctx);
}
