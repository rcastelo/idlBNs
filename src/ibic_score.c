#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

/* prototypes */

double
iBIC_node_score(const double* Sj, int p1, const int* pa, int lp, int node,
                double Nj, double n);

/*
 * C_iBIC_node_score
 *
 * Computes the iBIC score contribution for one node, replicating the inner
 * body of the iBIC() R loop using LAPACK Cholesky + triangular solve. This
 * is an R-facing API C wrapper for the iBIC_node_score() function, which does
 * the actual computation.
 *
 * Arguments
 * ---------
 * Sj_R    REALSXP  the (p+1)x(p+1) sufficient-statistics matrix S[[i]],
 *                  stored in column-major order (as R matrices always are)
 * pa_R    INTSXP   parent variable indices, 1-based (may be length 0)
 * node_R  INTSXP   scalar: the R loop variable i (1-based, range 1..p)
 *                  NOTE: the 0-based column of the response in S is node_C=i,
 *                  not i-1, because R accesses S at column i+1 (1-based) which
 *                  is column i in 0-based indexing
 * Nj_R    REALSXP  scalar: number of non-intervened observations for this node
 * n_R     REALSXP  scalar: total observation count n (used to compute lambda)
 *
 * Returns a length-1 REALSXP containing the node score s.
 */
SEXP
C_iBIC_node_score(SEXP Sj_R, SEXP pa_R, SEXP node_R, SEXP Nj_R, SEXP n_R) {
    /* dimensions */
    int p1     = (int)sqrt((double)LENGTH(Sj_R)); /* p+1, dim of square S     */
    int lp     = LENGTH(pa_R);                    /* number of parents        */
    int m      = lp + 1;                          /* intercept + parents      */
    int node_C = INTEGER(node_R)[0];              /* 0-based response column  */
    double Nj  = REAL(Nj_R)[0];
    double n   = REAL(n_R)[0];

    const double* Sj = REAL(Sj_R);
    const int*    pa = INTEGER(pa_R);

    double s = iBIC_node_score(Sj, p1, pa, lp, node_C, Nj, n);

    /* return as length-1 numeric vector */
    SEXP result = PROTECT(allocVector(REALSXP, 1));
    REAL(result)[0] = s;
    UNPROTECT(1);

    return result;
}

/*
 * iBIC_node_score
 *
 * Computes the iBIC score contribution for one node, replicating the inner
 * body of the iBIC() R loop using LAPACK Cholesky + triangular solve.
 *
 * Arguments
 * ---------
 * Sj      double* the (p+1)x(p+1) sufficient-statistics matrix S[[i]],
 *                  stored in column-major order (as R matrices always are)
 * p1      int      p+1, dim of square S
 * pa      INTSXP   parent variable indices, 1-based (may be length 0)
 * lp      int      number of parents
 * node    int      scalar: the R loop variable i (1-based, range 1..p)
 *                  NOTE: the 0-based column of the response in S is node=i,
 *                  not i-1, because R accesses S at column i+1 (1-based) which
 *                  is column i in 0-based indexing
 * Nj      REALSXP  scalar: number of non-intervened observations for this node
 * n       REALSXP  scalar: total observation count n (used to compute lambda)
 *
 * Returns a length-1 REALSXP containing the node score s.
 */
double
iBIC_node_score(const double* Sj, int p1, const int* pa, int lp, int node,
                double Nj, double n) {
    int m = lp + 1;                  /* intercept + parents */

    /* build 0-based index array: [0, pa[0], pa[1], ...] */
    int* idx = (int *) R_alloc(m, sizeof(int));
    idx[0] = 0;
    for (int k = 0; k < lp; k++)
        idx[k + 1] = pa[k];   /* pa[k] is 1-based R variable index;
                                  equals correct 0-based S column (see note) */

    /* extract ZtZ (mxm), column-major */
    double* ZtZ = (double *) R_alloc((size_t) m * m, sizeof(double));
    for (int c = 0; c < m; c++)
        for (int r = 0; r < m; r++)
            ZtZ[c * m + r] = Sj[(size_t)idx[c] * p1 + idx[r]];

    /* extract ZtY (mx1) */
    double* ZtY = (double *) R_alloc(m, sizeof(double));
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

/* ascending comparator for qsort(), avoiding the classic overflow-prone
   "return *a - *b" pattern */
static int
cmp_int(const void *a, const void *b) {
    int ia = *(const int *) a, ib = *(const int *) b;
    return (ia > ib) - (ia < ib);
}

/*
 * build_cache_key
 *
 * Builds the same cache key format as the R function .cached_scores_key(),
 * i.e. paste(sort.int(paset), collapse=":"), or the literal ":" when
 * paset is empty (matching .cached_scores_key()'s nchar(k)==0L special
 * case). Returns an R_alloc'd, NUL-terminated string, freed automatically
 * when the top-level .Call returns (or earlier, if the caller rewinds the
 * R_alloc "vmax" stack with vmaxset() after each node -- see
 * C_iBIC_score()).
 *
 * Assumes pa[] contains no NA_INTEGER values (guaranteed by construction
 * of 'pasets' on the R side; not defensively checked here).
 */
static char *
build_cache_key(const int *pa, int lp) {
    if (lp == 0)
        return ":";

    int *sorted = (int *) R_alloc(lp, sizeof(int));
    memcpy(sorted, pa, (size_t) lp * sizeof(int));
    qsort(sorted, lp, sizeof(int), cmp_int);

    /* worst case: 11 chars per int (sign + 10 digits) + 1 separator */
    size_t bufsize = (size_t) lp * 12 + 1;
    char *buf = (char *) R_alloc(bufsize, sizeof(char));
    size_t pos = 0;
    for (int k = 0; k < lp; k++) {
        int written = snprintf(buf + pos, bufsize - pos,
                               k == 0 ? "%d" : ":%d", sorted[k]);
        pos += (size_t) written;
    }
    return buf;
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
    int has_cache = (cached_scores_R != R_NilValue);
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
            /* R_existsVarInFrame()/R_getVar() (envir.c), not the
               legacy-only Rf_findVarInFrame() (declared in Rinternals.h
               only under #ifdef ENABLE_LEGACY_NONAPI_FUNS, not part of
               the default package-facing C API); inherits=FALSE matches
               cached.scores[[i]]'s single-frame (parent=emptyenv()) R
               "[[" lookup semantics */
            if (R_existsVarInFrame(env, sym)) {
                SEXP val = R_getVar(sym, env, FALSE);
                s = REAL(val)[0];
                found = 1;
            }
        }

        if (!found) {
            SEXP Sj_R = VECTOR_ELT(S_R, i);
            int p1 = (int) sqrt((double) LENGTH(Sj_R));
            double Nj = REAL(data_count_R)[i];
            s = iBIC_node_score(REAL(Sj_R), p1, pa, lp, i + 1, Nj, n);
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

        vmaxset(vmax); /* reclaim this node's scratch space (idx/ZtZ/ZtY
                           inside iBIC_node_score(), plus the key buffer) */
    }

    return Rf_ScalarReal(total);
}
