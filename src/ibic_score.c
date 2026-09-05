#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

/*
 * C_iBIC_node_score
 *
 * Computes the iBIC score contribution for one node, replicating the inner
 * body of the iBIC() R loop using LAPACK Cholesky + triangular solve.
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
        ZtY[r] = Sj[(size_t)node_C * p1 + idx[r]];

    /* extract YtY scalar */
    double YtY = Sj[(size_t)node_C * p1 + node_C];

    /* Cholesky factorisation of ZtZ (upper triangular in-place) */
    int info;
    F77_CALL(dpotrf)("U", &m, ZtZ, &m, &info FCONE);
    if (info != 0)
        error("C_iBIC_node_score: dpotrf failed (info=%d); "
              "ZtZ is not positive definite", info);

    /* triangular solve: R^T * cc = ZtY, result overwrites ZtY */
    /* R code: cc <- backsolve(R, ZtY, transpose=TRUE)
       backsolve with transpose=TRUE solves R^T x = b, matching dtrtrs "U","T","N" */
    int nrhs = 1;
    F77_CALL(dtrtrs)("U", "T", "N", &m, &nrhs, ZtZ, &m, ZtY, &m, &info
                     FCONE FCONE FCONE);
    if (info != 0)
        error("C_iBIC_node_score: dtrtrs failed (info=%d)", info);

    /* RSS = YtY - sum(cc^2) */
    double rss = YtY;
    for (int k = 0; k < m; k++)
        rss -= ZtY[k] * ZtY[k];

    /* iBIC node score */
    double lambda = 0.5 * log(n);
    double s = -0.5 * Nj * (1.0 + log(rss / Nj)) - lambda * (1.0 + lp);

    /* return as length-1 numeric vector */
    SEXP result = PROTECT(allocVector(REALSXP, 1));
    REAL(result)[0] = s;
    UNPROTECT(1);

    return result;
}
