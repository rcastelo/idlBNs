#include <R.h>
#include <Rinternals.h>
#include "cache_key.h"
#include "nh_scores.h"

/*
 * DELTA SCORING OF A WHOLE NEIGHBORHOOD
 *
 * Both the iBIC and the iBGe score of a DAG decompose as a sum over
 * vertices of a term that depends only on (vertex, its parent set):
 *
 *     score(G) = sum_i s(i, pa_G(i))
 *
 * A candidate move changes the parent set of exactly one vertex (an
 * addition or a removal both only touch the head v) or of exactly two (a
 * reversal of u -> v touches v, which loses u, and u, which gains v).
 * Every other vertex keeps the term it already had. So the total score of
 * a candidate is the current DAG's total with one or two terms swapped
 * out, and the whole neighborhood costs O(|NH|) node scores instead of
 * the O(p * |NH|) that re-summing all p vertices per candidate costs.
 *
 * The per-vertex terms of the *current* DAG are computed once here, up
 * front. What is saved is the *scoring* of vertices, not the summing of
 * their terms: each candidate's total is still accumulated over all p
 * vertices in vertex order, with the one or two changed terms substituted
 * in, so it comes out bit-identical to what re-scoring that candidate
 * from scratch would produce. p double additions per candidate are
 * nothing beside a cache lookup plus a possible Cholesky, and buying
 * exactness with them is worth it.
 *
 * Exactness here is not fussiness. Summing the same terms in a different
 * association order moves a total by an ulp or so, and two candidates
 * that a full re-sum separates by less than that -- which happens
 * constantly, since both iBIC and iBGe score Markov equivalent DAGs alike
 * on observational data -- can then swap places under which.max(). Each
 * such swap is harmless on its own, the two moves being equally good, but
 * a greedy search compounds it: one different choice leads to a different
 * neighborhood at the next step, and the run settles into a different
 * local optimum. Measured over 48 searches, that shifted final scores by
 * as much as 90 units in either direction. Summing over all p vertices
 * keeps every trajectory exactly where it was.
 *
 * The modified parent set handed to the score function is built to match,
 * element for element, what add.pasets()/remove.pasets()/reverse.pasets()
 * would have produced on the R side: an addition appends the new parent,
 * a removal keeps the surviving parents in their existing order. Cache
 * keys sort the parent set (see build_cache_key()), so this ordering does
 * not affect cache hits, but it does keep a cache-miss computation
 * bit-identical to the pre-delta-scoring code.
 */

/* score one node with the parent set 'pa', going through the per-node
   score cache when one was supplied */
static double
cached_node_score(node_score_fn fn, void *ctx, SEXP cached_scores_R,
                  int node0, const int *pa, int lp) {
    double s;

    if (cached_scores_R == R_NilValue)
        return fn(ctx, pa, lp, node0 + 1);

    SEXP env = VECTOR_ELT(cached_scores_R, node0);
    SEXP sym = R_NilValue;
    if (cache_lookup(env, pa, lp, &sym, &s))
        return s;

    s = fn(ctx, pa, lp, node0 + 1);
    cache_store(env, sym, s);

    return s;
}

/* copy pa[] into buf[] and append 'add1' (a 1-based vertex index),
   mirroring add.pasets()'s c(pasets[[v]], u); returns the new length */
static int
paset_with(const int *pa, int lp, int add1, int *buf) {
    for (int k = 0; k < lp; k++)
        buf[k] = pa[k];
    buf[lp] = add1;

    return lp + 1;
}

/* copy pa[] into buf[] dropping 'drop1' (a 1-based vertex index) while
   keeping the surviving parents in order, mirroring remove.pasets()'s
   setdiff(pasets[[v]], u); returns the new length */
static int
paset_without(const int *pa, int lp, int drop1, int *buf) {
    int m = 0;
    for (int k = 0; k < lp; k++)
        if (pa[k] != drop1)
            buf[m++] = pa[k];

    return m;
}

/*
 * total_with
 *
 * Sums the p per-vertex terms in 'base' in vertex order, with vertex c1's
 * term replaced by s1 and, when c2 >= 0, vertex c2's replaced by s2. The
 * replacements are written into 'base' and restored afterwards rather
 * than branched around inside the loop, so the accumulation is the same
 * straight-line "total += term[i]" over i = 0..p-1 that C_iBIC_score() and
 * C_iBGe_score() run, and therefore lands on the same double, bit for
 * bit. 'base' is restored on return, so it stays the current DAG's terms.
 */
static double
total_with(double *base, int p, int c1, double s1, int c2, double s2) {
    double sav1 = base[c1];
    double sav2 = 0.0;

    base[c1] = s1;
    if (c2 >= 0) {
        sav2 = base[c2];
        base[c2] = s2;
    }

    double total = 0.0;
    for (int i = 0; i < p; i++)
        total += base[i];

    base[c1] = sav1;
    if (c2 >= 0)
        base[c2] = sav2;

    return total;
}

SEXP
nh_scores_driver(SEXP pasets_R, SEXP cached_scores_R, SEXP op_R, SEXP u_R,
                 SEXP v_R, node_score_fn fn, void *ctx) {
    if (TYPEOF(pasets_R) != VECSXP)
        error("nh_scores_driver: 'pasets' must be a list");
    if (TYPEOF(op_R) != INTSXP || TYPEOF(u_R) != INTSXP ||
        TYPEOF(v_R) != INTSXP)
        error("nh_scores_driver: 'op', 'u' and 'v' must be integer vectors");

    int p = LENGTH(pasets_R);
    R_xlen_t k = XLENGTH(op_R);

    if (XLENGTH(u_R) != k || XLENGTH(v_R) != k)
        error("nh_scores_driver: 'op', 'u' and 'v' must have the same length");
    if (cached_scores_R != R_NilValue &&
        (TYPEOF(cached_scores_R) != VECSXP || LENGTH(cached_scores_R) != p))
        error("nh_scores_driver: 'cached_scores' must be either NULL or a list of length %d",
              p);

    const int *op = INTEGER(op_R);
    const int *uu = INTEGER(u_R);
    const int *vv = INTEGER(v_R);

    SEXP res_R = PROTECT(allocVector(REALSXP, k));
    double *res = REAL(res_R);

    /* per-vertex terms of the current DAG, and their sum in vertex order */
    void *vmax0 = vmaxget();
    double *base = (double *) R_alloc((size_t) p, sizeof(double));
    for (int i = 0; i < p; i++) {
        void *vmax = vmaxget();
        SEXP pa_R = VECTOR_ELT(pasets_R, i);
        base[i] = cached_node_score(fn, ctx, cached_scores_R, i,
                                    INTEGER(pa_R), LENGTH(pa_R));
        vmaxset(vmax);
    }
    /* scratch for one modified parent set; a vertex can gain at most one
       parent over its current set, so p ints always suffice */
    int *buf = (int *) R_alloc((size_t) p, sizeof(int));

    for (R_xlen_t m = 0; m < k; m++) {
        void *vmax = vmaxget();
        int u0 = uu[m] - 1;
        int v0 = vv[m] - 1;

        if (u0 < 0 || u0 >= p || v0 < 0 || v0 >= p || u0 == v0)
            error("nh_scores_driver: invalid move (u=%d, v=%d) for p=%d",
                  uu[m], vv[m], p);

        SEXP pav_R = VECTOR_ELT(pasets_R, v0);
        const int *pav = INTEGER(pav_R);
        int lpv = LENGTH(pav_R);
        int lp;
        double s;

        switch (op[m]) {
        case IDLBNS_OP_ADD:
            lp = paset_with(pav, lpv, uu[m], buf);
            s = cached_node_score(fn, ctx, cached_scores_R, v0, buf, lp);
            res[m] = total_with(base, p, v0, s, -1, 0.0);
            break;
        case IDLBNS_OP_REMOVE:
            lp = paset_without(pav, lpv, uu[m], buf);
            s = cached_node_score(fn, ctx, cached_scores_R, v0, buf, lp);
            res[m] = total_with(base, p, v0, s, -1, 0.0);
            break;
        case IDLBNS_OP_REVERSE: {
            lp = paset_without(pav, lpv, uu[m], buf);
            double sv = cached_node_score(fn, ctx, cached_scores_R, v0, buf,
                                          lp);
            SEXP pau_R = VECTOR_ELT(pasets_R, u0);
            lp = paset_with(INTEGER(pau_R), LENGTH(pau_R), vv[m], buf);
            double su = cached_node_score(fn, ctx, cached_scores_R, u0, buf,
                                          lp);
            res[m] = total_with(base, p, v0, sv, u0, su);
            break;
        }
        default:
            error("nh_scores_driver: unknown move operation code %d", op[m]);
        }

        vmaxset(vmax); /* reclaim this move's cache-key and score scratch */
    }

    vmaxset(vmax0);
    UNPROTECT(1);

    return res_R;
}
