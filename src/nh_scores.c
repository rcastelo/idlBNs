#include <R.h>
#include <Rinternals.h>
#include "cache_key.h"
#include "cache_ref.h"
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

/*
 * score one node, going through the per-vertex cache when there is one.
 *
 * TWO ARRAYS, and the difference between them is the whole subtlety of this
 * cache. 'pa' is the parent set in INSERTION order and is what gets scored:
 * the node score is a function of the parent SEQUENCE, because ZtZ is built
 * with the columns in that order and a symmetrically permuted ZtZ Choleskys
 * to a slightly different double. 'key' is the same set ASCENDING and is
 * what the cache is keyed by.
 *
 * So the cached value for a set is whichever order-variant of it missed
 * first, and the cache is part of the arithmetic rather than a transparent
 * memo of it. Both cache backends key identically, which is what makes them
 * store identical doubles given the same sequence of lookups.
 */
static double
cached_node_score(node_score_fn fn, void *ctx, idl_cache_ref *cr,
                  int node0, const int *pa, int lp,
                  const int *key, int klen) {
    double s;

    if (cr->kind == IDL_CACHE_NONE)
        return fn(ctx, pa, lp, node0 + 1);

    if (idl_cache_get(cr, node0, key, klen, &s))
        return s;

    s = fn(ctx, pa, lp, node0 + 1);
    idl_cache_put(cr, node0, key, klen, s);

    return s;
}

/* insert 'add1' into the ascending, duplicate-free 'key' -- the sorted
   counterpart of paset_with(). O(k): a binary-searched position and one
   copy, no sort. */
static int
key_with(const int *key, int klen, int add1, int *buf) {
    int lo = 0, hi = klen;
    while (lo < hi) {
        int mid = lo + ((hi - lo) >> 1);
        if (key[mid] < add1)
            lo = mid + 1;
        else
            hi = mid;
    }
    for (int k = 0; k < lo; k++)
        buf[k] = key[k];
    buf[lo] = add1;
    for (int k = lo; k < klen; k++)
        buf[k + 1] = key[k];

    return klen + 1;
}

/* drop 'drop1' from the ascending 'key' -- the sorted counterpart of
   paset_without() */
static int
key_without(const int *key, int klen, int drop1, int *buf) {
    int m = 0;
    for (int k = 0; k < klen; k++)
        if (key[k] != drop1)
            buf[m++] = key[k];

    return m;
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
    idl_cache_ref cr = idl_cache_resolve(cached_scores_R, p,
                                        "nh_scores_driver");

    const int *op = INTEGER(op_R);
    const int *uu = INTEGER(u_R);
    const int *vv = INTEGER(v_R);

    SEXP res_R = PROTECT(allocVector(REALSXP, k));
    double *res = REAL(res_R);

    void *vmax0 = vmaxget();

    /* the ascending form of every parent set, built once here into a flat
       arena, so that each candidate's cache key can be derived in O(k) by
       key_with()/key_without() instead of sorting per lookup -- which is
       what build_cache_key() had to do on every single call */
    size_t tot = 0;
    for (int i = 0; i < p; i++)
        tot += (size_t) LENGTH(VECTOR_ELT(pasets_R, i));
    int *sarena = (int *) R_alloc(tot + 1, sizeof(int));
    int *soff = (int *) R_alloc((size_t) p, sizeof(int));
    {
        size_t at = 0;
        for (int i = 0; i < p; i++) {
            SEXP pa_R = VECTOR_ELT(pasets_R, i);
            int lp = LENGTH(pa_R);
            soff[i] = (int) at;
            if (lp > 0) {
                memcpy(sarena + at, INTEGER(pa_R), (size_t) lp * sizeof(int));
                qsort(sarena + at, (size_t) lp, sizeof(int), cmp_int);
            }
            at += (size_t) lp;
        }
    }

    /* per-vertex terms of the current DAG, and their sum in vertex order */
    double *base = (double *) R_alloc((size_t) p, sizeof(double));
    for (int i = 0; i < p; i++) {
        void *vmax = vmaxget();
        SEXP pa_R = VECTOR_ELT(pasets_R, i);
        base[i] = cached_node_score(fn, ctx, &cr, i,
                                    INTEGER(pa_R), LENGTH(pa_R),
                                    sarena + soff[i], LENGTH(pa_R));
        vmaxset(vmax);
    }
    /* scratch for one modified parent set, and for its sorted key; a vertex
       can gain at most one parent over its current set, so p ints suffice */
    int *buf = (int *) R_alloc((size_t) p, sizeof(int));
    int *kbuf = (int *) R_alloc((size_t) p, sizeof(int));

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
            key_with(sarena + soff[v0], lpv, uu[m], kbuf);
            s = cached_node_score(fn, ctx, &cr, v0, buf, lp, kbuf, lp);
            res[m] = total_with(base, p, v0, s, -1, 0.0);
            break;
        case IDLBNS_OP_REMOVE:
            lp = paset_without(pav, lpv, uu[m], buf);
            key_without(sarena + soff[v0], lpv, uu[m], kbuf);
            s = cached_node_score(fn, ctx, &cr, v0, buf, lp, kbuf, lp);
            res[m] = total_with(base, p, v0, s, -1, 0.0);
            break;
        case IDLBNS_OP_REVERSE: {
            /* the head first, then the tail: the cache is order sensitive,
               so the sequence of lookups is part of the contract */
            lp = paset_without(pav, lpv, uu[m], buf);
            key_without(sarena + soff[v0], lpv, uu[m], kbuf);
            double sv = cached_node_score(fn, ctx, &cr, v0, buf, lp, kbuf, lp);
            SEXP pau_R = VECTOR_ELT(pasets_R, u0);
            int lpu = LENGTH(pau_R);
            lp = paset_with(INTEGER(pau_R), lpu, vv[m], buf);
            key_with(sarena + soff[u0], lpu, vv[m], kbuf);
            double su = cached_node_score(fn, ctx, &cr, u0, buf, lp, kbuf, lp);
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
