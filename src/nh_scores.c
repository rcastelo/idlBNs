#include <math.h>
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

/*
 * PER-CANDIDATE SCORING, SHARED BY THE TWO DRIVERS
 *
 * nh_scores_driver() returns every candidate's exact total;
 * nh_argmax_driver() returns only the winner, via the band below. Both must
 * issue exactly the same cache lookups in exactly the same order, because
 * the cache stores whichever order-variant of a parent set missed first --
 * so the per-candidate work lives here, once, rather than being written
 * twice and drifting.
 */
typedef struct {
    SEXP           pasets_R;
    idl_cache_ref *cr;
    node_score_fn  fn;
    void          *ctx;
    const int     *sarena;      /* ascending parent sets, flat            */
    const int     *soff;
    int           *buf;         /* scratch: modified set, scoring order   */
    int           *kbuf;        /* scratch: modified set, ascending        */
    int            p;
} idl_cand_ctx;

/* the vertices a move changes, and their new terms */
typedef struct {
    int    c1;                  /* head, always changed                    */
    double s1;
    int    c2;                  /* tail, changed only by a reversal, or -1 */
    double s2;
} idl_cand;

/* score the one or two vertices whose parent set move (op, u, v) changes.
   'u' and 'v' are 1-based. */
static idl_cand
score_candidate(const idl_cand_ctx *cc, int op, int u, int v) {
    int u0 = u - 1, v0 = v - 1;
    idl_cand out;
    out.c2 = -1;
    out.s2 = 0.0;

    if (u0 < 0 || u0 >= cc->p || v0 < 0 || v0 >= cc->p || u0 == v0)
        error("nh_scores: invalid move (u=%d, v=%d) for p=%d", u, v, cc->p);

    SEXP pav_R = VECTOR_ELT(cc->pasets_R, v0);
    const int *pav = INTEGER(pav_R);
    int lpv = LENGTH(pav_R);
    int lp;

    switch (op) {
    case IDLBNS_OP_ADD:
        lp = paset_with(pav, lpv, u, cc->buf);
        key_with(cc->sarena + cc->soff[v0], lpv, u, cc->kbuf);
        out.c1 = v0;
        out.s1 = cached_node_score(cc->fn, cc->ctx, cc->cr, v0, cc->buf, lp,
                                   cc->kbuf, lp);
        break;
    case IDLBNS_OP_REMOVE:
        lp = paset_without(pav, lpv, u, cc->buf);
        key_without(cc->sarena + cc->soff[v0], lpv, u, cc->kbuf);
        out.c1 = v0;
        out.s1 = cached_node_score(cc->fn, cc->ctx, cc->cr, v0, cc->buf, lp,
                                   cc->kbuf, lp);
        break;
    case IDLBNS_OP_REVERSE: {
        /* the head first, then the tail: the cache is order sensitive, so
           the sequence of lookups is part of the contract */
        lp = paset_without(pav, lpv, u, cc->buf);
        key_without(cc->sarena + cc->soff[v0], lpv, u, cc->kbuf);
        out.c1 = v0;
        out.s1 = cached_node_score(cc->fn, cc->ctx, cc->cr, v0, cc->buf, lp,
                                   cc->kbuf, lp);
        SEXP pau_R = VECTOR_ELT(cc->pasets_R, u0);
        int lpu = LENGTH(pau_R);
        lp = paset_with(INTEGER(pau_R), lpu, v, cc->buf);
        key_with(cc->sarena + cc->soff[u0], lpu, v, cc->kbuf);
        out.c2 = u0;
        out.s2 = cached_node_score(cc->fn, cc->ctx, cc->cr, u0, cc->buf, lp,
                                   cc->kbuf, lp);
        break;
    }
    default:
        error("nh_scores: unknown move operation code %d", op);
    }

    return out;
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

    idl_cand_ctx cc;
    cc.pasets_R = pasets_R; cc.cr = &cr; cc.fn = fn; cc.ctx = ctx;
    cc.sarena = sarena; cc.soff = soff; cc.buf = buf; cc.kbuf = kbuf;
    cc.p = p;

    for (R_xlen_t m = 0; m < k; m++) {
        void *vmax = vmaxget();

        idl_cand c = score_candidate(&cc, op[m], uu[m], vv[m]);
        res[m] = total_with(base, p, c.c1, c.s1, c.c2, c.s2);

        vmaxset(vmax); /* reclaim this move's cache-key and score scratch */
    }

    vmaxset(vmax0);
    UNPROTECT(1);

    return res_R;
}

/*
 * THE ARGMAX BAND: FINDING THE WINNER WITHOUT SUMMING EVERY CANDIDATE
 *
 * total_with() sums all p per-vertex terms for one candidate, and that is
 * what makes each candidate's total bit-identical to a full re-score. But
 * the neighbourhood holds O(p^2) candidates, so a step costs O(p^3)
 * additions: 1.0M at p = 100, 8.0M at p = 200, 125M at p = 500 -- about
 * 125 ms per step there, which would dominate everything else the port
 * saves.
 *
 * There is no exact shortcut. Recursive summation is not decomposable:
 * prefix sums let the fold start at the changed position (a factor of two)
 * but the tail depends on the running accumulator, so it stays O(p).
 *
 * So filter instead. For each candidate form the cheap estimate
 *
 *     est_m = fl(T_base + fl(delta1 + delta2))
 *
 * bound the error |T_m - est_m| rigorously, and compute the exact
 * total_with() only for the candidates whose error interval could still
 * contain the maximum. Everything else is PROVABLY worse, so the winner is
 * exactly the winner a full exact ranking would pick -- including its
 * position, since the band is scanned in ascending m with a strict '>',
 * which is what reproduces which.max()'s first-index tie-breaking.
 *
 * THE BOUND. With u = 2^-53 the unit roundoff and gamma_n = n*u/(1 - n*u),
 * recursive summation of n doubles satisfies |fl(sum) - sum| <= gamma_{n-1}
 * * sum|x_i|. Three error sources have to appear, and a bound that keeps
 * only the first is not a bound:
 *
 *   1. T_base and T_m are BOTH rounded sums, so both contribute; and the
 *      term-magnitude sum for candidate m is A_m, not A, because one or two
 *      terms have been swapped.
 *   2. the deltas are themselves rounded subtractions and one rounded
 *      addition.
 *   3. fl(T_base + d_m) rounds once more.
 *
 * giving
 *
 *   A'   = A * (1 + 2*p*u)                       covers A's own summation
 *   A_m  = A' - |e1| + |s1| - |e2| + |s2| + 4*u*A'
 *   B_m  = gamma_{p-1} * (A_m + A') + 2*u*(|d1| + |d2|) + u*|est_m|
 *
 * Per-candidate intervals are then strictly tighter than one global width,
 * and cost nothing extra:
 *
 *   L_m = nextafter(est_m - SAFETY*B_m, -inf)
 *   U_m = nextafter(est_m + SAFETY*B_m, +inf)
 *   M   = max_m L_m,   band = { m : U_m >= M }
 *
 * If m is outside the band then T_m <= U_m < M <= L_{m*} <= T_{m*} for the
 * m* attaining M, and m* is itself in the band since U >= L. So the true
 * maximum lies inside. SAFETY = 4 plus the outward nextafter covers the
 * rounding of the bound arithmetic itself; the band size is insensitive to
 * SAFETY over 1..16, so being generous costs nothing.
 *
 * At p = 500 with per-vertex terms of order 1e3 the band is about 5e-8
 * wide, against genuine inter-move score gaps of order 1e-1 to 1e2. What
 * lands inside it are the exactly-tied score-equivalence classes -- both
 * iBIC and iBGe score Markov equivalent DAGs alike on observational data --
 * so expect a handful, and the O(p^3) term collapses to O(p * |band|).
 *
 * TWO THINGS THIS MUST NOT DO. It must not change the cache-lookup sequence
 * -- every candidate still gets its node scores computed, in the same order,
 * through the same score_candidate() the exact driver uses; only the
 * p-term summation is skipped. And it must not be used where a caller wants
 * every candidate's total: tests/test_delta_scores.R asserts identical() on
 * the whole vector nh_scores_driver() returns, so that function is left
 * exactly as it was and this is a separate entry point.
 *
 * There is no correctness cliff. If the band degenerates to all k
 * candidates the answer is still exactly right and the cost is exactly what
 * it was before; only the speedup is lost. 'band' is returned so that can
 * be watched.
 */

#define IDL_UROUND  0x1p-53
#define IDL_SAFETY  4.0

static inline double
idl_gamma(double n) {
    double x = n * IDL_UROUND;

    return x / (1.0 - x);
}

SEXP
nh_argmax_driver(SEXP pasets_R, SEXP cached_scores_R, SEXP op_R, SEXP u_R,
                 SEXP v_R, SEXP stamp_R, int verify, node_score_fn fn,
                 void *ctx) {
    if (TYPEOF(pasets_R) != VECSXP)
        error("nh_argmax_driver: 'pasets' must be a list");
    if (TYPEOF(op_R) != INTSXP || TYPEOF(u_R) != INTSXP ||
        TYPEOF(v_R) != INTSXP)
        error("nh_argmax_driver: 'op', 'u' and 'v' must be integer vectors");

    int p = LENGTH(pasets_R);
    R_xlen_t k = XLENGTH(op_R);

    if (XLENGTH(u_R) != k || XLENGTH(v_R) != k)
        error("nh_argmax_driver: 'op', 'u' and 'v' must have the same length");
    if (k == 0)
        error("nh_argmax_driver: the neighbourhood is empty");
    idl_cache_ref cr = idl_cache_resolve(cached_scores_R, p,
                                        "nh_argmax_driver");

    /*
     * The addition memo is used only when the caller supplies the DAG's
     * per-vertex parent-set stamps AND there is a compiled cache to hold it.
     * Without stamps there is no way to tell a stale entry from a fresh one,
     * so the memo stays off and every candidate goes through the cache --
     * which is what direct callers and the differential tests do.
     */
    const int *stamp = NULL;
    idl_sc_cache *memo = NULL;
    if (stamp_R != R_NilValue) {
        if (TYPEOF(stamp_R) != INTSXP || LENGTH(stamp_R) != p)
            error("nh_argmax_driver: 'stamp' must be NULL or an integer vector of length %d",
                  p);
        if (cr.kind == IDL_CACHE_HASH) {
            stamp = INTEGER(stamp_R);
            memo = cr.hash;
            idl_sc_memo_enable(memo);
        }
    }

    const int *op = INTEGER(op_R);
    const int *uu = INTEGER(u_R);
    const int *vv = INTEGER(v_R);

    void *vmax0 = vmaxget();

    /* the ascending form of every parent set, exactly as the exact driver
       builds it, so the cache keys are identical */
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

    /* the current DAG's per-vertex terms, and their sum in vertex order */
    double *base = (double *) R_alloc((size_t) p, sizeof(double));
    for (int i = 0; i < p; i++) {
        void *vm = vmaxget();
        SEXP pa_R = VECTOR_ELT(pasets_R, i);
        base[i] = cached_node_score(fn, ctx, &cr, i, INTEGER(pa_R),
                                    LENGTH(pa_R), sarena + soff[i],
                                    LENGTH(pa_R));
        vmaxset(vm);
    }
    double total = 0.0, absA = 0.0;
    int finite_base = 1;
    for (int i = 0; i < p; i++) {
        total += base[i];
        absA += fabs(base[i]);
        if (!R_FINITE(base[i]))
            finite_base = 0;
    }

    int *buf = (int *) R_alloc((size_t) p, sizeof(int));
    int *kbuf = (int *) R_alloc((size_t) p, sizeof(int));
    idl_cand_ctx cc;
    cc.pasets_R = pasets_R; cc.cr = &cr; cc.fn = fn; cc.ctx = ctx;
    cc.sarena = sarena; cc.soff = soff; cc.buf = buf; cc.kbuf = kbuf;
    cc.p = p;

    /*
     * MEMORY TRAFFIC. An earlier version of this kept est, lo, hi, s1, s2,
     * c1 and c2 for every candidate -- 48 bytes each, 7.5 MB per call at
     * p = 400 -- when pass 2 needs the per-candidate terms only for the
     * BAND, which is normally a single candidate. So pass 1 now keeps just
     * one array, the upper interval end, and pass 2 re-scores the handful of
     * band members.
     *
     * Re-scoring is contents-neutral for the cache: the keys are the ones
     * pass 1 just looked up, so every one of them is a hit, and a hit
     * neither stores anything nor changes what is stored. It costs a few
     * extra lookups against 6 MB less memory written per call.
     */
    double *hi = (double *) R_alloc((size_t) k, sizeof(double));
    /* est is needed only to check the bound in verify mode */
    double *est = verify
        ? (double *) R_alloc((size_t) k, sizeof(double))
        : NULL;

    double Ap = absA * (1.0 + 2.0 * (double) p * IDL_UROUND);
    double gam = idl_gamma((double) (p > 1 ? p - 1 : 1));
    int all_finite = finite_base;
    double M = R_NegInf;

    /* pass 1: score every candidate -- the same lookups, in the same order,
       as the exact driver issues -- and bound its estimate */
    for (R_xlen_t m = 0; m < k; m++) {
        double d1, d2 = 0.0;

        if (memo != NULL && op[m] == IDLBNS_OP_ADD) {
            /* additions are ~99% of the candidates, and an addition's delta
               is a function of pa(v) alone -- through both s(v, pa(v) + u)
               and base[v] -- so one stamp compare decides validity */
            size_t slot = (size_t) (uu[m] - 1) * (size_t) p
                          + (size_t) (vv[m] - 1);
            if (memo->add_st[slot] == stamp[vv[m] - 1]) {
                d1 = memo->add_d[slot];
                memo->memo_hits++;
            } else {
                void *vm = vmaxget();
                idl_cand c = score_candidate(&cc, op[m], uu[m], vv[m]);
                vmaxset(vm);
                d1 = c.s1 - base[c.c1];
                memo->add_d[slot] = d1;
                memo->add_st[slot] = stamp[vv[m] - 1];
                memo->memo_misses++;
            }
        } else {
            void *vm = vmaxget();
            idl_cand c = score_candidate(&cc, op[m], uu[m], vv[m]);
            vmaxset(vm);
            d1 = c.s1 - base[c.c1];
            if (c.c2 >= 0)
                d2 = c.s2 - base[c.c2];
        }

        double em = total + (d1 + d2);

        /*
         * The magnitude bound uses |e| + |d| in place of |s|, which is valid
         * by the triangle inequality since s = e + d, and lets the bound be
         * computed from the delta alone -- a memo hit never reconstructs s.
         * It is marginally looser, so the band may be a candidate or two
         * wider; that cannot change the answer, because the true maximum is
         * inside either way and pass 2 takes the first candidate attaining
         * it in ascending order.
         */
        double Am = Ap + fabs(d1) + fabs(d2);
        Am += 4.0 * IDL_UROUND * Ap;

        double B = gam * (Am + Ap) + 2.0 * IDL_UROUND * (fabs(d1) + fabs(d2))
                   + IDL_UROUND * fabs(em);
        double w = IDL_SAFETY * B;

        /*
         * The interval ends are formed without nextafter(), which is a libm
         * call and was measurable at this candidate count. Rounding them to
         * nearest is safe because SAFETY = 4 leaves room: B already
         * contains the term u*|est|, so the rounding error in forming
         * est +/- w is at most about u*|est| <= B = w/4. Hence the computed
         * hi is at least est + 0.75w, comfortably above the est + B the
         * proof needs, and the computed lo is at most est - 0.75w, which
         * only ever makes M smaller and the band larger.
         */
        hi[m] = em + w;
        double lom = em - w;
        if (lom > M)
            M = lom;
        if (verify)
            est[m] = em;

        if (!R_FINITE(em) || !R_FINITE(B))
            all_finite = 0;
    }

    /* if anything is non-finite the interval arithmetic is meaningless, so
       score every candidate exactly instead. unreachable in practice -- the
       node score functions error on a failed Cholesky rather than returning
       a NaN -- but a guard against a silent wrong answer is worth having */
    R_xlen_t bidx = 0;
    double btot = 0.0;
    R_xlen_t nband = 0;
    double worst = 0.0;
    int first = 1;

    for (R_xlen_t m = 0; m < k; m++) {
        if (all_finite && hi[m] < M)
            continue;                           /* provably not the maximum */
        void *vm = vmaxget();
        idl_cand c = score_candidate(&cc, op[m], uu[m], vv[m]);
        vmaxset(vm);
        double t = total_with(base, p, c.c1, c.s1, c.c2, c.s2);
        nband++;
        /* ascending m with a strict '>' is what reproduces which.max() */
        if (first || t > btot) { btot = t; bidx = m; first = 0; }
    }
    if (first)                                  /* cannot happen: M = lo[j]
                                                   implies hi[j] >= M */
        error("nh_argmax_driver: the candidate band came out empty");

    /* verify mode: score every candidate exactly and check the band */
    if (verify) {
        double bestt = 0.0;
        R_xlen_t besti = 0;
        for (R_xlen_t m = 0; m < k; m++) {
            void *vm = vmaxget();
            idl_cand c = score_candidate(&cc, op[m], uu[m], vv[m]);
            vmaxset(vm);
            double t = total_with(base, p, c.c1, c.s1, c.c2, c.s2);
            if (m == 0 || t > bestt) { bestt = t; besti = m; }
            if (all_finite) {
                double B = (hi[m] - est[m]) / IDL_SAFETY;
                double dev = fabs(t - est[m]);
                if (B > 0.0 && dev / B > worst)
                    worst = dev / B;
                if (dev > IDL_SAFETY * B)
                    error("nh_argmax_driver: bound violated at candidate %lld: |T - est| = %g > %g",
                          (long long) (m + 1), dev, IDL_SAFETY * B);
            }
        }
        if (besti != bidx || bestt != btot)
            error("nh_argmax_driver: band picked candidate %lld (%.17g) but the exact maximum is %lld (%.17g)",
                  (long long) (bidx + 1), btot, (long long) (besti + 1), bestt);
    }

    SEXP ans = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, ScalarInteger((int) (bidx + 1)));   /* 1-based */
    SET_VECTOR_ELT(ans, 1, ScalarReal(btot));
    SET_VECTOR_ELT(ans, 2, ScalarReal((double) nband));
    SET_VECTOR_ELT(ans, 3, ScalarReal(worst));
    SEXP nms = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, mkChar("index"));
    SET_STRING_ELT(nms, 1, mkChar("total"));
    SET_STRING_ELT(nms, 2, mkChar("band"));
    SET_STRING_ELT(nms, 3, mkChar("worst"));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    vmaxset(vmax0);

    return ans;
}
