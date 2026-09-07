#include <R.h>
#include <Rinternals.h>
#include "dag.h"
#include "dag_R.h"
#include "nh_scores.h"          /* IDLBNS_OP_ADD / _REMOVE / _REVERSE */

/*
 * NEIGHBOURHOOD ENUMERATION
 *
 * A C port of nr.nh(), ar.nh() and ncr.nh() from R/search.R. The result is
 * the same shape those return: list(op=, u=, v=), three parallel integer
 * vectors of 1-based vertex indices, where for OP_REVERSE the pair names the
 * arc as it stands BEFORE the reversal.
 *
 * EMISSION ORDER IS PART OF THE SPECIFICATION, not an implementation detail.
 * hcmc()/hillclimbing() pick a move with which.max(), which breaks ties by
 * position, so two enumerations that produce the same move SET in a
 * different order send the search into different local optima. The order is:
 *
 *   for i = 1..p ascending:
 *       additions i -> w, for w ascending
 *       removals  i -> w, for w in the order the GRAPH stores i's children
 *   then, for the two reversal neighbourhoods, for i = 1..p ascending:
 *       reversals i -> w, again in stored child order
 *
 * The child order is edge-insertion order, because graph::addEdge appends
 * and graph::removeEdge compacts. It is routinely NOT ascending once a
 * search has done any removals or reversals, which is exactly why
 * tests/test_c_dag.R and tests/test_c_nh.R count how many non-ascending
 * child lists they saw: an implementation that emitted removals ascending
 * would agree on every move set and still change every trajectory.
 *
 * Legality:
 *   addition  i -> w is offered iff w is not i, not already a child of i,
 *       and not an ancestor of i (the last is nr.nh()'s !anc[na, i], and it
 *       subsumes "not a parent of i", since a parent is an ancestor).
 *   reversal  i -> w is offered iff no OTHER child of i is an ancestor of w
 *       (.reversible()), and, for ncr, iff the arc is not I-covered.
 *
 * An arc i -> w is covered iff pa(i) == pa(w) \ {i}, and I-covered iff it is
 * covered AND neither endpoint is an intervention target -- so a covered arc
 * touching a target stays in the neighbourhood. Because the DAG keeps an
 * ascending mirror of every parent set, the covered test is a linear merge
 * over two sorted int arrays, O(|pa(i)| + |pa(w)|), with no sorting and no
 * dependence on p.
 */

/* is the arc i -> w covered, i.e. pa(i) == pa(w) \ {i}?
   both lists are ascending and duplicate-free, and i is necessarily in
   pa(w), so the sizes must differ by exactly one. */
static int
arc_is_covered(const idl_dag *d, int i, int w) {
    const idl_ivec *pi = &d->pas[i];
    const idl_ivec *pw = &d->pas[w];

    if (pi->n != pw->n - 1)
        return 0;

    int a = 0;
    for (int b = 0; b < pw->n; b++) {
        int x = pw->v[b];
        if (x == i)
            continue;                       /* the setdiff */
        if (a >= pi->n || pi->v[a] != x)
            return 0;
        a++;
    }

    return a == pi->n;
}

/*
 * C_dag_nh
 *
 * st          externalptr  the DAG state
 * kind_R      INTSXP       1 = nr, 2 = ar, 3 = ncr
 * utargets_R  INTSXP       unique intervention target vertices, 1-based;
 *                          used only by kind 3. Values outside 1..p are
 *                          IGNORED rather than rejected, matching the R
 *                          code: utargets comes from unlist(targets), and
 *                          the documented examples include targets such as
 *                          list(0L, 2L), whose 0 simply never matches a
 *                          vertex index.
 *
 * Returns list(op=, u=, v=).
 */
SEXP
C_dag_nh(SEXP st, SEXP kind_R, SEXP utargets_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    int kind = asInteger(kind_R);

    if (kind != 1 && kind != 2 && kind != 3)
        error("C_dag_nh: 'kind' must be 1 (nr), 2 (ar) or 3 (ncr)");
    if (TYPEOF(utargets_R) != INTSXP)
        error("C_dag_nh: 'utargets' must be an integer vector");

    int p = d->p;
    size_t W = d->W;
    void *vmax = vmaxget();

    /* target flags, ignoring anything outside 1..p (see above) */
    char *istgt = NULL;
    if (kind == 3) {
        istgt = (char *) R_alloc((size_t) p, sizeof(char));
        memset(istgt, 0, (size_t) p);
        const int *ut = INTEGER(utargets_R);
        R_xlen_t nut = XLENGTH(utargets_R);
        for (R_xlen_t k = 0; k < nut; k++)
            if (ut[k] != NA_INTEGER && ut[k] >= 1 && ut[k] <= p)
                istgt[ut[k] - 1] = 1;
    }

    /* an exact upper bound: each vertex offers at most p-1 additions plus
       its children as removals, and the reversal block adds at most one
       entry per arc */
    size_t kmax = (size_t) p * (size_t) (p - 1) + (size_t) d->nedges;
    if (kmax == 0)
        kmax = 1;
    int *bop = (int *) R_alloc(kmax, sizeof(int));
    int *bu  = (int *) R_alloc(kmax, sizeof(int));
    int *bv  = (int *) R_alloc(kmax, sizeof(int));
    int *cand = (int *) R_alloc((size_t) p, sizeof(int));
    size_t n = 0;

    /* ---- the NR block: additions then removals, per vertex ---- */
    for (int i = 0; i < p; i++) {
        /* mask off everything i may not point at: itself, its children, and
           its ancestors. what remains, ascending, are the legal additions. */
        idl_bs_copy(d->sc_set, d->adj + (size_t) i * W, W);
        idl_bs_or(d->sc_set, d->anc + (size_t) i * W, W);
        idl_bs_set(d->sc_set, i);
        int nc = idl_bs_absent_to_ints(d->sc_set, W, p, cand);
        for (int k = 0; k < nc; k++) {
            bop[n] = IDLBNS_OP_ADD;
            bu[n] = i + 1;
            bv[n] = cand[k] + 1;
            n++;
        }
        /* removals, in stored child order */
        const idl_ivec *chi = &d->ch[i];
        for (int j = 0; j < chi->n; j++) {
            bop[n] = IDLBNS_OP_REMOVE;
            bu[n] = i + 1;
            bv[n] = chi->v[j] + 1;
            n++;
        }
    }

    /* ---- the reversal block, appended after the whole NR block ---- */
    if (kind != 1) {
        for (int i = 0; i < p; i++) {
            const idl_ivec *chi = &d->ch[i];
            for (int j = 0; j < chi->n; j++) {
                int w = chi->v[j];
                if (kind == 3) {
                    int touches = istgt[i] || istgt[w];
                    if (!touches && arc_is_covered(d, i, w))
                        continue;           /* I-covered: not in NCR */
                }
                if (!idl_dag_can_reverse(d, i, w))
                    continue;
                bop[n] = IDLBNS_OP_REVERSE;
                bu[n] = i + 1;
                bv[n] = w + 1;
                n++;
            }
        }
    }

    SEXP ans = PROTECT(allocVector(VECSXP, 3));
    SEXP op = PROTECT(allocVector(INTSXP, (R_xlen_t) n));
    SEXP u  = PROTECT(allocVector(INTSXP, (R_xlen_t) n));
    SEXP v  = PROTECT(allocVector(INTSXP, (R_xlen_t) n));
    memcpy(INTEGER(op), bop, n * sizeof(int));
    memcpy(INTEGER(u),  bu,  n * sizeof(int));
    memcpy(INTEGER(v),  bv,  n * sizeof(int));
    SET_VECTOR_ELT(ans, 0, op);
    SET_VECTOR_ELT(ans, 1, u);
    SET_VECTOR_ELT(ans, 2, v);
    SEXP nms = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nms, 0, mkChar("op"));
    SET_STRING_ELT(nms, 1, mkChar("u"));
    SET_STRING_ELT(nms, 2, mkChar("v"));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(5);
    vmaxset(vmax);

    return ans;
}

/* the covered-arc mask, in edgeMatrix() column order -- the C counterpart of
   cedges(). rcar() selects a covered arc by INDEX into this order, so the
   order is load-bearing. Needed in full for Step 4; exposed now because the
   same arc_is_covered() drives it and the differential test can pin both at
   once. */
SEXP
C_dag_cedges(SEXP st, SEXP utargets_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);

    if (TYPEOF(utargets_R) != INTSXP)
        error("C_dag_cedges: 'utargets' must be an integer vector");

    int p = d->p;
    void *vmax = vmaxget();
    char *istgt = (char *) R_alloc((size_t) p, sizeof(char));
    memset(istgt, 0, (size_t) p);
    const int *ut = INTEGER(utargets_R);
    R_xlen_t nut = XLENGTH(utargets_R);
    for (R_xlen_t k = 0; k < nut; k++)
        if (ut[k] != NA_INTEGER && ut[k] >= 1 && ut[k] <= p)
            istgt[ut[k] - 1] = 1;

    SEXP ans = PROTECT(allocVector(LGLSXP, d->nedges));
    int *a = LOGICAL(ans);
    int m = 0;
    for (int i = 0; i < p; i++) {
        const idl_ivec *chi = &d->ch[i];
        for (int j = 0; j < chi->n; j++) {
            int w = chi->v[j];
            a[m++] = arc_is_covered(d, i, w) && !(istgt[i] || istgt[w]);
        }
    }
    UNPROTECT(1);
    vmaxset(vmax);

    return ans;
}
