#ifndef IDLBNS_DAG_H
#define IDLBNS_DAG_H

#include <stddef.h>
#include "bitset.h"

/*
 * THE DAG STATE
 *
 * A C mirror of the four structures the R search maintained separately: the
 * graphNEL itself, the ancestor matrix ('anc' in R/search.R), the parent
 * sets ('pasets'), and implicitly the child lists that edgeL() served.
 *
 * Two representations, chosen per use:
 *
 *   reachability (anc, desc) -- dense, needs O(1) membership for the
 *       acyclicity test, so multi-word bitsets;
 *   parent and child lists   -- sparse, cost must follow the set size and
 *       not p, so int arrays.
 *
 * THREE ORDERING INVARIANTS, each of which is load-bearing for
 * bit-exactness. Getting any of them wrong changes search trajectories
 * while passing most tests:
 *
 *   pa[v]  parents of v in INSERTION order, exactly as add.pasets() and
 *          remove.pasets() maintain pasets[[v]]: append on add, and an
 *          order-preserving delete on remove. This order reaches the score
 *          functions, and the Cholesky of a symmetrically permuted ZtZ
 *          rounds differently -- measured 1.14e-13 on an iBIC total for
 *          parents c(1,2,3,4,5) against c(5,4,3,2,1).
 *
 *   pas[v] the same parents ASCENDING. This is the canonical form for the
 *          score cache key, and keeping it incrementally (binary-searched
 *          insertion point plus a memmove, O(k)) is what removes the qsort
 *          that build_cache_key() currently pays on every lookup.
 *
 *   ch[u]  children of u in graphNEL edgeL order, i.e. append on add and an
 *          order-preserving delete on remove. nr.nh() emits a vertex's
 *          removals in this order and edgeMatrix() enumerates edges in it,
 *          so it decides both which.max() tie-breaking and which covered
 *          arc rcar() picks for a given random draw.
 *
 * anc[v] (the v-th row of 'anc', W words) is the set of ANCESTORS of v,
 * which is column v of the R matrix -- R indexes anc[u, v] as "u is an
 * ancestor of v". desc[v] is the set of DESCENDANTS of v, i.e. row v of
 * the R matrix. Both are stored: desc could be recovered by scanning bit v
 * out of all p ancestor rows, but storing it costs 32 KB at p = 500, makes
 * the mutators O(|D| * W) instead of O(p * W), and gives a free
 * self-consistency invariant -- desc must be the exact transpose of anc --
 * which is precisely the class of bug that would silently change the legal
 * move set.
 */

/* growable int vector */
typedef struct {
    int *v;
    int  n;
    int  cap;
} idl_ivec;

typedef struct idl_dag {
    int       p;              /* number of vertices                        */
    size_t    W;              /* words per bitset row = IDL_NW(p)          */
    int       nedges;

    idl_ivec *pa;             /* p lists: parents, INSERTION order         */
    idl_ivec *pas;            /* p lists: parents, ASCENDING               */
    idl_ivec *ch;             /* p lists: children, edgeL order            */

    idl_word *pab;            /* p * W: pab + v * W = parents of v         */
    idl_word *adj;            /* p * W: adj + u * W = children of u        */
    idl_word *anc;            /* p * W: anc + v * W = ancestors of v       */
    idl_word *desc;           /* p * W: desc + v * W = descendants of v    */

    /* scratch, sized once at construction and reused by every mutation */
    idl_word *sc_delta;       /* W                                         */
    idl_word *sc_set;         /* W                                         */
    idl_word *sc_ancsave;     /* p * W: ancestor rows before a removal     */
    int      *sc_a;           /* p                                         */
    int      *sc_b;           /* p                                         */
    int      *sc_order;       /* p: topological order of the touched set   */
    int      *sc_indeg;       /* p                                         */
    int      *sc_queue;       /* p                                         */
    int      *sc_inD;         /* p: 0/1 membership flags                   */
} idl_dag;

/* an upper bound on p, purely to turn nonsense input into an error */
#define IDLBNS_PMAX 65536

/*
 * Lifecycle. Construction is deliberately in TWO steps, because the
 * allocators are R's (R_Calloc, per WRE 1.6.4 and 6.1.2) and those raise an
 * R error rather than returning NULL -- they longjmp.
 *
 *     idl_dag *d = idl_dag_alloc(p);        // the zeroed struct only
 *     ... wrap d in an external pointer ...
 *     ... register idl_dag_free as its finalizer ...
 *     idl_dag_init(d);                      // the interior; may longjmp
 *
 * Doing it in that order means a failure inside idl_dag_init() unwinds into
 * an object R already owns and will reclaim. Every interior pointer is NULL
 * until it is live, so the finalizer copes with a partially built DAG.
 * Calling idl_dag_init() on a struct that is not yet owned would leak
 * whatever it had taken.
 */
idl_dag *idl_dag_alloc(int p);
void     idl_dag_init(idl_dag *d);
void     idl_dag_free(idl_dag *d);

/* queries */
static inline int
idl_dag_is_anc(const idl_dag *d, int u, int v) {   /* R's anc[u, v] */
    return idl_bs_test(d->anc + (size_t) v * d->W, u);
}

static inline int
idl_dag_has_edge(const idl_dag *d, int u, int v) {
    return idl_bs_test(d->adj + (size_t) u * d->W, v);
}

/* would this move keep the graph acyclic? */
int idl_dag_can_add(const idl_dag *d, int u, int v);
int idl_dag_can_reverse(const idl_dag *d, int u, int v);

/* mutations. Preconditions are the caller's responsibility -- dag_R.c
   validates everything before calling, so that an error can never leave a
   half-applied move behind. */
void idl_dag_add_edge(idl_dag *d, int u, int v);
void idl_dag_remove_edge(idl_dag *d, int u, int v);
void idl_dag_reverse_edge(idl_dag *d, int u, int v);

/*
 * idl_dag_check
 *
 * Full internal self-consistency pass, for the debug options and the
 * differential tests. Returns 0 when consistent, or a non-zero code, and
 * writes a description into msg (of size msglen). O(p^2 / 64 + p * E).
 */
int idl_dag_check(const idl_dag *d, char *msg, size_t msglen);

#endif /* IDLBNS_DAG_H */
