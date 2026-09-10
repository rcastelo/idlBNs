#ifndef IDLBNS_CLIQUEPICK_H
#define IDLBNS_CLIQUEPICK_H

#include <stdint.h>

/* vertex sets inside one chain component; components are capped at 64 vertices */
typedef uint64_t cp_vset;

/*
 * One memoised subproblem: a vertex subset of the component, its AMO count,
 * and -- for sampling -- the per-clique weights, forbidden-prefix chains and
 * child subproblem ids that the count was assembled from.
 */
typedef struct {
    cp_vset  key;
    double   total;
    int      ncl;
    cp_vset *cl;        /* ncl clique vertex sets                          */
    double  *w;         /* ncl weights, summing to total                   */
    int     *nsp;       /* ncl subproblem counts                           */
    int    **sp;        /* ncl arrays of subproblem node ids               */
    int     *nfp;       /* ncl forbidden-prefix chain lengths              */
    cp_vset **fp;       /* ncl arrays of nested prefix sets, ascending     */
} cp_node;

typedef struct {
    const cp_vset *adj;
    int      m;
    int      nnode;
    int      cap;
    cp_node *node;
    int      ok;        /* 0 once anything overflowed or exceeded a bound  */
} cp_ctx;

/* Build the memo for the component `uni`; returns its node id, or -1 on
   failure. Allocation is R_alloc'd, so nothing needs freeing. */
int    cp_build(cp_ctx *ctx, const cp_vset *adj, int m, cp_vset uni);
double cp_count(cp_ctx *ctx, int id);

/*
 * Draw a uniform AMO of the component rooted at `id`, writing a topological
 * order into ord[] starting at position *tick. Consumes R's RNG; the caller
 * brackets it with Get/PutRNGstate.
 */
void   cp_sample(cp_ctx *ctx, int id, int *ord, int *tick);

/*
 * The k-th member of the class rooted at `id`, 0-based, as a topological order
 * written into ord[] from *tick. Enumerating a class of c members costs c
 * calls, i.e. time proportional to the output. Valid for k < cp_count(ctx, id).
 */
void   cp_member(cp_ctx *ctx, int id, double k, int *ord, int *tick);

#endif /* IDLBNS_CLIQUEPICK_H */
