#ifndef IDLBNS_NH_SCORES_H
#define IDLBNS_NH_SCORES_H

#include <R.h>
#include <Rinternals.h>

/* move operation codes, mirroring OP.ADD/OP.REMOVE/OP.REVERSE in
   R/search.R; a neighborhood is three parallel vectors (op, u, v) */
#define IDLBNS_OP_ADD     1
#define IDLBNS_OP_REMOVE  2
#define IDLBNS_OP_REVERSE 3

/*
 * node_score_fn
 *
 * Computes the score contribution of one node for a given parent set,
 * for whichever score function is driving the search. 'ctx' carries that
 * score function's global sufficient statistics (see the iBIC_ctx /
 * iBGe_ctx structs in ibic_score.c / ibge_score.c). 'pa' holds lp 1-based
 * parent variable indices and 'node' is a 1-based node index.
 */
typedef double (*node_score_fn)(void *ctx, const int *pa, int lp, int node);

/*
 * nh_scores_driver
 *
 * Scores a whole neighborhood of candidate moves in one call, returning a
 * REALSXP of length k = LENGTH(op_R) holding the total DAG score of each
 * candidate. See nh_scores.c for the delta-scoring rationale.
 */
SEXP
nh_scores_driver(SEXP pasets_R, SEXP cached_scores_R, SEXP op_R, SEXP u_R,
                 SEXP v_R, node_score_fn fn, void *ctx);

/*
 * nh_argmax_driver
 *
 * Returns only the WINNING candidate -- list(index=, total=, band=, worst=)
 * -- rather than all k totals, so that the O(p) exact summation per
 * candidate can be skipped for all but a provable handful. See nh_scores.c.
 *
 * verify != 0 additionally computes every candidate's exact total and checks
 * the band against it, which is how tests/test_c_band.R pins the error
 * bound. It is O(p * k) and so is only for testing.
 */
SEXP
nh_argmax_driver(SEXP pasets_R, SEXP cached_scores_R, SEXP op_R, SEXP u_R,
                 SEXP v_R, int verify, node_score_fn fn, void *ctx);

#endif /* IDLBNS_NH_SCORES_H */
