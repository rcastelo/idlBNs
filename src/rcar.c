#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include "dag.h"
#include "dag_R.h"

/*
 * REPEATED COVERED ARC REVERSAL, IN C
 *
 * A port of rcar() from R/search.R. It is the only part of the search that
 * consumes randomness, so it is the only part where a C port can silently
 * desynchronise the RNG stream and move every trajectory. Two R expressions
 * have to be reproduced draw for draw:
 *
 *     rr    <- sample(0:r, size = 1)
 *     rndce <- resample(which(cemask), size = 1)
 *
 * and tests/test_rng_equivalence.R pins the identities that let C do it:
 *
 *     sample(0:r, size = 1)  ==  (int) R_unif_index((double) length(0:r))
 *     resample(x, size = 1)  ==  x[(int) R_unif_index((double) length(x)) + 1]
 *
 * R_unif_index() is called rather than reimplemented: under the default
 * sample.kind = "Rejection" the number of unif_rand() calls per draw is data
 * dependent (1 to 3 for a size-1 draw), and R_unif_index() also dispatches
 * on R_sample_kind() so this stays correct under "Rounding".
 *
 * FOUR ORDERING RULES, each a correctness requirement:
 *
 *   1. Both early returns sit ABOVE GetRNGstate(). For stream VALUES it
 *      would not matter, since Get/Put with no draw in between round-trips
 *      the state. But in a session that has not yet used the RNG, an
 *      unconditional PutRNGstate() would CREATE .Random.seed from the clock
 *      -- and the first hcmc() iteration always takes the nedges == 0 branch,
 *      so this would fire on every single run. tests/test_rng_equivalence.R
 *      checks in a subprocess that a no-draw rcar leaves .Random.seed
 *      non-existent.
 *
 *   2. One Get/Put pair brackets the whole call rather than each draw. That
 *      is stream-identical, because Get/Put only synchronise the C generator
 *      with .Random.seed and nothing in between inspects it.
 *
 *   3. Nothing that can longjmp happens inside the region without closing it
 *      first. All allocation is done before GetRNGstate().
 *
 *   4. The covered set is recomputed from scratch on every iteration, as the
 *      R loop does (cemask <- cedges(tmp.g, utargets) inside the loop), and
 *      in edgeMatrix() column order -- vertex ascending, children in stored
 *      insertion order -- because the draw is an INDEX into that order.
 *      Enumerating the arcs any other way would pick a different arc for the
 *      same random number.
 */

/* collect the I-covered arcs in edgeMatrix() column order; returns how many */
static int
covered_arcs(const idl_dag *d, const char *istgt, int *from, int *to) {
    int n = 0;
    for (int i = 0; i < d->p; i++) {
        const idl_ivec *chi = &d->ch[i];
        for (int j = 0; j < chi->n; j++) {
            int w = chi->v[j];
            if (istgt != NULL && (istgt[i] || istgt[w]))
                continue;                   /* touches a target: not I-covered */
            if (!idl_dag_arc_is_covered(d, i, w))
                continue;
            from[n] = i;
            to[n] = w;
            n++;
        }
    }

    return n;
}

/*
 * C_dag_rcar
 *
 * st          externalptr  the DAG state, mutated in place
 * rlen_R      INTSXP       scalar length(0:r), computed on the R side so
 *                          that the coercion corner cases of 0:r for
 *                          non-integer or negative r stay where they already
 *                          behave correctly
 * utargets_R  INTSXP       unique target vertices, 1-based; values outside
 *                          1..p are ignored, matching R (see nbhd.c)
 *
 * Returns the number of reversals actually performed, as R's rr.
 */
SEXP
C_dag_rcar(SEXP st, SEXP rlen_R, SEXP utargets_R) {
    idl_dag *d = idlBNs_dag_from_extptr(st);
    int rlen = asInteger(rlen_R);

    if (rlen == NA_INTEGER || rlen < 1)
        error("C_dag_rcar: 'rlen' must be a positive integer (length(0:r))");
    if (TYPEOF(utargets_R) != INTSXP)
        error("C_dag_rcar: 'utargets' must be an integer vector");

    int p = d->p;
    void *vmax = vmaxget();

    /* --- everything that can allocate or fail, before any RNG contact --- */
    char *istgt = NULL;
    R_xlen_t nut = XLENGTH(utargets_R);
    if (nut > 0) {
        const int *ut = INTEGER(utargets_R);
        istgt = (char *) R_alloc((size_t) p, sizeof(char));
        memset(istgt, 0, (size_t) p);
        for (R_xlen_t k = 0; k < nut; k++)
            if (ut[k] != NA_INTEGER && ut[k] >= 1 && ut[k] <= p)
                istgt[ut[k] - 1] = 1;
    }

    /* rcar() never adds or removes an arc, only reverses, so the arc count
       is invariant and one allocation covers every iteration */
    int nmax = d->nedges > 0 ? d->nedges : 1;
    int *cfrom = (int *) R_alloc((size_t) nmax, sizeof(int));
    int *cto   = (int *) R_alloc((size_t) nmax, sizeof(int));

    /* --- the two early returns, still with zero RNG contact (rule 1) --- */
    if (d->nedges == 0) {
        vmaxset(vmax);
        return ScalarInteger(0);
    }
    if (covered_arcs(d, istgt, cfrom, cto) == 0) {
        vmaxset(vmax);
        return ScalarInteger(0);
    }

    GetRNGstate();

    int rr = (int) R_unif_index((double) rlen);      /* sample(0:r, size=1) */

    for (int k = 0; k < rr; k++) {
        /* recomputed every iteration, in edgeMatrix order (rule 4) */
        int nce = covered_arcs(d, istgt, cfrom, cto);
        if (nce == 0) {
            /* unreachable: reversing a covered arc u -> w leaves w -> u
               covered, since pa'(w) = pa(w) \ {u} = pa(u) = pa'(u) \ {w},
               and the reversal keeps both endpoints so the target filter is
               unchanged. R would error here out of sample.int(0, 1); guard
               it anyway, closing the RNG region first (rule 3). */
            PutRNGstate();
            vmaxset(vmax);
            error("C_dag_rcar: no covered arc remains after %d reversal(s)", k);
        }
        int pick = (int) R_unif_index((double) nce); /* resample(which(.), 1) */
        /* a covered arc cannot introduce a cycle, which is why the R code
           reverses without an acyclicity test and this does too */
        idl_dag_reverse_edge(d, cfrom[pick], cto[pick]);
    }

    PutRNGstate();
    vmaxset(vmax);

    return ScalarInteger(rr);
}
