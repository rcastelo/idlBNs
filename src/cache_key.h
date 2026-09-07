#ifndef IDLBNS_CACHE_KEY_H
#define IDLBNS_CACHE_KEY_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

/* Shared by ibic_score.c and ibge_score.c: both iBIC() and iBGe() cache
   node scores keyed by the same format as the R function
   .cached_scores_key(), i.e. paste(sort.int(paset), collapse=":"), or the
   literal ":" when paset is empty. Kept in one place so the key format
   has exactly one implementation. */

/* ascending comparator for qsort(), avoiding the classic overflow-prone
   "return *a - *b" pattern */
static inline int
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
 * R_alloc "vmax" stack with vmaxset() after each node).
 *
 * Assumes pa[] contains no NA_INTEGER values (guaranteed by construction
 * of 'pasets' on the R side; not defensively checked here).
 */
static inline char *
build_cache_key(const int *pa, int lp) {
    if (lp == 0) {
        char* buf = (char *) R_alloc(2, sizeof(char));
        buf[0] = ':';
        buf[1] = '\0';

        return buf; 
    }

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
 * cache_lookup
 *
 * Looks up the cached score of one node for the parent set pa[0..lp-1] in
 * the per-node environment 'env'. On a hit, writes the score to *out and
 * returns 1; on a miss, returns 0. Either way *sym receives the interned
 * symbol for the parent set's key, so a following cache_store() reuses it
 * instead of rebuilding and re-interning the key.
 *
 * R_existsVarInFrame()/R_getVar() (envir.c) are used rather than the
 * legacy-only Rf_findVarInFrame() (declared in Rinternals.h only under
 * ENABLE_LEGACY_NONAPI_FUNS, and not part of the default package-facing C
 * API); inherits=FALSE matches the single-frame (parent=emptyenv())
 * semantics of the R-level cached.scores[[i]][[k]] lookup.
 */
static inline int
cache_lookup(SEXP env, const int *pa, int lp, SEXP *sym, double *out) {
    char *key = build_cache_key(pa, lp);
    *sym = Rf_install(key); /* symbols are GC-safe unprotected */
    if (R_existsVarInFrame(env, *sym)) {
        *out = REAL(R_getVar(*sym, env, FALSE))[0];

        return 1;
    }

    return 0;
}

/*
 * cache_store
 *
 * Writes the score 's' into the per-node environment 'env' under the
 * symbol 'sym' returned by a preceding cache_lookup(). The Rf_ScalarReal()
 * result is protected before Rf_defineVar(): nesting an allocating call as
 * a sibling argument to another would leave it exposed to GC under C's
 * unspecified argument-evaluation order.
 */
static inline void
cache_store(SEXP env, SEXP sym, double s) {
    SEXP val = PROTECT(Rf_ScalarReal(s));
    Rf_defineVar(sym, val, env);
    UNPROTECT(1);
}

#endif /* IDLBNS_CACHE_KEY_H */
