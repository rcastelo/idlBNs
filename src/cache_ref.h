#ifndef IDLBNS_CACHE_REF_H
#define IDLBNS_CACHE_REF_H

#include <R.h>
#include <Rinternals.h>
#include "cache_key.h"
#include "sccache.h"

/*
 * ONE INTERFACE OVER TWO CACHE BACKENDS
 *
 * The 'cached.scores' argument of iBIC() and iBGe() is a DOCUMENTED public
 * argument, with a copy-pasteable constructor in man/iBIC.Rd telling users
 * to build a list of environments themselves, and two regression tests
 * compare the key sets of those environments across the C and R engines. So
 * it cannot be replaced.
 *
 * Instead it becomes polymorphic. The score drivers accept either form:
 *
 *   R_NilValue   no caching
 *   VECSXP       the documented list of p environments (cache_key.h)
 *   EXTPTRSXP    an idl_sc_cache, used only inside a search
 *
 * so the search can hand its own compiled cache straight through the
 * existing nh.scores.fun protocol without changing a signature, while a
 * direct iBIC(cached.scores = <list of environments>) call behaves exactly
 * as documented.
 *
 * The two backends must agree on the KEY, not merely produce correct scores:
 * the cached value for a parent set is whichever order-variant of that set
 * missed first, so an identical key sequence is what makes the two caches
 * store identical doubles. Both key on the parent set in ascending order --
 * the R backend by sorting it into a decimal string, this one by hashing the
 * sorted integers -- and tests/test_c_cache.R pins that they end a search
 * holding the same keys with the same values.
 */

#define IDL_CACHE_NONE 0
#define IDL_CACHE_ENVS 1
#define IDL_CACHE_HASH 2

typedef struct {
    int           kind;
    SEXP          envs;         /* IDL_CACHE_ENVS: the VECSXP             */
    idl_sc_cache *hash;         /* IDL_CACHE_HASH: the compiled cache     */
    /* carried between get and put, for whichever backend is in use */
    SEXP          sym;          /* ENVS: the interned key symbol          */
    idl_sc_cursor cur;          /* HASH: the probe position               */
} idl_cache_ref;

/* the tag installed on an idl_sc_cache external pointer; set in R_init */
extern SEXP idlBNs_sccache_tag;

/*
 * Classify 'cached_scores_R' and validate it against p. Errors on anything
 * that is neither NULL, nor a length-p list, nor a live cache pointer.
 */
static inline idl_cache_ref
idl_cache_resolve(SEXP cached_scores_R, int p, const char *who) {
    idl_cache_ref r;
    r.kind = IDL_CACHE_NONE;
    r.envs = R_NilValue;
    r.hash = NULL;
    r.sym = R_NilValue;

    if (cached_scores_R == R_NilValue)
        return r;

    if (TYPEOF(cached_scores_R) == EXTPTRSXP) {
        if (R_ExternalPtrTag(cached_scores_R) != idlBNs_sccache_tag)
            error("%s: 'cached_scores' external pointer has the wrong tag",
                  who);
        idl_sc_cache *c = (idl_sc_cache *) R_ExternalPtrAddr(cached_scores_R);
        if (c == NULL)
            error("%s: 'cached_scores' object is stale (finalized or deserialized)",
                  who);
        if (c->p != p)
            error("%s: 'cached_scores' is for %d vertices, not %d",
                  who, c->p, p);
        r.kind = IDL_CACHE_HASH;
        r.hash = c;

        return r;
    }

    if (TYPEOF(cached_scores_R) != VECSXP || LENGTH(cached_scores_R) != p)
        error("%s: 'cached_scores' must be NULL, a list of length %d, or a score cache object",
              who, p);
    r.kind = IDL_CACHE_ENVS;
    r.envs = cached_scores_R;

    return r;
}

/*
 * Look up the score of vertex node0 (0-based) for the parent set 'key'.
 * 'key' must be ASCENDING and duplicate-free for either backend: the hash
 * backend keys on it directly, and the environment backend sorts it again
 * into the documented ":"-joined string, which is a no-op on sorted input.
 * Returns 1 on a hit.
 */
static inline int
idl_cache_get(idl_cache_ref *r, int node0, const int *key, int klen,
              double *out) {
    switch (r->kind) {
    case IDL_CACHE_HASH:
        return idl_sc_get(r->hash, node0, key, klen, &r->cur, out);
    case IDL_CACHE_ENVS:
        return cache_lookup(VECTOR_ELT(r->envs, node0), key, klen, &r->sym,
                            out);
    default:
        return 0;
    }
}

static inline void
idl_cache_put(idl_cache_ref *r, int node0, const int *key, int klen,
              double s) {
    switch (r->kind) {
    case IDL_CACHE_HASH:
        idl_sc_put(r->hash, node0, key, klen, &r->cur, s);
        break;
    case IDL_CACHE_ENVS:
        cache_store(VECTOR_ELT(r->envs, node0), r->sym, s);
        break;
    default:
        break;
    }
}

#endif /* IDLBNS_CACHE_REF_H */
