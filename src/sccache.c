#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "sccache.h"
#include "cache_key.h"   /* build_cache_key, for the R-facing dump */

/* the tag carried by an idl_sc_cache external pointer, installed from
   R_init_idlBNs so that every entry point can verify what it was handed */
SEXP idlBNs_sccache_tag = NULL;

void
idl_sc_R_init(void) {
    idlBNs_sccache_tag = Rf_install("idlBNs_sccache");
}

/* initial slots per vertex. most vertices keep small parent sets, so start
   small and let the ones that matter grow. */
#define IDL_SC_CAP0      16u
/* grow at 0.6 rather than the usual 0.75: the search is miss-heavy early,
   and linear probing costs ~2.0 probes on a hit but ~3.6 on a miss at 0.6,
   against ~3.0 and ~8.5 at 0.75 */
#define IDL_SC_NUM       3u
#define IDL_SC_DEN       5u
#define IDL_SC_ARENA0    1024u

/*
 * idl_sc_hash
 *
 * FNV-1a over the key INTEGERS -- not over their decimal text, which is what
 * made the R key expensive -- with a splitmix64 finaliser. The finaliser is
 * not decoration: FNV-1a's low bits are weak, and linear probing consumes
 * exactly the low bits. klen is folded in so that the empty set does not
 * collide by construction with anything else.
 */
static inline uint64_t
idl_sc_hash(const int *key, int klen) {
    uint64_t h = 0xcbf29ce484222325ULL;              /* FNV-1a 64 basis */
    for (int i = 0; i < klen; i++) {
        h ^= (uint64_t) (uint32_t) key[i];
        h *= 0x00000100000001b3ULL;                  /* FNV-1a 64 prime */
    }
    h ^= (uint64_t) (uint32_t) klen;
    h *= 0x00000100000001b3ULL;
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL;        /* splitmix64 */
    h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;

    return h;
}

idl_sc_cache *
idl_sc_alloc(int p) {
    if (p < 1)
        error("idl_sc_alloc: 'p' must be positive");
    idl_sc_cache *c = R_Calloc(1, idl_sc_cache);
    c->p = p;

    return c;
}

void
idl_sc_init(idl_sc_cache *c) {
    c->tab = R_Calloc((size_t) c->p, idl_sc_tab);
    for (int i = 0; i < c->p; i++) {
        c->tab[i].cap = IDL_SC_CAP0;
        c->tab[i].mask = IDL_SC_CAP0 - 1u;
        c->tab[i].n = 0u;
        c->tab[i].slot = R_Calloc((size_t) IDL_SC_CAP0, idl_sc_slot);
        for (uint32_t k = 0; k < IDL_SC_CAP0; k++)
            c->tab[i].slot[k].klen = IDL_SC_EMPTY;
    }
    c->acap = IDL_SC_ARENA0;
    c->an = 0;
    c->arena = R_Calloc(c->acap, int);
}

void
idl_sc_memo_enable(idl_sc_cache *c) {
    if (c->add_d != NULL)
        return;
    size_t n = (size_t) c->p * (size_t) c->p;
    c->add_d = R_Calloc(n, double);
    c->add_st = R_Calloc(n, int);       /* 0 = never computed */
}

void
idl_sc_free(idl_sc_cache *c) {
    if (c == NULL)
        return;
    if (c->tab != NULL)
        for (int i = 0; i < c->p; i++)
            R_Free(c->tab[i].slot);
    R_Free(c->tab);
    R_Free(c->arena);
    R_Free(c->add_d);
    R_Free(c->add_st);
    R_Free(c);
}

/* double the table and re-insert every live slot. no key bytes are touched:
   the stored hash is enough to place an entry, and all keys are distinct so
   the first empty slot always wins. */
static void
tab_grow(idl_sc_tab *t) {
    uint32_t ocap = t->cap;
    idl_sc_slot *old = t->slot;
    uint32_t ncap = ocap << 1;

    if (ncap < ocap)                       /* uint32 overflow */
        error("idl_sc: score cache table too large");

    idl_sc_slot *fresh = R_Calloc((size_t) ncap, idl_sc_slot);
    for (uint32_t k = 0; k < ncap; k++)
        fresh[k].klen = IDL_SC_EMPTY;
    for (uint32_t k = 0; k < ocap; k++) {
        if (old[k].klen == IDL_SC_EMPTY)
            continue;
        uint32_t i = (uint32_t) (old[k].h & (ncap - 1u));
        while (fresh[i].klen != IDL_SC_EMPTY)
            i = (i + 1u) & (ncap - 1u);
        fresh[i] = old[k];
    }
    t->slot = fresh;
    t->cap = ncap;
    t->mask = ncap - 1u;
    R_Free(old);
}

static uint32_t
arena_push(idl_sc_cache *c, const int *key, int klen) {
    if (c->an + (size_t) klen > c->acap) {
        size_t ncap = c->acap;
        while (ncap < c->an + (size_t) klen)
            ncap <<= 1;
        if (ncap > 0xFFFFFFFFu)            /* koff is a uint32 offset */
            error("idl_sc: score cache key arena exceeded 4G entries");
        c->arena = R_Realloc(c->arena, ncap, int);
        c->acap = ncap;
    }
    uint32_t off = (uint32_t) c->an;
    if (klen > 0)
        memcpy(c->arena + c->an, key, (size_t) klen * sizeof(int));
    c->an += (size_t) klen;

    return off;
}

int
idl_sc_get(idl_sc_cache *c, int node0, const int *key, int klen,
           idl_sc_cursor *cur, double *out) {
    idl_sc_tab *t = &c->tab[node0];

    /* grow before probing, so that cur->slot stays valid for idl_sc_put */
    if ((uint64_t) (t->n + 1u) * IDL_SC_DEN > (uint64_t) t->cap * IDL_SC_NUM)
        tab_grow(t);

    uint64_t h = idl_sc_hash(key, klen);
    uint32_t i = (uint32_t) (h & t->mask);
    while (t->slot[i].klen != IDL_SC_EMPTY) {
        c->probes++;
        if (t->slot[i].h == h && t->slot[i].klen == klen &&
            (klen == 0 ||
             memcmp(c->arena + t->slot[i].koff, key,
                    (size_t) klen * sizeof(int)) == 0)) {
            *out = t->slot[i].s;
            c->hits++;

            return 1;
        }
        i = (i + 1u) & t->mask;
    }
    cur->h = h;
    cur->slot = i;
    c->misses++;

    return 0;
}

void
idl_sc_put(idl_sc_cache *c, int node0, const int *key, int klen,
           const idl_sc_cursor *cur, double s) {
    idl_sc_tab *t = &c->tab[node0];
    idl_sc_slot *sl = &t->slot[cur->slot];

    sl->h = cur->h;
    sl->koff = arena_push(c, key, klen);
    sl->klen = klen;
    sl->s = s;
    t->n++;
}

/* ------------------------------------------------------------------------ */
/* R interface                                                              */
/* ------------------------------------------------------------------------ */

static void
idl_sc_finalizer(SEXP ext) {
    idl_sc_cache *c = (idl_sc_cache *) R_ExternalPtrAddr(ext);
    if (c == NULL)
        return;
    idl_sc_free(c);
    R_ClearExternalPtr(ext);
}

/*
 * C_sccache_new
 *
 * Allocate an empty score cache for p vertices. Construction is in the same
 * two steps as the DAG's, and for the same reason: R_Calloc longjmps on
 * failure, so the object has to be owned by R -- external pointer wrapped
 * and finalizer registered -- before its interior is allocated.
 */
SEXP
C_sccache_new(SEXP p_R) {
    int p = asInteger(p_R);

    if (p == NA_INTEGER || p < 1)
        error("C_sccache_new: 'p' must be a positive integer");

    idl_sc_cache *c = idl_sc_alloc(p);
    SEXP ext = PROTECT(R_MakeExternalPtr(c, idlBNs_sccache_tag, R_NilValue));
    R_RegisterCFinalizerEx(ext, idl_sc_finalizer, TRUE);
    idl_sc_init(c);
    UNPROTECT(1);

    return ext;
}

static idl_sc_cache *
sccache_from_extptr(SEXP x) {
    if (TYPEOF(x) != EXTPTRSXP || R_ExternalPtrTag(x) != idlBNs_sccache_tag)
        error("not an idlBNs score cache object");
    idl_sc_cache *c = (idl_sc_cache *) R_ExternalPtrAddr(x);
    if (c == NULL)
        error("idlBNs score cache object is stale (finalized or deserialized)");

    return c;
}

/*
 * C_sccache_dump
 *
 * The whole cache as a list of p named numeric vectors, the names being the
 * same ":"-joined sorted parent-set keys the R-environment backend uses (and
 * ":" for the empty set). That is what lets tests/test_c_cache.R compare the
 * two backends' contents after a search directly against ls()/mget() on the
 * environments -- same keys, and bit-identical values.
 */
SEXP
C_sccache_dump(SEXP st) {
    idl_sc_cache *c = sccache_from_extptr(st);
    SEXP ans = PROTECT(allocVector(VECSXP, c->p));

    for (int i = 0; i < c->p; i++) {
        idl_sc_tab *t = &c->tab[i];
        SEXP v = PROTECT(allocVector(REALSXP, t->n));
        SEXP nm = PROTECT(allocVector(STRSXP, t->n));
        int m = 0;
        for (uint32_t k = 0; k < t->cap; k++) {
            if (t->slot[k].klen == IDL_SC_EMPTY)
                continue;
            void *vmax = vmaxget();
            REAL(v)[m] = t->slot[k].s;
            SET_STRING_ELT(nm, m,
                           mkChar(build_cache_key(c->arena + t->slot[k].koff,
                                                  t->slot[k].klen)));
            vmaxset(vmax);
            m++;
        }
        setAttrib(v, R_NamesSymbol, nm);
        SET_VECTOR_ELT(ans, i, v);
        UNPROTECT(2);
    }
    UNPROTECT(1);

    return ans;
}

/* hits / misses / probes, for the benchmark and for confirming the cache is
   actually being exercised rather than silently bypassed */
SEXP
C_sccache_stats(SEXP st) {
    idl_sc_cache *c = sccache_from_extptr(st);
    SEXP ans = PROTECT(allocVector(REALSXP, 6));
    SEXP nm = PROTECT(allocVector(STRSXP, 6));
    uint64_t entries = 0;
    for (int i = 0; i < c->p; i++)
        entries += c->tab[i].n;
    REAL(ans)[0] = (double) c->hits;
    REAL(ans)[1] = (double) c->misses;
    REAL(ans)[2] = (double) c->probes;
    REAL(ans)[3] = (double) entries;
    REAL(ans)[4] = (double) c->memo_hits;
    REAL(ans)[5] = (double) c->memo_misses;
    SET_STRING_ELT(nm, 0, mkChar("hits"));
    SET_STRING_ELT(nm, 1, mkChar("misses"));
    SET_STRING_ELT(nm, 2, mkChar("probes"));
    SET_STRING_ELT(nm, 3, mkChar("entries"));
    SET_STRING_ELT(nm, 4, mkChar("memo.hits"));
    SET_STRING_ELT(nm, 5, mkChar("memo.misses"));
    setAttrib(ans, R_NamesSymbol, nm);
    UNPROTECT(2);

    return ans;
}
