#ifndef IDLBNS_SCCACHE_H
#define IDLBNS_SCCACHE_H

#include <stdint.h>
#include <stddef.h>

/*
 * PER-VERTEX NODE-SCORE CACHE
 *
 * Replaces the R-environment cache for the duration of one search. What the
 * environments cost, per lookup, is: qsort the parent set, snprintf it into
 * a decimal string, Rf_install() that string, then R_existsVarInFrame() plus
 * R_getVar(). Two of those are the problem --
 *
 *   snprintf     builds a string whose length grows with the parent set, at
 *                a large constant, purely to be thrown away; and
 *   Rf_install   interns the string in R's GLOBAL symbol table, which is a
 *                fixed-bucket hash that is never garbage collected. A long
 *                search at p = 500 interns millions of symbols that then
 *                survive for the rest of the session and degrade every
 *                later Rf_install in the process.
 *
 * This keys on the parent vector's INTEGERS instead: hashed in O(|pa|),
 * compared in O(|pa|), stored in O(|pa|). Nothing about it scales with p,
 * which was the whole point of not using a p-bit mask as the key. And it is
 * freed when the search ends.
 *
 * ONE INVARIANT MATTERS ABOVE PERFORMANCE. A node score is a function of the
 * parent SEQUENCE -- iBIC_node_score() builds ZtZ with the columns in the
 * given order and Choleskys it, and a symmetrically permuted ZtZ rounds
 * differently (measured: 1.14e-13 on an iBIC total for parents c(1,2,3,4,5)
 * against c(5,4,3,2,1)). But the key is the parent SET. So whichever
 * order-variant of a set misses first is the value every later hit returns,
 * and the cache is therefore part of the arithmetic, not a transparent
 * memoisation of it.
 *
 * Consequences, all of which this implementation obeys:
 *
 *   - no eviction, no capacity cap, no lazy insertion. Any of those would
 *     change which variant is cached and move the search trajectory.
 *   - created empty at the start of a search and destroyed at the end, the
 *     same lifetime the R environments had.
 *   - the caller must issue the same lookups, for the same keys, in the same
 *     order as before: the p base terms for i = 0..p-1 first, then the
 *     candidates in neighbourhood order, and within a reversal the head
 *     before the tail.
 *
 * Given the same key sequence and no eviction, the hit/miss sequence is
 * identical to the R cache's, so the stored doubles are too -- which
 * tests/test_c_cache.R checks by comparing the two caches' full contents
 * after a search, keys and values.
 */

#define IDL_SC_EMPTY (-1)

typedef struct {
    uint64_t h;          /* full 64-bit hash: a colliding probe costs one
                            integer compare, not a memcmp                  */
    uint32_t koff;       /* key offset into the shared arena, in ints      */
    int32_t  klen;       /* |pa|, or IDL_SC_EMPTY for a free slot          */
    double   s;          /* the cached node score                          */
} idl_sc_slot;

typedef struct {
    idl_sc_slot *slot;
    uint32_t     cap;    /* power of two */
    uint32_t     mask;   /* cap - 1      */
    uint32_t     n;      /* live entries */
} idl_sc_tab;

typedef struct idl_sc_cache {
    int         p;
    idl_sc_tab *tab;     /* p tables, mirroring cached.scores one for one  */
    int        *arena;   /* bump arena of key ints, shared by all p tables */
    size_t      an, acap;
    uint64_t    hits, misses, probes;      /* diagnostics */

    /*
     * THE ADDITION MEMO.
     *
     * Even with the hash above, the dominant per-step cost at large p is the
     * cache lookup itself: measured 26.6 ns per candidate at p = 400 against
     * a 2.1 ns floor for a bare array read, i.e. the lookup is about 92% of
     * the work. And at p = 400 with 800 arcs, 154,000 of the 155,764
     * candidates are ADDITIONS -- removals and reversals number one per arc
     * each. So memoising additions alone captures 99% of the candidates for
     * p^2 entries rather than 3 p^2.
     *
     * add_d[u * p + v] is the score delta of adding u -> v, and add_st the
     * value of the DAG's pav_stamp[v] when it was computed. The delta
     * depends only on pa(v) -- through both s(v, pa(v) + u) and base[v] --
     * so a single stamp compare decides validity, with no invalidation
     * sweep. 0 means "never computed", and pav_stamp starts at 1.
     *
     * Allocated lazily, on the first neighbourhood scored with stamps
     * supplied, so a plain iBIC(cached.scores = ...) call never pays for it.
     * 12 bytes per (u, v): 3 MB at p = 500.
     */
    double     *add_d;
    int        *add_st;
    uint64_t    memo_hits, memo_misses;
} idl_sc_cache;

/* allocate the addition memo if it is not there yet */
void idl_sc_memo_enable(idl_sc_cache *c);

/* a probe position handed from get() to put(), so put() need not re-hash */
typedef struct {
    uint64_t h;
    uint32_t slot;
} idl_sc_cursor;

/* lifecycle. idl_sc_new() uses R_Calloc and so raises an R error rather than
   returning NULL; the caller must already own the result through an external
   pointer with a finalizer before it allocates anything further. */
idl_sc_cache *idl_sc_alloc(int p);
void          idl_sc_init(idl_sc_cache *c);
void          idl_sc_free(idl_sc_cache *c);

/*
 * idl_sc_get / idl_sc_put
 *
 * 'key' must be the parent set ASCENDING and duplicate-free -- the same
 * canonical form the R cache's sorted string key encodes. get() returns 1 on
 * a hit (writing *out) and 0 on a miss, and either way leaves a cursor that
 * put() consumes.
 *
 * get() also handles growth, BEFORE probing, so that the cursor it returns
 * is still valid when put() runs. Growing inside put() would invalidate it.
 */
int  idl_sc_get(idl_sc_cache *c, int node0, const int *key, int klen,
                idl_sc_cursor *cur, double *out);
void idl_sc_put(idl_sc_cache *c, int node0, const int *key, int klen,
                const idl_sc_cursor *cur, double s);

#endif /* IDLBNS_SCCACHE_H */
