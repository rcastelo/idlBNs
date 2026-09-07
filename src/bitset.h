#ifndef IDLBNS_BITSET_H
#define IDLBNS_BITSET_H

#include <stdint.h>
#include <string.h>

/*
 * MULTI-WORD BITSETS
 *
 * The DAG's reachability relation is dense -- it is a transitive closure --
 * so a bitset is the right representation for it, and it is a large memory
 * win over the p x p R logical matrix the search used to carry (4 bytes per
 * entry): at p = 500 the ancestor and descendant matrices together take
 * 64 KB against 1 MB for the single R matrix.
 *
 * Every row is W = ceil(p / 64) words and W is a runtime value held in the
 * DAG struct, so there is NO 64-vertex limit anywhere. That limit was an
 * artifact of sketching this with a single uint64_t per row.
 *
 * Parent and child lists are deliberately NOT bitsets -- they are sparse,
 * and their cost has to follow the parent-set size rather than p. They stay
 * int arrays (see idl_ivec in dag.h).
 */

typedef uint64_t idl_word;

#define IDL_WBITS      64
#define IDL_NW(p)      (((size_t) (p) + 63u) >> 6)
#define IDL_WBYTES(W)  ((W) * sizeof(idl_word))

/* count trailing zeros of a non-zero word */
static inline int
idl_ctz64(idl_word x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll((unsigned long long) x);
#else
    /* portable fallback: binary search on the low bits */
    int n = 0;
    if ((x & 0xFFFFFFFFu) == 0) { n += 32; x >>= 32; }
    if ((x & 0x0000FFFFu) == 0) { n += 16; x >>= 16; }
    if ((x & 0x000000FFu) == 0) { n +=  8; x >>=  8; }
    if ((x & 0x0000000Fu) == 0) { n +=  4; x >>=  4; }
    if ((x & 0x00000003u) == 0) { n +=  2; x >>=  2; }
    if ((x & 0x00000001u) == 0) { n +=  1; }
    return n;
#endif
}

static inline int
idl_popcount64(idl_word x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll((unsigned long long) x);
#else
    int n = 0;
    while (x) { x &= x - 1; n++; }
    return n;
#endif
}

static inline int
idl_bs_test(const idl_word *b, int i) {
    return (int) ((b[(unsigned) i >> 6] >> ((unsigned) i & 63u)) & 1u);
}

static inline void
idl_bs_set(idl_word *b, int i) {
    b[(unsigned) i >> 6] |= (idl_word) 1 << ((unsigned) i & 63u);
}

static inline void
idl_bs_clear(idl_word *b, int i) {
    b[(unsigned) i >> 6] &= ~((idl_word) 1 << ((unsigned) i & 63u));
}

static inline void
idl_bs_zero(idl_word *b, size_t W) {
    memset(b, 0, IDL_WBYTES(W));
}

static inline void
idl_bs_copy(idl_word *dst, const idl_word *src, size_t W) {
    memcpy(dst, src, IDL_WBYTES(W));
}

static inline void
idl_bs_or(idl_word *dst, const idl_word *src, size_t W) {
    for (size_t k = 0; k < W; k++)
        dst[k] |= src[k];
}

static inline int
idl_bs_eq(const idl_word *a, const idl_word *b, size_t W) {
    return memcmp(a, b, IDL_WBYTES(W)) == 0;
}

static inline int
idl_bs_any(const idl_word *b, size_t W) {
    for (size_t k = 0; k < W; k++)
        if (b[k])
            return 1;

    return 0;
}

static inline int
idl_bs_count(const idl_word *b, size_t W) {
    int n = 0;
    for (size_t k = 0; k < W; k++)
        n += idl_popcount64(b[k]);

    return n;
}

/*
 * idl_bs_to_ints
 *
 * Writes the indices of the set bits of 'b', in ASCENDING order, into
 * out[] (which must have room for p entries), and returns how many there
 * were. Used to walk a descendant or ancestor set; the ascending order is
 * incidental to the callers, which all treat these as sets.
 */
static inline int
idl_bs_to_ints(const idl_word *b, size_t W, int *out) {
    int n = 0;
    for (size_t k = 0; k < W; k++) {
        idl_word x = b[k];
        while (x) {
            out[n++] = (int) (k * IDL_WBITS) + idl_ctz64(x);
            x &= x - 1;          /* clear the lowest set bit */
        }
    }

    return n;
}

#endif /* IDLBNS_BITSET_H */
